from sys import stderr
from urllib.request import urlopen
from json import loads
from os.path import join, isdir, isfile
from os import makedirs, remove
from requests import get
from requests.exceptions import InvalidSchema
from urllib.error import URLError
from re import sub, search
from sqlite3 import connect
from hashlib import sha512
from pandas import read_csv
from tempfile import TemporaryDirectory


GENELAB_ROOT = "https://genelab-data.ndc.nasa.gov"
API_ROOT = "https://genelab-data.ndc.nasa.gov/genelab"
DELIM_AS_IS = "as.is"
DELIM_DEFAULT = "-"
STORAGE_PREFIX = ".genelab"
DB_NAME = join(STORAGE_PREFIX, "genefab.db")


def get_json(url, verbose=False):
    """HTTP get, decode, parse"""
    if verbose:
        print("Parsing url:", url, file=stderr)
    with urlopen(url) as response:
        return loads(response.read().decode())


def fetch_file(file_name, url, target_directory, update=False, verbose=False, http_fallback=True):
    """Soon to be deprecated: perform checks, download file"""
    if not isdir(target_directory):
        if isfile(target_directory):
            raise OSError("Local storage exists and is not a directory")
        makedirs(target_directory)
    target_file = join(target_directory, file_name)
    if not update:
        if isdir(target_file):
            raise OSError("Directory with target name exists: " + target_file)
        if isfile(target_file):
            if verbose:
                print("Reusing", file_name, file=stderr)
            return target_file
    try:
        stream = get(url, stream=True)
    except InvalidSchema:
        if http_fallback:
            stream = get(sub(r'^ftp:\/\/', "http://", url), stream=True)
        else:
            raise
    if stream.status_code != 200:
        raise URLError("{}: status code {}".format(url, stream.status_code))
    total_bytes = int(stream.headers.get("content-length", 0))
    with open(target_file, "wb") as output_handle:
        written_bytes = 0
        for block in stream.iter_content(1024):
            output_handle.write(block)
            written_bytes += len(block)
    if total_bytes != written_bytes:
        remove(target_file)
        raise URLError("Failed to download the correct number of bytes")
    return target_file


def guess_format(filemask):
    """Guess whether the file is a CSV or a TSV and whether it is compressed"""
    comp_ext_match = search(r'\.(gz|bz2|xz)$', filemask)
    if comp_ext_match:
        compression = comp_ext_match.group(1)
    else:
        compression = None
    fmt_ext_match = search(r'\.csv', filemask)
    if fmt_ext_match:
        sep = ","
    else:
        sep = "\t"
    return sep, compression


def update_table(filemask, url, table_prefix, verbose=False, http_fallback=True):
    """Check if table already in database, if not, download and store; return table name"""
    db = connect(DB_NAME)
    filemask_hash = sha512(filemask.encode("utf-8")).hexdigest()
    table_name = "{}/{}".format(table_prefix, filemask_hash)
    query_mask = (
        "SELECT count(name) FROM sqlite_master WHERE type='table' AND name='{}'"
    )
    cursor = db.execute(query_mask.format(table_name))
    if cursor.fetchone()[0] > 0: # table exists
        return table_name
    # otherwise, table doesn't exist and we need to download data:
    try:
        stream = get(url, stream=True)
    except InvalidSchema:
        if http_fallback:
            stream = get(sub(r'^ftp:\/\/', "http://", url), stream=True)
        else:
            raise
    if stream.status_code != 200:
        raise URLError("{}: status code {}".format(url, stream.status_code))
    total_bytes = int(stream.headers.get("content-length", 0))
    with TemporaryDirectory() as tempdir:
        sep, compression = guess_format(filemask)
        target_file = join(tempdir, filemask_hash)
        if compression:
            target_file = target_file + "." + compression
        with open(target_file, "wb") as output_handle:
            written_bytes = 0
            for block in stream.iter_content(1024):
                output_handle.write(block)
                written_bytes += len(block)
        if total_bytes != written_bytes:
            remove(target_file)
            raise URLError("Failed to download the correct number of bytes")
        read_csv(target_file, sep=sep).to_sql(table_name, db)
        db.commit()
        db.close()
    return table_name


def to_cls(dataframe, target, continuous="infer", space_sub=lambda s: sub(r'\s', "", s)):
    """Convert a presumed annotation/factor dataframe to CLS format"""
    sample_count = dataframe.shape[0]
    if continuous == "infer":
        try:
            _ = dataframe[target].astype(float)
            continuous = True
        except ValueError:
            continuous = False
    elif not isinstance(continuous, bool):
        if continuous == "0":
            continuous = False
        elif continuous == "1":
            continuous = True
        else:
            error_message = "`continuous` can be either boolean-like or 'infer'"
            raise TypeError(error_message)
    if continuous:
        cls_data = [
            ["#numeric"], ["#" + target],
            dataframe[target].astype(float)
        ]
    else:
        if space_sub is None:
            space_sub = lambda s: s
        classes = dataframe[target].unique()
        cls_data = [
            [sample_count, len(classes), 1],
            ["# "+space_sub(classes[0])] + [space_sub(c) for c in classes[1:]],
            [space_sub(v) for v in dataframe[target]]
        ]
    return "\n".join([
        "\t".join([str(f) for f in fields]) for fields in cls_data
    ])


FFIELD_VALUES = {
    "Project+Type": [
        "Spaceflight Study", "Spaceflight Project", "Spaceflight",
        "Flight Study", "Flight", "ground", "parabolic"
    ],
    "Study+Factor+Name": [
        "Absorbed Radiation Dose", "Age", "Altitude", "animal housing",
        "Antibiotic concentration", "Atmospheric Pressure", "Bed Rest",
        "Bleomycin Treatment", "cage", "CANONT:Part", "cell culture",
        "Cell cycle phase", "Cell Line", "clinical treatment", "collection set",
        "condition", "control group", "culture", "Culture Condition",
        "culture media", "development", "developmental condition",
        "developmental stage", "Diet", "dissection condition",
        "Dissection Timeline", "Donor", "dose", "ecotype", "EFO:light",
        "Electromagnetic Fields", "environment exposure",
        "Environmental Stress", "environmentalstress", "Exercise",
        "exposure duration", "food deprivation", "Fractionated Dose",
        "Freezing", "freezing profile", "Gender", "generation",
        "generation number", "genotype", "gravitation", "gravity",
        "gravity type", "Gravity, Altered", "growth environment",
        "hindlimb unloading", "hypergravity", "Individual", "infection",
        "Injection", "Ionizing Radiation", "Ionzing Radiation", "irradiate",
        "Irradiated", "Irradiation", "light cycle", "location",
        "Magnetic field", "MESH:Atmospheric Pressure", "MESH:Gravitation",
        "Microgravity", "mouse strain", "Muscle, Skeletal", "Neoplasm",
        "Nutrition", "O2", "organism part", "organism_part", "osteo-induced",
        "post radiation timepoint", "Preservation method", "pressure",
        "protocol host", "protocoltype", "Radiation", "Radiation Distance",
        "Radiation dosage", "radiation dose", "radiation type",
        "Radiation, Ionzing", "RNA sequencing", "sample collection protocol",
        "sample type", "Sampling time", "Sex", "Simulated microgravity",
        "Smoking Status", "Space Flight", "spaceflight",
        "Stimulated gravity (g level)", "stimulus", "strain", "strain/genotype",
        "Stress", "target gene specification", "temperature", "time",
        "time after treatment", "timepoint", "tissue", "tissue storage time",
        "treated with", "Treatment", "Treatment Duration", "treatment group",
        "treatment time", "variation", "viral load", "water deprivation",
        "weightlessness", "Weightlessness Simulation", "zygosity"
    ],
    "organism": [
        "Acinetobacter pittii", "Arabidopsis thaliana", "Aspergillus fumigatus",
        "Aspergillus niger", "Aspergillus terreus", "Aureobasidium pullulans",
        "Bacillus", "Bacillus subtilis", "Beauveria bassiana", "Brassica rapa",
        "Caenorhabditis elegans", "Candida albicans", "cellular organisms",
        "Ceratopteris richardii", "Cladosporium cladosporioides",
        "Cladosporium sphaerospermum", "Danio rerio", "Daphnia magna",
        "Drosophila melanogaster", "Enterobacter",
        "Enterobacteria phage lambda", "environmental samples",
        "Escherichia coli", "Euprymna scolopes", "Fusarium solani",
        "Helix lucorum", "Homo sapiens", "Klebsiella", "metagenomic data",
        "Microbiota", "Mus musculus", "Mycobacterium marinum",
        "Oryzias latipes", "Pantoea conspicua", "Pseudomonas aeruginosa",
        "Rattus norvegicus", "Rhodospirillum rubrum",
        "Saccharomyces cerevisiae", "Staphylococcus", "Staphylococcus aureus",
        "Streptococcus mutans", "Trichoderma virens"
    ],
    "Study+Assay+Measurement+Type": [
        "deletion pool profiling", "DNA methylation profiling",
        "environmental gene survey", "genome sequencing",
        "metabolite profiling", "protein expression profiling",
        "RNA methylation profiling", "transcription profiling"
    ]
}


FFIELD_ALIASES = {
    "ptype": "Project+Type",
    "factor": "Study+Factor+Name",
    "organism": "organism",
    "assay": "Study+Assay+Measurement+Type"
}
