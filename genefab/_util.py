from sys import stderr
from urllib.request import urlopen
from json import loads
from os.path import join, isdir, isfile
from os import makedirs, remove
from requests import get
from requests.exceptions import InvalidSchema
from urllib.error import URLError
from re import sub
from ._exceptions import GeneLabException

GENELAB_ROOT = "https://genelab-data.ndc.nasa.gov"
API_ROOT = "https://genelab-data.ndc.nasa.gov/genelab"
DELIM_AS_IS = "as.is"
DELIM_DEFAULT = "-"


def get_json(url, verbose=False):
    """HTTP get, decode, parse"""
    if verbose:
        print("Parsing url:", url, file=stderr)
    with urlopen(url) as response:
        return loads(response.read().decode())


def fetch_file(file_name, url, target_directory, update=False, verbose=False, http_fallback=True):
    """Perform checks, download file"""
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


def to_cls(dataframe, space_sub=lambda s: sub(r'\s', "", s)):
    """Convert a presumed annotation/factor dataframe to CLS format"""
    sample_count, factor_count = dataframe.shape
    if factor_count != 1:
        raise GeneLabException("CLS invalid with multiple factors")
    factor_name = dataframe.columns[0]
    classes = dataframe[factor_name].unique()
    if space_sub is None:
        space_sub = lambda s: s
    cls_data = [
        [sample_count, len(classes), 1],
        ["# " + space_sub(classes[0])] + [space_sub(c) for c in classes[1:]],
        [space_sub(v) for v in dataframe[factor_name]]
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
