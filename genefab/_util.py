from argparse import Namespace
from csv import Sniffer
from copy import deepcopy
from re import sub
from hashlib import sha512
from datetime import datetime
from os import path
from contextlib import closing
from sqlite3 import connect


GENELAB_ROOT = "https://genelab-data.ndc.nasa.gov"
API_ROOT = "https://genelab-data.ndc.nasa.gov/genelab"
DELIM_AS_IS = "as.is"
DELIM_DEFAULT = "-"
STORAGE_PREFIX = ".genelab"
LOG_SCHEMA = "('time' INTEGER, 'url' TEXT, 'ip' TEXT, 'exception' TEXT, 'comment' TEXT)"


DEFAULT_RARGS = Namespace(
    data_rargs = {
        "fields": None,
        "index": None,
        "file_filter": ".*",
        "name_delim": DELIM_DEFAULT,
        "melted": False, # TODO: 'descriptive' conflictable
        "descriptive": False,
        "any_below": None,
    },
    data_filter_rargs = {
        "filter": None,
        "sort_by": None,
        "ascending": True,
    },
    display_rargs = {
        "fmt": "tsv", # TODO: 'raw' conflictable
        "header": False, # TODO: 'top' conflictable
        "top": None,
        "showcol": None,
        "hidecol": None,
    },
    non_data_rargs = {
        "diff": True,
        "named_only": True,
        "cls": None,
        "continuous": "infer",
    }
)


def parse_rargs(request_args):
    """Get all common arguments from request.args"""
    rargs = deepcopy(DEFAULT_RARGS)
    for rarg_type, rargs_of_type in DEFAULT_RARGS.__dict__.items():
        for rarg, rarg_default_value in rargs_of_type.items():
            if rarg in request_args:
                if not isinstance(rarg_default_value, bool):
                    if len(request_args.getlist(rarg)) == 1:
                        getattr(rargs, rarg_type)[rarg] = request_args[rarg]
                    else:
                        getattr(rargs, rarg_type)[rarg] = request_args.getlist(rarg)
                elif request_args[rarg] == "0":
                    getattr(rargs, rarg_type)[rarg] = False
                else:
                    getattr(rargs, rarg_type)[rarg] = True
    return rargs


def data_rargs_digest(accession, assay_name, data_rargs):
    """Convert data_rargs to a string digest"""
    raw_digest = [accession, assay_name]
    for key, value in sorted(data_rargs.items()):
        raw_digest.append(key+"="+str(value))
    raw_string_digest = ",".join(raw_digest)
    string_digest = sub(r'[^0-9A-Za-z_.=-]', ".", raw_string_digest)
    hexdigest = sha512(raw_string_digest.encode("utf-8")).hexdigest()[:32]
    return string_digest + "/" + hexdigest


def guess_format(target_file):
    """Guess whether the file is a CSV or a TSV and whether it is compressed"""
    with open(target_file, mode="rb") as handle:
        magic = handle.read(3)
        if magic == b"\x1f\x8b\x08":
            compression = "gzip"
            from gzip import open as _open
        elif magic == b"\x42\x5a\x68":
            compression = "bz2"
            from bz2 import open as _open
        else:
            compression = None
            _open = open
    with _open(target_file, newline="", mode="rt") as handle:
        dialect = Sniffer().sniff(handle.read(2**20))
        sep = dialect.delimiter
    return sep, compression


def date2stamp(fd, key="date_modified", fallback_key="date_created", fallback_value=-1):
    """Convert date like 'Fri Oct 11 22:02:48 EDT 2019' to timestamp"""
    strdate = fd.get(key)
    if strdate is None:
        strdate = fd.get(fallback_key)
    if strdate is None:
        return fallback_value
    else:
        try:
            dt = datetime.strptime(strdate, "%a %b %d %H:%M:%S %Z %Y")
        except ValueError:
            return fallback_value
        else:
            return int(dt.timestamp())


def log(request, exception):
    """Save exception context to sqlite3 log database"""
    db_name = path.join(STORAGE_PREFIX, "log.sqlite3")
    with closing(connect(db_name)) as db:
        db.cursor().execute("CREATE TABLE IF NOT EXISTS 'log' " + LOG_SCHEMA)
        db.cursor().execute(
            "INSERT INTO 'log' ('time', 'url', 'ip', 'exception', 'comment') " +
            "VALUES({}, '{}', '{}', '{}', '{}')".format(
                int(datetime.timestamp(datetime.now())),
                request.url, request.remote_addr,
                type(exception).__name__,
                sub(r'[^0-9A-Za-z_ ]', "_", str(exception))
            )
        )
        db.commit()


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
