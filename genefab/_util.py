from urllib.request import urlopen
from json import loads
from os.path import join, isfile, isdir, exists
from requests import get
from math import ceil
from tqdm import tqdm
from os import mkdir, remove
from re import sub
from shutil import rmtree, copyfileobj
from gzip import open as gzopen
from urllib.error import URLError
from sys import stderr

URL_ROOT = "https://genelab-data.ndc.nasa.gov/genelab"
DEFAULT_STORAGE = ".genefab"

def get_json(url):
    """HTTP get, decode, parse"""
    print("Parsing url: ", url, file=stderr)
    with urlopen(url) as response:
        return loads(response.read().decode())

def fetch_file(file_name, url, target_directory=DEFAULT_STORAGE, update=False):
    """Perform checks, download file"""
    if not isdir(target_directory):
        if isfile(target_directory):
            raise OSError("Local storage exists and is not a directory")
        mkdir(target_directory)
    target_file = join(target_directory, file_name)
    if not update:
        if isdir(target_file):
            raise OSError("Directory with target name exists: " + target_file)
        if isfile(target_file):
            print("Reusing", file_name, file=stderr)
            return target_file
    stream = get(url, stream=True)
    if stream.status_code != 200:
        raise URLError("{}: status code {}".format(url, stream.status_code))
    total_bytes = int(stream.headers.get("content-length", 0))
    total_kb = ceil(total_bytes / 1024)
    with open(target_file, "wb") as output_handle:
        written_bytes = 0
        stream_iterator = tqdm(
            stream.iter_content(1024), desc="Downloading "+file_name,
            total=total_kb, unit="KB"
        )
        for block in stream_iterator:
            output_handle.write(block)
            written_bytes += len(block)
    if total_bytes != written_bytes:
        remove(target_file)
        raise URLError("Failed to download the correct number of bytes")
    return target_file

def ensure_dir(target_dir, force_new_dir):
    if exists(target_dir):
        if not isdir(target_dir):
            raise OSError("Target name exists and is not a directory")
        elif force_new_dir:
            rmtree(target_dir)
        else:
            raise OSError("Target directory exists")
    mkdir(target_dir)

def gunzip(source_file, target_file=None, target_dir=None, keep_original=False):
    """Alternative to call(["gunzip", ...]), because the latter fails on Windows"""
    if target_file is None:
        target_file = sub(r'\.gz$', "", source_file)
    with gzopen(source_file, "rb") as compressed:
        with open(target_file, "wb") as uncompressed:
            copyfileobj(compressed, uncompressed)
    if target_dir is not None:
        call(["mv", target_file, target_dir])
    if not keep_original:
        remove(source_file)

FFIELD_VALUES = {
    "Project+Type": [
        "Spaceflight Study", "Spaceflight Project", "Spaceflight", "Flight Study", "Flight", "ground", "parabolic"
    ],
    "Study+Factor+Name": [
        "Absorbed Radiation Dose", "Age", "Altitude", "animal housing", "Antibiotic concentration", "Atmospheric Pressure",
        "Bed Rest", "Bleomycin Treatment", "cage", "CANONT:Part", "cell culture", "Cell cycle phase", "Cell Line",
        "clinical treatment", "collection set", "condition", "control group", "culture", "Culture Condition", "culture media",
        "development", "developmental condition", "developmental stage", "Diet", "dissection condition", "Dissection Timeline",
        "Donor", "dose", "ecotype", "EFO:light", "Electromagnetic Fields", "environment exposure", "Environmental Stress",
        "environmentalstress", "Exercise", "exposure duration", "food deprivation", "Fractionated Dose", "Freezing",
        "freezing profile", "Gender", "generation", "generation number", "genotype", "gravitation", "gravity",
        "gravity type", "Gravity, Altered", "growth environment", "hindlimb unloading", "hypergravity", "Individual",
        "infection", "Injection", "Ionizing Radiation", "Ionzing Radiation", "irradiate", "Irradiated", "Irradiation",
        "light cycle", "location", "Magnetic field", "MESH:Atmospheric Pressure", "MESH:Gravitation", "Microgravity",
        "mouse strain", "Muscle, Skeletal", "Neoplasm", "Nutrition", "O2", "organism part", "organism_part", "osteo-induced",
        "post radiation timepoint", "Preservation method", "pressure", "protocol host", "protocoltype", "Radiation", "Radiation Distance",
        "Radiation dosage", "radiation dose", "radiation type", "Radiation, Ionzing", "RNA sequencing", "sample collection protocol",
        "sample type", "Sampling time", "Sex", "Simulated microgravity", "Smoking Status", "Space Flight", "spaceflight",
        "Stimulated gravity (g level)", "stimulus", "strain", "strain/genotype", "Stress", "target gene specification",
        "temperature", "time", "time after treatment", "timepoint", "tissue", "tissue storage time", "treated with", "Treatment",
        "Treatment Duration", "treatment group", "treatment time", "variation", "viral load", "water deprivation", "weightlessness",
        "Weightlessness Simulation", "zygosity"
    ],
    "organism": [
        "Acinetobacter pittii", "Arabidopsis thaliana", "Aspergillus fumigatus", "Aspergillus niger", "Aspergillus terreus",
        "Aureobasidium pullulans", "Bacillus", "Bacillus subtilis", "Beauveria bassiana", "Brassica rapa", "Caenorhabditis elegans",
        "Candida albicans", "cellular organisms", "Ceratopteris richardii", "Cladosporium cladosporioides", "Cladosporium sphaerospermum",
        "Danio rerio", "Daphnia magna", "Drosophila melanogaster", "Enterobacter", "Enterobacteria phage lambda", "environmental samples",
        "Escherichia coli", "Euprymna scolopes", "Fusarium solani", "Helix lucorum", "Homo sapiens", "Klebsiella", "metagenomic data",
        "Microbiota", "Mus musculus", "Mycobacterium marinum", "Oryzias latipes", "Pantoea conspicua", "Pseudomonas aeruginosa",
        "Rattus norvegicus", "Rhodospirillum rubrum", "Saccharomyces cerevisiae", "Staphylococcus", "Staphylococcus aureus",
        "Streptococcus mutans", "Trichoderma virens"
    ],
    "Study+Assay+Measurement+Type": [
        "deletion pool profiling", "DNA methylation profiling", "environmental gene survey", "genome sequencing",
        "metabolite profiling", "protein expression profiling", "RNA methylation profiling", "transcription profiling"
    ]
}

FFIELD_ALIASES = {
    "ptype": "Project+Type",
    "factor": "Study+Factor+Name",
    "organism": "organism",
    "assay": "Study+Assay+Measurement+Type"
}
