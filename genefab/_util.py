from sys import stderr
from urllib.request import urlopen
from json import loads

URL_ROOT = "https://genelab-data.ndc.nasa.gov/genelab"

def get_json(url):
    """HTTP get, decode, parse"""
    print("Parsing url: ", url, file=stderr)
    with urlopen(url) as response:
        return loads(response.read().decode())

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
