#!/usr/bin/env python3
from sys import stderr
from genefab import GLDS

easy_accessions = [
    "GLDS-182",  "GLDS-157",  "GLDS-7",    "GLDS-43",   "GLDS-9",    "GLDS-144",
    "GLDS-15",   "GLDS-208",  "GLDS-158",  "GLDS-42",   "GLDS-40",   "GLDS-160",
    "GLDS-109",  "GLDS-195",  "GLDS-155",  "GLDS-190",  "GLDS-45",   "GLDS-156",
    "GLDS-205",  "GLDS-8",    "GLDS-29",   "GLDS-149",  "GLDS-123",  "GLDS-11",
    "GLDS-17",   "GLDS-174",  "GLDS-176",  "GLDS-80",   "GLDS-188",  "GLDS-134",
    "GLDS-116",  "GLDS-28",   "GLDS-79",   "GLDS-33",   "GLDS-140",  "GLDS-128",
    "GLDS-172",  "GLDS-35",   "GLDS-183",  "GLDS-22",   "GLDS-41",   "GLDS-92",
    "GLDS-63",   "GLDS-71",   "GLDS-175",  "GLDS-131",  "GLDS-130",  "GLDS-165",
    "GLDS-13",   "GLDS-166",  "GLDS-3",    "GLDS-39",   "GLDS-50",   "GLDS-151",
    "GLDS-189",  "GLDS-31",   "GLDS-78",   "GLDS-121",  "GLDS-88",   "GLDS-94",
    "GLDS-44",   "GLDS-129",  "GLDS-32",   "GLDS-117",  "GLDS-76",   "GLDS-89",
    "GLDS-125",  "GLDS-167",  "GLDS-51",   "GLDS-159"
]

# success/failure: GLDS-116

easy_accessions = [
    "GLDS-50"
]

output_mask = "data/{}-{}.tsv.gz"

for accession in easy_accessions:
    glds = GLDS(accession)
    for assay in glds.assays:
        print(assay.metadata.T.to_csv(sep="\t").replace("\t \t", "\tNA\t"))
        #print(assay.fields.keys())
        try:
            if "Comment:  Processed Data Archive file" in assay.fields:
                matrix = assay.get_combined_matrix("Derived Array Data File", "processed")
            elif "Comment:  Derived Array Data File Name" in assay.fields:
                matrix = assay.get_combined_matrix("Comment:  Derived Array Data File Name" )
            elif "Derived Array Data File Name" in assay.fields:
                matrix = assay.get_combined_matrix("Derived Array Data File Name")
            else:
                matrix = assay.get_combined_matrix()
        except:
            print("Something broke for:", accession, assay.name, file=stderr, flush=True)
            continue
        matrix.to_csv(
            output_mask.format(accession, assay.name),
            sep="\t", compression="gzip", index=False
        )
        print("Success:", accession, assay.name, flush=True)
    break
