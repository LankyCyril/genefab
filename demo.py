#!/usr/bin/env python3
from genefab import get_datasets

datasets = get_datasets(maxcount=1000, verbose=True)
with_derived, successful = 0, 0
failed_accessions = set()

for dataset in datasets:
    for assay in dataset.assays:
        if assay.available_derived_file_types:
            with_derived += 1
            try:
                matrix = assay.get_combined_matrix()
            except Exception as e:
                print("FAIL:", dataset.accession, assay.name, e)
                failed_accessions.add(dataset.accession)
            else:
                print("SUCCESS:", dataset.accession, assay.name, matrix.shape)
                successful += 1

print("# total datasets found:", len(datasets))
print("# assays with derived data:", with_derived)
print("# successful matrices: {}/{}".format(successful, with_derived))
print("# failed matrices: {}/{}".format(with_derived-successful, with_derived))

print("Problematic datasets:")
print(*sorted(failed_accessions), sep=", ")
