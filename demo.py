#!/usr/bin/env python3
from genefab import get_datasets

datasets = get_datasets(assay="transcription", maxcount=1000, verbose=True)

header = [
    "#Accession", "AssayName", "Success", "MatrixRows", "MatrixCols", "Error"
]
print(*header, sep="\t")

for dataset in datasets:
    for assay in dataset.assays:
        if assay.has_arrays and assay.available_derived_file_types:
            try:
                matrix = assay.get_combined_matrix()
            except Exception as e:
                success = False
                error = repr(e)
                shape = (-1, -1)
            else:
                success = True
                error = None
                shape = matrix.shape
            fields = [
                dataset.accession, assay.name, success, *shape, error
            ]
            print(*fields, sep="\t", flush=True)
