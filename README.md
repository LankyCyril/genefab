# GeneFab

## Quick start (demo)

```
from genefab import GeneLabDataSet

glds = GeneLabDataSet("GLDS-158")
assay = glds.assays[0]
matrix = assay.get_combined_matrix()
filename = "GLDS-158-" + assay.name + ".tsv"
matrix.to_csv(filename, sep="\t", index=False)
```

## Description

Implements classes and methods:
* `GeneLabDataSet(accession)` which stores metadata and data files associated
with a GLDS accession number
* `GLDS(accession)`: alias to `GeneLabDataSet`
* `get_datasets(ptype, organism, factor, assay, maxcount=25)` which looks up and
stores all GLDS datasets matching passed search terms

## Requirements

* Python 3.5+

Python packages:
* pandas
* requests
* tqdm

It is recommended to use a user installed
[Conda](https://www.anaconda.com/download/) environment.  
Full requirements can be fulfilled by using the \*.yaml files provided in the
top-level directory of the repo:

```
$ conda env create --name genefab --file environment-linux.yaml
$ source activate genefab
$ python
>>> from genefab import get_datasets, GLDS
```
