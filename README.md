# GeneFab

Work in progress with an admittedly silly name.  
There are **a lot** of corner cases and something is probably very broken.

## Quick start (demo)

(coming up)

## Description

Implements classes and methods:
* `GeneLabDataSet(accession)` which stores metadata and data files associated
with a GLDS accession number
* `GLDS(accession)`: alias to `GeneLabDataSet`
* `get_datasets(ptype, organism, factor, assay, maxcount=25)` which looks up and
stores all GLDS datasets matching passed search terms

## Requirements

* Python 3.5+
* R (must be user-installed so as to allow installation of packages on the fly)

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
