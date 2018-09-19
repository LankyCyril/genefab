# GeneFab

Work in progress with an admittedly silly name.

## Quick start (demo)

...

## Description

Implements classes:
* `GLDS(accession)` which stores metadata and data files associated with a
GLDS accession number
* `GLDSCollection(ptype, organism, factor, assay, maxcount=25)` which looks up
and stores all GLDS datasets matching passed search terms
* `MicroarrayExperiment(glds)` which converts a `GLDS` instance into an object
capable of performing a differential expression analysis based on microarray
data.

## Requirements

* Python 3.5+
* R (must be user-installed so as to allow installation of packages on the fly)

Python packages:
* pandas
* requests
* tqdm

It is recommended to use a user installed [Conda]() environment.  
Full requirements can be fulfilled by using the \*.yaml files provided in the
top-level directory of the repo:

```
$ conda env create --name genefab --file environment-linux.yaml
$ source activate genefab
$ python
>>> from genefab import GLDSCollection, MicroarrayExperiment
```
