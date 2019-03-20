# GeneFab

## Quick start (demo)

```
from genefab import GeneLabDataSet

```

Under construction lol

## Description

### GeneLabDataSet(accession), and its alias GLDS(accession)

### Assay()

Instances of this class are generated and contained by `GeneLabDataSet`.

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
