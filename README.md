# GeneFab

## Installation

GeneFab has been tested under Python 3.5+.  
It is recommended to use a user installed
[Conda](https://www.anaconda.com/download/) environment.  
Currently GeneFab can be installed via **pip** directly from github:  
`pip install -e git://github.com/LankyCyril/genefab.git#egg=genefab`

For development purposes, it can also be cloned and used like this:

```
$ git clone https://github.com/LankyCyril/genefab
$ cd genefab
$ conda env create --name genefab --file environment-linux.yaml
$ conda activate genefab
$ python
>>> from genefab import get_datasets, GLDS
    ...
```

## Demo

https://github.com/LankyCyril/genefab/blob/master/genefab-demo.ipynb

## Description

The current iteration of GeneFab supports GeneLab datasets processed according
to API version 2.1. As of April 2, 2019, these are only the datasets
GLDS-4, GLDS-30, and GLDS-42.

#### GeneLabDataSet(accession, storage_prefix=".genelab", index_by="Sample Name", verbose=False)

Initializes a GeneLabDataSet instance corresponding to `accession`.  
Each GeneLabDataSet contains references to assays performed in the study; they
are stored under the field `assays` as a dictionary of `Assay` instances.  
The dictionary can be enumerated by keys (assay names) and values
(Assay objects).

Argument `storage_prefix` controls where the assay data is stored. For example,
if `storage_prefix == ".genelab"`, the dataset accession is "GLDS-30", and the
assay name is "a_GLDS-30_microarray_metadata-txt", the data will be stored under
".genelab/GLDS-30/a_GLDS-30_microarray_metadata-txt/", relative to the current
path.  
This path must be relative to the current working directory and does not use
shell expansion (e.g., starting the path with "~" will fail).

Argument `index_by` is inherited by the dataset's assays as well and controls
which field becomes the index for their metadata. Usually, the default setting
("Sample Name") need not be changed. This field accepts regular expressions and
treats them as case-insensitive, full-length matchers. In other words,
".*mple name" would still match "Sample Name" (and maybe something else...), but
"mple name" wouldn't match it.

`GLDS` is an alias to `GeneLabDataSet`.

```python
glds = GeneLabDataSet("GLDS-30", storage="data/genelab")
assay = glds.assays[0]
```

#### Assay()

Objects of this class are created implicitly after initializing a GeneLabDataSet
object and are stored in the `assays` field of the GeneLabDataSet object.  
Data in each assay is indexed by the same field as requested when calling
the GeneLabDataSet constructor (`index_by`, default "Sample Name").

`Assay.samples`: list of samples in the assay.

`Assay.factors`: a dataframe of factor values associated with each sample.

`Assay.normalized_data`: a dataframe of normalized array data for each sample.

`Assay.processed_data`: a dataframe of processed (normalized and annotated) data
for each sample.

`Assay.normalized_annotated_data`: an alias to `Assay.processed_data`.

`Assay.metadata`: a DataFrame-like object that can be interrogated by
human-readable field names. It maps these names to internal field names and
returns DataFrame slices according to requests (see `Assay.metadata.fields` for
more information).  
As human-readable names can map to more than one internal field/column, this
object can only be indexed by double-bracket syntax and can only return
DataFrames (not Series).  
Indexing is done by regular expressions that are treated as case-insensitive,
full-length matchers. For example, `assay.metadata[["extract name"]]` matches
fields exposed as "Extract Name", while `assay.metadata[[".*extract name"]]` can
match both "Extract Name" and "Labeled Extract Name".  
Indexing by sample name is possible: `assay.metadata.loc[sample_name]`  
Simple logical indexing with `.loc` is also possible:  
`space_labels = assay.metadata.loc[assay.factors=="Space Flight", ["Labeled Extract Name"]]`

`Assay.metadata.fields`: a dictionary mapping human-readable field names to
internal field names in the API. For example:
```
>>> assay.metadata.fields["Protocol REF"]
{'a100017protocolref', 'a100009protocolref'}
```

`Assay.metadata.columns`: an alias to `list(Assay.metadata.fields.keys())`.

`Assay.metadata.index`: sample names in the assay.

```python
annotation = assay.factors
data = assay.processed_data

space_samples = assay.metadata[assay.factors=="Space Flight"].index
space_data = assay.processed_data[space_samples]
```
