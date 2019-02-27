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

Additionally, see `demo.py`, which also has code for collecting statistics on
the currently supported and unsupported transcription datasets.

## Description

`get_datasets(maxcount=25, storage=".genelab", **search_kwargs, **additional_kwargs)`:  
Finds up to `maxcount` datasets and sets up assay data storage under `storage`.
Returns a list of `GeneLabDataSet` instances.  
`**search_kwargs` can be `ptype, organism, factor, assay` and allow regexes,
a simple example being `organism="mus", assay="transcription"`.  
`**additional_kwargs` are `verbose`, `onerror`, and `assay_strict_indexing`.  
`verbose` is self-explanatory;  
`onerror` controls what happens if an error occurs when accessing a dataset;
when `onerror="warn"` (default), prints a warning to stderr and continues; if
it is `"ignore"`, continues silently; otherwise raises the error and fails;  
`assay_strict_indexing` (default `False`) is inherited by the instances of
`Assay` inside each `GeneLabDataSet` (see below).

`GeneLabDataSet(accession)` (and its alias `GLDS(accession)`):  
Initializes a GeneLabDataSet instance corresponding to `accession`.  
Additional arguments are `assay_strict_indexing` (see notes for `get_datasets`
and notes for `Assay`), `verbose`, and `storage_prefix` (same as `storage`
in `get_datasets`).  
After successful initialization, contains field `assays`, which is a list
of `Assay` instances corresponding to assays that are referenced by this GLDS.  
Also contains field `factors` that lists experimental factors in the GLDS.

`Assay()`:  
Instances of this class are generated and contained by `GeneLabDataSet`.  
Stores a dataframe with assay metadata in raw form (`Assay.metadata`), but also
has property `fields`, which maps human-readable field names to internal field
names. Using these human-readable field names, metadata can be accessed with
pandas-like syntax: `assay[["Protocol REF"]]`.  
Note the double square brackets: because one human-readable field name can map
to multiple raw fields, this ensures that it will return a DataFrame and not a
Series (which can be ambiguous).  
This behavior can be overridden by setting `Assay.strict_indexing` to `True`
(and, upstream, by calling `get_datasets` and `GeneLabDataSet` with
`assay_strict_indexing=True`), but may lead to conflicts and errors if you're
not extra careful.  
Instances of `Assay` can also be indexed with `.loc[]` and with boolean queries,
just like in pandas.  
`Assay.available_protocols` returns a set with all "Protocol REF" values
referenced by the assay metadata; this is useful for checks like
`if "nucleic acid hybridization" in assay.available_protocols: ...`.  
`Assay.available_file_types` returns a set with all available file types
referenced by the assay metadata.  
`Assay.available_derived_file_types` returns a set with all available derived
file types (such as "Processed Array Data File" and so on).  
`Assay.get_combined_matrix()` attempts to download, parse and combine derived
data into a single DataFrame. If derived data is referenced in the metadata
individually for each sample, returns an instance of `PerSampleMatrix`, which
inherits from `pandas.DataFrame` and acts like one. If derived data is
referenced in combined, complex files, will in the future return an instance
of `CompoundMatrix`, but currently fails with a `NotImplementedError`.

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
