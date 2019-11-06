# GeneFab

The anatomy of a request:

`{base_URL}/{dataset_accession}/{assay_name}/{data_category}/{data_type}/{transform}/?{get_arguments}`  
Here we assume that `{base_URL}` is "http://genelab-visualization.usra.edu/genefab/".

The request can be routed to any depth (for example, only including
`dataset_accession`, `assay_name`, and `data_category`).  
`get_arguments` are optional (see section **GET arguments**).

### /{dataset_accession}/

Returns information about the dataset and the assays it contains, for example:

type    | name                             | factors          | normalized annotated data file | differential expression analysis data transformation | normalized counts data file | genelab microarray data processing protocol | genelab rnaseq data processing protocol
--------|----------------------------------|------------------|--------------------------------|------------------------------------------------------|-----------------------------|---------------------------------------------|----------------------------------------
dataset | GLDS-4                           | Spaceflight      | NA                             | NA                                                   | NA                          | NA                                          | NA
dataset | GLDS-4                           | Cosmic Radiation | NA                             | NA                                                   | NA                          | NA                                          | NA
assay   | a_GLDS-4_microarray_metadata-txt | Spaceflight      | True                           | True                                                 | False                       | False                                       | False

### /{dataset_accession}/{assay_name}/

Returns assay metadata. If the dataset only has one assay, "assay" can be used
in place of the actual assay name.  
The header of the metadata represents human readable field names (`field`) and
internally used field names (`internal_field`).  
The row names of the metadata correspond to sample names in the assay.

### /{dataset_accession}/{assay_name}/factors/

Returns a table with factor values per sample.  
With `?cls={factor_name}` returns values of one factor for all samples in CLS
format.  
Whether the factor is continuous or discrete is inferred automatically; to
explicitly choose, pass `&continuous=1` for continuous and `&continuous=0` for
discrete.

### /{dataset_accession}/{assay_name}/annotation/

Returns all named annotation information that is different between the samples.
Unnamed fields can be included by passing `?named_only=0`; the fields that are
the same for all samples can be included by passing `?diff=0`.  
Just like with **/factors/**, can be converted to the CLS format with `?cls=`
and controlled with `&continuous=` if necessary.

### /{dataset_accession}/{assay_name}/data/{data_type}/{transform}/

Possible values for `data_type`:  
"processed" returns normalized and annotated counts/array files;  
"deg" returns the analysis table for differentially expressed genes;  
"viz-table" returns the expanded analysis table;  
"pca" returns the results of the principal component analysis.

Possible values for `transform`:  
"gct" (only for "processed") returns processed data in GCT format;  
"melted" returns data melted by the sample name;  
"descriptive" returns data melted by the sample name and described with the
information from **/annotation/**.

## GET arguments

**fmt**: "tsv", "json" (supported everywhere); "html", "raw" (partial support)  
*controls presentation of tables*.

**header**: "0" or "1" (boolean)  
*when set to "1", only outputs the header of the table*.

**showcol**, **hidecol**: 'column_name'  
*only show / hide the specified columns*  
Both showcol and hidecol can be specified multiple times. If showcol is not
present, it is presumed to be set to all column values (show all).

**top**: positive integer value  
*only print the first `top` rows of the table*.

**cls**, **continuous**, **diff**, **named_only**: see above (sections for
**/factors/** and **/annotation/**).

**sort_by**: column name  
*sorts the output table by the values in `sort_by`*.

**ascending**: "0" or "1" (boolean, default "1")  
*when set to "0", reverses the order of sorting (sort order becomes descending)*

**filter**: 'column_name<value' (supported comparison operators are `<`, `<=`,
`>=`, `>`, `==`, `!=`);  
The comparison must be enclosed in single quotes, for example:  
`filter='Adj-p-value-(Space Flight)v(Ground Control)<.05'`.  
filter can be specified multiple times, and all comparisons will be carried
out (with a logical AND between all comparisons).  
Values can be strings, integers, floats, or booleans
(capitalized as True/False).  
*only print the rows that pass the comparison*.

**any_below**: float value between 0 and 1 (supported only for **deg** and
**viz-table**)  
*only print the rows where at least one of the adjusted p-values is below the
specified threshold*.
