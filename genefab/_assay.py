from os.path import join
from genefab._exceptions import GeneLabJSONException
from genefab._exceptions import GeneLabException
from collections import defaultdict
from pandas import concat, Series, Index, DataFrame, merge
from re import search, fullmatch, IGNORECASE, sub
from genefab._util import DELIM_AS_IS
from genefab._display import to_cls

ASSAY_CHARACTERISTICS = [
    "normalized annotated data file",
    "differential expression analysis data transformation",
    "normalized counts data file"
]
ASSAY_PROTOCOL_LIST = [
    "genelab microarray data processing protocol",
    "genelab rnaseq data processing protocol"
]


class AssayMetadataLocator():
    """Emulate behavior of Pandas `.loc` for class AssayMetadata()"""

    def __init__(self, parent):
        """Point to parent"""
        self.parent = parent

    def __getitem__(self, key):
        """Query parent.raw_metadata with .loc, using field titles instead of internal field ids"""
        if isinstance(key, tuple): # called with .loc[x, y]
            try:
                index_patterns, title_patterns = key
            except ValueError:
                raise IndexError("Incorrect index for assay metadata")
            else: # assume both indices and titles are collections of regexes
                indices = set.union(*({
                        ix for ix in self.parent.parent.raw_metadata.index
                        if fullmatch(pattern, ix, flags=IGNORECASE)
                    } for pattern in index_patterns
                ))
                row_subset = self.parent.loc[list(indices)]
                field_titles = set.union(*(
                    self.parent.parent._match_field_titles(t, method=fullmatch)
                    for t in title_patterns
                ))
                fields = set.union(*(
                    self.parent.parent._fields[title]
                    for title in field_titles
                ))
                return row_subset[list(fields)]
        else: # assume called with .loc[x] and interpret `x` the best we can
            if isinstance(key, DataFrame) and (key.shape[1] == 1):
                # assume being indexed by boolean column, delegate to parent[]:
                return self.parent[key]
            elif isinstance(key, (tuple, list, set)):
                # assume it is a collection of regexes:
                indices = set.union(*({
                        ix for ix in self.parent.parent.raw_metadata.index
                        if fullmatch(pattern, ix, flags=IGNORECASE)
                    } for pattern in key
                ))
                return self.parent.parent.raw_metadata.loc[list(indices)]
            else: # last resort; just pass it to raw_metadata directly
                return self.parent.parent.raw_metadata.loc[key]


class AssayMetadata():
    """Makes individual assay metadata accessible with Pandas-like indexing"""

    def __init__(self, parent):
        """Point to parent and initialize children"""
        self.parent = parent
        self.loc = AssayMetadataLocator(self)

    def to_frame(self):
        """Raw metadata with multiindex columns (human-readable -> internal)"""
        multicols = ["field", "internal_field"]
        fields_df = DataFrame(
            data=[[k, v] for k, vv in self.parent._fields.items() for v in vv],
            columns=multicols
        )
        columns_df = DataFrame(
            data=self.parent.raw_metadata.columns, columns=["internal_field"]
        )
        multiindex_df = merge(columns_df, fields_df, sort=False, how="outer") \
            .fillna("Unknown")
        mdv = multiindex_df["internal_field"].values
        rmv = self.parent.raw_metadata.columns.values
        if (mdv != rmv).any():
            raise GeneLabException("Could not generate extended raw metadata")
        as_frame = self.parent.raw_metadata.copy()
        as_frame.columns = multiindex_df.set_index(multicols).index
        return as_frame.sort_index(by="field", axis="columns")

    def __getitem__(self, patterns):
        """Get metadata by field title (rather than internal field id)"""
        if isinstance(patterns, DataFrame) and (patterns.shape[1] == 1):
            # assume being indexed by boolean column, check if able to coerce:
            indexer = patterns.iloc[:,0]
            if indexer.dtype == bool:
                return self.parent.raw_metadata.loc[indexer]
            else:
                raise IndexError("Cannot index by arbitrary DataFrame")
        if isinstance(patterns, (tuple, list, set, Series, Index)):
            # assume being indexed by column name regex:
            titles = set.union(set(), *(
                self.parent._match_field_titles(p, method=fullmatch)
                for p in patterns
            ))
            if titles:
                return self.parent.raw_metadata[
                    list(set.union(*(self.parent._fields[t] for t in titles)))
                ]
            else:
                return DataFrame()
        else:
            raise IndexError("AssayMetadata: column indexer must be list-like")

    @property
    def index(self):
        """List of samples"""
        return self.parent.raw_metadata.index

    @property
    def fields(self):
        """Alias to self.parent._fields"""
        return self.parent._fields

    @property
    def columns(self):
        """List of full-length column indexing options"""
        return list(self.parent._fields.keys())


class Assay():
    """Stores individual assay information and metadata in raw form"""
    name = None
    _fields, raw_metadata, metadata = None, None, None
    parent, glds_file_urls, glds_file_dates = None, None, None
    storage = None
    _normalized_data, _processed_data = None, None
    _indexed_by, _name_delim, _field_indexed_by = None, True, None

    def __init__(self, parent, name, json, glds_file_urls, glds_file_dates, storage_prefix, index_by, name_delim):
        """Parse JSON into assay metadata"""
        self.parent, self.name, self._json = parent, name, json
        self.glds_file_urls = glds_file_urls
        self.glds_file_dates = glds_file_dates
        self.storage = join(storage_prefix, name)
        self._raw, self._header = self._json["raw"], self._json["header"]
        # populate and freeze self._fields (this can be refactored...):
        self._field2title = {e["field"]: e["title"] for e in self._header}
        if len(self._field2title) != len(self._header):
            raise GeneLabJSONException("Conflicting IDs of data fields")
        self._fields = defaultdict(set)
        for field, title in self._field2title.items():
            self._fields[title].add(field)
        self._fields = dict(self._fields)
        # populate metadata and index with `index_by`:
        self.raw_metadata = concat(map(Series, self._raw), axis=1).T
        self._field_indexed_by = self._get_unique_field_from_title(index_by)
        maybe_indexed_by = self._match_field_titles(index_by, method=fullmatch)
        if len(maybe_indexed_by) != 1:
            raise IndexError(
                "Nonexistent or ambiguous index_by value: '{}'".format(index_by)
            )
        self._indexed_by = maybe_indexed_by.pop()
        self.raw_metadata = self.raw_metadata.set_index(self._field_indexed_by)
        self._name_delim = name_delim
        if name_delim != DELIM_AS_IS:
            self.raw_metadata.index = self.raw_metadata.index.map(
                lambda f: sub(r'[._-]', name_delim, f)
            )
        del self._fields[self._indexed_by]
        # initialize indexing functions:
        self.metadata = AssayMetadata(self)

    def _match_field_titles(self, pattern, flags=IGNORECASE, method=search):
        """Find fields matching pattern"""
        if self._indexed_by:
            field_pool = set(self._fields) | {self._indexed_by}
        else:
            field_pool = self._fields
        return {
            title for title in field_pool
            if method(pattern, title, flags=flags)
        }

    def _get_unique_field_from_title(self, title):
        """Get unique raw metadata column name; fail if anything is ambiguous"""
        matching_titles = self._match_field_titles(title)
        if len(matching_titles) == 0:
            raise IndexError("Nonexistent '{}'".format(title))
        elif len(matching_titles) > 1:
            raise IndexError("Ambiguous '{}'".format(title))
        else:
            matching_title = matching_titles.pop()
            if matching_title == self._indexed_by:
                matching_fields = {self._field_indexed_by}
            else:
                matching_fields = self._fields[matching_title]
        if len(matching_fields) == 0:
            raise IndexError("Nonexistent '{}'".format(title))
        elif len(matching_fields) > 1:
            raise IndexError("Ambiguous '{}'".format(title))
        else:
            return list(matching_fields)[0]

    @property
    def factor_values(self):
        """Get factor names and their values"""
        return {
            field_title: set(self.metadata[[field_title]].values.flatten())
            for field_title in self._match_field_titles(r'^factor value:  ')
        }

    def annotation(self, differential_annotation=True, named_only=True, index_by="Sample Name", cls=None, continuous="infer"):
        """Get annotation of samples: entries that differ (default) or all entries"""
        samples_keys = set(self.parent.samples.keys())
        if len(samples_keys) == 1:
            samples_key = samples_keys.pop()
        else:
            samples_key = sub(r'^a', "s", self.name)
        if samples_key not in self.parent.samples:
            error_message = "Could not find an unambiguous samples key"
            raise GeneLabJSONException(error_message)
        annotation_dataframe = concat([
            Series(raw_sample_annotation)
            for raw_sample_annotation in self.parent.samples[samples_key]["raw"]
        ], axis=1)
        samples_field2title = {
            entry["field"]: entry["title"]
            for entry in self.parent.samples[samples_key]["header"]
        }
        if named_only:
            index_subset = [
                field for field in annotation_dataframe.index
                if field in samples_field2title
            ]
            annotation_dataframe = annotation_dataframe.loc[index_subset]
        annotation_dataframe.index = annotation_dataframe.index.map(
            lambda field: samples_field2title.get(field, field)
        )
        if differential_annotation:
            differential_rows = annotation_dataframe.apply(
                lambda r: len(set(r.values))>1, axis=1
            )
            annotation_dataframe = annotation_dataframe[differential_rows]
        annotation_dataframe = annotation_dataframe.T.set_index(index_by).T
        if self._name_delim != DELIM_AS_IS:
            annotation_dataframe.columns = annotation_dataframe.columns.map(
                lambda f: sub(r'[._-]', self._name_delim, f)
            )
        annotation_dataframe.columns.name = index_by
        if cls:
            return to_cls(
                annotation_dataframe.T, target=cls, continuous=continuous
            )
        else:
            return annotation_dataframe.T

    def factors(self, cls=None, continuous="infer"):
        """Get DataFrame of samples and factors in human-readable form"""
        annotation = self.annotation()
        factor_fields = [
            field for field in annotation.columns
            if search(r'^factor value', field, flags=IGNORECASE)
        ]
        factors_dataframe = annotation[factor_fields]
        factors_dataframe.index.name, factors_dataframe.columns.name = (
            self._indexed_by, "Factor"
        )
        if cls == "*":
            if factors_dataframe.shape[1] != 1:
                raise KeyError("one of multiple factors needs to be specified")
            else:
                cls = str(factors_dataframe.columns[0])
            return to_cls(factors_dataframe, target=cls, continuous=continuous)
        elif cls is not None:
            return to_cls(factors_dataframe, target=cls, continuous=continuous)
        else:
            return factors_dataframe

    def _get_file_url(self, filemask):
        """Get URL of file defined by file mask (such as *SRR1781971_*)"""
        regex_filemask = filemask.split("/")[0].replace("*", ".*")
        matching_names = {
            filename for filename in self.glds_file_urls.keys()
            if search(regex_filemask, filename)
        }
        if len(matching_names) == 0:
            return None
        elif len(matching_names) > 1:
            raise GeneLabJSONException("Multiple file URLs match name")
        else:
            return self.glds_file_urls[matching_names.pop()]


class AssayDispatcher(dict):
    """Contains Assay objects, indexable by name or by attributes"""

    def __init__(self, parent, json, glds_file_urls, glds_file_dates, storage_prefix, index_by, name_delim):
        """Populate dictionary of assay_name -> Assay()"""
        try:
            for assay_name, assay_json in json.items():
                super().__setitem__(
                    assay_name,
                    Assay(
                        parent, assay_name, assay_json, index_by=index_by,
                        name_delim=name_delim, storage_prefix=storage_prefix,
                        glds_file_urls=glds_file_urls,
                        glds_file_dates=glds_file_dates
                    )
                )
        except KeyError:
            raise GeneLabJSONException(
                "Malformed assay JSON ({})".format(self.accession)
            )

    @property
    def _summary_dataframe(self):
        """List assay names and types"""
        repr_rows = []
        for assay_name, assay in self.items():
            protocol_set = set(map(
                str.lower,
                assay.metadata[["Protocol REF"]].values.flatten()
            ))
            factors = set(
                sub('^factor value:\s+', "", f, flags=IGNORECASE)
                for f in assay.factors().columns
            )
            for factor in factors:
                repr_row = [assay_name, factor]
                for characteristic in ASSAY_CHARACTERISTICS:
                    repr_row.append(
                        len(assay._match_field_titles(characteristic)) > 0
                    )
                for protocol in ASSAY_PROTOCOL_LIST:
                    repr_row.append(protocol in protocol_set)
                repr_rows.append(repr_row)
        repr_dataframe = DataFrame(
            data=repr_rows,
            columns=(
                ["name", "factors"] +
                ASSAY_CHARACTERISTICS +
                ASSAY_PROTOCOL_LIST
            )
        )
        return repr_dataframe.copy()
