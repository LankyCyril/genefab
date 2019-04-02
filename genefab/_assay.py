from os.path import join
from ._exceptions import GeneLabJSONException, GeneLabFileException
from collections import defaultdict
from pandas import concat, Series, Index, DataFrame, read_csv
from re import search, fullmatch, split, IGNORECASE
from numpy import nan
from ._util import fetch_file


class MetadataRow():
    """Implements a slice of assay metadata for one sample, Series-like"""

    def __init__(self, parent, sample, raw_row):
        """Inherit from parent(s)"""
        self.parent = parent
        self.sample = sample
        self.raw_row = raw_row.copy()

    def __getitem__(self, key):
        """Reuse parent methods"""
        if isinstance(key, str):
            return self.parent.metadata.loc[self.sample, [key]].iloc[0]
        else:
            return self.parent.metadata.loc[self.sample, key]

    def __repr__(self):
        """Short description of fields and samples"""
        return "\n".join([
            "Sample: " + self.sample,
            "Fields: [" + ", ".join(
                repr(k) for k in self.parent.fields.keys()
            ) + "]"
        ])


class AssayMetadataLocator():
    """Emulate behavior of Pandas `.loc` for class AssayMetadata()"""

    def __init__(self, parent):
        """Point to parent"""
        self.parent = parent

    def __getitem__(self, key):
        """Query parent.raw_metadata with .loc, using field titles instead of internal field ids"""
        if isinstance(key, tuple): # called with .loc[x, y]
            try:
                indices, titles = key
            except ValueError:
                raise IndexError("Incorrect index for assay metadata")
            else:
                row_subset = self.parent.loc[indices]
                field_titles = set.union(*(
                    self.parent.parent._match_field_titles(t, method=fullmatch)
                    for t in titles
                ))
                fields = set.union(*(
                    self.parent.parent.fields[title]
                    for title in field_titles
                ))
                return row_subset[list(fields)]
        else: # assume called with .loc[x] and interpret `x` the best we can
            if isinstance(key, DataFrame) and (key.shape[1] == 1):
                # assume being indexed by boolean column, delegate to parent[]:
                return self.parent[key]
            else:
                return self.parent.parent.raw_metadata.loc[key]


class AssayMetadata():
    """Makes individual assay metadata accessible with Pandas-like indexing"""

    def __init__(self, parent):
        """Point to parent and initialize children"""
        self.parent = parent
        self.loc = AssayMetadataLocator(self)

    def __repr__(self):
        """Short description of fields and samples"""
        return "\n".join([
            "Samples: [" + ", ".join(
                repr(ix) for ix in self.parent.raw_metadata.index
            ) + "]",
            "Fields: [" + ", ".join(
                repr(k) for k in self.parent.fields.keys()
            ) + "]",
            "Factor values: " + repr(self.parent.factor_values)
        ])

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
            titles = set.union(set(), *(
                self.parent._match_field_titles(p, method=fullmatch)
                for p in patterns
            ))
            if titles:
                return self.parent.raw_metadata[
                    list(set.union(*(self.parent.fields[t] for t in titles)))
                ]
            else:
                return DataFrame()
        else:
            raise IndexError("AssayMetadata: column indexer must be list-like")

    def iterrows(self):
        """Iterate over metadata slices for each sample"""
        for sample, raw_row in self.parent.raw_metadata.iterrows():
            yield sample, MetadataRow(self.parent, sample, raw_row)


class Assay():
    """Stores individual assay information and metadata in raw form"""
    name = None
    fields, raw_metadata, metadata = None, None, None
    parent, glds_file_urls = None, None
    storage = None

    _normalized_data, _processed_data = None, None
    _indexed_by, _field_indexed_by = None, None

    def __init__(self, parent, name, json, glds_file_urls, storage_prefix, index_by="Sample Name"):
        """Parse JSON into assay metadata"""
        self.parent = parent
        self.name = name
        self.glds_file_urls = glds_file_urls
        self.storage = join(storage_prefix, name)
        self._json = json
        self._raw, self._header = self._json["raw"], self._json["header"]
        # populate and freeze self.fields (this can be refactored...):
        self._field2title = {
            entry["field"]: entry["title"] for entry in self._header
        }
        if len(self._field2title) != len(self._header):
            raise GeneLabJSONException("Conflicting IDs of data fields")
        self.fields = defaultdict(set)
        for field, title in self._field2title.items():
            self.fields[title].add(field)
        self.fields = dict(self.fields)
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
        del self.fields[self._indexed_by]
        # initialize indexing functions:
        self.metadata = AssayMetadata(self)

    def _match_field_titles(self, pattern, flags=IGNORECASE, method=search):
        """Find fields matching pattern"""
        if self._indexed_by:
            field_pool = set(self.fields) | {self._indexed_by}
        else:
            field_pool = self.fields
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
                matching_fields = self.fields[matching_title]
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

    @property
    def factors(self):
        """Get DataFrame of samples and factors in human-readable form"""
        factor_field2title = {}
        for factor in self.factor_values:
            factor_titles = self.fields[factor]
            if len(factor_titles) != 1:
                raise GeneLabJSONException(
                    "Nonexistent or ambiguous factor fields: '{}'".format(
                        factor
                    )
                )
            else:
                factor_field2title[list(factor_titles)[0]] = factor
        raw_factors_dataframe = self.metadata[list(self.factor_values.keys())]
        factors_dataframe = raw_factors_dataframe.copy()
        factors_dataframe.columns = [
            factor_field2title[field] for field in raw_factors_dataframe.columns
        ]
        factors_dataframe.index.name, factors_dataframe.columns.name = (
            self._indexed_by, "Factor"
        )
        return factors_dataframe

    @property
    def samples(self):
        """Get sample names"""
        return list(self.raw_metadata.index)

    @property
    def has_arrays(self):
        return "Array Design REF" in self.fields

    @property
    def has_normalized_data(self):
        return len(self._match_field_titles("normalized data files")) > 0

    @property
    def has_processed_data(self):
        return (
            len(self._match_field_titles("normalized annotated data files")) > 0
        )

    # alias:
    @property
    def has_normalized_annotated_data(self):
        return self.has_processed_data

    @property
    def available_file_types(self):
        """List file types referenced in metadata"""
        file_types = set()
        for title in self._match_field_titles(r'\bfile\b'):
            available_files = self.metadata[[title]].values.flatten()
            if not set(available_files) <= {"", None, nan}:
                file_types.add(title)
        return file_types

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

    def _translate_data_sample_names(self, data, data_columns="hybridization assay name"):
        """Convert data header to match metadata index"""
        field_from = self._get_unique_field_from_title(data_columns)
        field_to = self._field_indexed_by
        column_from, column_to = (
            self.raw_metadata.reset_index()[field_from],
            self.raw_metadata.reset_index()[field_to]
        )
        if len(column_from) == len(set(column_from)) == len(set(column_to)):
            column_translator = dict(zip(column_from, column_to))
        else:
            raise IndexError("Cannot reindex '{}' to ambiguous '{}'".format(
                self._indexed_by, data_columns
            ))
        translated_data = data.copy()
        translated_columns = []
        for column in data.columns:
            matching_keys = {
                k for k in column_translator.keys() if search(k, column)
            }
            if len(matching_keys) == 1:
                translated_columns.append(
                    column_translator[matching_keys.pop()]
                )
            else:
                raise IndexError("Cannot reindex '{}' to ambiguous '{}'".format(
                    self._indexed_by, data_columns
                ))
        translated_data.columns = Index(
            translated_columns, name=data.columns.name
        )
        return translated_data

    def _read_data_from(self, field_title, blacklist_regex, force_redownload, translate_sample_names, data_columns, sep="\t"):
        """Download (if necessary) and parse data contained in a single target file linked to by target field"""
        meta_files = self.metadata[[field_title]]
        if len(meta_files):
            filenames = set.union(*(
                set(split(r'\s*,\s*', e)) for e in meta_files.values.flatten()
            ))
            target_filenames = {
                fn for fn in filenames if not search(blacklist_regex, fn)
            }
            if len(target_filenames) == 0:
                raise GeneLabFileException(
                    "No suitable normalized annotated data files found"
                )
            elif len(target_filenames) > 1:
                raise GeneLabFileException(
                    "Multiple normalized annotated data files found"
                )
            else:
                filename = target_filenames.pop()
                url = self._get_file_url(filename)
                fetch_file(filename, url, self.storage, update=force_redownload)
                csv = join(self.storage, filename)
                data = read_csv(csv, sep=sep, index_col=0)
                data.columns.name = self._indexed_by
                if translate_sample_names:
                    return self._translate_data_sample_names(
                        data, data_columns=data_columns
                    )
                else:
                    return data
        else:
            return None

    def get_normalized_data(self, force_redownload=False, translate_sample_names=True, data_columns="hybridization assay name"):
        """Get normalized data from file(s) listed under 'normalized data files'"""
        self._normalized_data = self._read_data_from(
            ".*normalized data files.*",
            blacklist_regex=r'\.rda(ta)?(\.gz)?$',
            force_redownload=force_redownload,
            translate_sample_names=translate_sample_names,
            data_columns=data_columns
        )
        return self._normalized_data

    @property
    def normalized_data(self):
        if self._normalized_data is None:
            return self.get_normalized_data()
        else:
            return self._normalized_data

    def get_processed_data(self, force_redownload=False, translate_sample_names=True, data_columns="hybridization assay name"):
        """Get processed data from file(s) listed under 'normalized annotated data files'"""
        self._processed_data = self._read_data_from(
            ".*normalized annotated data files.*",
            blacklist_regex=r'\.rda(ta)?(\.gz)?$',
            force_redownload=force_redownload,
            translate_sample_names=translate_sample_names,
            data_columns=data_columns
        )
        return self._processed_data

    @property
    def processed_data(self):
        if self._processed_data is None:
            return self.get_processed_data()
        else:
            return self._processed_data

    # alias:
    def get_normalized_annotated_data(self, force_redownload=False):
        return self.get_processed_data(force_redownload=force_redownload)

    # alias:
    @property
    def normalized_annotated_data(self):
        return self.processed_data
