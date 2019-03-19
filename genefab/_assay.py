from os.path import join
from ._exceptions import GeneLabJSONException, GeneLabFileException
from collections import defaultdict
from pandas import concat, Series, Index, DataFrame, read_csv
from re import search, fullmatch, split, IGNORECASE
from numpy import nan
from ._util import fetch_file


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
                return self.parent[titles].loc[indices]
        else: # assume called with .loc[x] and interpret `x` the best we can
            return self.parent.parent.raw_metadata.loc[key]


class AssayMetadata():
    """Makes individual assay metadata accessible with Pandas-like indexing"""

    def __init__(self, parent):
        """Point to parent and initialize children"""
        self.parent = parent
        self.loc = AssayMetadataLocator(self)

    def __getitem__(self, patterns):
        """Get metadata by field title (rather than internal field id)"""
        if isinstance(patterns, Series) and (patterns.dtype == bool):
            return self.parent.raw_metadata.loc[patterns]
        if isinstance(patterns, dict):
            _patterns = list(patterns.keys())
        else:
            _patterns = patterns
        if isinstance(_patterns, (tuple, list, set, Series, Index)):
            titles = set.union(*[
                self.parent._match_field_titles(p, method=fullmatch)
                for p in _patterns
            ])
            if titles:
                return self.parent.raw_metadata[
                    list(set.union(*[self.parent.fields[t] for t in titles]))
                ]
            else:
                return DataFrame()
        else:
            raise IndexError("AssayMetadata: column indexer must be list-like")


class Assay():
    """Stores individual assay information and metadata in raw form"""
    name = None
    parent = None
    raw_metadata = None
    metadata = None
    glds_file_urls = None
    fields = None
    storage = None

    _normalized_data = None
    _processed_data = None
    _indexed_by = None

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
        # populate metadata and index with Sample Name:
        self.raw_metadata = concat(map(Series, self._raw), axis=1).T
        if index_by not in self.fields:
            raise GeneLabJSONException(
                "Cannot index by nonexistent field: '{}'".format(index_by)
            )
        elif len(self.fields[index_by]) != 1:
            raise GeneLabJSONException(
                "Cannot index by ambiguous field: '{}'".format(index_by)
            )
        else:
            self._indexed_by = index_by
            index_field = list(self.fields[index_by])[0]
            self.raw_metadata = self.raw_metadata.set_index(index_field)
        # initialize indexing functions:
        self.metadata = AssayMetadata(self)

    def _match_field_titles(self, pattern, flags=IGNORECASE, method=search):
        """Find fields matching pattern"""
        return {
            title for title in self.fields
            if method(pattern, title, flags=flags)
        }

    @property
    def factors(self):
        """Get factor names and their values"""
        return {
            field_title: set(self.metadata[[field_title]].values.flatten())
            for field_title in self._match_field_titles(r'^factor value:  ')
        }

    @property
    def has_arrays(self):
        return "Array Design REF" in self.fields

    @property
    def has_normalized_data(self):
        return len(self._match_field_titles("normalized data files")) > 0

    @property
    def has_normalized_annotated_data(self):
        return (
            len(self._match_field_titles("normalized annotated data files")) > 0
        )

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
        if not search(data_columns, self._indexed_by, flags=IGNORECASE):
            column_translator = self.metadata[
                self._match_field_titles(data_columns)
            ]
            _h, _w = column_translator.shape
            if _w == 1:
                if len(set(column_translator.index)) == _h:
                    if len(set(column_translator.values.flatten())) == _h:
                        translated_data = data.copy()
                        column_translator = column_translator \
                            .reset_index() \
                            .set_index(column_translator.columns[0])
                        translated_data.columns = [
                            column_translator.loc[colname].values[0]
                            for colname in data.columns
                        ]
                        return translated_data
            else:
                raise IndexError("Cannot reindex '{}' to ambiguous '{}".format(
                    self._indexed_by, data_columns
                ))
        else:
            return data

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
                if translate_sample_names:
                    data = self._translate_data_sample_names(
                        data, data_columns=data_columns
                    )
        else:
            data = DataFrame()
        return data

    def get_normalized_data(self, force_redownload=False, translate_sample_names=True, data_columns="hybridization assay name"):
        """Get normalized data from file(s) listed under 'normalized data files'"""
        self._normalized_data = self._read_data_from(
            ".*normalized data files.*",
            blacklist_regex=r'\.rda(ta)?(\.gz)?$',
            force_redownload=force_redownload,
            translate_sample_names=translate_sample_names,
            data_columns=data_columns
        )

    @property
    def normalized_data(self):
        if self._normalized_data is None:
            self.get_normalized_data()
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

    @property
    def processed_data(self):
        if self._processed_data is None:
            self.get_processed_data()
        return self._processed_data

    # alias:
    def get_normalized_annotated_data(self, force_redownload=False):
        self.get_processed_data(force_redownload=force_redownload)

    # alias:
    @property
    def normalized_annotated_data(self):
        return self.processed_data
