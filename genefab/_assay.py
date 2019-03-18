from os.path import join
from ._exceptions import GeneLabJSONException
from collections import defaultdict
from pandas import concat, Series, Index, DataFrame
from re import search, fullmatch, IGNORECASE
from numpy import nan


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
        """Inherit values from parent Assay()"""
        self.parent = parent
        self.loc = AssayMetadataLocator(self)

    def __getitem__(self, patterns):
        """Get metadata by field title (rather than internal field id)"""
        if isinstance(patterns, Series) and (patterns.dtype == bool):
            return self.parent.raw_metadata.loc[patterns]
        elif isinstance(patterns, (tuple, list, set, Series, Index)):
            titles = set.union(*[
                self.parent._match_field_titles(p, method=fullmatch)
                for p in patterns
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

    def __init__(self, parent, name, json, glds_file_urls, storage_prefix):
        """Prase JSON into assay metadata"""
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
        if len(self.fields["Sample Name"]) != 1:
            raise GeneLabJSONException("Number of 'Sample Name' fields != 1")
        else:
            sample_name_field = list(self.fields["Sample Name"])[0]
            self.raw_metadata = self.raw_metadata.set_index(sample_name_field)
        # initialize indexing functions:
        self.metadata = AssayMetadata(self)

    def _match_field_titles(self, pattern, flags=IGNORECASE, method=search):
        """Find fields matching pattern"""
        return {
            title for title in self.fields
            if method(pattern, title, flags=flags)
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
