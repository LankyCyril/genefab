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

    _processed_data = None

    def __init__(self, parent, name, json, glds_file_urls, storage_prefix):
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

    def get_processed_data(self, force_redownload=False):
        """Get processed data from file(s) listed under 'normalized annotated data files'"""
        meta_files = self.metadata[[".*normalized annotated data files.*"]]
        if len(meta_files):
            filenames = set.union(*(
                set(split(r'\s*,\s*', entry))
                for entry in meta_files.values.flatten()
            ))
            target_filenames = {
                filename for filename in filenames
                if not search(r'\.rda(ta)?(\.gz)?$', filename)
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
                self._processed_data = read_csv(
                    join(self.storage, filename), sep="\t", index_col=0
                )
        else:
            self._processed_data = DataFrame()

    @property
    def processed_data(self):
        if self._processed_data is None:
            self.get_processed_data()
        return self._processed_data

    # aliases:
    @property
    def normalized_annotated_data(self): return self.processed_data
