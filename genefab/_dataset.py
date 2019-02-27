from sys import stderr
from re import search, IGNORECASE, compile
from urllib.parse import quote_plus
from ._util import get_json, fetch_file, flat_extract, permissive_search_group
from ._util import FFIELD_ALIASES, FFIELD_VALUES, API_ROOT, GENELAB_ROOT
from ._exceptions import GeneLabJSONException
from ._checks import safe_file_name
from pandas import concat, Series, Index, read_csv
from collections import defaultdict
from numpy import nan
from os.path import join
from os import walk


class AssayMetadataLocator():
    """Emulate behavior of Pandas `.loc` for class Assay()"""

    def __init__(self, parent):
        """Point to parent"""
        self.parent = parent

    def __getitem__(self, key):
        """Query parent.metadata with .loc, using field titles instead of internal field ids"""
        if isinstance(key, tuple): # called with .loc[x, y]
            try:
                indices, titles = key
            except ValueError:
                raise IndexError("Incorrect index for assay metadata")
            subset = self.parent[titles].loc[indices]
        else: # assume called with .loc[x] and interpret `x` the best we can
            subset = self.parent.metadata.loc[key]
        if (subset.shape == (1,)) and (not self.parent.strict_indexing):
            return subset.iloc[0]
        else:
            return subset


class Assay():
    """Stores individual assay metadata"""
    name = None
    parent = None
    metadata = None
    glds_file_urls = None
    fields = None
    strict_indexing = True
    storage = None

    def __init__(self, parent, name, json, glds_file_urls, storage_prefix, strict_indexing=True):
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
        self.metadata = concat(map(Series, self._raw), axis=1).T
        if len(self.fields["Sample Name"]) != 1:
            raise GeneLabJSONException("Number of 'Sample Name' fields != 1")
        else:
            sample_name_field = list(self.fields["Sample Name"])[0]
            self.metadata = self.metadata.set_index(sample_name_field)
        # initialize indexing functions:
        self.strict_indexing = strict_indexing
        self.loc = AssayMetadataLocator(self)

    def __getitem__(self, titles):
        """Get metadata by field title (rather than internal field id)"""
        if isinstance(titles, (tuple, list, Series, Index)):
            if isinstance(titles, Series) and (titles.dtype == bool):
                return self.metadata.loc[titles]
            else:
                return self.metadata[
                    list(set.union(*[self.fields[t] for t in titles]))
                ]
        elif self.strict_indexing:
            raise IndexError("Assay: column indexer must be list-like")
        else:
            subset = self.metadata[list(self.fields[titles])]
            if subset.shape[1] > 1:
                return subset
            else:
                return subset.iloc[:,0]

    @property
    def available_file_types(self):
        """List file types referenced in metadata"""
        file_types = set()
        for title in self.fields:
            if search(r'\bfile\b', title, flags=IGNORECASE):
                available_files = self[[title]].values.flatten()
                if not set(available_files) <= {"", None, nan}:
                    file_types.add(title)
        return file_types

    @property
    def available_derived_file_types(self):
        """List file types with derived data referenced in metadata"""
        return {
            permissive_search_group(r'^.*(processed|derived).+file.*$', ft)
            for ft in self.available_file_types
        } - {None}

    @property
    def available_protocols(self):
        """List protocol REFs referenced in metadata"""
        return set(self[["Protocol REF"]].values.flatten()) - {"", None, nan}

    def _get_file_url(self, filemask):
        """Get URL of file defined by file mask (such as *SRR1781971_*)"""
        regex_filemask = filemask.split("/")[0].replace("*", ".*")
        matching_names = [
            filename for filename in self.glds_file_urls.keys()
            if search(regex_filemask, filename)
        ]
        if len(matching_names) == 0:
            return None
        elif len(matching_names) > 1:
            raise GeneLabJSONException("Multiple file URLs match name")
        else:
            return self.glds_file_urls[matching_names[0]]

    def _download_archive_files(self, force_reload=False):
        """Download files/archives etc that contain files targeted by get_combined_matrix()"""
        derived_entries = self[list(self.available_derived_file_types)]
        for field in derived_entries.columns:
            if derived_entries[field].str.endswith(".zip").all():
                internal_file_urls = set(map(
                    self._get_file_url, derived_entries[field].values
                ))
                break
        else:
            internal_file_urls = set()
        external_file_urls = set(filter(
            compile(r'^(ftp|http|https):\/\/').search,
            derived_entries.values.flatten()
        ))
        file_urls = (internal_file_urls | external_file_urls) - {None}
        for file_url in file_urls:
            file_name = safe_file_name(file_url)
            fetch_file(file_name, file_url, self.storage, update=force_reload)

    def get_combined_matrix(self, force_reload=False):
        """Download (if necessary), parse and combine derived files"""
        self._download_archive_files(force_reload=force_reload)
        filenames = next(walk(self.storage))[2]
        zips = [filename for filename in filenames if filename.endswith(".zip")]
        for zip_filename in zips:
            flat_extract(
                join(self.storage, zip_filename),
                target_directory=self.storage
            )
        derived_entries = self[list(self.available_derived_file_types)]
        for field in derived_entries.columns:
            if len(set(derived_entries[field])) == derived_entries.shape[0]:
                derived_files = derived_entries[field]
                break
        else:
            raise NotImplementedError(
                "Derived files not referenced individually"
            )
        sample_dataframes = []
        for sample_name, derived_filename in derived_files.iteritems():
            filename = safe_file_name(
                derived_filename, as_mask=True, directory=self.storage
            )
            sample_dataframe = read_csv(
                join(self.storage, filename), sep="\t"
            )
            sample_dataframe["Sample Name"] = sample_name
            sample_dataframes.append(sample_dataframe)
        return concat(sample_dataframes, axis=0, ignore_index=True)


class GeneLabDataSet():
    """Stores GLDS metadata associated with an accession number"""
    accession = None
    assays = None
    file_urls = None
    verbose = False
    storage = None

    def __init__(self, accession, assay_strict_indexing=True, verbose=False, storage_prefix=".genelab"):
        """Request JSON representation of ISA metadata and store fields"""
        self.accession = accession
        self.verbose = verbose
        self.storage = join(storage_prefix, accession)
        data_json = get_json(
            "{}/data/study/data/{}/".format(API_ROOT, accession), self.verbose
        )
        if len(data_json) > 1:
            raise GeneLabJSONException("Too many results returned")
        self._json = data_json[0]
        try:
            self.internal_id = self._json["_id"]
            self.metadata_id = self._json["metadata_id"]
            if len(self._json["foreignFields"]) != 1:
                raise NotImplementedError("Multiple foreignFields")
            self._isa2json = self._json["foreignFields"][0]["isa2json"]
            self._info = self._isa2json["additionalInformation"]
            for field in "description", "samples", "ontologies", "organisms":
                setattr(self, field, self._info[field])
        except KeyError:
            raise GeneLabJSONException(
                "Malformed JSON ({})".format(self.accession)
            )
        try:
            self.assays = [
                Assay(
                    self, assay_name, assay_json, storage_prefix=self.storage,
                    glds_file_urls=self._get_file_urls(),
                    strict_indexing=assay_strict_indexing
                )
                for assay_name, assay_json in self._info["assays"].items()
            ]
        except KeyError:
            raise GeneLabJSONException(
                "Malformed assay JSON ({})".format(self.accession)
            )

    @property
    def factors(self):
        """List factors"""
        return [
            factor_info["factor"] for factor_info in self.description["factors"]
        ]
 
    def __repr__(self):
        """Simple description, for now"""
        return (
            self.accession +
            " (number of assays: {}".format(len(self.assays)) +
            "; factors: " + ", ".join(self.factors) + ")"
        )
 
    def _get_file_urls(self, force_reload=False):
        """Get filenames and associated URLs"""
        if self.accession is None:
            raise ValueError("Uninitialized GLDS instance")
        elif (not force_reload) and (self.file_urls is not None):
            return self.file_urls
        else:
            getter_url = "{}/data/glds/files/{}"
            acc_nr = search(r'\d+$', self.accession).group()
            files_json = get_json(
                getter_url.format(API_ROOT, acc_nr), self.verbose
            )
            try:
                filedata = files_json["studies"][self.accession]["study_files"]
            except KeyError:
                raise GeneLabJSONException("Malformed JSON")
            self.file_urls = {
                fd["file_name"]: GENELAB_ROOT+fd["remote_url"]
                for fd in filedata
            }
            return self.file_urls


def get_ffield_matches(verbose=False, **ffield_kwargs):
    """Expand passed regexes to all matching ffield values"""
    for ffield_alias, ffregex in ffield_kwargs.items():
        if verbose:
            print("looking up", ffield_alias, end="(s): ", file=stderr)
        if ffield_alias in FFIELD_ALIASES:
            ffield = FFIELD_ALIASES[ffield_alias]
        else:
            raise IndexError("Unrecognized field: " + ffield_alias)
        for ffvalue in FFIELD_VALUES[ffield]:
            if search(ffregex, ffvalue, IGNORECASE):
                if verbose:
                    print('"{}"'.format(ffvalue), end=", ", file=stderr)
                yield ffield, ffvalue
        if verbose:
            print("\b", file=stderr)


def get_datasets(maxcount="25", assay_strict_indexing=True, storage=".genelab", verbose=False, onerror="warn", **ffield_kwargs):
    """Match passed regexes and combine into search URL, get JSON and parse for accessions"""
    url_lead_components = [
        API_ROOT+"/data/search/?term=GLDS", "type=cgene", "size="+str(maxcount)
    ]
    url_ffield_components = [
        "ffield={}&fvalue={}".format(ffield, quote_plus(ffvalue))
        for ffield, ffvalue
        in get_ffield_matches(verbose=verbose, **ffield_kwargs)
    ]
    url = "&".join(url_lead_components + url_ffield_components)
    try:
        json = get_json(url, verbose=verbose)["hits"]["hits"]
    except:
        raise GeneLabJSONException("Unrecognized JSON structure")
    datasets = []
    for hit in json:
        try:
            datasets.append(
                GeneLabDataSet(
                    hit["_id"], assay_strict_indexing=assay_strict_indexing,
                    storage_prefix=storage, verbose=verbose
                )
            )
        except GeneLabJSONException as e:
            if onerror == "ignore":
                pass
            elif onerror == "warn":
                print("Warning:", e, file=stderr)
            else:
                raise
    return datasets
