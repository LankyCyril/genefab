from sys import stderr
from re import search, IGNORECASE
from urllib.parse import quote_plus
from ._util import get_json
from ._util import FFIELD_ALIASES, FFIELD_VALUES, API_ROOT, GENELAB_ROOT
from pandas import concat, Series
from pandas.core.indexes.base import Index
from collections import defaultdict


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
                raise ValueError("Incorrect index for assay metadata")
            return self.parent[titles].loc[indices]
        else: # assume called with .loc[x] and interpret `x` the best we can
            return self.parent.metadata.loc[key]


class Assay():
    """Stores individual assay metadata"""
    name = None
    metadata = None
    glds_file_urls = {}
    fields = defaultdict(set)
 
    def __init__(self, assay_name, assay_json, glds_file_urls):
        """Prase JSON into assay metadata"""
        self.name = assay_name
        self.glds_file_urls = glds_file_urls
        self._json = assay_json
        self._raw, self._header = self._json["raw"], self._json["header"]
        self._field2title = {
            entry["field"]: entry["title"] for entry in self._header
        }
        if len(self._field2title) != len(self._header):
            raise ValueError("Conflicting IDs of data fields")
        for field, title in self._field2title.items():
            self.fields[title].add(field)
        self.metadata = concat(map(Series, self._raw), axis=1).T
        self.loc = AssayMetadataLocator(self)
 
    def __getitem__(self, titles):
        """Get metadata by field title (rather than internal field id)"""
        if not isinstance(titles, (tuple, list, Series, Index)):
            raise IndexError("Assay: column indexer must be list-like")
        return self.metadata[
            list(set.union(*[self.fields[t] for t in titles]))
        ]
 
    def get_file_url(self, filemask):
        """Get URL of file defined by file mask (such as *SRR1781971_*)"""
        regex_filemask = filemask.replace("*", ".*")
        matching_names = [
            filename for filename in self.glds_file_urls.keys()
            if search(regex_filemask, filename)
        ]
        if len(matching_names) == 0:
            raise ValueError("No URL found")
        elif len(matching_names) > 1:
            raise ValueError("Multiple file URLs match name")
        else:
            return self.glds_file_urls[matching_names[0]]


class GeneLabDataSet():
    """Stores GLDS metadata associated with an accession number"""
    accession = None
    assays = []
    file_urls = None
    verbose = False
 
    def __init__(self, accession, verbose=False):
        """Request JSON representation of ISA metadata and store fields"""
        self.accession = accession
        self.verbose = verbose
        getter_url = "{}/data/study/data/{}/"
        data_json = get_json(
            getter_url.format(API_ROOT, accession), self.verbose
        )
        if len(data_json) > 1:
            raise ValueError("Too many results returned, unexpected behavior")
        else:
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
            raise ValueError("Malformed JSON")
        self.assays = [
            Assay(assay_name, assay_json, glds_file_urls=self._get_file_urls())
            for assay_name, assay_json in self._info["assays"].items()
        ]
 
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
                raise ValueError("Malformed JSON")
            self.file_urls = {
                fd["file_name"]: GENELAB_ROOT+fd["remote_url"]
                for fd in filedata
            }
            return self.file_urls


def get_ffield_matches(**kwargs):
    """Expand passed regexes to all matching ffield values"""
    for ffield_alias, ffregex in kwargs.items():
        print("looking up", ffield_alias, end="(s): ", file=stderr)
        if ffield_alias in FFIELD_ALIASES:
            ffield = FFIELD_ALIASES[ffield_alias]
        else:
            raise ValueError("Unrecognized field: " + ffield_alias)
        for ffvalue in FFIELD_VALUES[ffield]:
            if search(ffregex, ffvalue, IGNORECASE):
                print('"{}"'.format(ffvalue), end=", ", file=stderr)
                yield ffield, ffvalue
        print("\b", file=stderr)


def get_datasets(**kwargs):
    """Match passed regexes and combine into search URL, get JSON and parse for accessions"""
    if "maxcount" in kwargs:
        maxcount = str(kwargs["maxcount"])
        del kwargs["maxcount"]
    else:
        maxcount = "25"
    if "verbose" in kwargs:
        verbose = bool(kwargs["verbose"])
        del kwargs["verbose"]
    else:
        verbose = False
    term_pairs = [
        "ffield={}&fvalue={}".format(ffield, quote_plus(ffvalue))
        for ffield, ffvalue in get_ffield_matches(**kwargs)
    ]
    url = "&".join(
        [API_ROOT+"/data/search/?term=GLDS", "type=cgene", "size="+maxcount]
        + term_pairs
    )
    try:
        json = get_json(url, verbose=verbose)["hits"]["hits"]
    except:
        raise ValueError("Unrecognized JSON structure")
    return [
        GeneLabDataSet(hit["_id"], verbose=verbose) for hit in json
    ]
