from sys import stderr
from re import search, IGNORECASE
from urllib.parse import quote_plus
from ._util import get_json
from ._util import FFIELD_ALIASES, FFIELD_VALUES, API_ROOT, GENELAB_ROOT
from ._exceptions import GeneLabJSONException
from ._assay import AssayDispatcher
from os.path import join


class GeneLabDataSet():
    """Stores GLDS metadata associated with an accession number"""
    accession = None
    assays = None
    file_urls = None
    verbose = False
    storage = None

    def __init__(self, accession, verbose=False, storage_prefix=".genelab", index_by="Sample Name"):
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
        self.assays = AssayDispatcher(
            parent=self, json=self._info["assays"], index_by=index_by,
            storage_prefix=self.storage, glds_file_urls=self._get_file_urls()
        )

    @property
    def factors(self):
        """List factors"""
        return [
            factor_info["factor"] for factor_info in self.description["factors"]
        ]
 
    def __repr__(self):
        """Simple description, for now"""
        return "\n".join([
            "name: " + self.accession,
            "assays: [" + ", ".join(
                repr(assay_name) for assay_name in self.assays.keys()
            ) + "]",
            "factors: [" + ", ".join(
                repr(factor) for factor in self.factors
            ) + "]"
        ])
 
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


def get_datasets(maxcount="25", storage=".genelab", verbose=False, onerror="warn", **ffield_kwargs):
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
                    hit["_id"], storage_prefix=storage, verbose=verbose
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
