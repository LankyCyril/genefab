from sys import stderr
from re import search, IGNORECASE
from urllib.parse import quote_plus
from genefab._util import get_json
from genefab._util import FFIELD_ALIASES, FFIELD_VALUES, API_ROOT, GENELAB_ROOT
from genefab._util import DELIM_DEFAULT, STORAGE_PREFIX
from genefab._exceptions import GeneLabJSONException
from genefab._assay import AssayDispatcher
from pandas import DataFrame, concat
from os.path import join


class GeneLabDataSet():
    """Stores GLDS metadata associated with an accession number"""
    accession, assays, storage = None, None, None
    verbose = False

    def __init__(self, accession, verbose=False, storage_prefix=STORAGE_PREFIX, index_by="Sample Name", name_delim=DELIM_DEFAULT):
        """Request JSON representation of ISA metadata and store fields"""
        self.accession = accession
        self.verbose = verbose
        self.storage = join(storage_prefix, accession)
        data_json = get_json(
            "{}/data/study/data/{}/".format(API_ROOT, accession), self.verbose
        )
        if len(data_json) == 0:
            raise GeneLabJSONException("Invalid JSON (GLDS does not exist?)")
        if len(data_json) > 1:
            raise GeneLabJSONException("Invalid JSON, too many sections")
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
            parent=self, json=self._info["assays"], storage_prefix=self.storage,
            name_delim=name_delim, index_by=index_by,
            glds_file_urls=self.get_files_info("urls"),
            glds_file_dates=self.get_files_info("dates")
        )

    @property
    def factors(self):
        """List factors"""
        return [fi["factor"] for fi in self.description["factors"]]

    @property
    def summary_dataframe(self):
        """List factors, assay names and types"""
        assays_df = self.assays.summary_dataframe.copy()
        assays_df["type"] = "assay"
        factors_df = DataFrame(
            columns=["type", "name", "factors"],
            data=[
                ["dataset", self.accession, factor]
                for factor in self.factors
            ]
        )
        return concat([factors_df, assays_df], axis=0, sort=False)

    def get_files_info(self, kind="urls"):
        """Get filenames and associated URLs"""
        if self.accession is None:
            raise ValueError("Uninitialized GLDS instance")
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
            if kind == "urls":
                return {
                    fd["file_name"]: GENELAB_ROOT+fd["remote_url"]
                    for fd in filedata
                }
            elif kind == "dates":
                return {fd["file_name"]: fd["date_created"] for fd in filedata}
            else:
                raise ValueError("Unrecognized parameter: '{}'".format(kind))


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


def get_datasets(maxcount="25", storage=STORAGE_PREFIX, verbose=False, onerror="warn", **ffield_kwargs):
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
        except Exception as e:
            if onerror == "ignore":
                pass
            elif onerror == "warn":
                msgmask = "Warning: Could not process {} due to error:"
                print(msgmask.format(hit["_id"]), e, file=stderr)
            else:
                raise
    return datasets
