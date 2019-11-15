from sys import stderr
from re import search, IGNORECASE
from urllib.parse import quote_plus
from genefab._util import get_json, date2stamp
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
            error_message = "Malformed JSON ({})".format(self.accession)
            raise GeneLabJSONException(error_message)
        self.assays = AssayDispatcher(
            parent=self, json=self._info["assays"], storage_prefix=self.storage,
            name_delim=name_delim, glds_file_urls=self.get_files_info("urls"),
            index_by=index_by, glds_file_dates=self.get_files_info("dates")
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
        elif kind == "urls":
            getter_url = "{}/data/glds/files/{}"
            acc_nr = search(r'\d+$', self.accession).group()
            files_json = get_json(
                getter_url.format(API_ROOT, acc_nr), self.verbose
            )
            try:
                filedata = files_json["studies"][self.accession]["study_files"]
            except KeyError:
                raise GeneLabJSONException("Malformed JSON")
            return {
                fd["file_name"]: GENELAB_ROOT+fd["remote_url"]
                for fd in filedata
            }
        elif kind == "dates":
            getter_url = "{}/data/study/filelistings/{}"
            filedata = get_json(
                getter_url.format(API_ROOT, self.internal_id), self.verbose
            )
            return {fd["file_name"]: date2stamp(fd) for fd in filedata}
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
