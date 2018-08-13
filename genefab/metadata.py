from urllib.request import urlopen
from json import loads
from copy import deepcopy
from numpy import nan
from re import sub, search
from pandas import DataFrame, concat
from argparse import Namespace
from warnings import warn
from os.path import join, exists, isdir
from subprocess import call

URL_ROOT = "https://genelab-data.ndc.nasa.gov/genelab"

def get_json(url):
    """HTTP get, decode, parse"""
    with urlopen(url) as response:
        return loads(response.read().decode())

def get_hits_json(max_hit_count=10000):
    """Get all possible 'hits' to projects from the GLDS database"""
    poller_url = "{}/data/search?term=GLDS&type=cgene&size={}"
    url = poller_url.format(URL_ROOT, max_hit_count)
    meta = get_json(url)
    if meta["timed_out"]:
        raise ConnectionError("Request timed out")
    else:
        try:
            return meta["hits"]["hits"]
        except KeyError:
            raise ValueError("Unexpected JSON format")

def humanize_hit_source(hit_source, na_values=(), sort_subfields=False):
    """Flatten fields, sprinkle separators where necessary"""
    hhs = deepcopy(hit_source)
    hhs["Study Person"] = " ".join(hit_source["Study Person"].values())
    hhs["Mission Start Date"] = hit_source["Mission"]["Start Date"]
    hhs["Mission End Date"] = hit_source["Mission"]["End Date"]
    hhs["Mission Name"] = hit_source["Mission"]["Name"]
    del hhs["Mission"]
    for key, value in hhs.items():
        if not isinstance(value, str):
            err = "{} field {} is not a string".format(hhs["Accession"], key)
            raise ValueError(err)
        else:
            if value in na_values:
                value = nan
            elif "  " in value:
                # NB! This conversion is crude and does not yet account for
                # many edge cases, or 'parallel' subfields in different fields
                value = sub(r'\s\s+', ";", value)
                if sort_subfields:
                    subfields = sorted(value.split(";"))
                    value = ";".join(subfields)
            hhs[key] = value
    return hhs

def hits_to_dataframe(hits, na_values=("NA", "NaN", ""), sort_subfields=True):
    """Find all fields that exist in sources, combine all data into one dataframe"""
    fields = set()
    humanized_sources = []
    for hit in hits:
        hs = humanize_hit_source(hit["_source"], na_values, sort_subfields)
        fields |= set(hs.keys())
        humanized_sources.append(hs)
    meta_table = DataFrame(
        columns=sorted(fields), index=range(len(humanized_sources))
    )
    for i, hs in enumerate(humanized_sources):
        for key, value in hs.items():
            meta_table.loc[i, key] = value
    if "Accession" in meta_table.columns:
        meta_table.index = meta_table["Accession"]
        return meta_table
    else:
        raise KeyError("Key 'Accession' missing from metadata")

class GeneLabDataSet():
    accession = None
    _frame = None
    _storage = None
    """"""
    def __init__(self, accession, storage=".genefab"):
        """Request JSON representation of ISA metadata and store fields"""
        self.accession = accession
        self._storage = join(storage, accession)
        if not exists(self._storage):
            call(["mkdir", "-p", self._storage])
        elif not isdir(self._storage):
            err_mask = "Storage directory {} cannot be created (file exists)"
            raise OSError(err_mask.format(self._storage))
        getter_url = "{}/data/study/data/{}/"
        data_json = get_json(getter_url.format(URL_ROOT, accession))
        if len(data_json) > 1:
            raise ValueError("Too many results returned, unexpected behavior")
        try:
            self._isa2json = data_json[0]["foreignFields"][0]["isa2json"]
            assay_data = self._isa2json["additionalInformation"]["assays"]
            assay_json = list(assay_data.values())[0]
            self._header = assay_json["header"]
            self._raw = assay_json["raw"]
        except KeyError:
            raise ValueError("Malformed JSON")
    """"""
    def _field_name_to_id(self, field_name):
        """Convert external field name to internal field id"""
        for record in self._header:
            if record["title"] == field_name:
                # trust that field is unique (write checker!)
                if "field" in record:
                    return record["field"]
                elif "columns" in record:
                    # trust that field is unique (write checker!)
                    if "field" in record["columns"][0]:
                        return record["columns"][0]["field"]
        else:
            raise KeyError("Field {} missing from metadata".format(field_name))
    """"""
    def _is_consistent(self, raw):
        """Check if keys are the same for each record in _raw"""
        key_set = set(raw[0].keys())
        for record in raw:
            if set(record.keys()) != key_set:
                return False
        else:
            return True
    """"""
    def frame(self):
        """Convert _raw field of _isa2json to pandas DataFrame"""
        if self._frame is None:
            if not self._is_consistent(self._raw):
                raise KeyError("_raw field keys are inconsistent")
            else:
                self._frame = DataFrame(
                    columns=sorted(self._raw[0].keys()),
                    index=range(len(self._raw))
                )
                for i, record in enumerate(self._raw):
                    for key, value in record.items():
                        self._frame.loc[i, key] = value
                sample_names_field_id = self._field_name_to_id("Sample Name")
                self._frame.set_index(sample_names_field_id, inplace=True)
        return self._frame
    """"""
    def get_factors(self, as_fields=False):
        """Get factor type from _header"""
        factors = {}
        for record in self._header:
            if record["title"] == "Factor Value":
                for column in record["columns"]:
                    factor = column["title"]
                    if as_fields:
                        values = column["field"]
                    else:
                        values = set(self.frame()[column["field"]].values)
                    factors[factor] = values
        if not factors:
            raise KeyError("No factor associated with dataset")
        return factors
    """"""
    def datafiles(self, as_relpath=False):
        """Return DataFrame subset to filenames and factor values"""
        factor_fields = self.get_factors(as_fields=True)
        factor_dataframe = self.frame()[list(factor_fields.values())]
        factor_dataframe.columns = list(factor_fields.keys())
        datafiles_field_id = self._field_name_to_id("Array Data File")
        names_only = self.frame()[datafiles_field_id].apply(
            lambda fn: "{}_microarray_{}.gz".format(
                self.accession, sub(r'^\*', "", fn)
            )
        )
        names_only.name = "filename"
        return concat([factor_dataframe, names_only], axis=1)
    """"""
    def datafile_urls(self):
        """Return DataFrame subset to filenames (converted to URLs) and factor values"""
        datafiles = self.datafiles()
        datafiles_urls = datafiles["filename"].apply(
            lambda fn: "{}/static/media/dataset/{}".format(
                URL_ROOT, fn
            )
        )
        datafiles_urls.name = "URL"
        return concat([datafiles, datafiles_urls], axis=1)
    """"""
    def retrieve_datafiles(self):
        """Actually download datafiles"""
        meta = self.datafile_urls()[["filename", "URL"]]
        for _, filename, url in meta.itertuples():
            call(["wget", "-O", filename, url], cwd=self._storage)
