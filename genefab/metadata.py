from urllib.request import urlopen
from json import loads
from re import sub
from pandas import DataFrame, concat

URL_ROOT = "https://genelab-data.ndc.nasa.gov/genelab"

def get_json(url):
    """HTTP get, decode, parse"""
    with urlopen(url) as response:
        return loads(response.read().decode())

class GeneLabDataSet():
    accession = None
    _frame = None
    """"""
    def __init__(self, accession):
        """Request JSON representation of ISA metadata and store fields"""
        self.accession = accession
        getter_url = "{}/data/study/data/{}/"
        data_json = get_json(getter_url.format(URL_ROOT, accession))
        if len(data_json) > 1:
            raise ValueError("Too many results returned, unexpected behavior")
        self._json = data_json[0]
        try:
            self._internal_id = self._json["_id"]
            self._isa2json = self._json["foreignFields"][0]["isa2json"]
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
    def get_file_table(self):
        """Return DataFrame subset to filenames and factor values"""
        factor_fields = self.get_factors(as_fields=True)
        factor_dataframe = self.frame()[list(factor_fields.values())]
        factor_dataframe.columns = list(factor_fields.keys())
        datafiles_field_id = self._field_name_to_id("Array Data File")
        datafiles_names = self.frame()[datafiles_field_id]
        datafiles_names = datafiles_names.apply(
            lambda fn: "{}_microarray_{}.gz".format(
                self.accession, sub(r'^\*', "", fn)
            )
        )
        datafiles_names.name = "filename"
        datafiles = concat([factor_dataframe, datafiles_names], axis=1)
        return datafiles
    """"""
    def get_data_urls(self):
        """Return URLs of data files stored on GeneLab servers"""
        listing_url = "{}/data/study/filelistings/{}".format(
            URL_ROOT, self._internal_id
        )
        return {
            record["file_name"]: "{}/static/media/dataset/{}".format(
                URL_ROOT, record["file_name"]
            )
            for record in get_json(listing_url)
        }
