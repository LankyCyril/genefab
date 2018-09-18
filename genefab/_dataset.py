from re import sub
from pandas import DataFrame, concat
from functools import lru_cache
from ._util import get_json, fetch_file, URL_ROOT, DEFAULT_STORAGE

class GLDS():
    """Implements single dataset interface (generated from accession id)"""
    accession = None
 
    def __init__(self, accession, storage=DEFAULT_STORAGE):
        """Request JSON representation of ISA metadata and store fields"""
        self.accession = accession
        self._storage = storage
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
 
    def _is_consistent(self, raw):
        """Check if keys are the same for each record in _raw"""
        key_set = set(raw[0].keys())
        for record in raw:
            if set(record.keys()) != key_set:
                return False
        else:
            return True
 
    @lru_cache(maxsize=None)
    def field_ids(self, field_name, column_name=None):
        """Convert external field name to internal field id"""
        fields = []
        for record in self._header:
            if record["title"] == field_name:
                if column_name is not None:
                    if "columns" in record:
                        for column in record["columns"]:
                            if "field" in column:
                                if column.get("title", None) == column_name:
                                    fields.append(column["field"])
                elif "field" in record:
                    fields.append(record["field"])
        return fields
 
    @property
    @lru_cache(maxsize=None)
    def frame(self):
        """Convert _raw field of _isa2json to pandas DataFrame"""
        if not self._is_consistent(self._raw):
            raise KeyError("_raw field keys are inconsistent")
        _frame = DataFrame(
            columns=sorted(self._raw[0].keys()),
            index=range(len(self._raw))
        )
        for i, record in enumerate(self._raw):
            for key, value in record.items():
                _frame.loc[i, key] = value
        sample_names_field_ids = self.field_ids("Sample Name")
        if len(sample_names_field_ids) != 1:
            raise ValueError("Number of 'Sample Name' fields is not 1")
        return _frame.set_index(sample_names_field_ids[0])
 
    @lru_cache(maxsize=None)
    def field_values(self, field_name, column_name=None):
        """Convert external field name to internal field id, return set of possible values for the field"""
        return set.union(*(
            set(self.frame[field_id])
            for field_id in self.field_ids(field_name, column_name)
        ), set())
 
    @lru_cache(maxsize=None)
    def factors(self, *, as_fields):
        """Get factor type from _header"""
        _factors = {}
        for record in self._header:
            if record["title"] == "Factor Value":
                for column in record["columns"]:
                    factor = column["title"]
                    if as_fields:
                        values = column["field"]
                    else:
                        values = set(self.frame[column["field"]].values)
                    _factors[factor] = values
        return _factors
 
    @lru_cache(maxsize=None)
    def property_table(self, field_name):
        """Return DataFrame subset to filenames and factor values"""
        factor_fields = self.factors(as_fields=True)
        if not factor_fields:
            return None
        factor_dataframe = self.frame[list(factor_fields.values())]
        factor_dataframe.columns = list(factor_fields.keys())
        property_field_ids = self.field_ids(field_name)
        if len(property_field_ids) != 1:
            raise ValueError("Number of '{}' ids is not 1".format(field_name))
        property_names = self.frame[property_field_ids[0]]
        property_names.name = field_name
        return concat([factor_dataframe, property_names], axis=1)
 
    @property
    def has_raw_arrays(self):
        return (len(self.field_values("Array Data File")) != 0)
    @property
    def has_derived_arrays(self):
        return (len(self.field_values("Derived Array Data File")) != 0)
    @property
    def is_microarray(self):
        return self.has_raw_arrays or self.has_derived_arrays
 
    @property
    @lru_cache(maxsize=None)
    def file_list(self):
        """Return names and URLs of data files stored on GeneLab servers"""
        listing_url = "{}/data/study/filelistings/{}".format(
            URL_ROOT, self._internal_id
        )
        return {
            record["file_name"]: "{}/static/media/dataset/{}".format(
                URL_ROOT, record["file_name"]
            )
            for record in get_json(listing_url)
        }
 
    def fetch_files(self, update=False):
        """Alias for file_list(self, fetch=True, update=update)"""
        for file_name, url in self.file_list.items():
            fetch_file(file_name, url, self._storage, update=update)
