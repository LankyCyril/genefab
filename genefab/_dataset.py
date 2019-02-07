from sys import stderr
from re import search, IGNORECASE
from urllib.parse import quote_plus
from ._util import get_json, FFIELD_ALIASES, FFIELD_VALUES, URL_ROOT
from pandas import DataFrame, concat


class Assay():
    """Stores individual assay metadata"""
    dataframe = None
 
    def __init__(self, assay_json):
        """Prase JSON into assay metadata"""
        self._json = assay_json
        self._raw, self._header = self._json["raw"], self._json["header"]
        field2title = {entry["field"]: entry["title"] for entry in self._header}
        if len(field2title) != len(self._header):
            raise ValueError("Conflicting IDs of data fields")
        converted_entries = []
        for entry in self._raw:
            converted_entry_data = [
                [field2title.get(field, field), value]
                for field, value in entry.items()
                if value
            ]
            converted_entry = DataFrame(converted_entry_data).set_index(0)
            converted_entry.index.name = None
            converted_entry = converted_entry[1].rename(None)
            converted_entries.append(converted_entry)
        self.dataframe = concat(converted_entries, axis=1).T
 
    def __repr__(self):
        return repr(self.dataframe)


class GeneLabDataSet():
    """Stores GLDS metadata associated with an accession number"""
    accession = None
    assays = []
 
    def __init__(self, accession):
        """Request JSON representation of ISA metadata and store fields"""
        self.accession = accession
        getter_url = "{}/data/study/data/{}/"
        data_json = get_json(getter_url.format(URL_ROOT, accession))
        if len(data_json) > 1:
            raise ValueError("Too many results returned, unexpected behavior")
        else:
            self._json = data_json[0]
        try:
            self.internal_id = self._json["_id"]
            self.metadata_id = self._json["metadata_id"]
            self._isa2json = self._json["foreignFields"][0]["isa2json"]
            self._info = self._isa2json["additionalInformation"]
            for field in "description", "samples", "ontologies", "organisms":
                setattr(self, field, self._info[field])
        except KeyError:
            raise ValueError("Malformed JSON")
        self.assays = [
            Assay(assay_json)
            for assay_json in self._info["assays"].values()
        ]


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
    term_pairs = [
        "ffield={}&fvalue={}".format(ffield, quote_plus(ffvalue))
        for ffield, ffvalue in get_ffield_matches(**kwargs)
    ]
    url = "&".join(
        [URL_ROOT+"/data/search/?term=GLDS", "type=cgene", "size="+maxcount]
        + term_pairs
    )
    try:
        json = get_json(url)["hits"]["hits"]
    except:
        raise ValueError("Unrecognized JSON structure")
    return [
        GeneLabDataSet(hit["_id"]) for hit in json
    ]
