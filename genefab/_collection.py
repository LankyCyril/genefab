from re import search, IGNORECASE
from urllib.parse import quote_plus
from sys import stderr
from ._dataset import GLDS
from ._util import get_json, URL_ROOT, FFIELD_ALIASES, FFIELD_VALUES

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

class GLDSCollection():
    """Implements collection of GLDS instances (generated from search terms)"""
    accessions = None
    _datasets = {}
 
    def __init__(self, **kwargs):
        """Match passed regexes and combine into search URL, store JSON"""
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
            self._json = get_json(url)["hits"]["hits"]
        except:
            raise ValueError("Unrecognized JSON structure")
        self.accessions = [
            hit["_id"] for hit in self._json
        ]
 
    def __getitem__(self, accession):
        if accession not in self.accessions:
            raise KeyError("No such dataset in collection")
        if accession not in self._datasets:
            self._datasets[accession] = GLDS(accession)
        return self._datasets[accession]
 
    def __iter__(self):
        for accession in self.accessions:
            yield self[accession]
 
    def __len__(self):
        return len(self.accessions or [])
