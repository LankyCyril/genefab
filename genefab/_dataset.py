from sys import stderr
from re import search, IGNORECASE
from urllib.parse import quote_plus
from ._util import get_json, FFIELD_ALIASES, FFIELD_VALUES, URL_ROOT

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
        hit["_id"] for hit in json
    ]

GeneLabDataSet = None
