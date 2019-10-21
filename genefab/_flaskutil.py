from argparse import Namespace
from genefab._util import DELIM_DEFAULT
from copy import deepcopy
from flask import Response
from json import JSONEncoder, dumps
from pandas import DataFrame, option_context
from re import sub


DEFAULT_RARGS = Namespace(
    data_rargs = {
        "file_filter": ".*",
        "name_delim": DELIM_DEFAULT,
        "melted": False,
        "descriptive": False,
        "any_below": None,
    },
    data_filter_rargs = {
        "filter": None,
        "sort_by": None,
        "ascending": True,
    },
    display_rargs = {
        "fmt": "tsv",
        "header": False,
        "top": None,
        "cols": None,
        "excludecols": None,
    },
    non_data_rargs = {
        "diff": True,
        "named_only": True,
        "cls": None,
        "continuous": "infer",
    }
)


def parse_rargs(request_args):
    """Get all common arguments from request.args"""
    rargs = deepcopy(DEFAULT_RARGS)
    for rarg_type, rargs_of_type in DEFAULT_RARGS.__dict__.items():
        for rarg, rarg_default_value in rargs_of_type.items():
            if rarg in request_args:
                if not isinstance(rarg_default_value, bool):
                    getattr(rargs, rarg_type)[rarg] = request_args[rarg]
                elif request_args[rarg] == "0":
                    getattr(rargs, rarg_type)[rarg] = False
                else:
                    getattr(rargs, rarg_type)[rarg] = True
    return rargs


class SetEnc(JSONEncoder):
    """Allow dumps to convert sets to serializable lists"""
    def default(self, entry):
        if isinstance(entry, set):
            return list(entry)
        else:
            return JSONEncoder.default(self, entry)


def to_dataframe(obj):
    """Convert simple structured objects (dicts, list, tuples) to a DataFrame representation"""
    if isinstance(obj, dict):
        return DataFrame(
            columns=["key", "value"],
            data=[[k, v] for k, vv in obj.items() for v in vv]
        )
    else:
        return DataFrame(columns=["value"], data=obj)


def display_object(obj, fmt, index="auto"):
    """Select appropriate converter and mimetype for fmt"""
    if isinstance(obj, (dict, tuple, list)):
        if fmt == "json":
            return Response(dumps(obj, cls=SetEnc), mimetype="text/json")
        else:
            obj, index = to_dataframe(obj), False
    if index == "auto":
        index = (fmt == "json")
    if isinstance(obj, DataFrame):
        if fmt == "list":
            if obj.shape == (1, 1):
                obj_repr = sub(r'\s*,\s*', "\n", str(obj.iloc[0, 0]))
                return Response(obj_repr, mimetype="text/plain")
            else:
                raise ValueError("multiple cells selected")
        elif fmt == "tsv":
            obj_repr = obj.to_csv(sep="\t", index=index, na_rep="")
            return Response(obj_repr, mimetype="text/plain")
        elif fmt == "html":
            with option_context("display.max_colwidth", -1):
                return obj.to_html(index=index, na_rep="", justify="left")
        elif fmt == "json":
            try:
                obj_repr = obj.to_json(index=index, orient="index")
                return Response(obj_repr, mimetype="text/json")
            except ValueError:
                obj_repr = obj.to_json(index=index, orient="records")
                return Response(obj_repr, mimetype="text/json")
        else:
            raise ValueError("wrong extension or type?")
    elif fmt == "raw":
        return Response(obj, mimetype="application")
    else:
        return NotImplementedError(
            "{} cannot be displayed".format(type(obj).__name__)
        )
