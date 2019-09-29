from genefab._util import DELIM_DEFAULT
from flask import Response
from json import JSONEncoder, dumps
from pandas import DataFrame, option_context
from re import sub


def parse_rargs(rargs):
    """Get all common arguments from request.args"""
    return {
        "fmt": rargs.get("fmt", "tsv"),
        "name_delim": rargs.get("name_delim", DELIM_DEFAULT),
        "header": rargs.get("header", "0"),
        "melted": rargs.get("melted", "0"),
        "descriptive": rargs.get("descriptive", "0"),
        "file_filter": rargs.get("file_filter", ".*"),
        "filter": rargs.get("filter", None),
        "cls": rargs.get("cls", None),
        "continuous": rargs.get("continuous", "infer")
    }


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
