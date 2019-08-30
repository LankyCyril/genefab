#!/usr/bin/env python
from flask import Flask, Response, request
from genefab import GLDS, GeneLabJSONException
from genefab._util import fetch_file, DELIM_AS_IS, DELIM_DEFAULT
from pandas import DataFrame, option_context, read_csv, merge
from json import dumps, JSONEncoder
from html import escape
from re import sub, split, search

app = Flask("genefab")

@app.route("/", methods=["GET"])
def hello_space():
    """Hello, Space!"""
    return "Hello, {}!".format(request.args.get("name", "Space"))


def ResponseError(mask, code, *args):
    """Generate HTTP error code and message"""
    if code == 400: explanation = "; bad request"
    elif code == 404: explanation = "; not found"
    elif code == 501: explanation = "; not implemented"
    else: explanation = ""
    if args:
        converted_args = [
            escape(str(arg)) if type(arg) == type(type) else arg for arg in args
        ]
        message = str(code) + explanation + ": " + mask.format(*converted_args)
    else:
        message = str(code) + explanation + ": " + mask
    return message, code


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
                return ResponseError("multiple cells selected", 400)
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
            return ResponseError("wrong extension or type?", 400)
    else:
        return ResponseError("{} cannot be displayed", 501, type(obj))


def parse_rargs(rargs):
    """Get all common arguments from request.args"""
    return {
        "fmt": rargs.get("fmt", "tsv"),
        "name_delim": rargs.get("name_delim", DELIM_DEFAULT)
    }


def get_assay(accession, assay_name, rargs):
    """Get assay object via GLDS accession and assay name"""
    rargdict = parse_rargs(rargs)
    try:
        glds = GLDS(accession, name_delim=rargdict["name_delim"])
    except GeneLabJSONException as e:
        return None, "404; not found: {}".format(e), 404
    if assay_name == "assay":
        if len(glds.assays) == 1:
            assay_name = list(glds.assays.keys())[0]
            return glds.assays[assay_name], None, 200
        else:
            mask = "400; bad request: ambiguous assay name {} for {}"
            return None, mask.format(assay_name, accession), 400
    elif assay_name in glds.assays:
        return glds.assays[assay_name], None, 200
    else:
        mask = "404; not found: assay {} does not exist under {}"
        return None, mask.format(assay_name, accession), 404


def subset_metadata(metadata, rargs):
    """Subset metadata by fields and index"""
    fields, index = rargs.get("fields", None), rargs.get("index", None)
    is_subset = True
    if fields and index:
        repr_df = metadata.loc[[index], [fields]]
    elif fields:
        repr_df = metadata[[fields]]
    elif index:
        repr_df = metadata.loc[[index]]
    else:
        repr_df = metadata.to_frame()
        is_subset = False
    return repr_df, is_subset


def filter_cells(subset, filename_filter):
    """Filter values in cells based on filename filter"""
    filtered_values = set()
    for cell in map(str, subset.values.flatten()):
        filenames = split(r'[,\s]+', cell)
        for filename in filenames:
            if search(filename_filter, filename):
                filtered_values.add(filename)
    return filtered_values


def melt_file_data(repr_df, melting):
    """Melt dataframe with the use of sample annotation"""
    foundry = repr_df.T
    foundry.index.name = "Sample Name"
    foundry = foundry.reset_index().melt(id_vars="Sample Name")
    if melting is True:
        return foundry
    else:
        return merge(melting.T.reset_index(), foundry, how="outer")


def serve_file_data(assay, filemask, rargs, melting=False):
    """Find file URL that matches filemask, redirect to download or interpret"""
    try:
        url = assay._get_file_url(filemask)
    except GeneLabJSONException:
        return ResponseError("multiple files match mask", 400)
    if url is None:
        return ResponseError("file not found", 404)
    local_filepath = fetch_file(filemask, url, assay.storage)
    rargdict = parse_rargs(rargs)
    if rargdict["fmt"] == "raw":
        if melting is not False:
            return ResponseError("cannot melt/annotate raw object", 400)
        else:
            with open(local_filepath, mode="rb") as handle:
                return Response(handle.read(), mimetype="application")
    elif rargdict["fmt"] in {"tsv", "json"}:
        try:
            repr_df = read_csv(local_filepath, sep="\t", index_col=0)
        except Exception as e:
            return ResponseError(format(e), 400)
        if rargdict["name_delim"] != DELIM_AS_IS:
            repr_df.columns = repr_df.columns.map(
                lambda f: sub(r'[._-]', rargdict["name_delim"], f)
            )
        if melting is not False:
            repr_df = melt_file_data(repr_df, melting=melting)
        return display_object(repr_df, rargdict["fmt"], index=True)
    else:
        return ResponseError("fmt={}".format(rargdict["fmt"]), 501)


@app.route("/<accession>/", methods=["GET"])
def glds_summary(accession):
    """Report factors, assays, and/or raw JSON"""
    try:
        glds = GLDS(accession)
    except GeneLabJSONException as e:
        return ResponseError("{}", 404, format(e))
    rargdict = parse_rargs(request.args)
    if rargdict["fmt"] == "raw":
        return display_object([glds._json], "json")
    else:
        return display_object(glds._summary_dataframe, rargdict["fmt"])


@app.route("/<accession>/<assay_name>/factors/", methods=["GET"])
def assay_factors(accession, assay_name):
    """DataFrame of samples and factors in human-readable form"""
    assay, message, status = get_assay(accession, assay_name, request.args)
    if assay is None:
        return message, status
    else:
        rargdict = parse_rargs(request.args)
        return display_object(assay.factors, rargdict["fmt"], index=True)


@app.route("/<accession>/<assay_name>/annotation/", methods=["GET"])
def assay_annotation(accession, assay_name):
    """DataFrame of samples and factors in human-readable form"""
    assay, message, status = get_assay(accession, assay_name, request.args)
    if assay is None:
        return message, status
    else:
        rargdict = parse_rargs(request.args)
        annotation = assay.annotation(
            differential_annotation=rargdict.get("diff", True),
            named_only=rargdict.get("named_only", True)
        )
        return display_object(annotation, rargdict["fmt"], index=True)


@app.route("/<accession>/<assay_name>/", methods=["GET", "POST"])
def assay_metadata(accession, assay_name):
    """DataFrame view of metadata, optionally queried"""
    assay, message, status = get_assay(accession, assay_name, request.args)
    if assay is None:
        return message, status
    subset, _ = subset_metadata(assay.metadata, request.args)
    rargdict = parse_rargs(request.args)
    return display_object(subset, rargdict["fmt"], index=True)


@app.route("/<accession>/<assay_name>/<prop>/", methods=["GET"])
def assay_summary(accession, assay_name, prop):
    """Provide overview of samples, fields, factors in metadata"""
    assay, message, status = get_assay(accession, assay_name, request.args)
    if assay is None:
        return message, status
    rargdict = parse_rargs(request.args)
    if prop == "fields":
        return display_object(assay._fields, rargdict["fmt"])
    elif prop == "index":
        return display_object(list(assay.raw_metadata.index), rargdict["fmt"])
    else:
        return ResponseError("{} is not a valid property", 400, prop)


@app.route("/<accession>/<assay_name>/data/", methods=["GET"])
def get_data(accession, assay_name):
    """Serve any kind of data"""
    assay, message, status = get_assay(accession, assay_name, request.args)
    if assay is None:
        return message, status
    subset, is_subset = subset_metadata(assay.metadata, request.args)
    if not is_subset:
        return ResponseError("no entries selected", 400)
    filename_filter = request.args.get("filter", r'.*')
    filtered_values = filter_cells(subset, filename_filter)
    if len(filtered_values) == 0:
        return ResponseError("no data", 404)
    elif len(filtered_values) > 1:
        return ResponseError("multiple data files match search criteria", 400)
    elif request.args.get("annotated", "0") == "1":
        return serve_file_data(
            assay, filtered_values.pop(), request.args,
            melting=assay.annotation()
        )
    elif request.args.get("melted", "0") == "1":
        return serve_file_data(
            assay, filtered_values.pop(), request.args,
            melting=True
        )
    else:
        return serve_file_data(assay, filtered_values.pop(), request.args)
