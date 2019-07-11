#!/usr/bin/env python
from flask import Flask, Response, request
from genefab import GLDS, GeneLabJSONException
from genefab._util import fetch_file
from pandas import DataFrame, concat, option_context
from json import dumps, JSONEncoder
from html import escape, unescape
from re import sub

app = Flask("genefab")

@app.route("/", methods=["GET"])
def hello_space():
    """Hello, Space!"""
    return "Hello, {}!".format(request.args.get("name", "Space"))


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


def display_object(obj, rettype, index="auto"):
    """Select appropriate converter and mimetype for rettype"""
    if isinstance(obj, (dict, tuple, list)):
        if rettype == "json":
            return Response(dumps(obj, cls=SetEnc), mimetype="text/json")
        else:
            obj, index = to_dataframe(obj), False
    if index == "auto":
        index = (rettype == "json")
    if isinstance(obj, DataFrame):
        if rettype == "list":
            if obj.shape == (1, 1):
                obj_repr = sub(r'\s*,\s*', "\n", str(obj.iloc[0, 0]))
                return Response(obj_repr, mimetype="text/plain")
            else:
                return "400; bad request (multiple cells selected)", 400
        elif rettype == "tsv":
            obj_repr = obj.to_csv(sep="\t", index=index, na_rep="")
            return Response(obj_repr, mimetype="text/plain")
        elif rettype == "html":
            with option_context("display.max_colwidth", -1):
                return obj.to_html(index=index, na_rep="", justify="left")
        elif rettype == "json":
            try:
                obj_repr = obj.to_json(index=index, orient="index")
                return Response(obj_repr, mimetype="text/json")
            except ValueError:
                obj_repr = obj.to_json(index=index, orient="records")
                return Response(obj_repr, mimetype="text/json")
        else:
            return "400; bad request (wrong extension/type?)", 400
    else:
        mask = "501; not implemented: {} cannot be displayed"
        return mask.format(escape(str(type(obj)))), 501


def get_bool(rargs, key, default_value):
    """Get boolean value from GET request args"""
    raw_value = rargs.get(key, default_value)
    return (raw_value not in {"False", "0"}) and bool(raw_value)


def get_assay(accession, assay_name, rargs):
    """Get assay object via GLDS accession and assay name"""
    spaces_in_sample_names = get_bool(rargs, "spaces_in_sample_names", True)
    try:
        glds = GLDS(accession, spaces_in_sample_names=spaces_in_sample_names)
    except GeneLabJSONException as e:
        return None, "404; not found: {}".format(e), 404
    if assay_name in glds.assays:
        return glds.assays[assay_name], None, 200
    else:
        mask = "404; not found: assay {} does not exist under {}"
        return None, mask.format(assay_name, accession), 404


@app.route("/<accession>.<rettype>")
def glds_summary(accession, rettype):
    """Report factors, assays, and/or raw JSON"""
    try:
        glds = GLDS(accession)
    except GeneLabJSONException as e:
        return "404; not found: {}".format(e), 404
    if rettype == "rawjson":
        return display_object([glds._json], "json")
    else:
        assays_df = glds.assays._as_dataframe.copy()
        assays_df.index.name = "name"
        assays_df["type"] = "assay"
        factors_df = DataFrame(
            columns=["type", "name", "factors"],
            data=[["dataset", accession, factor] for factor in glds.factors]
        )
        repr_df = concat(
            [factors_df, assays_df.reset_index()],
            axis=0, sort=False
        )
        return display_object(repr_df, rettype)


@app.route("/<accession>/<assay_name>/factors.<rettype>", methods=["GET"])
def assay_factors(accession, assay_name, rettype):
    """DataFrame of samples and factors in human-readable form"""
    assay, message, status = get_assay(accession, assay_name, request.args)
    if assay is None:
        return message, status
    else:
        return display_object(assay.factors, rettype, index=True)


@app.route("/<accession>/<assay_name>.<rettype>", methods=["GET", "POST"])
def assay_metadata(accession, assay_name, rettype):
    """DataFrame view of metadata, optionally queried"""
    assay, message, status = get_assay(accession, assay_name, request.args)
    if assay is None:
        return message, status
    elif request.args:
        fields = set(unescape(request.args.get("fields", "")).split()) - {""}
        index = set(unescape(request.args.get("index", "")).split()) - {""}
        if len(index) and len(fields):
            repr_df = assay.metadata.loc[list(index), list(fields)]
        elif index:
            repr_df = assay.metadata.loc[list(index)]
        elif fields:
            repr_df = assay.metadata[list(fields)]
        else:
            repr_df = assay.metadata.to_frame()
    else:
        repr_df = assay.metadata.to_frame()
    return display_object(repr_df, rettype, index=True)


@app.route("/<accession>/<assay_name>/<prop>.<rettype>", methods=["GET"])
def assay_summary(accession, assay_name, prop, rettype):
    """Provide overview of samples, fields, factors in metadata"""
    assay, message, status = get_assay(accession, assay_name, request.args)
    if assay is None:
        return message, status
    if prop == "fields":
        return display_object(assay._fields, rettype)
    elif prop == "index":
        return display_object(list(assay.raw_metadata.index), rettype)
    elif prop == "factors":
        return display_object(assay.factor_values, rettype)
    else:
        mask = "400; bad request: {} is not a valid property"
        return mask.format(prop), 400


@app.route("/<accession>/<assay_name>/file/<filemask>", methods=["GET"])
def locate_file(accession, assay_name, filemask):
    """Find file URL that matches filemask, redirect to download"""
    assay, message, status = get_assay(accession, assay_name, request.args)
    if assay is None:
        return message, status
    try:
        url = assay._get_file_url(filemask)
    except GeneLabJSONException:
        return "400; bad request: multiple files match mask", 400
    if url is None:
        return "404; not found: file not found", 404
    else:
        local_filepath = fetch_file(filemask, url, assay.storage)
        with open(local_filepath, mode="rb") as handle:
            return Response(handle.read(), mimetype="application/octet-stream")


@app.route("/<accession>/<assay_name>/<kind>_data.<rettype>", methods=["GET"])
def get_data(accession, assay_name, kind, rettype):
    """Serve up normalized/processed data"""
    assay, message, status = get_assay(accession, assay_name, request.args)
    if assay is None:
        return message, status
    translate_sample_names, data_columns = (
        get_bool(request.args, "translate_sample_names", False),
        request.args.get("data_columns", None)
    )
    if kind == "normalized":
        if translate_sample_names or data_columns:
            repr_df = assay.get_normalized_data(
                translate_sample_names=translate_sample_names,
                data_columns=data_columns
            )
        else:
            repr_df = assay.normalized_data
    elif kind in {"processed", "normalized_annotated"}:
        if translate_sample_names or data_columns:
            repr_df = assay.get_processed_data(
                translate_sample_names=translate_sample_names,
                data_columns=data_columns
            )
        else:
            repr_df = assay.processed_data
    else:
        return "400; bad request: unknown data request", 400
    return display_object(repr_df, rettype, index=True)
