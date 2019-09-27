#!/usr/bin/env python
from flask import Flask, request
from genefab import GLDS, GeneLabJSONException, GeneLabException
from pandas import DataFrame
from genefab._flaskutil import parse_rargs, ResponseError, display_object
from genefab._flaskbridge import get_assay, subset_metadata, filter_cells
from genefab._flaskbridge import serve_file_data
from re import sub

app = Flask("genefab")

@app.route("/", methods=["GET"])
def hello_space():
    """Hello, Space!"""
    return "Hello, {}!".format(request.args.get("name", "Space"))


@app.errorhandler(Exception)
def exception_catcher(e):
    e_repr = sub(r'[<>]', "", repr(type(e)))
    return "400; <b>{}</b>: {}".format(e_repr, str(e)), 400


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
        if rargdict["cls"]:
            if rargdict["fmt"] != "tsv":
                error_mask = "{} format is unsuitable for CLS (use tsv)"
                return ResponseError(error_mask.format(rargdict["fmt"]), 400)
            else:
                try:
                    obj = assay.factors(
                        cls=rargdict["cls"], continuous=rargdict["continuous"]
                    )
                except Exception as e:
                    return ResponseError(format(e), 400)
                else:
                    return display_object(obj, "raw")
        else:
            return display_object(assay.factors(), rargdict["fmt"], index=True)


@app.route("/<accession>/<assay_name>/annotation/", methods=["GET"])
def assay_annotation(accession, assay_name):
    """DataFrame of samples and factors in human-readable form"""
    assay, message, status = get_assay(accession, assay_name, request.args)
    if assay is None:
        return message, status
    else:
        rargdict = parse_rargs(request.args)
        if rargdict["cls"]:
            if rargdict["fmt"] != "tsv":
                error_mask = "{} format is unsuitable for CLS (use tsv)"
                return ResponseError(error_mask.format(rargdict["fmt"]), 400)
            else:
                try:
                    annotation = assay.annotation(
                        differential_annotation=rargdict["diff"],
                        named_only=rargdict["named_only"],
                        cls=rargdict["cls"], continuous=rargdict["continuous"]
                    )
                except GeneLabException as e:
                    return ResponseError(format(e), 400)
                return display_object(annotation, "raw")
        else:
            annotation = assay.annotation(
                differential_annotation=rargdict["diff"],
                named_only=rargdict["named_only"]
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
def get_data(accession, assay_name, rargs=None):
    """Serve any kind of data"""
    if rargs is None:
        rargs = request.args
    assay, message, status = get_assay(accession, assay_name, rargs)
    if assay is None:
        return message, status
    subset, is_subset = subset_metadata(assay.metadata, rargs)
    if not is_subset: # no specific cells selected
        if "file_filter" not in rargs: # no filenames selected either
            return ResponseError("no entries selected", 400)
        else: # filenames selected, should match just one
            filtered_values = filter_cells(
                DataFrame(assay.glds_file_urls.keys()), rargs["file_filter"]
            )
    else:
        filename_filter = rargs.get("file_filter", r'.*')
        filtered_values = filter_cells(subset, filename_filter)
    if len(filtered_values) == 0:
        return ResponseError("no data", 404)
    elif len(filtered_values) > 1:
        return ResponseError("multiple data files match search criteria", 400)
    else:
        fv = filtered_values.pop()
        try:
            if rargs.get("descriptive", "0") == "1":
                return serve_file_data(
                    assay, fv, rargs, melting=assay.annotation()
                )
            elif rargs.get("melted", "0") == "1":
                return serve_file_data(
                    assay, fv, rargs, melting=list(assay.annotation().columns)
                )
            else:
                return serve_file_data(assay, fv, rargs)
        except Exception as e:
            return ResponseError(format(e), 400)


def get_gct(accession, assay_name, rargs):
    """Get GCT formatted processed data"""
    assay, message, status = get_assay(accession, assay_name, rargs)
    if assay is None:
        return message, status
    rargdict = parse_rargs(rargs)
    are_rargs_sane = (
        (rargdict["fmt"] == "tsv") and (rargdict["header"] == "0") and
        (rargdict["melted"] == "0") and (rargdict["descriptive"] == "0")
    )
    if not are_rargs_sane:
        error_message = (
            "None of the 'fmt', 'header', 'melted', 'descriptive' " +
            "arguments make sense with the GCT format"
        )
        return ResponseError(error_message, 400)
    return display_object(assay.gct, fmt="raw")


def get_data_alias_helper(accession, assay_name, data_type, rargs, transformation_type=None):
    """Dispatch data for URL aliases"""
    if data_type == "processed":
        data_fields = {
            "fields": ".*normalized.*annotated.*",
            "file_filter": "txt"
        }
    elif data_type == "deg":
        data_fields = {
            "fields": ".*differential.*expression.*",
            "file_filter": "expression.csv"
        }
    elif data_type == "viz-table":
        data_fields = {
            "file_filter": ".*visualization_output_table.csv"
        }
    else:
        error_mask = "Unknown data alias: '{}'"
        return ResponseError(error_mask.format(data_type), 400)
    if transformation_type is None:
        query_fields = {**data_fields, **rargs}
    elif transformation_type == "melted":
        query_fields = {**data_fields, **{"melted": "1"}, **rargs}
    elif transformation_type == "descriptive":
        query_fields = {**data_fields, **{"descriptive": "1"}, **rargs}
    elif transformation_type == "gct":
        if data_type == "processed":
            return get_gct(accession, assay_name, rargs=rargs)
        else:
            return ResponseError("GCT only available for processed data", 400)
    else:
        error_mask = "Unknown transformation alias: '{}'"
        return ResponseError(error_mask.format(transformation_type), 400)
    return get_data(accession, assay_name, rargs=query_fields)


@app.route("/<accession>/<assay_name>/data/<data_type>/", methods=["GET"])
def get_data_plain_alias(accession, assay_name, data_type):
    """Alias 'processed', 'deg', and 'viz-table' endpoints"""
    return get_data_alias_helper(accession, assay_name, data_type, request.args)


@app.route("/<accession>/<assay_name>/data/<data_type>/<transformation_type>/", methods=["GET"])
def get_data_transformed_alias(accession, assay_name, data_type, transformation_type):
    """Alias 'melted', 'descriptive', and 'gct' endpoints"""
    return get_data_alias_helper(
        accession, assay_name, data_type, request.args, transformation_type
    )
