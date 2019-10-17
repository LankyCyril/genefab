#!/usr/bin/env python
from sys import stderr
from flask import Flask, request
from genefab import GLDS, GeneLabJSONException, GeneLabException
from pandas import DataFrame
from genefab._flaskutil import parse_rargs, display_object
from genefab._flaskbridge import get_assay, subset_metadata, filter_cells
from genefab._flaskbridge import serve_file_data, try_sqlite, dump_to_sqlite
from os import environ
from re import sub
from io import BytesIO


FLASK_DEBUG_MARKERS = {"development", "staging", "stage", "debug", "debugging"}
app = Flask("genefab")

try:
    from flask_compress import Compress
    COMPRESS_MIMETYPES = [
        "text/plain", "text/html", "text/css", "text/xml",
        "application/json", "application/javascript"
    ]
    Compress(app)
except Exception as e: # I'm sorry
    print("Warning: Could not apply auto-compression", file=stderr)
    print("The error was:", e, file=stderr)


@app.route("/", methods=["GET"])
def hello_space():
    """Hello, Space!"""
    return "Hello, {}!".format(request.args.get("name", "Space"))


def exception_catcher(e):
    if isinstance(e, FileNotFoundError):
        code, explanation = 404, "Not Found"
    elif isinstance(e, NotImplementedError):
        code, explanation = 501, "Not Implemented"
    else:
        code, explanation = 400, "Bad Request"
    error_mask = "<b>HTTP error</b>: {} ({})<br><b>{}</b>: {}"
    return error_mask.format(code, explanation, type(e).__name__, str(e)), code

if environ.get("FLASK_ENV", None) not in FLASK_DEBUG_MARKERS:
    exception_catcher = app.errorhandler(Exception)(exception_catcher)
else:
    try:
        from flask_cors import CORS
        CORS(app)
    except ModuleNotFoundError:
        print("No module flask_cors (in FLASK_ENV==development)", file=stderr)
        pass


@app.route("/<accession>/", methods=["GET"])
def glds_summary(accession):
    """Report factors, assays, and/or raw JSON"""
    try:
        glds = GLDS(accession)
    except GeneLabJSONException as e:
        raise FileNotFoundError(e)
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
                return GeneLabException(error_mask.format(rargdict["fmt"]), 400)
            else:
                obj = assay.factors(
                    cls=rargdict["cls"], continuous=rargdict["continuous"]
                )
                return display_object(obj, "raw")
        else:
            return display_object(assay.factors(), rargdict["fmt"], index=True)


@app.route("/<accession>/<assay_name>/factors/cls/", methods=["GET"])
def assay_factors_cls(accession, assay_name):
    """DataFrame of samples and factors in CLS format"""
    assay, message, status = get_assay(accession, assay_name, request.args)
    if assay is None:
        return message, status
    else:
        rargdict = parse_rargs(request.args)
        if rargdict["fmt"] != "tsv":
            error_mask = "{} format is unsuitable for CLS (use tsv)"
            return GeneLabException(error_mask.format(rargdict["fmt"]), 400)
        else:
            if rargdict["cls"]:
                obj = assay.factors(
                    cls=rargdict["cls"],
                    continuous=rargdict["continuous"]
                )
            else:
                obj = assay.factors(cls="*", continuous=rargdict["continuous"])
            return display_object(obj, "raw")


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
                raise GeneLabException(error_mask.format(rargdict["fmt"]))
            else:
                annotation = assay.annotation(
                    differential_annotation=rargdict["diff"],
                    named_only=rargdict["named_only"],
                    cls=rargdict["cls"], continuous=rargdict["continuous"]
                )
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
        raise AttributeError("{} is not a valid property".format(prop))


@app.route("/<accession>/<assay_name>/data/", methods=["GET"])
def get_data(accession, assay_name, rargs=None):
    """Serve any kind of data"""
    assay, message, status = get_assay(accession, assay_name, rargs)
    if assay is None:
        return message, status
    if rargs is None:
        rargs = request.args
    file_data = try_sqlite(accession, assay.name, request.url, rargs)
    if file_data is not None:
        return file_data
    subset, is_subset = subset_metadata(assay.metadata, rargs)
    if not is_subset: # no specific cells selected
        if "file_filter" not in rargs: # no filenames selected either
            raise ValueError("no entries selected")
        else: # filenames selected, should match just one
            filtered_values = filter_cells(
                DataFrame(assay.glds_file_urls.keys()), rargs["file_filter"]
            )
    else:
        filename_filter = rargs.get("file_filter", r'.*')
        filtered_values = filter_cells(subset, filename_filter)
    if len(filtered_values) == 0:
        raise FileNotFoundError("no data")
    elif len(filtered_values) > 1:
        raise ValueError("multiple data files match search criteria")
    else:
        fv = filtered_values.pop()
        if rargs.get("descriptive", "0") == "1":
            file_data = serve_file_data(
                assay, fv, rargs, melting=assay.annotation().T
            )
        elif rargs.get("with_factors", "0") == "1":
            file_data = serve_file_data(
                assay, fv, rargs, melting=assay.factors().T
            )
        elif rargs.get("melted", "0") == "1":
            file_data = serve_file_data(
                assay, fv, rargs, melting=list(assay.annotation().T.columns)
            )
        else:
            file_data = serve_file_data(assay, fv, rargs)
    dump_to_sqlite(accession, assay.name, file_data, request.url)
    return file_data


def get_gct(accession, assay_name, rargs):
    """Get GCT formatted processed data; needs to be refactored!"""
    assay, message, status = get_assay(accession, assay_name, rargs)
    if assay is None:
        return message, status
    pdata_request_url = sub(r'/gct/$', "/", request.url)
    processed_file_data = try_sqlite(
        accession, assay.name, pdata_request_url, rargs
    )
    if processed_file_data is None:
        rargdict = parse_rargs(rargs)
        are_rargs_sane = (
            (rargdict["fmt"] == "tsv") and (rargdict["header"] == "0") and
            (rargdict["melted"] == "0") and (rargdict["descriptive"] == "0")
        )
        if not are_rargs_sane:
            raise AttributeError(
                "None of the 'fmt', 'header', 'melted', 'descriptive' "
                "arguments make sense with the GCT format"
            )
        processed_file_data = get_data(accession, assay.name, rargs=rargs)
    pdata_bytes = processed_file_data.data
    converted_lines, nrows, ncols = [], None, None
    for byteline in map(bytes.strip, BytesIO(pdata_bytes)):
        if nrows is None:
            converted_lines.append(
                sub(br'^[^\t]+', b'Name\tDescription', byteline).decode()
            )
            nrows, ncols = 0, byteline.count(b'\t')
        elif byteline:
            converted_lines.append(
                sub(br'^([^\t]+)', br'\1\t\1', byteline).decode()
            )
            nrows += 1
    gct_header = "#1.2\n{}\t{}\n".format(nrows, ncols)
    file_data = gct_header + "\n".join(converted_lines)
    return display_object(file_data, fmt="raw")


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
        raise GeneLabException("Unknown data alias: '{}'".format(data_type))
    if transformation_type is None:
        query_fields = {**data_fields, **rargs}
    elif transformation_type == "melted":
        query_fields = {**data_fields, **{"melted": "1"}, **rargs}
    elif transformation_type == "with-factors":
        query_fields = {**data_fields, **{"with_factors": "1"}, **rargs}
    elif transformation_type == "descriptive":
        query_fields = {**data_fields, **{"descriptive": "1"}, **rargs}
    elif transformation_type == "gct":
        if data_type == "processed":
            return get_gct(accession, assay_name, {**data_fields, **rargs})
        else:
            raise ValueError("GCT only available for processed data")
    else:
        error_mask = "Unknown transformation alias: '{}'"
        raise GeneLabException(error_mask.format(transformation_type))
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
