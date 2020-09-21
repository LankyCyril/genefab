#!/usr/bin/env python
from sys import stderr
from flask import Flask, request
from flask_caching import Cache
from genefab import GLDS, GeneLabJSONException, GeneLabException
from genefab._readme import html
from genefab._display import display_object, traceback_printer, exception_catcher
from genefab._util import parse_rargs
from genefab._bridge import get_assay, subset_metadata, resolve_file_name, filter_table_data
from genefab._sqlite import retrieve_table_data, try_sqlite, dump_to_sqlite
from os import environ
from copy import deepcopy
from urllib.request import urlopen
from json import loads
from pandas import DataFrame


FLASK_DEBUG_MARKERS = {"development", "staging", "stage", "debug", "debugging"}
CACHE_CONFIG = {"CACHE_TYPE": "filesystem", "CACHE_DIR": ".genelab-ttl-cache"}
DEG_CSV_REGEX = r'^GLDS-[0-9]+_(array|rna_seq)(_all-samples)?_differential_expression.csv$'
VIZ_CSV_REGEX = r'^GLDS-[0-9]+_(array|rna_seq)(_all-samples)?_visualization_output_table.csv$'
PCA_CSV_REGEX = r'^GLDS-[0-9]+_(array|rna_seq)(_all-samples)?_visualization_PCA_table.csv$'


app = Flask("genefab")
cache = Cache(config=CACHE_CONFIG)
cache.init_app(app)


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


if environ.get("FLASK_ENV", None) in FLASK_DEBUG_MARKERS:
    traceback_printer = app.errorhandler(Exception)(traceback_printer)
else:
    exception_catcher = app.errorhandler(Exception)(exception_catcher)


@cache.memoize(timeout=60)
def get_json(url):
    """HTTP get, decode, parse"""
    with urlopen(url) as response:
        return loads(response.read().decode())


@app.route("/", methods=["GET"])
def hello_space():
    """Hello, Space!"""
    return html.format(url_root=request.url_root.rstrip("/"))


@app.route("/favicon.<imgtype>")
def favicon(imgtype):
    """Catch request for favicons"""
    return ""


@app.route("/<accession>/", methods=["GET"])
def glds_summary(accession):
    """Report factors, assays, and/or raw JSON"""
    rargs = parse_rargs(request.args)
    try:
        glds = GLDS(accession, get_json=get_json)
    except GeneLabJSONException as e:
        raise FileNotFoundError(e)
    if rargs.display_rargs["fmt"] == "raw":
        return display_object([glds._json], {"fmt": "json"})
    else:
        return display_object(
            glds.summary_dataframe, rargs.display_rargs
        )


@app.route("/<accession>/<assay_name>/", methods=["GET", "POST"])
def assay_metadata(accession, assay_name):
    """DataFrame view of metadata, optionally queried"""
    rargs = parse_rargs(request.args)
    assay, message, status = get_assay(accession, assay_name, rargs, get_json)
    if assay is None:
        return message, status
    subset, _ = subset_metadata(assay.metadata, rargs)
    return display_object(subset, rargs.display_rargs, index=True)


@app.route("/<accession>/<assay_name>/factors/", methods=["GET"])
def assay_factors(accession, assay_name):
    """DataFrame of samples and factors in human-readable form"""
    rargs = parse_rargs(request.args)
    assay, message, status = get_assay(accession, assay_name, rargs, get_json)
    if assay is None:
        return message, status
    else:
        if rargs.non_data_rargs["cls"]:
            if rargs.display_rargs["fmt"] != "tsv":
                error_mask = "{} format is unsuitable for CLS (use tsv)"
                raise GeneLabException(
                    error_mask.format(rargs.display_rargs["fmt"])
                )
            else:
                obj = assay.factors(
                    cls=rargs.non_data_rargs["cls"],
                    continuous=rargs.non_data_rargs["continuous"]
                )
                return display_object(obj, {"fmt": "raw"})
        else:
            return display_object(
                assay.factors(), rargs.display_rargs, index=True
            )


@app.route("/<accession>/<assay_name>/annotation/", methods=["GET"])
def assay_annotation(accession, assay_name):
    """DataFrame of samples and factors in human-readable form"""
    rargs = parse_rargs(request.args)
    assay, message, status = get_assay(accession, assay_name, rargs, get_json)
    if assay is None:
        return message, status
    else:
        if rargs.non_data_rargs["cls"]:
            if rargs.display_rargs["fmt"] != "tsv":
                error_mask = "{} format is unsuitable for CLS (use tsv)"
                raise GeneLabException(
                    error_mask.format(rargs.display_rargs["fmt"])
                )
            else:
                annotation = assay.annotation(
                    differential_annotation=rargs.non_data_rargs["diff"],
                    named_only=rargs.non_data_rargs["named_only"],
                    cls=rargs.non_data_rargs["cls"],
                    continuous=rargs.non_data_rargs["continuous"]
                )
                return display_object(annotation, {"fmt": "raw"})
        else:
            annotation = assay.annotation(
                differential_annotation=rargs.non_data_rargs["diff"],
                named_only=rargs.non_data_rargs["named_only"]
            )
            return display_object(
                annotation, rargs.display_rargs, index=True
            )


@app.route("/<accession>/<assay_name>/data/", methods=["GET"])
def get_data(accession, assay_name, rargs=None, return_raw=False):
    """Serve any kind of data"""
    if rargs is None:
        rargs = parse_rargs(request.args)
    assay, message, status = get_assay(accession, assay_name, rargs, get_json)
    if assay is None:
        return message, status
    filename = resolve_file_name(assay, rargs)
    table_data = try_sqlite(
        accession, assay.name, rargs.data_rargs,
        expect_date=assay.glds_file_dates.get(filename, -1)
    )
    if table_data is None:
        table_data = retrieve_table_data(assay, filename, rargs.data_rargs)
        dump_to_sqlite(
            accession, assay.name, rargs.data_rargs, table_data,
            set_date=assay.glds_file_dates.get(filename, -1)
        )
    filtered_table_data = filter_table_data(table_data, rargs.data_filter_rargs)
    if return_raw:
        return filtered_table_data
    elif rargs.display_rargs["fmt"] in {"tsv", "json"}:
        return display_object(
            filtered_table_data, rargs.display_rargs, index="auto"
        )
    else:
        raise NotImplementedError("fmt={}".format(rargs.display_rargs["fmt"]))


def get_gct(accession, assay_name, rargs):
    """Get GCT formatted processed data; needs to be refactored!"""
    assay, message, status = get_assay(accession, assay_name, rargs, get_json)
    if assay is None:
        return message, status
    pdata = get_data(accession, assay.name, rargs=rargs, return_raw=True)
    if isinstance(pdata, DataFrame):
        gct_header = "#1.2\n{}\t{}\n".format(pdata.shape[0], pdata.shape[1]-1)
        pdata.columns = ["Description"] + list(pdata.columns)[1:]
        pdata.insert(loc=0, column="Name", value=pdata["Description"])
        gct_data = gct_header + pdata.to_csv(sep="\t", index=False)
        return display_object(gct_data, {"fmt": "raw"})
    else:
        raise TypeError("Unexpected type: expected `DataFrame`")


def assess_data_alias(data_type, rargs, transform):
    """Checks if the URL alias is resolvable"""
    if data_type in {"processed", "deg", "viz-table", "pca"}:
        are_fields_dirty = (rargs.data_rargs["fields"] is not None)
        is_file_filter_dirty = (rargs.data_rargs["file_filter"] != ".*")
        if are_fields_dirty or is_file_filter_dirty:
            error_mask = "'{}' already sets 'fields' and 'file_filter', {}"
            return GeneLabException(error_mask.format(
                data_type, "cannot overwrite them with GET arguments"
            ))
    else:
        return GeneLabException("Unknown data alias: '{}'".format(data_type))
    if transform in {"melted", "descriptive"}:
        for rarg in {"melted", "descriptive"}:
            if rargs.data_rargs[rarg]:
                error_mask = "'{}' already sets '{}', {}"
                return GeneLabException(error_mask.format(
                    transform, rarg, "cannot overwrite it with GET arguments"
                ))
    elif transform == "gct":
        if data_type != "processed":
            return ValueError("GCT only available for processed data")
        are_rargs_sane = (
            (rargs.display_rargs["fmt"] == "tsv") and
            (rargs.display_rargs["header"] == False) and
            (rargs.data_rargs["melted"] == False) and
            (rargs.data_rargs["descriptive"] == False)
        )
        if not are_rargs_sane:
            return AttributeError(
                "None of the 'fmt', 'header', 'melted', 'descriptive' "
                "arguments make sense with the GCT format"
            )
    elif transform is not None:
        error_mask = "Unknown transformation alias: '{}'"
        return GeneLabException(error_mask.format(transform))


def get_data_alias_helper(accession, assay_name, data_type, rargs, transform=None):
    """Dispatch data for URL aliases"""
    AssessmentError = assess_data_alias(data_type, rargs, transform)
    if AssessmentError is not None:
        raise AssessmentError
    modified_rargs = deepcopy(rargs)
    if data_type == "processed":
        modified_rargs.data_rargs["fields"] = ".*normalized.*annotated.*"
        modified_rargs.data_rargs["file_filter"] = "txt"
    elif data_type == "deg":
        modified_rargs.data_rargs["fields"] = ".*differential.*expression.*"
        modified_rargs.data_rargs["file_filter"] = DEG_CSV_REGEX
    elif data_type == "viz-table":
        modified_rargs.data_rargs["file_filter"] = VIZ_CSV_REGEX
    elif data_type == "pca":
        if transform == "melted":
            transform = None
        modified_rargs.data_rargs["file_filter"] = PCA_CSV_REGEX
    if transform == "gct":
        return get_gct(accession, assay_name, modified_rargs)
    elif transform is not None:
        modified_rargs.data_rargs[transform] = True
    return get_data(accession, assay_name, rargs=modified_rargs)


@app.route("/<accession>/<assay_name>/data/<data_type>/", methods=["GET"])
def get_data_plain_alias(accession, assay_name, data_type):
    """Alias 'processed', 'deg', and 'viz-table' endpoints"""
    rargs = parse_rargs(request.args)
    return get_data_alias_helper(accession, assay_name, data_type, rargs)


@app.route("/<accession>/<assay_name>/data/<data_type>/<transform>/", methods=["GET"])
def get_data_transformed_alias(accession, assay_name, data_type, transform):
    """Alias 'melted', 'descriptive', and 'gct' endpoints"""
    rargs = parse_rargs(request.args)
    return get_data_alias_helper(
        accession, assay_name, data_type, rargs, transform
    )
