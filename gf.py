#!/usr/bin/env python
from sys import stderr
from flask import Flask, request
from genefab import GLDS, GeneLabJSONException, GeneLabException
from genefab import GeneLabDataManagerException
from pandas import DataFrame
from genefab._readme import html
from genefab._display import display_object
from genefab._util import parse_rargs
from genefab._bridge import get_assay, subset_metadata, filter_metadata_cells
from genefab._sqlite import retrieve_table_data, try_sqlite, dump_to_sqlite
from genefab._bridge import filter_table_data
from os import environ
from copy import deepcopy


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
    return html.format(url_root=request.url_root.rstrip("/"))


def exception_catcher(e):
    if isinstance(e, FileNotFoundError):
        code, explanation = 404, "Not Found"
    elif isinstance(e, NotImplementedError):
        code, explanation = 501, "Not Implemented"
    elif isinstance(e, GeneLabDataManagerException):
        code, explanation = 500, "GeneLab Data Manager Internal Server Error"
    else:
        code, explanation = 400, "Bad Request"
    error_mask = "<b>HTTP error</b>: {} ({})<br><b>{}</b>: {}"
    return error_mask.format(code, explanation, type(e).__name__, str(e)), code

if environ.get("FLASK_ENV", None) not in FLASK_DEBUG_MARKERS:
    exception_catcher = app.errorhandler(Exception)(exception_catcher)


@app.route("/favicon.<imgtype>")
def favicon(imgtype):
    """Catch request for favicons"""
    return ""


@app.route("/<accession>/", methods=["GET"])
def glds_summary(accession):
    """Report factors, assays, and/or raw JSON"""
    rargs = parse_rargs(request.args)
    try:
        glds = GLDS(accession)
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
    assay, message, status = get_assay(accession, assay_name, rargs)
    if assay is None:
        return message, status
    subset, _ = subset_metadata(assay.metadata, rargs)
    return display_object(subset, rargs.display_rargs, index=True)


@app.route("/<accession>/<assay_name>/factors/", methods=["GET"])
def assay_factors(accession, assay_name):
    """DataFrame of samples and factors in human-readable form"""
    rargs = parse_rargs(request.args)
    assay, message, status = get_assay(accession, assay_name, rargs)
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
    assay, message, status = get_assay(accession, assay_name, rargs)
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
def get_data(accession, assay_name, rargs=None):
    """Serve any kind of data"""
    if rargs is None:
        rargs = parse_rargs(request.args)
    assay, message, status = get_assay(accession, assay_name, rargs)
    if assay is None:
        return message, status
    subset, is_subset = subset_metadata(assay.metadata, rargs)
    if not is_subset: # no specific cells selected
        if "file_filter" not in rargs.data_rargs: # no filenames selected either
            raise ValueError("no entries selected")
        else: # filenames selected, should match just one
            filtered_values = filter_metadata_cells(
                DataFrame(assay.glds_file_urls.keys()),
                rargs.data_rargs["file_filter"]
            )
    else:
        filtered_values = filter_metadata_cells(
            subset, rargs.data_rargs["file_filter"]
        )
    if len(filtered_values) == 0:
        raise FileNotFoundError("no data")
    elif len(filtered_values) > 1:
        raise ValueError("multiple data files match search criteria")
    else:
        filename = filtered_values.pop()
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
    if rargs.display_rargs["fmt"] == "raw":
        raise NotImplementedError("fmt=raw with SQLite3")
    elif rargs.display_rargs["fmt"] in {"tsv", "json"}:
        return display_object(
            filter_table_data(table_data, rargs.data_filter_rargs),
            rargs.display_rargs, index="auto"
        )
    else:
        raise NotImplementedError("fmt={}".format(rargs.display_rargs["fmt"]))


def get_gct(accession, assay_name, rargs):
    """Get GCT formatted processed data; needs to be refactored!"""
    assay, message, status = get_assay(accession, assay_name, rargs)
    if assay is None:
        return message, status
    pdata = try_sqlite(accession, assay.name, rargs.data_rargs)
    if pdata is None:
        are_rargs_sane = (
            (rargs.display_rargs["fmt"] == "tsv") and
            (rargs.display_rargs["header"] == False) and
            (rargs.data_rargs["melted"] == False) and
            (rargs.data_rargs["descriptive"] == False)
        )
        if not are_rargs_sane:
            raise AttributeError(
                "None of the 'fmt', 'header', 'melted', 'descriptive' "
                "arguments make sense with the GCT format"
            )
        pdata = get_data(accession, assay.name, rargs=rargs)
    gct_header = "#1.2\n{}\t{}\n".format(*pdata.shape)
    pdata.columns = ["Description"] + list(pdata.columns)[1:]
    pdata.insert(loc=0, column="Name", value=pdata["Description"])
    gct_data = gct_header + pdata.to_csv(sep="\t", index=False)
    return display_object(gct_data, {"fmt": "raw"})


def get_data_alias_helper(accession, assay_name, data_type, rargs, transform=None):
    """Dispatch data for URL aliases"""
    if data_type in {"processed", "deg", "viz-table", "pca"}:
        are_fields_dirty = (rargs.data_rargs["fields"] is not None)
        is_file_filter_dirty = (rargs.data_rargs["file_filter"] != ".*")
        if are_fields_dirty or is_file_filter_dirty:
            error_mask = "'{}' already sets 'fields' and 'file_filter', {}"
            raise GeneLabException(error_mask.format(
                data_type, "cannot overwrite them with GET arguments"
            ))
    else:
        raise GeneLabException("Unknown data alias: '{}'".format(data_type))
    if transform in {"melted", "descriptive"}:
        for rarg in {"melted", "descriptive"}:
            if rargs.data_rargs[rarg]:
                error_mask = "'{}' already sets '{}', {}"
                raise GeneLabException(error_mask.format(
                    transform, rarg, "cannot overwrite it with GET arguments"
                ))
    elif transform == "gct":
        if data_type != "processed":
            raise ValueError("GCT only available for processed data")
    elif transform is not None:
        error_mask = "Unknown transformation alias: '{}'"
        raise GeneLabException(error_mask.format(transform))
    modified_rargs = deepcopy(rargs)
    if data_type == "processed":
        modified_rargs.data_rargs["fields"] = ".*normalized.*annotated.*"
        modified_rargs.data_rargs["file_filter"] = "txt"
    elif data_type == "deg":
        modified_rargs.data_rargs["fields"] = ".*differential.*expression.*"
        modified_rargs.data_rargs["file_filter"] = "expression.csv"
    elif data_type == "viz-table":
        modified_rargs.data_rargs["file_filter"] = ".*vis.*_output_table.csv"
    elif data_type == "pca":
        if transform == "melted":
            transform = None
        modified_rargs.data_rargs["file_filter"] = ".*vis.*_PCA_table.csv"
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
