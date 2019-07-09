#!/usr/bin/env python
from flask import Flask, Response
from genefab import GLDS, GeneLabJSONException
from pandas import DataFrame, concat, option_context
from json import dumps, JSONEncoder
from html import escape

app = Flask("genefab")

@app.route("/")
def hello_space():
    """Hello, Space!"""
    return "Hello, Space!"


class SetEnc(JSONEncoder):
    """Allow dumps to convert sets to serializable lists"""
    def default(self, entry):
        if isinstance(entry, set):
            return list(entry)
        else:
            return JSONEncoder.default(self, entry)


def display_object(obj, rettype, index=False):
    """Select appropriate converter and mimetype for rettype"""
    if isinstance(obj, (dict, tuple, list)):
        if rettype == "json":
            return Response(dumps(obj, cls=SetEnc), mimetype="text/json")
        else:
            obj = DataFrame(
                columns=["key", "value"],
                data=[[k, v] for k, vv in obj.items() for v in vv]
            )
            index = False
    if isinstance(obj, DataFrame):
        if rettype == "tsv":
            return Response(
                obj.to_csv(sep="\t", index=index, na_rep=""),
                mimetype="text/plain"
            )
        elif rettype == "html":
            with option_context("display.max_colwidth", -1):
                return obj.to_html(index=index, na_rep="", justify="left")
        elif rettype == "json":
            return Response(
                obj.to_json(index=True, orient="records"),
                mimetype="text/json"
            )
        else:
            return "400; bad request (wrong extension/type?)", 400
    else:
        mask = "501; not implemented: {} cannot be displayed"
        return mask.format(escape(str(type(obj)))), 501


def get_assay(accession, assay_name):
    """Get assay object via GLDS accession and assay name"""
    try:
        glds = GLDS(accession)
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


@app.route("/<accession>/<assay_name>/metadata.<rettype>")
def assay_metadata(accession, assay_name, rettype):
    assay, message, status = get_assay(accession, assay_name)
    if assay is None:
        return message, status
    if rettype == "raw":
        return display_object(assay.raw_metadata, "tsv", index=True)
    elif rettype == "pkl":
        return "501; not implemented: pickles coming soon", 501
    else:
        return display_object(assay.extended_raw_metadata, rettype, index=True)


@app.route("/<accession>/<assay_name>/metadata/<prop>.<rettype>")
def assay_summary(accession, assay_name, prop, rettype):
    """Provide overview of samples, fields, factors in metadata"""
    assay, message, status = get_assay(accession, assay_name)
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
