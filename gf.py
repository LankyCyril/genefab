#!/usr/bin/env python
from flask import Flask, Response
from genefab import GLDS, GeneLabJSONException
from pandas import DataFrame, concat
from json import dumps, JSONEncoder
from numpy import nan

app = Flask("genefab")

@app.route("/")
def hello_space():
    """Hello, Space!"""
    return "Hello, Space!"


def display_dataframe(df, rettype):
    """Select appropriate converter and mimetype for rettype"""
    if rettype == "tsv":
        return Response(
            df.to_csv(sep="\t", index=False, na_rep=""), mimetype="text/plain"
        )
    elif rettype == "html":
        return df.to_html(index=False, na_rep="")
    else:
        return "400; bad request (wrong extension/type?)", 400


class SetEnc(JSONEncoder):
    """Allow dumps to convert sets to serializable lists"""
    def default(self, entry):
        if isinstance(entry, set):
            return list(entry)
        else:
            return JSONEncoder.default(self, entry)


def get_assay(accession, assay_name):
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
    if rettype == "json":
        return Response(dumps([glds._json]), mimetype="text/json")
    assays_df = glds.assays._as_dataframe.copy()
    assays_df.index.name = "name"
    assays_df["type"] = "assay"
    factors_df = DataFrame(
        columns=["type", "name", "factors"],
        data=[["dataset", accession, factor] for factor in glds.factors]
    )
    repr_df = concat([factors_df, assays_df.reset_index()], axis=0, sort=False)
    return display_dataframe(repr_df, rettype)


@app.route("/<accession>/<assay_name>/<selection>")
def assay_summary(accession, assay_name, selection):
    """Provide overview of samples, fields, factors in metadata"""
    try:
        glds = GLDS(accession)
    except GeneLabJSONException as e:
        return "404; not found: {}".format(e), 404
    if assay_name in glds.assays:
        assay = glds.assays[assay_name]
    else:
        mask = "404; not found: assay {} does not exist under {}"
        return mask.format(assay_name, accession), 404
    if selection == "fields":
        return Response(dumps(assay._fields, cls=SetEnc), mimetype="text/json")
    elif selection == "index":
        return Response(
            dumps(list(assay.raw_metadata.index)), mimetype="text/json"
        )
    elif selection == "factors":
        return Response(
            dumps(assay.factor_values, cls=SetEnc), mimetype="text/json"
        )
    else:
        mask = "400; bad request: {} is not a valid selection"
        return mask.format(selection), 400


@app.route("/<accession>/<assay_name>/metadata.<selection>")
def assay_metadata(accession, assay_name, selection):
    assay, message, status = get_assay(accession, assay_name)
    if assay is None:
        return message, status
    if selection == "raw":
        return Response(
            assay.raw_metadata.to_csv(sep="\t"), mimetype="text/plain"
        )
    elif selection == "pkl":
        return "501; not implemented: Pickles coming soon", 501
    else:
        mask = "400; bad request: {} is not a valid selection"
        return mask.format(selection), 400
