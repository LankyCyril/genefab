#!/usr/bin/env python
from flask import Flask, Response
from genefab import GLDS, GeneLabJSONException
from pandas import DataFrame, concat
from json import dumps

app = Flask("genefab")

@app.route("/")
def hello_space():
    """Hello, Space!"""
    return "Hello, Space!"

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
    if rettype == "tsv":
        return Response(
            repr_df.to_csv(sep="\t", index=False, na_rep=""),
            mimetype="text/plain"
        )
    elif rettype == "html":
        return repr_df.to_html(index=False, na_rep="")
    else:
        return "400; bad request (wrong extension/type?)", 400
