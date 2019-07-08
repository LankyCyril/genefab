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


@app.route("/<accession>/<assay_name>.<rettype>")
def assay_summary(accession, assay_name, rettype):
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
    if rettype == "fields":
        return Response(dumps(assay._fields, cls=SetEnc), mimetype="text/json")
    if rettype == "index":
        return Response(
            dumps(list(assay.raw_metadata.index)), mimetype="text/json"
        )
    repr_list = (
        [["index", assay._indexed_by, ix] for ix in assay.raw_metadata.index] +
        [["field", nan, f] for f in assay._fields.keys()] +
        [["factor", k, v] for k, vs in assay.factor_values.items() for v in vs]
    )
    repr_df = DataFrame(data=repr_list, columns=["type", "name", "value"])
    return display_dataframe(repr_df, rettype)
