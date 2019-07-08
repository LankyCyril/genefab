#!/usr/bin/env python
from flask import Flask, Response
from genefab import GLDS
from json import dumps

app = Flask("genefab")

@app.route("/")
def hello_space():
    return "Hello, space!"

@app.route("/<accession>.json")
def glds_json(accession):
    glds = GLDS(accession)
    return Response(dumps([glds._json]), mimetype="text/json")

@app.route("/<accession>.tsv")
def glds_tsv(accession):
    glds = GLDS(accession)
    return Response(
        glds.assays._as_dataframe.to_csv(sep="\t"),
        mimetype="text/plain"
    )
