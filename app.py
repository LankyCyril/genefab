#!/usr/bin/env python
from flask import Flask
app = Flask("genefab")

@app.route("/")
def hello_world():
    return "Hello, world!"
