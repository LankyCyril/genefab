#!/usr/bin/env python
from flask import Flask

def hello_space():
    return "Hello, space!"

app = Flask("genefab")
app.route("/")(hello_space)
