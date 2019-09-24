from genefab import GLDS, GeneLabJSONException
from genefab._util import fetch_file, DELIM_AS_IS
from genefab._flaskutil import parse_rargs, display_object, ResponseError
from re import sub, split, search
from pandas import DataFrame, read_csv, Index, merge
from flask import Response


def get_assay(accession, assay_name, rargs):
    """Get assay object via GLDS accession and assay name"""
    rargdict = parse_rargs(rargs)
    try:
        glds = GLDS(accession, name_delim=rargdict["name_delim"])
    except GeneLabJSONException as e:
        return None, "404; not found: {}".format(e), 404
    if assay_name == "assay":
        if len(glds.assays) == 1:
            assay_name = list(glds.assays.keys())[0]
            return glds.assays[assay_name], None, 200
        else:
            mask = "400; bad request: ambiguous assay name {} for {}"
            return None, mask.format(assay_name, accession), 400
    elif assay_name in glds.assays:
        return glds.assays[assay_name], None, 200
    else:
        mask = "404; not found: assay {} does not exist under {}"
        return None, mask.format(assay_name, accession), 404


def subset_metadata(metadata, rargs):
    """Subset metadata by fields and index"""
    fields, index = rargs.get("fields", None), rargs.get("index", None)
    is_subset = True
    if fields and index:
        repr_df = metadata.loc[[index], [fields]]
    elif fields:
        repr_df = metadata[[fields]]
    elif index:
        repr_df = metadata.loc[[index]]
    else:
        repr_df = metadata.to_frame()
        is_subset = False
    return repr_df, is_subset


def filter_cells(subset, filename_filter):
    """Filter values in cells based on filename filter"""
    filtered_values = set()
    for cell in map(str, subset.values.flatten()):
        filenames = split(r'[,\s]+', cell)
        for filename in filenames:
            if search(filename_filter, filename):
                filtered_values.add(filename)
    return filtered_values


def melt_file_data(repr_df, melting):
    """Melt dataframe with the use of sample annotation"""
    foundry = repr_df.reset_index().copy()
    if isinstance(melting, (list, Index)):
        id_vars = [c for c in foundry.columns if c not in melting]
        return foundry.melt(
            id_vars=id_vars, value_vars=melting, var_name="Sample Name"
        )
    elif isinstance(melting, DataFrame):
        id_vars = [c for c in foundry.columns if c not in melting.columns]
        foundry = foundry.melt(
            id_vars=id_vars, value_vars=melting, var_name="Sample Name"
        )
        return merge(melting.T.reset_index(), foundry, how="outer")
    else:
        raise TypeError("cannot melt/describe with a non-dataframe object")


def serve_file_data(assay, filemask, rargs, melting=False):
    """Find file URL that matches filemask, redirect to download or interpret"""
    try:
        url = assay._get_file_url(filemask)
    except GeneLabJSONException:
        return ResponseError("multiple files match mask", 400)
    if url is None:
        return ResponseError("file not found", 404)
    local_filepath = fetch_file(filemask, url, assay.storage)
    rargdict = parse_rargs(rargs)
    if rargdict["fmt"] == "raw":
        if melting is not False:
            return ResponseError("cannot melt/describe raw object", 400)
        else:
            with open(local_filepath, mode="rb") as handle:
                return Response(handle.read(), mimetype="application")
    elif rargdict["fmt"] in {"tsv", "json"}:
        try:
            if local_filepath.endswith(".csv"):
                repr_df = read_csv(local_filepath, index_col=0)
            else:
                repr_df = read_csv(local_filepath, sep="\t", index_col=0)
        except Exception as e:
            return ResponseError(format(e), 400)
        if rargdict["name_delim"] != DELIM_AS_IS:
            repr_df.columns = repr_df.columns.map(
                lambda f: sub(r'[._-]', rargdict["name_delim"], f)
            )
        if melting is not False:
            try:
                repr_df = melt_file_data(repr_df, melting=melting)
            except Exception as e:
                return ResponseError(format(e), 400)
        else:
            repr_df = repr_df.reset_index()
        return display_object(repr_df, rargdict["fmt"], index="auto")
    else:
        return ResponseError("fmt={}".format(rargdict["fmt"]), 501)
