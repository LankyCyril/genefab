from genefab import GLDS, GeneLabJSONException, GeneLabException
from genefab._util import fetch_file, DELIM_AS_IS
from genefab._flaskutil import parse_rargs, display_object
from re import sub, split, search
from pandas import DataFrame, read_csv, Index, merge
from flask import Response
from csv import reader
from operator import __lt__, __le__, __eq__, __ne__, __ge__, __gt__


OPERATOR_MAPPER = {
    "<": __lt__, "<=": __le__, ">=": __ge__, ">": __gt__,
    "==": __eq__, "!=": __ne__
}


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


def get_filtered_repr_df(repr_df, value_filter_raw):
    """Interpret the filter request argument and subset the repr dataframe"""
    value_filter = sub(r'(^\')|(\'$)', "", value_filter_raw)
    indexer = None
    for field_filter in next(reader([value_filter])):
        match = search(r'(^[^<>=]+)([<>=]+)(.+)$', field_filter)
        if not match:
            raise ValueError("Malformed `filter`")
        field, comparison, value = match.groups()
        if comparison not in OPERATOR_MAPPER:
            error_mask = "Bad comparison: '{}'"
            raise ValueError(error_mask.format(comparison))
        else:
            compare = OPERATOR_MAPPER[comparison]
        if field not in repr_df.columns:
            error_mask = "Unknown field (column): '{}'"
            raise ValueError(error_mask.format(field))
        try:
            value = float(value)
        except ValueError:
            if value == "True":
                value = True
            elif value == "False":
                value = False
        try:
            field_indexer = compare(repr_df[field], value)
        except TypeError:
            emsk = "Invalid comparison (TypeError): {} {} {}"
            raise TypeError(emsk.format(field, comparison, value))
        if indexer is None:
            indexer = field_indexer
        else:
            indexer &= field_indexer
    return repr_df[indexer]


def serve_formatted_file_data(local_filepath, rargdict, melting):
    """Format file data accoring to rargdict and melting"""
    if rargdict["header"] == "1":
        read_kwargs = {"index_col": 0, "nrows": 1}
    else:
        read_kwargs = {"index_col": 0}
    if local_filepath.endswith(".csv"):
        repr_df = read_csv(local_filepath, **read_kwargs)
    else:
        repr_df = read_csv(local_filepath, sep="\t", **read_kwargs)
    if rargdict["name_delim"] != DELIM_AS_IS:
        convert_delim = lambda f: sub(r'[._-]', rargdict["name_delim"], f)
        repr_df.columns = repr_df.columns.map(convert_delim)
    if melting is not False:
        repr_df = melt_file_data(repr_df, melting=melting)
    else:
        repr_df = repr_df.reset_index()
    if rargdict["filter"] is not None:
        repr_df = get_filtered_repr_df(repr_df, rargdict["filter"])
    if rargdict["header"] == "1":
        repr_df = DataFrame(columns=repr_df.columns, index=["header"])
    return display_object(repr_df, rargdict["fmt"], index="auto")


def serve_file_data(assay, filemask, rargs, melting=False):
    """Find file URL that matches filemask, redirect to download or interpret"""
    try:
        url = assay._get_file_url(filemask)
    except GeneLabJSONException:
        raise ValueError("multiple files match mask")
    if url is None:
        raise FileNotFoundError
    local_filepath = fetch_file(filemask, url, assay.storage)
    rargdict = parse_rargs(rargs)
    if rargdict["fmt"] == "raw":
        if melting is not False:
            raise GeneLabException("cannot melt/describe raw object")
        elif rargdict["filter"] is not None:
            raise GeneLabException("cannot filter raw object")
        else:
            with open(local_filepath, mode="rb") as handle:
                return Response(handle.read(), mimetype="application")
    elif rargdict["fmt"] in {"tsv", "json"}:
        return serve_formatted_file_data(local_filepath, rargdict, melting)
    else:
        raise NotImplementedError("fmt={}".format(rargdict["fmt"]))
