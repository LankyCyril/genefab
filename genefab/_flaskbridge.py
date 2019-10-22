from genefab import GLDS, GeneLabJSONException
from genefab._util import update_table, DELIM_AS_IS, STORAGE_PREFIX
from re import sub, split, search, IGNORECASE
from pandas import DataFrame, Index, merge, read_sql_query
from csv import reader
from operator import __lt__, __le__, __eq__, __ne__, __ge__, __gt__
from sqlite3 import connect
from hashlib import sha512
from pandas.io.sql import DatabaseError as PandasDatabaseError
from os import path


OPERATOR_MAPPER = {
    "<": __lt__, "<=": __le__, ">=": __ge__, ">": __gt__,
    "==": __eq__, "!=": __ne__
}


def get_assay(accession, assay_name, rargs):
    """Get assay object via GLDS accession and assay name"""
    try:
        glds = GLDS(accession, name_delim=rargs.data_rargs["name_delim"])
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
    fields, index = rargs.data_rargs["fields"], rargs.data_rargs["index"]
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


def filter_metadata_cells(subset, filename_filter):
    """Filter values in metadata cells based on filename filter"""
    filtered_values = set()
    for cell in map(str, subset.values.flatten()):
        filenames = split(r'[,\s]+', cell)
        for filename in filenames:
            if search(filename_filter, filename):
                filtered_values.add(filename)
    return filtered_values


def melt_table_data(repr_df, melting):
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


def get_padj_filtered_repr_df(repr_df, any_below):
    """Only pass entries where at least one Adj-p-value field is significant"""
    cutoff = float(any_below)
    filterable_fields = [
        field for field in repr_df.columns
        if search(r'^adj-p-value', field, flags=IGNORECASE)
    ]
    indexer = None
    for field in filterable_fields:
        if indexer is None:
            indexer = (repr_df[field] < cutoff)
        else:
            indexer |= (repr_df[field] < cutoff)
    return repr_df[indexer]


def retrieve_formatted_table_data(accession, assay_name, table_name, data_rargs, melting):
    """Format file data accoring to rargdict and melting"""
    db = connect(path.join( # FIXME
        STORAGE_PREFIX, accession + "-" + assay_name + ".sqlite3"
    ))
    query = "SELECT * FROM '{}'".format(table_name)
    try:
        repr_df = read_sql_query(query, db, index_col="index")
        repr_df = repr_df.set_index(repr_df.columns[0])
    except PandasDatabaseError:
        raise PandasDatabaseError("Expected a table but none found")
    if data_rargs["name_delim"] != DELIM_AS_IS:
        conv_delim = lambda f: sub(r'[._-]', data_rargs["name_delim"], f)
        repr_df.columns = repr_df.columns.map(conv_delim)
    if data_rargs["any_below"] is not None:
        repr_df = get_padj_filtered_repr_df(repr_df, data_rargs["any_below"])
    if melting is not False:
        repr_df = melt_table_data(repr_df, melting=melting)
    else:
        repr_df = repr_df.reset_index()
    return repr_df


def filter_table_data(repr_df, data_filter_rargs):
    """Filter dataframe"""
    if data_filter_rargs["filter"] is not None:
        repr_df = get_filtered_repr_df(repr_df, data_filter_rargs["filter"])
    if data_filter_rargs["sort_by"] is not None:
        if data_filter_rargs["sort_by"] in repr_df.columns:
            repr_df = repr_df.sort_values(
                by=data_filter_rargs["sort_by"],
                ascending=data_filter_rargs["ascending"]
            )
        else:
            error_mask = "Unknown field (column) '{}'"
            raise IndexError(error_mask.format(data_filter_rargs["sort_by"]))
    return repr_df


def retrieve_table_data(assay, filemask, data_rargs, melting=False):
    """Find file URL that matches filemask, redirect to download or interpret"""
    try:
        url = assay._get_file_url(filemask)
    except GeneLabJSONException:
        raise ValueError("multiple files match mask")
    if url is None:
        raise FileNotFoundError
    table_name = update_table( # FIXME
        assay.parent.accession, assay.name, filemask, url, assay.storage
    )
    return retrieve_formatted_table_data(
        assay.parent.accession, assay.name, table_name, data_rargs, melting
    )


def try_sqlite(accession, assay_name, url):
    """Try to load dataframe from DB_NAME"""
    return None # TODO FREEME FIXME
    db = connect(path.join(
        STORAGE_PREFIX, accession + "-" + assay_name + ".sqlite3"
    ))
    table_name = "flaskbridge-" + sha512(url.encode("utf-8")).hexdigest()
    query = "SELECT * FROM '{}'".format(table_name)
    try:
        return read_sql_query(query, db, index_col="index")
    except PandasDatabaseError:
        return None


def dump_to_sqlite(accession, assay_name, table_data, url):
    """Save transformed dataframe to DB_NAME"""
    return None # TODO FREEME FIXME
    table_name = "flaskbridge-" + sha512(url.encode("utf-8")).hexdigest()
    db = connect(path.join(
        STORAGE_PREFIX, accession + "-" + assay_name + ".sqlite3"
    ))
    try:
        table_data.to_sql(table_name, db)
    except ValueError:
        pass
    else:
        db.commit()
    db.close()
