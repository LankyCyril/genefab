from genefab import GeneLabJSONException
from os.path import join
from os import remove, path
from requests import get
from requests.exceptions import InvalidSchema
from urllib.error import URLError
from sqlite3 import connect
from hashlib import sha512
from pandas import read_csv, read_sql_query, DataFrame, Index, merge
from pandas.io.sql import DatabaseError as PandasDatabaseError
from tempfile import TemporaryDirectory
from genefab._util import STORAGE_PREFIX, guess_format, DELIM_AS_IS
from re import sub, search, IGNORECASE


def update_table(accession, assay_name, filemask, url, table_prefix, verbose=False, http_fallback=True):
    """Check if table already in database, if not, download and store; return table name"""
    db = connect(join(
        STORAGE_PREFIX, accession + "-" + assay_name + ".sqlite3"
    ))
    filemask_hash = sha512(filemask.encode("utf-8")).hexdigest()
    table_name = "{}/{}".format(table_prefix, filemask_hash)
    query_mask = (
        "SELECT count(name) FROM sqlite_master WHERE type='table' AND name='{}'"
    )
    cursor = db.execute(query_mask.format(table_name))
    if cursor.fetchone()[0] > 0: # table exists
        return table_name
    # otherwise, table doesn't exist and we need to download data:
    try:
        stream = get(url, stream=True)
    except InvalidSchema:
        if http_fallback:
            stream = get(sub(r'^ftp:\/\/', "http://", url), stream=True)
        else:
            raise
    if stream.status_code != 200:
        raise URLError("{}: status code {}".format(url, stream.status_code))
    total_bytes = int(stream.headers.get("content-length", 0))
    with TemporaryDirectory() as tempdir:
        sep, compression = guess_format(filemask)
        target_file = join(tempdir, filemask_hash)
        if compression:
            target_file = target_file + "." + compression
        with open(target_file, "wb") as output_handle:
            written_bytes = 0
            for block in stream.iter_content(1024):
                output_handle.write(block)
                written_bytes += len(block)
        if total_bytes != written_bytes:
            remove(target_file)
            raise URLError("Failed to download the correct number of bytes")
        read_csv(target_file, sep=sep).to_sql(table_name, db)
        db.commit()
        db.close()
    return table_name


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
    if data_rargs["any_below"] is not None: # FIXME MOVEME
        repr_df = get_padj_filtered_repr_df(repr_df, data_rargs["any_below"])
    if melting is not False: # FIXME MOVEME
        repr_df = melt_table_data(repr_df, melting=melting)
    else: # FIXME MOVEME
        repr_df = repr_df.reset_index()
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
    return None # FIXME
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
    return None # FIXME
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
