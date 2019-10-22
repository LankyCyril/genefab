from genefab import GeneLabJSONException, GeneLabDataManagerException
from os import remove, path
from requests import get
from requests.exceptions import InvalidSchema
from urllib.error import URLError
from sqlite3 import connect
from hashlib import sha512
from pandas import read_csv, read_sql_query, DataFrame, Index, merge
from pandas.io.sql import DatabaseError as PandasDatabaseError
from tempfile import TemporaryDirectory
from genefab._util import STORAGE_PREFIX, DELIM_AS_IS
from genefab._util import guess_format, data_rargs_digest
from re import sub, search, IGNORECASE


def download_table(accession, assay_name, filemask, url, verbose=False, http_fallback=True):
    """Download and interpret table file"""
    try:
        stream = get(url, stream=True)
    except InvalidSchema:
        if http_fallback:
            stream = get(sub(r'^ftp:\/\/', "http://", url), stream=True)
        else:
            raise
    if stream.status_code != 200:
        raise GeneLabDataManagerException("HTTP Error {}: {}".format(
            stream.status_code, url
        ))
    total_bytes = int(stream.headers.get("content-length", 0))
    with TemporaryDirectory() as tempdir:
        sep, compression = guess_format(filemask)
        filemask_hash = sha512(filemask.encode("utf-8")).hexdigest()
        target_file = path.join(tempdir, filemask_hash)
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
        return read_csv(target_file, sep=sep)


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


def format_table_data(repr_df, assay, data_rargs):
    """Format file data accoring to rargdict and melting"""
    if data_rargs["name_delim"] != DELIM_AS_IS:
        conv_delim = lambda f: sub(r'[._-]', data_rargs["name_delim"], f)
        repr_df.columns = repr_df.columns.map(conv_delim)
    if data_rargs["any_below"] is not None: # FIXME MOVEME
        repr_df = get_padj_filtered_repr_df(repr_df, data_rargs["any_below"])
    if data_rargs["descriptive"]:
        repr_df = melt_table_data(repr_df, melting=assay.annotation().T)
    elif data_rargs["melted"]:
        repr_df = melt_table_data(repr_df, melting=list(assay.annotation().T.columns))
    else:
        repr_df = repr_df.reset_index()
    return repr_df


def retrieve_table_data(assay, filemask, data_rargs):
    """Find file URL that matches filemask, redirect to download or interpret"""
    try:
        url = assay._get_file_url(filemask)
    except GeneLabJSONException:
        raise ValueError("multiple files match mask")
    if url is None:
        raise FileNotFoundError
    repr_df = download_table(
        assay.parent.accession, assay.name, filemask, url
    )
    return format_table_data(repr_df, assay, data_rargs)


def try_sqlite(accession, assay_name, data_rargs):
    """Try to load dataframe from DB_NAME"""
    table_name = data_rargs_digest(data_rargs)
    db = connect(path.join(
        STORAGE_PREFIX, accession + "-" + assay_name + ".sqlite3"
    ))
    query = "SELECT * FROM '{}'".format(table_name)
    try:
        return read_sql_query(query, db, index_col="index")
    except PandasDatabaseError:
        return None


def dump_to_sqlite(accession, assay_name, data_rargs, table_data):
    """Save transformed dataframe to DB_NAME"""
    table_name = data_rargs_digest(data_rargs)
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
