from genefab import GeneLabJSONException, GeneLabDataManagerException
from os import remove, path
from requests import get
from requests.exceptions import InvalidSchema
from urllib.error import URLError
from contextlib import closing
from sqlite3 import connect, OperationalError, Binary
from hashlib import sha512
from pandas import read_csv, read_json, DataFrame, Index, merge
from pandas.io.sql import DatabaseError as PandasDatabaseError
from tempfile import TemporaryDirectory
from genefab._util import STORAGE_PREFIX, DELIM_AS_IS
from genefab._util import guess_format, data_rargs_digest
from genefab._display import fix_cols
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
        filemask_hash = sha512(filemask.encode("utf-8")).hexdigest()
        target_file = path.join(tempdir, filemask_hash)
        with open(target_file, "wb") as output_handle:
            written_bytes = 0
            for block in stream.iter_content(1024):
                output_handle.write(block)
                written_bytes += len(block)
        if total_bytes != written_bytes:
            remove(target_file)
            raise URLError("Failed to download the correct number of bytes")
        sep, compression = guess_format(target_file)
        return read_csv(target_file, sep=sep, compression=compression)


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


def melt_table_data(repr_df, melting, cols_to_fix={"Unnamed: 0": "Sample Name"}):
    """Melt dataframe with the use of sample annotation"""
    if cols_to_fix:
        repr_df = fix_cols(repr_df, cols_to_fix)
    foundry = repr_df.reset_index().copy()
    if isinstance(melting, (list, Index)):
        id_vars = [c for c in foundry.columns if c not in melting]
        melted_data = foundry.melt(
            id_vars=id_vars, value_vars=melting, var_name="Sample Name"
        )
    elif isinstance(melting, DataFrame):
        try:
            id_vars = [c for c in foundry.columns if c not in melting.columns]
            foundry = foundry.melt(
                id_vars=id_vars, value_vars=melting, var_name="Sample Name"
            )
            melted_data = merge(melting.T.reset_index(), foundry, how="outer")
        except KeyError:
            melted_data = merge(
                melting.T.reset_index(), foundry, how="outer", on="Sample Name"
            )
    else:
        raise TypeError("cannot melt/describe with a non-dataframe object")
    if "index" in melted_data.columns:
        return melted_data.drop(columns=["index"])
    else:
        return melted_data


def format_table_data(repr_df, assay, data_rargs, cols_to_fix={"Unnamed: 0"}):
    """Format file data accoring to rargdict and melting"""
    if data_rargs["name_delim"] != DELIM_AS_IS:
        conv_delim = lambda f: sub(r'[._-]', data_rargs["name_delim"], f)
        repr_df.columns = repr_df.columns.map(conv_delim)
        for col_to_fix in cols_to_fix:
            if col_to_fix in repr_df.columns:
                repr_df[col_to_fix] = repr_df[col_to_fix].apply(conv_delim)
    if data_rargs["any_below"] is not None:
        repr_df = get_padj_filtered_repr_df(repr_df, data_rargs["any_below"])
    if data_rargs["descriptive"]:
        repr_df = melt_table_data(repr_df, melting=assay.annotation().T)
    elif data_rargs["melted"]:
        repr_df = melt_table_data(
            repr_df,
            melting=list(assay.annotation().T.columns)
        )
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


def read_blob_as_table(table_name, db):
    query_mask = "SELECT blob FROM 'table_blobs' WHERE name = '{}' LIMIT 1"
    row = db.cursor().execute(query_mask.format(table_name)).fetchone()
    if isinstance(row, tuple) and (len(row) == 1):
        return read_json(row[0].decode(), orient="split")
    else:
        return PandasDatabaseError


def drop_blob(table_name, db, drop_date=True):
    try:
        db.cursor().execute(
            "DELETE FROM 'table_blobs' WHERE name = '{}'".format(table_name)
        )
    except OperationalError:
        pass
    if drop_date:
        try:
            db.cursor().execute(
                "DELETE FROM 'table_dates' WHERE name = '{}'".format(table_name)
            )
        except OperationalError:
            pass


def try_sqlite(accession, assay_name, data_rargs, expect_date):
    """Try to load dataframe from DB_NAME"""
    table_name = data_rargs_digest(accession, assay_name, data_rargs)
    db_name = path.join(STORAGE_PREFIX, "tables.sqlite3")
    with closing(connect(db_name)) as db:
        date_query_mask = "SELECT date FROM 'table_dates' WHERE name = '{}'"
        date_query = date_query_mask.format(table_name)
        try:
            stored_dates = db.cursor().execute(date_query).fetchall()
        except OperationalError:
            stored_dates = []
        is_stored_date_expected = (
            isinstance(stored_dates, list) and (len(stored_dates) == 1) and
            isinstance(stored_dates[0], tuple) and (len(stored_dates[0]) == 1)
            and (stored_dates[0][0] == expect_date)
        )
        if is_stored_date_expected:
            try:
                return read_blob_as_table(table_name, db)
            except (PandasDatabaseError, OperationalError):
                pass
        # otherwise, the table is too old and needs to be destroyed:
        drop_blob(table_name, db)
        return None


def write_table_as_blob(table_data, table_name, db):
    db.cursor().execute(
        "CREATE TABLE IF NOT EXISTS 'table_blobs' ('name' TEXT, 'blob' BLOB)"
    )
    table_data.to_json("tmp.json", orient="split")
    db.cursor().execute(
        "INSERT INTO 'table_blobs' (name, blob) VALUES (?, ?)",
        [table_name, Binary(table_data.to_json(orient="split").encode("utf-8"))]
    )


def dump_to_sqlite(accession, assay_name, data_rargs, table_data, set_date):
    """Save transformed dataframe to DB_NAME"""
    table_name = data_rargs_digest(accession, assay_name, data_rargs)
    db_name = path.join(STORAGE_PREFIX, "tables.sqlite3")
    with closing(connect(db_name)) as db:
        drop_blob(table_name, db, drop_date=False)
        write_table_as_blob(table_data, table_name, db)
        test_query_mask = "SELECT date FROM '{}' WHERE name = '{}'"
        test_query = test_query_mask.format("table_dates", table_name)
        try:
            date_test = db.cursor().execute(test_query).fetchall()
        except OperationalError:
            new_date_table = DataFrame(
                data=[[table_name, set_date]],
                index=[0], columns=["name", "date"]
            )
            new_date_table.to_sql("table_dates", db, index=False)
        else:
            if len(date_test):
                qmask = "UPDATE '{}' SET date = {} WHERE name = '{}'"
                query = qmask.format("table_dates", set_date, table_name)
            else:
                qmask = "INSERT INTO '{}' (name, date) VALUES ('{}', {})"
                query = qmask.format("table_dates", table_name, set_date)
            db.cursor().execute(query)
        db.commit()
