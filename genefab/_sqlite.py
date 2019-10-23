from genefab import GeneLabJSONException, GeneLabDataManagerException
from os import remove, path
from requests import get
from requests.exceptions import InvalidSchema
from urllib.error import URLError
from contextlib import closing
from sqlite3 import connect, OperationalError
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
        melted_data = foundry.melt(
            id_vars=id_vars, value_vars=melting, var_name="Sample Name"
        )
    elif isinstance(melting, DataFrame):
        id_vars = [c for c in foundry.columns if c not in melting.columns]
        foundry = foundry.melt(
            id_vars=id_vars, value_vars=melting, var_name="Sample Name"
        )
        melted_data = merge(melting.T.reset_index(), foundry, how="outer")
    else:
        raise TypeError("cannot melt/describe with a non-dataframe object")
    if "index" in melted_data.columns:
        return melted_data.drop(columns=["index"])
    else:
        return melted_data


def format_table_data(repr_df, assay, data_rargs):
    """Format file data accoring to rargdict and melting"""
    if data_rargs["name_delim"] != DELIM_AS_IS:
        conv_delim = lambda f: sub(r'[._-]', data_rargs["name_delim"], f)
        repr_df.columns = repr_df.columns.map(conv_delim)
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


def try_sqlite(accession, assay_name, data_rargs, expect_date):
    """Try to load dataframe from DB_NAME"""
    table_name = data_rargs_digest(data_rargs)
    db_name = path.join(
        STORAGE_PREFIX, accession + "-" + assay_name + ".sqlite3"
    )
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
            query = "SELECT * FROM '{}'".format(table_name)
            try:
                return read_sql_query(query, db, index_col="index")
            except PandasDatabaseError:
                pass
        query = "DELETE FROM 'table_dates' WHERE name = '{}'".format(table_name)
        try:
            db.cursor().execute(query)
            db.commit()
        except OperationalError:
            pass
        db.cursor().execute("DROP TABLE IF EXISTS '{}'".format(table_name))
        db.commit()
        return None


def dump_to_sqlite(accession, assay_name, data_rargs, table_data, set_date):
    """Save transformed dataframe to DB_NAME"""
    table_name = data_rargs_digest(data_rargs)
    db_name = path.join(
        STORAGE_PREFIX, accession + "-" + assay_name + ".sqlite3"
    )
    with closing(connect(db_name)) as db:
        db.cursor().execute("DROP TABLE IF EXISTS '{}'".format(table_name))
        db.commit()
        table_data.to_sql(table_name, db)
        db.commit()
        test_query_mask = "SELECT date FROM 'table_dates' WHERE name = '{}'"
        test_query = test_query_mask.format(table_name)
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
                qmask = "UPDATE 'table_dates' SET date = {} WHERE name = '{}'"
                query = qmask.format(set_date, table_name)
            else:
                qmask = "INSERT INTO 'table_dates' (name, date) VALUES ('{}', {})"
                query = qmask.format(table_name, set_date)
            db.cursor().execute(query)
        db.commit()
