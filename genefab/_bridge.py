from genefab import GLDS, GeneLabJSONException
from re import sub, split, search
from operator import __lt__, __le__, __eq__, __ne__, __ge__, __gt__


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


def get_filtered_repr_df(repr_df, field_filters_raw):
    """Interpret the filter request argument and subset the repr dataframe"""
    if not isinstance(field_filters_raw, (list, tuple, set)):
        field_filters = [field_filters_raw]
    else:
        field_filters = field_filters_raw
    indexer = None
    for field_filter in field_filters:
        field_filter_stripped = sub(r'(^\')|(\'$)', "", field_filter)
        match = search(r'(^[^<>=]+)([<>=]+)(.+)$', field_filter_stripped)
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
            emsk = "Invalid comparison (TypeError): '{}' is `{}` and {} is `{}`"
            raise TypeError(emsk.format(
                field, repr_df[field].dtype, value, type(value).__name__
            ))
        if indexer is None:
            indexer = field_indexer
        else:
            indexer &= field_indexer
    return repr_df[indexer]


def filter_table_data(repr_df, data_filter_rargs):
    """Filter dataframe"""
    if data_filter_rargs["filter"] is not None:
        repr_df = get_filtered_repr_df(repr_df, data_filter_rargs["filter"])
    if data_filter_rargs["sort_by"] is not None:
        sort_by = sub(r'(^\')|(\'$)', "", data_filter_rargs["sort_by"])
        if sort_by in repr_df.columns:
            repr_df = repr_df.sort_values(
                by=sort_by, ascending=data_filter_rargs["ascending"]
            )
        else:
            error_mask = "Unknown field (column) '{}'"
            raise IndexError(error_mask.format(sort_by))
    return repr_df
