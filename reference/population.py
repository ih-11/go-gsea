"""
filters/population.py

Stage A of the two-stage labeling architecture: determines which genes/rows
are ELIGIBLE for a given analysis before any class label (High/Low/etc.) is
assigned. This is the "population" that enrichment gets compared against --
getting this right is what the Yamasaki lab's methodology notes flagged as
the single most important correctness issue in GO enrichment.

Generic: operates on any dataframe + column name, knows nothing about PR,
TPM, or Chlamydomonas specifically.
"""


def read_depth_filter(df, col, thresh):
    """Keep rows where col >= thresh. E.g. N_READS_THRESH=100 in the precedent."""
    return df[df[col] >= thresh].copy()


def usage_filter(df, col, thresh):
    """Keep rows where a usage fraction col >= thresh. E.g. USAGE_THRESH=0.05."""
    return df[df[col] >= thresh].copy()


def min_group_size_filter(df, group_col, min_n):
    """
    Keep only groups (e.g. transcripts) that have >= min_n surviving rows
    (e.g. sites). E.g. N_SITES_THRESH=2 in the precedent -- a transcript
    needs at least 2 eligible sites to be usable at all.
    """
    counts = df.groupby(group_col).size()
    valid_groups = counts[counts >= min_n].index
    return df[df[group_col].isin(valid_groups)].copy()


def chain_filters(df, *filter_calls):
    """
    Applies a sequence of (filter_fn, kwargs) pairs in order, so a caller
    can compose e.g. read_depth_filter -> usage_filter -> min_group_size_filter
    exactly like Yamasaki's fltr_data() does, without hardcoding that
    sequence into this module.
    """
    result = df
    for fn, kwargs in filter_calls:
        result = fn(result, **kwargs)
    return result