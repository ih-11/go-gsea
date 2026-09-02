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


def drop_uninformative_for_clustering(df, cols, thresh=0):
    """
    Drops rows where EVERY one of cols is <= thresh (default 0).

    WHEN TO USE THIS: before labelers.cluster() specifically, and ONLY
    before cluster(). Do not use this before rank_tail, explicit_threshold,
    or boolean_flag -- those strategies already handle zero-value genes
    correctly (a gene at the floor simply sorts into the bottom tail, or
    fails/passes a threshold, exactly as it should), dropping them there
    would only needlessly shrink the population and discard real,
    meaningful Low-value signal for no benefit.

    WHY THIS EXISTS: it is a tool-mismatch fix, not a claim that
    zero-expression genes are biologically uninformative. A gene with
    zero expression is telling you something real (the gene is off) --
    that is a legitimate result, not noise. The actual problem is that
    Ward's-method (or any continuous-distance) clustering mechanically
    merges every gene tied at exactly the same floor value into one
    trivially "perfect" cluster (mean/min silhouette suspiciously at or
    near 1.0), regardless of whether any real relationship exists among
    them, simply because they sit at identical distance from each other
    after the Yeo-Johnson + Z-score transform. This is a known,
    real-data-confirmed distance-clustering artifact (see the k=13
    rejection in docs/examples/chlamydomonas.md demo C, and the
    artifact checks in demos F/G), first and most clearly seen in a
    SINGLE-column (1D) case, not only multi-column.

    IMPORTANT: dropping these rows removes them from the clustering
    pass's population entirely, changing every GO term's expected count
    and p-value for that run. Do not treat a --drop-uninformative run as
    the complete picture -- the recommended practice is a companion run
    on the SAME population using boolean_flag (or explicit_threshold at
    the same thresh) to test the excluded, zero-value genes' own GO
    enrichment on their own terms, so nothing is silently lost from the
    overall analysis, it is just analyzed with a tool that actually fits
    its distribution instead of one that produces an artifact on it.

    Rows failing this filter are dropped from the returned dataframe
    entirely, not silently folded into the background population under
    a different label.

    cols: list of column names. A row is dropped only if ALL of them are
    <= thresh, not if just one is. Works identically for a single column
    (the original, clearest real case) or several.
    """
    if not cols:
        raise ValueError("cols must contain at least one column name")
    mask = (df[cols] <= thresh).all(axis=1)
    return df[~mask].copy()


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