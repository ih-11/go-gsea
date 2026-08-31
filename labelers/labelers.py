"""
labelers/labelers.py

Stage B of the two-stage labeling architecture: given an ELIGIBLE population
(already filtered by filters/population.py), assigns a class label to a
subset of rows for GO enrichment testing.

Every labeler returns a dataframe with an added 'class' column -- rows
NOT selected for any class get class=None, meaning they remain part of
the background population but aren't in any tested group. This matches
the real -c input format confirmed from stats_test_based_on_go.py.

Generic: operates on any dataframe + column name(s), knows nothing about
PR, TPM, or Chlamydomonas specifically.
"""
import pandas as pd
from scipy import stats as scipy_stats
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import pdist


def rank_tail(df, col, pct=10, id_col=None):
    """
    Labels the top pct% and bottom pct% of rows (by col) as 'High'/'Low'.
    Matches the confirmed real-data convention (exactly 10%/10% in the
    AT-T87-HS-PR.tsv etc. precedent).
    """
    df = df.copy()
    df = df.sort_values(col, ascending=False)
    n = len(df)
    n_tail = int(n * pct / 100)

    df["class"] = None
    df.iloc[:n_tail, df.columns.get_loc("class")] = "High"
    df.iloc[-n_tail:, df.columns.get_loc("class")] = "Low"
    return df


def explicit_threshold(df, col, high_thresh=None, low_thresh=None):
    """
    Labels rows above high_thresh as 'High' and below low_thresh as 'Low'.
    Generic form of the real N_READS_THRESH/USAGE_THRESH-style fixed cutoffs.
    Either bound can be omitted to only label one side.
    """
    df = df.copy()
    df["class"] = None
    if high_thresh is not None:
        df.loc[df[col] >= high_thresh, "class"] = "High"
    if low_thresh is not None:
        df.loc[df[col] <= low_thresh, "class"] = "Low"
    return df


def boolean_flag(df, bool_col, true_label="Present", false_label=None):
    """
    Labels rows by an existing boolean/categorical column -- e.g. G4
    presence, motif presence, model-predicted-vs-not. false_label=None
    (the default) leaves non-matching rows as background-only, matching
    a typical presence/absence enrichment question; pass a string to
    make it a real two-class comparison instead.
    """
    df = df.copy()
    df["class"] = None
    df.loc[df[bool_col] == True, "class"] = true_label
    if false_label is not None:
        df.loc[df[bool_col] == False, "class"] = false_label
    return df


def cluster(df, cols, method="ward", transform="yeo-johnson", n_clusters=None,
            t=None):
    """
    Ward's-method hierarchical clustering, matching the Yamasaki lab's
    confirmed translatome analysis (Yeo-Johnson -> Z-score -> Euclidean ->
    Ward's method). class = cluster ID as a string ("1", "2", ...).

    n_clusters: cut the tree to exactly this many clusters.
    t: alternative to n_clusters -- cut at this distance threshold instead.
    Exactly one of n_clusters/t must be given.
    """
    if (n_clusters is None) == (t is None):
        raise ValueError("Provide exactly one of n_clusters or t")

    df = df.copy()
    data = df[cols].copy()

    if transform == "yeo-johnson":
        # yeojohnson() returns (transformed_array, lambda) -- keep only [0]
        data = data.apply(lambda c: scipy_stats.yeojohnson(c)[0])

    z = (data - data.mean()) / data.std()

    dist = pdist(z.values, metric="euclidean")
    link = linkage(dist, method=method)

    if n_clusters is not None:
        labels = fcluster(link, t=n_clusters, criterion="maxclust")
    else:
        labels = fcluster(link, t=t, criterion="distance")

    df["class"] = [str(l) for l in labels]
    return df