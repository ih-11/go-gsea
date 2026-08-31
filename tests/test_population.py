"""
tests/test_population.py

Verifies filters/population.py against small, hand-checkable synthetic data --
not real Chlamydomonas data, since this module must stay species-agnostic.
"""
import pandas as pd
import pytest

from filters.population import (
    read_depth_filter,
    usage_filter,
    min_group_size_filter,
    chain_filters,
)


def test_read_depth_filter_keeps_only_rows_at_or_above_thresh():
    df = pd.DataFrame({"gene_id": ["a", "b", "c"], "n_reads": [50, 100, 150]})
    result = read_depth_filter(df, "n_reads", thresh=100)
    assert set(result["gene_id"]) == {"b", "c"}


def test_read_depth_filter_boundary_is_inclusive():
    # thresh itself must pass, not just strictly-greater-than
    df = pd.DataFrame({"gene_id": ["a"], "n_reads": [100]})
    result = read_depth_filter(df, "n_reads", thresh=100)
    assert len(result) == 1


def test_usage_filter_keeps_only_rows_at_or_above_thresh():
    df = pd.DataFrame({"site_id": ["x", "y", "z"], "usage": [0.01, 0.05, 0.20]})
    result = usage_filter(df, "usage", thresh=0.05)
    assert set(result["site_id"]) == {"y", "z"}


def test_min_group_size_filter_drops_undersized_groups():
    # transcript t1 has 3 sites, t2 has 1 site -- with min_n=2, only t1 survives
    df = pd.DataFrame({
        "transcript_id": ["t1", "t1", "t1", "t2"],
        "site_id": ["s1", "s2", "s3", "s4"],
    })
    result = min_group_size_filter(df, "transcript_id", min_n=2)
    assert set(result["transcript_id"]) == {"t1"}
    assert len(result) == 3  # all 3 of t1's rows kept, not just 2


def test_min_group_size_filter_boundary_is_inclusive():
    df = pd.DataFrame({"transcript_id": ["t1", "t1"], "site_id": ["s1", "s2"]})
    result = min_group_size_filter(df, "transcript_id", min_n=2)
    assert len(result) == 2


def test_chain_filters_applies_in_order_and_composes_correctly():
    # mimics the real precedent: read-depth filter, then usage filter,
    # then require >=2 surviving sites per transcript
    df = pd.DataFrame({
        "transcript_id": ["t1", "t1", "t1", "t2", "t2"],
        "site_id":       ["s1", "s2", "s3", "s4", "s5"],
        "n_reads":       [100, 100, 50, 100, 100],   # s3 fails read depth
        "usage":         [0.05, 0.01, 0.05, 0.05, 0.05],  # s2 fails usage
    })
    # t1 loses s3 (reads) and s2 (usage) -> only s1 survives -> t1 has 1 site, dropped by min_group_size
    # t2 keeps both s4 and s5 -> 2 sites -> survives
    result = chain_filters(
        df,
        (read_depth_filter, {"col": "n_reads", "thresh": 100}),
        (usage_filter, {"col": "usage", "thresh": 0.05}),
        (min_group_size_filter, {"group_col": "transcript_id", "min_n": 2}),
    )
    assert set(result["transcript_id"].unique()) == {"t2"}
    assert len(result) == 2


def test_filters_never_mutate_the_input_dataframe():
    # .copy() must actually be happening, or a caller's original df gets corrupted
    df = pd.DataFrame({"gene_id": ["a", "b"], "n_reads": [50, 150]})
    original = df.copy()
    _ = read_depth_filter(df, "n_reads", thresh=100)
    pd.testing.assert_frame_equal(df, original)