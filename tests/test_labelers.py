"""
tests/test_labelers.py
"""
import pandas as pd
import numpy as np
import pytest

from labelers.labelers import (
    rank_tail, explicit_threshold, boolean_flag, cluster,
    elbow_curve, silhouette_scores, sample_silhouette,
)

def test_rank_tail_labels_exactly_top_and_bottom_pct():
    df = pd.DataFrame({"gene_id": [f"g{i}" for i in range(100)], "value": range(100)})
    result = rank_tail(df, "value", pct=10)
    assert (result["class"] == "High").sum() == 10
    assert (result["class"] == "Low").sum() == 10
    assert result[result["class"] == "High"]["value"].min() == 90
    assert result[result["class"] == "Low"]["value"].max() == 9


def test_rank_tail_middle_rows_are_unlabeled():
    df = pd.DataFrame({"gene_id": [f"g{i}" for i in range(100)], "value": range(100)})
    result = rank_tail(df, "value", pct=10)
    assert result["class"].isna().sum() == 80


def test_explicit_threshold_labels_both_sides():
    df = pd.DataFrame({"gene_id": ["a", "b", "c", "d"], "value": [1, 5, 10, 15]})
    result = explicit_threshold(df, "value", high_thresh=10, low_thresh=5)
    assert result.set_index("gene_id")["class"].to_dict() == {
        "a": "Low", "b": "Low", "c": "High", "d": "High"
    }


def test_explicit_threshold_one_sided_only():
    df = pd.DataFrame({"gene_id": ["a", "b"], "value": [1, 100]})
    result = explicit_threshold(df, "value", high_thresh=50)
    assert result.set_index("gene_id")["class"].to_dict() == {"a": None, "b": "High"}


def test_boolean_flag_default_leaves_false_unlabeled():
    df = pd.DataFrame({"gene_id": ["a", "b"], "has_motif": [True, False]})
    result = boolean_flag(df, "has_motif")
    assert result.set_index("gene_id")["class"].to_dict() == {"a": "Present", "b": None}


def test_boolean_flag_two_class_mode():
    df = pd.DataFrame({"gene_id": ["a", "b"], "has_motif": [True, False]})
    result = boolean_flag(df, "has_motif", true_label="WithG4", false_label="NoG4")
    assert result.set_index("gene_id")["class"].to_dict() == {"a": "WithG4", "b": "NoG4"}


def test_cluster_requires_exactly_one_of_n_clusters_or_t():
    df = pd.DataFrame({"gene_id": ["a", "b"], "val": [1.0, 2.0]})
    with pytest.raises(ValueError):
        cluster(df, ["val"], n_clusters=2, t=1.0)
    with pytest.raises(ValueError):
        cluster(df, ["val"])


def test_cluster_separates_two_obviously_distinct_groups():
    np.random.seed(0)
    clump_a = np.random.normal(0, 0.1, 20)
    clump_b = np.random.normal(50, 0.1, 20)
    df = pd.DataFrame({
        "gene_id": [f"g{i}" for i in range(40)],
        "val": np.concatenate([clump_a, clump_b]),
    })
    result = cluster(df, ["val"], n_clusters=2)
    labels_a = set(result.iloc[:20]["class"])
    labels_b = set(result.iloc[20:]["class"])
    assert len(labels_a) == 1
    assert len(labels_b) == 1
    assert labels_a != labels_b

def test_elbow_curve_wss_decreases_from_k1_to_k2():
    np.random.seed(1)
    clump_a = np.random.normal(0, 0.1, 20)
    clump_b = np.random.normal(50, 0.1, 20)
    df = pd.DataFrame({
        "gene_id": [f"g{i}" for i in range(40)],
        "val": np.concatenate([clump_a, clump_b]),
    })
    k_list, wss_list = elbow_curve(df, ["val"], max_k=5)
    assert k_list == [1, 2, 3, 4, 5]
    assert wss_list[0] > wss_list[1]


def test_elbow_curve_wss_is_monotonically_non_increasing():
    np.random.seed(2)
    df = pd.DataFrame({
        "gene_id": [f"g{i}" for i in range(30)],
        "val": np.random.normal(0, 1, 30),
    })
    k_list, wss_list = elbow_curve(df, ["val"], max_k=10)
    assert all(wss_list[i] >= wss_list[i + 1] - 1e-9 for i in range(len(wss_list) - 1))


def test_silhouette_scores_prefers_the_true_number_of_clumps():
    np.random.seed(3)
    clump_a = np.random.normal(0, 0.1, 20)
    clump_b = np.random.normal(50, 0.1, 20)
    df = pd.DataFrame({
        "gene_id": [f"g{i}" for i in range(40)],
        "val": np.concatenate([clump_a, clump_b]),
    })
    k_list, score_list = silhouette_scores(df, ["val"], k_range=range(2, 6))
    best_k = k_list[int(np.argmax(score_list))]
    assert best_k == 2


def test_sample_silhouette_returns_one_row_per_input_row():
    np.random.seed(4)
    clump_a = np.random.normal(0, 0.1, 10)
    clump_b = np.random.normal(50, 0.1, 10)
    df = pd.DataFrame({
        "gene_id": [f"g{i}" for i in range(20)],
        "val": np.concatenate([clump_a, clump_b]),
    })
    result = sample_silhouette(df, ["val"], n_clusters=2)
    assert len(result) == 20
    assert "cluster" in result.columns
    assert "silhouette" in result.columns
    assert result["silhouette"].mean() > 0.8