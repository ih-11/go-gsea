"""
tests/test_ora.py
"""
import pandas as pd
from enrichment.ora import fisher_exact_for_term, fold_enrichment, bh_correct
from enrichment.ora import count_term_occurrences, run_ora

def test_fisher_exact_detects_strong_overrepresentation():
    p_value = fisher_exact_for_term(
        n_pop_with_term=100, n_pop_total=1000,
        n_class_with_term=40, n_class_total=50,
    )
    assert p_value < 0.001


def test_fisher_exact_no_enrichment_when_proportions_match():
    p_value = fisher_exact_for_term(
        n_pop_with_term=100, n_pop_total=1000,
        n_class_with_term=5, n_class_total=50,
    )
    assert p_value > 0.5


def test_fisher_exact_detects_underrepresentation():
    p_value = fisher_exact_for_term(
        n_pop_with_term=500, n_pop_total=1000,
        n_class_with_term=0, n_class_total=50,
    )
    assert p_value < 0.001


def test_fold_enrichment_positive_when_overrepresented():
    fe = fold_enrichment(
        n_pop_with_term=100, n_pop_total=1000,
        n_class_with_term=40, n_class_total=50,
    )
    assert fe > 0


def test_fold_enrichment_negative_when_underrepresented():
    fe = fold_enrichment(
        n_pop_with_term=500, n_pop_total=1000,
        n_class_with_term=0, n_class_total=50,
    )
    assert fe < 0


def test_fold_enrichment_near_zero_when_matched():
    fe = fold_enrichment(
        n_pop_with_term=100, n_pop_total=1000,
        n_class_with_term=5, n_class_total=50,
    )
    assert abs(fe) < 0.2

def test_bh_correct_returns_same_length_as_input():
    q_values = bh_correct([0.001, 0.01, 0.05, 0.5, 0.9])
    assert len(q_values) == 5


def test_bh_correct_q_values_are_never_smaller_than_p_values():
    # BH correction can only inflate (or keep equal), never shrink a p-value
    p_values = [0.001, 0.01, 0.05, 0.5, 0.9]
    q_values = bh_correct(p_values)
    assert all(q >= p for p, q in zip(p_values, q_values))


def test_bh_correct_preserves_relative_order():
    # smaller p-values should still map to smaller-or-equal q-values
    p_values = [0.5, 0.001, 0.9, 0.01]
    q_values = bh_correct(p_values)
    order_by_p = sorted(range(4), key=lambda i: p_values[i])
    q_at_order = [q_values[i] for i in order_by_p]
    assert q_at_order == sorted(q_at_order)


def test_bh_correct_empty_input():
    assert bh_correct([]) == []

def test_count_term_occurrences_basic():
    gene_go = {
        "g1": {"GO:0001", "GO:0002"},
        "g2": {"GO:0001"},
        "g3": {"GO:0003"},
    }
    counts = count_term_occurrences(gene_go, ["g1", "g2", "g3"])
    assert counts == {"GO:0001": 2, "GO:0002": 1, "GO:0003": 1}


def test_count_term_occurrences_only_counts_requested_genes():
    gene_go = {
        "g1": {"GO:0001"},
        "g2": {"GO:0001"},
    }
    # only ask about g1 -- g2 should not contribute even though it's in gene_go
    counts = count_term_occurrences(gene_go, ["g1"])
    assert counts == {"GO:0001": 1}


def test_count_term_occurrences_handles_gene_not_in_godb():
    gene_go = {"g1": {"GO:0001"}}
    # g2 isn't in gene_go at all -- should be silently skipped, not error
    counts = count_term_occurrences(gene_go, ["g1", "g2"])
    assert counts == {"GO:0001": 1}


def test_count_term_occurrences_empty_gene_list():
    gene_go = {"g1": {"GO:0001"}}
    assert count_term_occurrences(gene_go, []) == {}

def test_run_ora_end_to_end_toy_example():
    # 4 genes, term GO:0001 is strongly enriched in class "High"
    gene_go = {
        "g1": {"GO:0001"},
        "g2": {"GO:0001"},
        "g3": {"GO:0002"},
        "g4": {"GO:0002"},
    }
    labeled_df = pd.DataFrame({
        "gene_id": ["g1", "g2", "g3", "g4"],
        "class": ["High", "High", None, None],
    })
    result = run_ora(gene_go, labeled_df)

    assert set(result["class"].unique()) == {"High"}
    row_0001 = result[result["go_id"] == "GO:0001"].iloc[0]
    assert row_0001["observed"] == 2
    assert row_0001["population"] == 2
    assert row_0001["fold_enrichment"] > 0


def test_run_ora_handles_multiple_classes_independently():
    gene_go = {
        "g1": {"GO:0001"}, "g2": {"GO:0001"},
        "g3": {"GO:0002"}, "g4": {"GO:0002"},
    }
    labeled_df = pd.DataFrame({
        "gene_id": ["g1", "g2", "g3", "g4"],
        "class": ["High", "High", "Low", "Low"],
    })
    result = run_ora(gene_go, labeled_df)
    assert set(result["class"].unique()) == {"High", "Low"}
    # each class should have its own row for both terms tested independently
    assert len(result) == 4  # 2 classes x 2 terms


def test_run_ora_significance_uses_p_threshold_by_default():
    gene_go = {"g1": {"GO:0001"}, "g2": {"GO:0001"}, "g3": {"GO:0002"}}
    labeled_df = pd.DataFrame({
        "gene_id": ["g1", "g2", "g3"],
        "class": ["High", "High", None],
    })
    result = run_ora(gene_go, labeled_df, thresh_type="p", thresh=0.01)
    assert "significance" in result.columns
    assert (result["significance"] == (result["p_value"] < 0.01)).all()


def test_run_ora_empty_class_returns_empty_dataframe():
    gene_go = {"g1": {"GO:0001"}}
    labeled_df = pd.DataFrame({"gene_id": ["g1"], "class": [None]})
    result = run_ora(gene_go, labeled_df)
    assert len(result) == 0