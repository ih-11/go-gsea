"""
tests/test_ora.py
"""
from enrichment.ora import fisher_exact_for_term, fold_enrichment


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