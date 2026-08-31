"""
tests/test_ora.py
"""
from enrichment.ora import fisher_exact_for_term


def test_fisher_exact_detects_strong_overrepresentation():
    # population: 1000 genes, 100 have the term (10%)
    # class: 50 genes, 40 have the term (80%) -- obviously over-represented
    odds_ratio, p_value = fisher_exact_for_term(
        n_pop_with_term=100, n_pop_total=1000,
        n_class_with_term=40, n_class_total=50,
    )
    assert p_value < 0.001
    assert odds_ratio > 1


def test_fisher_exact_no_enrichment_when_proportions_match():
    # class proportion (10%) matches population proportion (10%) exactly
    odds_ratio, p_value = fisher_exact_for_term(
        n_pop_with_term=100, n_pop_total=1000,
        n_class_with_term=5, n_class_total=50,
    )
    assert p_value > 0.5


def test_fisher_exact_detects_underrepresentation():
    # population: 50% have the term, class: 0% have it -- clearly under-represented
    odds_ratio, p_value = fisher_exact_for_term(
        n_pop_with_term=500, n_pop_total=1000,
        n_class_with_term=0, n_class_total=50,
    )
    assert p_value < 0.001
    assert odds_ratio < 1