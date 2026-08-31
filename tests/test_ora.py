"""
tests/test_ora.py
"""
from enrichment.ora import fisher_exact_for_term, fold_enrichment, bh_correct


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