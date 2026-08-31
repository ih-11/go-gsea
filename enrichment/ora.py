"""
enrichment/ora.py

Over-representation analysis (ORA): Fisher's exact test per GO term,
comparing a labeled class against the full population.

Built incrementally -- this first piece is just the per-term test.
BH correction, fold-enrichment direction, and the full godb-driven
multi-term loop come next.
"""
from scipy import stats as scipy_stats


def fisher_exact_for_term(n_pop_with_term, n_pop_total, n_class_with_term, n_class_total):
    """
    Runs Fisher's exact test for one GO term, matching the exact 2x2 table
    construction confirmed from stats_test_based_on_go.py:

        [[n_pop_with_term,   n_pop_total - n_pop_with_term],
         [n_class_with_term, n_class_total - n_class_with_term]]

    Returns (odds_ratio, p_value).
    """
    table = [
        [n_pop_with_term, n_pop_total - n_pop_with_term],
        [n_class_with_term, n_class_total - n_class_with_term],
    ]
    return scipy_stats.fisher_exact(table)