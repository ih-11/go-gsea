"""
enrichment/ora.py

Over-representation analysis (ORA): Fisher's exact test per GO term,
comparing a labeled class against the full population.
"""
import math
from scipy import stats as scipy_stats


def fisher_exact_for_term(n_pop_with_term, n_pop_total, n_class_with_term, n_class_total):
    """
    Runs Fisher's exact test for one GO term. Table layout matches
    stats_test_based_on_go.py exactly: population row first, class row second.

    Returns p_value only -- direction is NOT taken from scipy's odds_ratio
    (its sign/magnitude depends on row orientation, which is confusing and
    not how the real precedent determines direction -- see fold_enrichment
    below instead).
    """
    table = [
        [n_pop_with_term, n_pop_total - n_pop_with_term],
        [n_class_with_term, n_class_total - n_class_with_term],
    ]
    _, p_value = scipy_stats.fisher_exact(table)
    return p_value


def fold_enrichment(n_pop_with_term, n_pop_total, n_class_with_term, n_class_total):
    """
    log2 fold-enrichment, matching stats_test_based_on_go.py exactly:
    log2((observed+1)/(expected+1)), where expected is the class-size-scaled
    count you'd see if the class matched the population's proportion.

    >0 means over-represented in the class, <0 means under-represented.
    """
    pop_ratio = n_pop_with_term / n_pop_total
    n_expected = pop_ratio * n_class_total
    return math.log2((n_class_with_term + 1) / (n_expected + 1))