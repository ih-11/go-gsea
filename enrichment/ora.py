"""
enrichment/ora.py

Over-representation analysis (ORA): Fisher's exact test per GO term,
comparing a labeled class against the full population.
"""
import math
import pandas as pd
from scipy import stats as scipy_stats
from statsmodels.stats.multitest import multipletests


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


def bh_correct(p_values):
    """
    Benjamini-Hochberg FDR correction across a list of p-values.
    Matches stats_test_based_on_go.py exactly: multipletests(..., method='fdr_bh').

    Returns a list of q-values in the same order as the input.
    """
    if len(p_values) == 0:
        return []
    _, q_values, _, _ = multipletests(p_values, method="fdr_bh")
    return list(q_values)


def count_term_occurrences(gene_go, gene_ids):
    """
    Given gene_go (dict[gene_id] -> set(GO_id), from a .godb) and a list of
    gene_ids (e.g. the population, or one class), returns dict[GO_id] -> count
    of how many of those genes carry that term.
    """
    counts = {}
    for gene_id in gene_ids:
        for term_id in gene_go.get(gene_id, set()):
            counts[term_id] = counts.get(term_id, 0) + 1
    return counts


def restrict_to_annotated_genes(labeled_df, gene_go, id_col="gene_id"):
    """
    Matches stats_test_based_on_go.py's population_list construction: keeps
    only rows whose gene_id has at least one entry in gene_go, since genes
    with zero GO annotation can never contribute to any term's count and
    must not inflate the population/class denominators.

    Returns (filtered_df, n_unknown, unknown_ratio).
    """
    is_annotated = labeled_df[id_col].isin(gene_go.keys())
    n_unknown = int((~is_annotated).sum())
    unknown_ratio = n_unknown / len(labeled_df) if len(labeled_df) else 0.0
    return labeled_df[is_annotated].copy(), n_unknown, unknown_ratio


def run_ora(gene_go, labeled_df, id_col="gene_id", class_col="class",
            thresh_type="p", thresh=0.01, unknown_ratio_thresh=0.9):
    """
    Runs GO over-representation analysis for every class present in
    labeled_df's class_col, against the full labeled_df as population
    (matching the -c input format: every row is population, non-null
    class rows are also tested as a foreground group).

    Genes with no entry in gene_go are dropped from the population BEFORE
    any counting happens (matches stats_test_based_on_go.py's population_list
    construction) -- otherwise every expected-count and Fisher's-exact
    denominator would be inflated by genes that can never contribute to any
    term's count, silently understating every real enrichment signal.

    Raises ValueError if the unknown-gene ratio exceeds unknown_ratio_thresh.

    Returns a DataFrame with one row per (class, GO_term) tested, columns
    matching the real stats_test_based_on_go.py output: go_id, class,
    population, observed, expected, fold_enrichment, p_value, q_value,
    significance.
    """
    labeled_df, n_unknown, unknown_ratio = restrict_to_annotated_genes(
        labeled_df, gene_go, id_col
    )
    if unknown_ratio >= unknown_ratio_thresh:
        raise ValueError(
            f"Too many genes ({n_unknown}, {unknown_ratio:.1%}) have no GO "
            f"annotation -- exceeds unknown_ratio_thresh={unknown_ratio_thresh:.1%}."
        )

    population_ids = labeled_df[id_col].tolist()
    pop_counts = count_term_occurrences(gene_go, population_ids)
    n_pop_total = len(population_ids)

    classes = [c for c in labeled_df[class_col].dropna().unique()]

    rows = []
    for cls in classes:
        class_ids = labeled_df.loc[labeled_df[class_col] == cls, id_col].tolist()
        n_class_total = len(class_ids)
        class_counts = count_term_occurrences(gene_go, class_ids)

        class_rows = []
        for term_id, n_pop_with_term in pop_counts.items():
            n_class_with_term = class_counts.get(term_id, 0)
            p_value = fisher_exact_for_term(
                n_pop_with_term, n_pop_total, n_class_with_term, n_class_total
            )
            fe = fold_enrichment(
                n_pop_with_term, n_pop_total, n_class_with_term, n_class_total
            )
            pop_ratio = n_pop_with_term / n_pop_total
            n_expected = pop_ratio * n_class_total
            class_rows.append({
                "go_id": term_id,
                "class": cls,
                "population": n_pop_with_term,
                "observed": n_class_with_term,
                "expected": n_expected,
                "fold_enrichment": fe,
                "p_value": p_value,
            })

        q_values = bh_correct([r["p_value"] for r in class_rows])
        for row, q in zip(class_rows, q_values):
            row["q_value"] = q

        rows.extend(class_rows)

    result = pd.DataFrame(rows)
    if len(result) == 0:
        return result

    sig_col = "p_value" if thresh_type == "p" else "q_value"
    result["significance"] = result[sig_col] < thresh
    return result