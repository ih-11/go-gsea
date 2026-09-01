"""
tests/test_merge_tables.py
"""
import pandas as pd
import pytest

from dataprep.merge_tables import merge_gene_tables


def write_table(path, gene_ids, values, value_col="PR_gene"):
    pd.DataFrame({"gene_id": gene_ids, value_col: values}).to_csv(path, sep="\t", index=False)


def test_merge_two_sources_inner_join_keeps_only_shared_genes(tmp_path):
    path_a = tmp_path / "a.tsv"
    path_b = tmp_path / "b.tsv"
    write_table(path_a, ["g1", "g2", "g3"], [1.0, 2.0, 3.0])
    write_table(path_b, ["g2", "g3", "g4"], [20.0, 30.0, 40.0])

    result = merge_gene_tables([
        (str(path_a), "PR_gene", "cond_a"),
        (str(path_b), "PR_gene", "cond_b"),
    ])

    assert set(result["gene_id"]) == {"g2", "g3"}
    assert list(result.columns) == ["gene_id", "cond_a", "cond_b"]


def test_merge_preserves_correct_values_per_gene(tmp_path):
    path_a = tmp_path / "a.tsv"
    path_b = tmp_path / "b.tsv"
    write_table(path_a, ["g1", "g2"], [1.0, 2.0])
    write_table(path_b, ["g1", "g2"], [10.0, 20.0])

    result = merge_gene_tables([
        (str(path_a), "PR_gene", "cond_a"),
        (str(path_b), "PR_gene", "cond_b"),
    ])
    result = result.set_index("gene_id")

    assert result.loc["g1", "cond_a"] == 1.0
    assert result.loc["g1", "cond_b"] == 10.0
    assert result.loc["g2", "cond_a"] == 2.0
    assert result.loc["g2", "cond_b"] == 20.0


def test_merge_three_sources(tmp_path):
    path_a = tmp_path / "a.tsv"
    path_b = tmp_path / "b.tsv"
    path_c = tmp_path / "c.tsv"
    write_table(path_a, ["g1", "g2", "g3"], [1.0, 2.0, 3.0])
    write_table(path_b, ["g1", "g2", "g3"], [10.0, 20.0, 30.0])
    write_table(path_c, ["g1", "g2"], [100.0, 200.0])  # g3 missing here

    result = merge_gene_tables([
        (str(path_a), "PR_gene", "cond_a"),
        (str(path_b), "PR_gene", "cond_b"),
        (str(path_c), "PR_gene", "cond_c"),
    ])

    # inner join across all three -- g3 must drop out since path_c lacks it
    assert set(result["gene_id"]) == {"g1", "g2"}
    assert list(result.columns) == ["gene_id", "cond_a", "cond_b", "cond_c"]


def test_merge_different_value_col_names_across_sources(tmp_path):
    path_a = tmp_path / "a.tsv"
    path_b = tmp_path / "b.tsv"
    write_table(path_a, ["g1", "g2"], [1.0, 2.0], value_col="TPM")
    write_table(path_b, ["g1", "g2"], [10.0, 20.0], value_col="PTPM")

    result = merge_gene_tables([
        (str(path_a), "TPM", "tpm_cond_a"),
        (str(path_b), "PTPM", "ptpm_cond_a"),
    ])
    assert list(result.columns) == ["gene_id", "tpm_cond_a", "ptpm_cond_a"]


def test_merge_raises_on_empty_sources():
    with pytest.raises(ValueError):
        merge_gene_tables([])


def test_merge_raises_on_duplicate_output_col_names(tmp_path):
    path_a = tmp_path / "a.tsv"
    path_b = tmp_path / "b.tsv"
    write_table(path_a, ["g1"], [1.0])
    write_table(path_b, ["g1"], [10.0])

    with pytest.raises(ValueError):
        merge_gene_tables([
            (str(path_a), "PR_gene", "same_name"),
            (str(path_b), "PR_gene", "same_name"),
        ])


def test_merge_raises_on_missing_value_col(tmp_path):
    path_a = tmp_path / "a.tsv"
    write_table(path_a, ["g1"], [1.0], value_col="PR_gene")

    with pytest.raises(ValueError):
        merge_gene_tables([(str(path_a), "TPM", "cond_a")])  # TPM doesn't exist in this file


def test_merge_supports_gzipped_input(tmp_path):
    path_a = tmp_path / "a.tsv.gz"
    path_b = tmp_path / "b.tsv.gz"
    pd.DataFrame({"gene_id": ["g1", "g2"], "PR_gene": [1.0, 2.0]}).to_csv(
        path_a, sep="\t", index=False, compression="gzip"
    )
    pd.DataFrame({"gene_id": ["g1", "g2"], "PR_gene": [10.0, 20.0]}).to_csv(
        path_b, sep="\t", index=False, compression="gzip"
    )

    result = merge_gene_tables([
        (str(path_a), "PR_gene", "cond_a"),
        (str(path_b), "PR_gene", "cond_b"),
    ])
    assert set(result["gene_id"]) == {"g1", "g2"}


def test_merge_outer_join_keeps_union_when_requested(tmp_path):
    path_a = tmp_path / "a.tsv"
    path_b = tmp_path / "b.tsv"
    write_table(path_a, ["g1", "g2"], [1.0, 2.0])
    write_table(path_b, ["g2", "g3"], [20.0, 30.0])

    result = merge_gene_tables([
        (str(path_a), "PR_gene", "cond_a"),
        (str(path_b), "PR_gene", "cond_b"),
    ], how="outer")

    assert set(result["gene_id"]) == {"g1", "g2", "g3"}