"""
tests/test_merge_manifest_cli.py

Tests for --merge-manifest support in run_pipeline.py: the new
multi-source path, kept separate from the existing single-source
--input-table tests since it exercises genuinely new code
(read_merge_manifest, load_merged_table) on top of everything already
tested via --input-table.
"""
import pickle
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_pipeline.py"


def make_synthetic_godb(path):
    gene_go = {
        f"gene{i}": {"GO:0000001"} if i < 10 else {"GO:0000002"}
        for i in range(20)
    }
    godb = {"gene_go": gene_go, "provenance": {}, "source": "synthetic", "obo": "synthetic"}
    with open(path, "wb") as f:
        pickle.dump(godb, f)


def run_cli(args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH)] + args,
        capture_output=True, text=True,
    )


def test_cli_input_table_and_merge_manifest_are_mutually_exclusive(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    make_synthetic_godb(godb_path)
    input_path = tmp_path / "input.tsv"
    pd.DataFrame({"gene_id": ["gene0"], "metric": [1.0]}).to_csv(input_path, sep="\t", index=False)
    manifest_path = tmp_path / "merge.txt"
    manifest_path.write_text("[cond_a]\npath = x.tsv\nvalue_col = metric\n")

    result = run_cli([
        "--input-table", str(input_path),
        "--merge-manifest", str(manifest_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(tmp_path / "results"),
        "--dataset-name", "synthtest",
    ])
    assert result.returncode != 0


def test_cli_neither_input_table_nor_merge_manifest_fails(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    make_synthetic_godb(godb_path)

    result = run_cli([
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(tmp_path / "results"),
        "--dataset-name", "synthtest",
    ])
    assert result.returncode != 0


def test_cli_merge_manifest_two_conditions_cluster(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    make_synthetic_godb(godb_path)

    cond_a_path = tmp_path / "cond_a.tsv"
    cond_b_path = tmp_path / "cond_b.tsv"
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(20)],
        "PR_gene": list(range(20)),
    }).to_csv(cond_a_path, sep="\t", index=False)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(20)],
        "PR_gene": list(range(20, 40)),
    }).to_csv(cond_b_path, sep="\t", index=False)

    manifest_path = tmp_path / "merge.txt"
    manifest_path.write_text(
        f"[cond_a]\n"
        f"path = {cond_a_path}\n"
        f"value_col = PR_gene\n"
        f"\n"
        f"[cond_b]\n"
        f"path = {cond_b_path}\n"
        f"value_col = PR_gene\n"
    )

    output_dir = tmp_path / "results"
    result = run_cli([
        "--merge-manifest", str(manifest_path),
        "--godb", str(godb_path),
        "--metric-col", "cond_a", "cond_b",
        "--label-strategy", "cluster", "--n-clusters", "2",
        "--output-dir", str(output_dir),
        "--dataset-name", "merged_cluster_test",
    ])
    assert result.returncode == 0, result.stderr
    assert "Loaded 20 genes (merged from" in result.stdout
    assert (output_dir / "summary.merged_cluster_test.all.tsv").exists()


def test_cli_merge_manifest_only_keeps_genes_present_in_both_sources(tmp_path):
    """Inner join is the default -- a gene missing from one source must
    not survive into the merged table."""
    godb_path = tmp_path / "test.godb.pkl"
    make_synthetic_godb(godb_path)

    cond_a_path = tmp_path / "cond_a.tsv"
    cond_b_path = tmp_path / "cond_b.tsv"
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(20)],
        "PR_gene": list(range(20)),
    }).to_csv(cond_a_path, sep="\t", index=False)
    # cond_b is missing gene19 entirely
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(19)],
        "PR_gene": list(range(20, 39)),
    }).to_csv(cond_b_path, sep="\t", index=False)

    manifest_path = tmp_path / "merge.txt"
    manifest_path.write_text(
        f"[cond_a]\n"
        f"path = {cond_a_path}\n"
        f"value_col = PR_gene\n"
        f"\n"
        f"[cond_b]\n"
        f"path = {cond_b_path}\n"
        f"value_col = PR_gene\n"
    )

    output_dir = tmp_path / "results"
    result = run_cli([
        "--merge-manifest", str(manifest_path),
        "--godb", str(godb_path),
        "--metric-col", "cond_a", "cond_b",
        "--label-strategy", "cluster", "--n-clusters", "2",
        "--output-dir", str(output_dir),
        "--dataset-name", "inner_join_test",
    ])
    assert result.returncode == 0, result.stderr
    assert "Loaded 19 genes (merged from" in result.stdout  # gene19 correctly dropped


def test_cli_merge_manifest_missing_file_fails_clearly(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    make_synthetic_godb(godb_path)

    result = run_cli([
        "--merge-manifest", str(tmp_path / "does_not_exist.txt"),
        "--godb", str(godb_path),
        "--metric-col", "cond_a", "cond_b",
        "--label-strategy", "cluster", "--n-clusters", "2",
        "--output-dir", str(tmp_path / "results"),
        "--dataset-name", "test",
    ])
    assert result.returncode != 0