"""
tests/test_drop_uninformative_cli.py

Tests for the new --drop-uninformative CLI flag in run_pipeline.py, which
formalizes the real by-hand zero-floor artifact check documented in
docs/examples/chlamydomonas.md.
"""
import pickle
import subprocess
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_pipeline.py"


def make_synthetic_godb(path, n=20):
    gene_go = {
        f"gene{i}": {"GO:0000001"} if i < n // 2 else {"GO:0000002"}
        for i in range(n)
    }
    godb = {"gene_go": gene_go, "provenance": {}, "source": "synthetic", "obo": "synthetic"}
    with open(path, "wb") as f:
        pickle.dump(godb, f)


def run_cli(args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH)] + args,
        capture_output=True, text=True,
    )


def test_cli_without_drop_uninformative_flag_behaves_unchanged(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    make_synthetic_godb(godb_path)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(20)],
        "metric": [0] * 5 + list(range(1, 16)),
    }).to_csv(input_path, sep="\t", index=False)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(tmp_path / "results"),
        "--dataset-name", "no_drop_test",
    ])
    assert result.returncode == 0, result.stderr
    assert "Loaded 20 genes" in result.stdout
    assert "Dropped uninformative rows" not in result.stdout


def test_cli_drop_uninformative_removes_zero_rows_single_column(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    make_synthetic_godb(godb_path)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(20)],
        "metric": [0] * 6 + list(range(1, 15)),
    }).to_csv(input_path, sep="\t", index=False)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--drop-uninformative",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(tmp_path / "results"),
        "--dataset-name", "drop_test",
    ])
    assert result.returncode == 0, result.stderr
    assert "Dropped uninformative rows (all metric columns <= 0.0): 20 -> 14 genes" in result.stdout


def test_cli_drop_uninformative_requires_ALL_columns_zero_for_multi_column(tmp_path):
    """The multi-column case: cluster on two columns, only drop genes that
    are zero in BOTH, not either. Note: the original real-data artifact
    this filter formalizes (see docs/examples/chlamydomonas.md demo C's
    k=13 rejection) was actually SINGLE-column (1D), the case covered by
    test_cli_drop_uninformative_removes_zero_rows_single_column above --
    this test just confirms the same logic generalizes correctly to more
    than one column, not that multi-column is the primary motivating case."""
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    make_synthetic_godb(godb_path)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(20)],
        "col_a": [0] * 5 + [0] * 5 + list(range(1, 11)),
        "col_b": [0] * 5 + list(range(1, 6)) + list(range(1, 11)),
        # genes 0-4: 0,0 (uninformative). genes 5-9: 0, nonzero (informative).
        # genes 10-19: nonzero, nonzero (informative).
    }).to_csv(input_path, sep="\t", index=False)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "col_a", "col_b",
        "--drop-uninformative",
        "--label-strategy", "cluster", "--n-clusters", "2",
        "--output-dir", str(tmp_path / "results"),
        "--dataset-name", "drop_2col_test",
    ])
    assert result.returncode == 0, result.stderr
    assert "Dropped uninformative rows (all metric columns <= 0.0): 20 -> 15 genes" in result.stdout


def test_cli_drop_uninformative_respects_custom_threshold(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    make_synthetic_godb(godb_path)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(20)],
        "metric": [0.5] * 8 + list(range(2, 14)),
    }).to_csv(input_path, sep="\t", index=False)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--drop-uninformative", "--drop-uninformative-thresh", "1.0",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(tmp_path / "results"),
        "--dataset-name", "drop_custom_thresh_test",
    ])
    assert result.returncode == 0, result.stderr
    assert "Dropped uninformative rows (all metric columns <= 1.0): 20 -> 12 genes" in result.stdout