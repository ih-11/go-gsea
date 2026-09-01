"""
tests/test_multi_metric_col.py

Tests specifically for the --metric-col nargs='+' widening: validation
that single-column strategies reject multiple values, and that cluster
genuinely accepts and uses more than one.
"""
import pickle
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_pipeline.py"

sys.path.insert(0, str(REPO_ROOT))
from scripts.run_pipeline import validate_metric_cols


def make_synthetic_godb(path):
    gene_go = {
        f"gene{i}": {"GO:0000001"} if i < 5 else {"GO:0000002"}
        for i in range(10)
    }
    godb = {"gene_go": gene_go, "provenance": {}, "source": "synthetic", "obo": "synthetic"}
    with open(path, "wb") as f:
        pickle.dump(godb, f)


def run_cli(args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH)] + args,
        capture_output=True, text=True,
    )


def test_validate_metric_cols_rejects_multiple_for_rank_tail():
    with pytest.raises(ValueError):
        validate_metric_cols(["col_a", "col_b"], "rank_tail")


def test_validate_metric_cols_rejects_multiple_for_explicit_threshold():
    with pytest.raises(ValueError):
        validate_metric_cols(["col_a", "col_b"], "explicit_threshold")


def test_validate_metric_cols_rejects_multiple_for_boolean_flag():
    with pytest.raises(ValueError):
        validate_metric_cols(["col_a", "col_b"], "boolean_flag")


def test_validate_metric_cols_accepts_one_for_single_column_strategies():
    validate_metric_cols(["col_a"], "rank_tail")  # should not raise


def test_validate_metric_cols_accepts_multiple_for_cluster():
    validate_metric_cols(["col_a", "col_b", "col_c"], "cluster")  # should not raise


def test_validate_metric_cols_accepts_one_for_cluster_too():
    validate_metric_cols(["col_a"], "cluster")  # should not raise


def test_cli_rejects_multiple_metric_cols_for_rank_tail(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    make_synthetic_godb(godb_path)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(10)],
        "metric_a": list(range(10)),
        "metric_b": list(range(10, 20)),
    }).to_csv(input_path, sep="\t", index=False)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric_a", "metric_b",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(tmp_path / "results"),
        "--dataset-name", "synthtest",
    ])
    assert result.returncode != 0
    assert "exactly one" in result.stderr


def test_cli_cluster_accepts_multiple_metric_cols(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    output_dir = tmp_path / "results"
    make_synthetic_godb(godb_path)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(10)],
        "metric_a": list(range(10)),
        "metric_b": list(range(10, 20)),
    }).to_csv(input_path, sep="\t", index=False)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric_a", "metric_b",
        "--label-strategy", "cluster", "--n-clusters", "2",
        "--output-dir", str(output_dir),
        "--dataset-name", "synthtest_cluster",
    ])
    assert result.returncode == 0, result.stderr
    assert (output_dir / "summary.synthtest_cluster.all.tsv").exists()


def test_cli_cluster_still_works_with_single_metric_col(tmp_path):
    """Backward compatibility: cluster with exactly one --metric-col
    (the way every real run so far, including the real k=4 CR_3D run,
    has actually used it) must still work unchanged."""
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    output_dir = tmp_path / "results"
    make_synthetic_godb(godb_path)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(10)],
        "metric": list(range(10)),
    }).to_csv(input_path, sep="\t", index=False)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--label-strategy", "cluster", "--n-clusters", "2",
        "--output-dir", str(output_dir),
        "--dataset-name", "synthtest_single",
    ])
    assert result.returncode == 0, result.stderr
    assert (output_dir / "summary.synthtest_single.all.tsv").exists()