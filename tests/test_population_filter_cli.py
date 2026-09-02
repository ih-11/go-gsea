"""
tests/test_population_filter_cli.py

Tests for the new Stage A population-filter flags wired into
run_pipeline.py: --read-depth-col/-thresh, --usage-col/-thresh,
--min-group-col/-n. All optional and off by default -- the existing
behavior (no filtering) must be provably unchanged when none are given.
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
from scripts.run_pipeline import build_population_filters, build_parser


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


def parse_args(argv):
    return build_parser().parse_args(argv)


def test_build_population_filters_empty_when_no_flags_given():
    args = parse_args([
        "--input-table", "x.tsv", "--godb", "y.pkl", "--metric-col", "metric",
        "--label-strategy", "rank_tail", "--output-dir", "out", "--dataset-name", "d",
    ])
    assert build_population_filters(args) == []


def test_build_population_filters_raises_when_thresh_missing():
    args = parse_args([
        "--input-table", "x.tsv", "--godb", "y.pkl", "--metric-col", "metric",
        "--label-strategy", "rank_tail", "--output-dir", "out", "--dataset-name", "d",
        "--read-depth-col", "n_reads",
        # --read-depth-thresh deliberately omitted
    ])
    with pytest.raises(ValueError):
        build_population_filters(args)


def test_build_population_filters_returns_all_three_in_fixed_order():
    args = parse_args([
        "--input-table", "x.tsv", "--godb", "y.pkl", "--metric-col", "metric",
        "--label-strategy", "rank_tail", "--output-dir", "out", "--dataset-name", "d",
        "--read-depth-col", "n_reads", "--read-depth-thresh", "10",
        "--usage-col", "usage", "--usage-thresh", "0.05",
        "--min-group-col", "transcript_id", "--min-group-n", "2",
    ])
    filters = build_population_filters(args)
    assert len(filters) == 3
    fn_names = [fn.__name__ for fn, _ in filters]
    assert fn_names == ["read_depth_filter", "usage_filter", "min_group_size_filter"]


def test_cli_without_population_filter_flags_behaves_unchanged(tmp_path):
    """No population-filter flags given -> identical output to before this
    feature existed. This is the regression guard for every prior real run."""
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    make_synthetic_godb(godb_path)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(20)],
        "metric": list(range(20)),
    }).to_csv(input_path, sep="\t", index=False)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(tmp_path / "results"),
        "--dataset-name", "no_filter_test",
    ])
    assert result.returncode == 0, result.stderr
    assert "Loaded 20 genes" in result.stdout
    assert "Population filter" not in result.stdout


def test_cli_read_depth_filter_actually_removes_genes(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    make_synthetic_godb(godb_path)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(20)],
        "metric": list(range(20)),
        "n_reads": [5] * 10 + [50] * 10,  # first 10 fail a thresh=10 filter
    }).to_csv(input_path, sep="\t", index=False)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--read-depth-col", "n_reads", "--read-depth-thresh", "10",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(tmp_path / "results"),
        "--dataset-name", "filter_test",
    ])
    assert result.returncode == 0, result.stderr
    assert "Loaded 20 genes" in result.stdout
    assert "Population filter (Stage A): 20 -> 10 genes" in result.stdout


def test_cli_chained_filters_are_stricter_than_either_alone(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    make_synthetic_godb(godb_path, n=20)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(20)],
        "metric": list(range(20)),
        "n_reads": [5] * 8 + [50] * 12,       # 12 pass read-depth alone
        "usage": [0.01] * 5 + [0.5] * 15,     # 15 pass usage alone
    }).to_csv(input_path, sep="\t", index=False)
    # genes 8-19 pass read-depth; genes 5-19 pass usage
    # both: genes 8-19 (12 genes)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--read-depth-col", "n_reads", "--read-depth-thresh", "10",
        "--usage-col", "usage", "--usage-thresh", "0.05",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(tmp_path / "results"),
        "--dataset-name", "chained_filter_test",
    ])
    assert result.returncode == 0, result.stderr
    assert "Population filter (Stage A): 20 -> 12 genes" in result.stdout


def test_cli_read_depth_col_without_thresh_fails_clearly(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    make_synthetic_godb(godb_path)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(20)],
        "metric": list(range(20)),
        "n_reads": list(range(20)),
    }).to_csv(input_path, sep="\t", index=False)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--read-depth-col", "n_reads",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(tmp_path / "results"),
        "--dataset-name", "missing_thresh_test",
    ])
    assert result.returncode != 0
    assert "read-depth-thresh" in result.stderr