"""
tests/test_run_pipeline.py

Tests scripts/run_pipeline.py as a real CLI, via subprocess, against small
synthetic data built fresh in tmp_path for each test. No dependency on
real Chlamydomonas files, this stays fast and fully self-contained.
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
    """A tiny, hand-built .godb: 10 genes, 2 GO terms, no propagation needed."""
    gene_go = {
        f"gene{i}": {"GO:0000001"} if i < 5 else {"GO:0000002"}
        for i in range(10)
    }
    godb = {"gene_go": gene_go, "provenance": {}, "source": "synthetic", "obo": "synthetic"}
    with open(path, "wb") as f:
        pickle.dump(godb, f)


def make_synthetic_input_table(path, extra_unmatched_ids=0):
    """10 genes with a numeric metric column, optionally padded with IDs
    that won't match the synthetic godb, to exercise unknown_ratio_thresh."""
    rows = {
        "gene_id": [f"gene{i}" for i in range(10)],
        "metric": list(range(10)),
    }
    for j in range(extra_unmatched_ids):
        rows["gene_id"].append(f"unmatched{j}")
        rows["metric"].append(999)
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def run_cli(args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH)] + args,
        capture_output=True, text=True,
    )


def test_cli_writes_expected_output_files(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    output_dir = tmp_path / "results"

    make_synthetic_godb(godb_path)
    make_synthetic_input_table(input_path)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(output_dir),
        "--dataset-name", "synthtest",
    ])

    assert result.returncode == 0, result.stderr
    all_file = output_dir / "summary.synthtest.all.tsv"
    assert all_file.exists()
    df = pd.read_csv(all_file, sep="\t")
    assert set(df["class"].unique()) <= {"High", "Low"}


def test_cli_slim_godb_is_optional_and_writes_separately(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    slim_godb_path = tmp_path / "test.slim.godb.pkl"
    input_path = tmp_path / "input.tsv"
    output_dir = tmp_path / "results" / "GO"
    slim_output_dir = tmp_path / "results" / "GO_slim"

    make_synthetic_godb(godb_path)
    make_synthetic_godb(slim_godb_path)  # reuse the same builder, contents don't matter here
    make_synthetic_input_table(input_path)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--slim-godb", str(slim_godb_path),
        "--metric-col", "metric",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(output_dir),
        "--slim-output-dir", str(slim_output_dir),
        "--dataset-name", "synthtest",
    ])

    assert result.returncode == 0, result.stderr
    assert (output_dir / "summary.synthtest.all.tsv").exists()
    assert (slim_output_dir / "summary.synthtest.all.tsv").exists()


def test_cli_without_slim_godb_does_not_create_slim_output(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    output_dir = tmp_path / "results"

    make_synthetic_godb(godb_path)
    make_synthetic_input_table(input_path)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(output_dir),
        "--dataset-name", "synthtest",
    ])

    assert result.returncode == 0, result.stderr
    assert "[GO-slim]" not in result.stdout


def test_cli_missing_required_flag_fails_with_clear_error(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    make_synthetic_godb(godb_path)
    make_synthetic_input_table(input_path)

    # --output-dir deliberately omitted
    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--dataset-name", "synthtest",
    ])

    assert result.returncode != 0
    assert "output-dir" in result.stderr or "output_dir" in result.stderr


def test_cli_slim_godb_without_slim_output_dir_fails_clearly(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    slim_godb_path = tmp_path / "test.slim.godb.pkl"
    input_path = tmp_path / "input.tsv"
    output_dir = tmp_path / "results"

    make_synthetic_godb(godb_path)
    make_synthetic_godb(slim_godb_path)
    make_synthetic_input_table(input_path)

    # --slim-godb given, --slim-output-dir deliberately omitted
    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--slim-godb", str(slim_godb_path),
        "--metric-col", "metric",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(output_dir),
        "--dataset-name", "synthtest",
    ])

    assert result.returncode != 0
    assert "slim-output-dir" in result.stderr or "slim_output_dir" in result.stderr


def test_cli_strip_id_suffix_and_exclude_id(tmp_path):
    """
    --strip-id-suffix is applied BEFORE --exclude-id (see run_pipeline.py's
    load_input_table docstring), so --exclude-id must be given in its
    already-stripped form here, not the raw suffixed ID, or the exclusion
    will silently fail to match.
    """
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    output_dir = tmp_path / "results"

    make_synthetic_godb(godb_path)
    # gene IDs carry a version suffix, plus one row to be excluded outright
    df = pd.DataFrame({
        "gene_id": [f"gene{i}.v1.0" for i in range(10)] + ["spike_in_control.v1.0"],
        "metric": list(range(10)) + [999],
    })
    df.to_csv(input_path, sep="\t", index=False)

    result = run_cli([
        "--input-table", str(input_path),
        "--godb", str(godb_path),
        "--metric-col", "metric",
        "--id-col", "gene_id",
        "--strip-id-suffix", r"\.v\d+\.\d+$",
        "--exclude-id", "spike_in_control",  # post-strip form, not the raw .v1.0-suffixed ID
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(output_dir),
        "--dataset-name", "synthtest",
    ])

    assert result.returncode == 0, result.stderr
    assert "Loaded 10 genes" in result.stdout