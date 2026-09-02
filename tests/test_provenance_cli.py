"""
tests/test_provenance_cli.py

Confirms scripts/run_pipeline.py actually calls enrichment.output's
provenance support with real run parameters, end to end through the CLI,
not just that enrichment/output.py's own functions work in isolation.
"""
import json
import pickle
import subprocess
import sys
from pathlib import Path

import pandas as pd


REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_PATH = REPO_ROOT / "scripts" / "run_pipeline.py"


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


def test_cli_output_files_carry_a_provenance_line(tmp_path):
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
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(output_dir),
        "--dataset-name", "provenance_test",
        "--thresh-type", "p", "--thresh", "0.01",
    ])
    assert result.returncode == 0, result.stderr

    all_file = output_dir / "summary.provenance_test.all.tsv"
    with open(all_file) as f:
        first_line = f.readline()

    assert first_line.startswith("#"), "expected a provenance line, output looks unstamped"
    prov = json.loads(first_line[1:])

    assert prov["tool"] == "go-gsea"
    assert prov["run_parameters"]["label_strategy"] == "rank_tail"
    assert prov["run_parameters"]["metric_col"] == ["metric"]
    assert prov["run_parameters"]["thresh_type"] == "p"
    assert prov["run_parameters"]["thresh"] == 0.01
    assert prov["run_parameters"]["dataset_name"] == "provenance_test"

    # the file must still parse correctly with the provenance line present
    loaded = pd.read_csv(all_file, sep="\t", comment="#")
    assert len(loaded) > 0


def test_cli_provenance_run_parameters_reflect_merge_manifest_when_used(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    output_dir = tmp_path / "results"
    make_synthetic_godb(godb_path)
    pd.DataFrame({
        "gene_id": [f"gene{i}" for i in range(10)],
        "metric": list(range(10)),
    }).to_csv(input_path, sep="\t", index=False)

    manifest_path = tmp_path / "manifest.txt"
    manifest_path.write_text(
        f"[cond_a]\n"
        f"path = {input_path}\n"
        f"value_col = metric\n"
    )

    result = run_cli([
        "--merge-manifest", str(manifest_path),
        "--godb", str(godb_path),
        "--metric-col", "cond_a",
        "--label-strategy", "rank_tail", "--pct", "20",
        "--output-dir", str(output_dir),
        "--dataset-name", "merge_provenance_test",
    ])
    assert result.returncode == 0, result.stderr

    all_file = output_dir / "summary.merge_provenance_test.all.tsv"
    with open(all_file) as f:
        first_line = f.readline()
    prov = json.loads(first_line[1:])

    assert prov["run_parameters"]["merge_manifest"] == str(manifest_path)
    assert prov["run_parameters"]["input_table"] is None