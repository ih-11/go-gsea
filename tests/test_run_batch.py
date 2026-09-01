"""
tests/test_run_batch.py
"""
import pickle
import pandas as pd
import pytest

from scripts.run_batch import build_argv_from_row, run_batch


def make_synthetic_godb(path):
    gene_go = {
        f"gene{i}": {"GO:0000001"} if i < 5 else {"GO:0000002"}
        for i in range(10)
    }
    godb = {"gene_go": gene_go, "provenance": {}, "source": "synthetic", "obo": "synthetic"}
    with open(path, "wb") as f:
        pickle.dump(godb, f)


def make_synthetic_input_table(path):
    rows = {"gene_id": [f"gene{i}" for i in range(10)], "metric": list(range(10))}
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False)


def test_build_argv_from_row_auto_derives_dataset_name():
    row = pd.Series({
        "species": "AT", "condition": "HS", "input_table": "in.tsv",
        "godb": "at.godb.pkl", "metric_col": "TPM", "label_strategy": "rank_tail",
    })
    argv = build_argv_from_row(row)
    assert "--dataset-name" in argv
    assert argv[argv.index("--dataset-name") + 1] == "AT_HS.TPM"


def test_build_argv_from_row_respects_explicit_dataset_name():
    row = pd.Series({
        "species": "AT", "condition": "HS", "input_table": "in.tsv",
        "godb": "at.godb.pkl", "metric_col": "TPM", "label_strategy": "rank_tail",
        "dataset_name": "custom_name",
    })
    argv = build_argv_from_row(row)
    assert argv[argv.index("--dataset-name") + 1] == "custom_name"


def test_build_argv_from_row_raises_on_missing_required_column():
    row = pd.Series({
        "species": "AT", "condition": "HS", "input_table": "in.tsv",
        "godb": "at.godb.pkl", "metric_col": "TPM",
        # label_strategy deliberately missing
    })
    with pytest.raises(ValueError):
        build_argv_from_row(row)


def test_build_argv_from_row_includes_slim_godb_when_given():
    row = pd.Series({
        "species": "AT", "condition": "HS", "input_table": "in.tsv",
        "godb": "at.godb.pkl", "metric_col": "TPM", "label_strategy": "rank_tail",
        "slim_godb": "at.slim.godb.pkl",
    })
    argv = build_argv_from_row(row)
    assert "--slim-godb" in argv
    assert "--slim-output-dir" in argv


def test_build_argv_from_row_omits_slim_godb_when_not_given():
    row = pd.Series({
        "species": "AT", "condition": "HS", "input_table": "in.tsv",
        "godb": "at.godb.pkl", "metric_col": "TPM", "label_strategy": "rank_tail",
        "slim_godb": None,
    })
    argv = build_argv_from_row(row)
    assert "--slim-godb" not in argv


def test_build_argv_from_row_handles_multiple_exclude_ids():
    row = pd.Series({
        "species": "AT", "condition": "HS", "input_table": "in.tsv",
        "godb": "at.godb.pkl", "metric_col": "TPM", "label_strategy": "rank_tail",
        "exclude_id": "id1;id2;id3",
    })
    argv = build_argv_from_row(row)
    assert argv.count("--exclude-id") == 3


def test_run_batch_dry_run_does_not_execute_anything(tmp_path):
    manifest_path = tmp_path / "manifest.csv"
    pd.DataFrame([{
        "species": "AT", "condition": "HS",
        "input_table": "nonexistent.tsv", "godb": "nonexistent.godb.pkl",
        "metric_col": "TPM", "label_strategy": "rank_tail",
    }]).to_csv(manifest_path, index=False)

    summary = run_batch(str(manifest_path), dry_run=True)
    assert summary["succeeded"] == ["AT_HS.TPM"]
    assert summary["failed"] == []


def test_run_batch_real_run_two_rows(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    make_synthetic_godb(godb_path)
    make_synthetic_input_table(input_path)

    manifest_path = tmp_path / "manifest.csv"
    output_dir = tmp_path / "results"
    pd.DataFrame([
        {
            "species": "SP1", "condition": "C1", "input_table": str(input_path),
            "godb": str(godb_path), "metric_col": "metric", "label_strategy": "rank_tail",
            "pct": 20, "output_dir": str(output_dir),
        },
        {
            "species": "SP2", "condition": "C1", "input_table": str(input_path),
            "godb": str(godb_path), "metric_col": "metric", "label_strategy": "rank_tail",
            "pct": 20, "output_dir": str(output_dir),
        },
    ]).to_csv(manifest_path, index=False)

    summary = run_batch(str(manifest_path))
    assert set(summary["succeeded"]) == {"SP1_C1.metric", "SP2_C1.metric"}
    assert summary["failed"] == []
    assert (output_dir / "summary.SP1_C1.metric.all.tsv").exists()
    assert (output_dir / "summary.SP2_C1.metric.all.tsv").exists()


def test_run_batch_one_bad_row_does_not_stop_the_others(tmp_path):
    godb_path = tmp_path / "test.godb.pkl"
    input_path = tmp_path / "input.tsv"
    make_synthetic_godb(godb_path)
    make_synthetic_input_table(input_path)

    manifest_path = tmp_path / "manifest.csv"
    output_dir = tmp_path / "results"
    pd.DataFrame([
        {
            "species": "BAD", "condition": "C1", "input_table": "does_not_exist.tsv",
            "godb": str(godb_path), "metric_col": "metric", "label_strategy": "rank_tail",
            "pct": 20, "output_dir": str(output_dir),
        },
        {
            "species": "GOOD", "condition": "C1", "input_table": str(input_path),
            "godb": str(godb_path), "metric_col": "metric", "label_strategy": "rank_tail",
            "pct": 20, "output_dir": str(output_dir),
        },
    ]).to_csv(manifest_path, index=False)

    summary = run_batch(str(manifest_path))
    assert summary["succeeded"] == ["GOOD_C1.metric"]
    assert len(summary["failed"]) == 1
    assert summary["failed"][0][0] == "BAD_C1.metric"