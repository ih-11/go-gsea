"""
tests/test_output.py
"""
import os
import pandas as pd
import pytest

from enrichment.output import write_results


@pytest.fixture
def sample_result_df():
    return pd.DataFrame({
        "go_id": ["GO:0001", "GO:0002", "GO:0003", "GO:0004"],
        "class": ["High", "High", "Low", "Low"],
        "population": [10, 10, 10, 10],
        "observed": [8, 1, 8, 1],
        "expected": [5, 5, 5, 5],
        "fold_enrichment": [0.7, -1.3, 0.7, -1.3],
        "p_value": [0.001, 0.001, 0.001, 0.001],
        "q_value": [0.01, 0.01, 0.01, 0.01],
        "significance": [True, True, True, True],
    })


def test_write_results_always_writes_all_file(tmp_path, sample_result_df):
    written = write_results(sample_result_df, str(tmp_path), "testrun")
    assert os.path.exists(written["all"])
    loaded = pd.read_csv(written["all"], sep="\t")
    assert len(loaded) == 4


def test_write_results_splits_over_and_under_per_class(tmp_path, sample_result_df):
    written = write_results(sample_result_df, str(tmp_path), "testrun")
    assert "High.over" in written
    assert "High.under" in written
    assert "Low.over" in written
    assert "Low.under" in written

    high_over = pd.read_csv(written["High.over"], sep="\t")
    assert list(high_over["go_id"]) == ["GO:0001"]

    high_under = pd.read_csv(written["High.under"], sep="\t")
    assert list(high_under["go_id"]) == ["GO:0002"]


def test_write_results_skips_empty_over_under_files():
    df = pd.DataFrame({
        "go_id": ["GO:0001"], "class": ["High"], "population": [10],
        "observed": [5], "expected": [5], "fold_enrichment": [0.0],
        "p_value": [0.9], "q_value": [0.9], "significance": [False],
    })
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        written = write_results(df, tmp_dir, "testrun")
        assert "all" in written
        assert "High.over" not in written
        assert "High.under" not in written


def test_write_results_handles_empty_dataframe():
    empty_df = pd.DataFrame()
    import tempfile
    with tempfile.TemporaryDirectory() as tmp_dir:
        written = write_results(empty_df, tmp_dir, "testrun")
        assert os.path.exists(written["all"])
        assert len(written) == 1  # only "all", nothing else


def test_write_results_creates_output_dir_if_missing(tmp_path, sample_result_df):
    nested_dir = tmp_path / "does" / "not" / "exist" / "yet"
    written = write_results(sample_result_df, str(nested_dir), "testrun")
    assert os.path.exists(written["all"])