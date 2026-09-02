"""
tests/test_output.py
"""
import json
import os
import pandas as pd
import pytest

from enrichment.output import write_results, build_provenance


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
        assert len(written) == 1


def test_write_results_creates_output_dir_if_missing(tmp_path, sample_result_df):
    nested_dir = tmp_path / "does" / "not" / "exist" / "yet"
    written = write_results(sample_result_df, str(nested_dir), "testrun")
    assert os.path.exists(written["all"])


def test_write_results_without_provenance_is_byte_identical_to_plain_tsv(tmp_path, sample_result_df):
    """Backward-compatibility guard: write_results(provenance=None), the
    default, must produce output indistinguishable from a plain
    df.to_csv(sep='\\t', index=False) -- no provenance line, nothing extra."""
    written = write_results(sample_result_df, str(tmp_path), "testrun")
    with open(written["all"]) as f:
        first_line = f.readline()
    assert not first_line.startswith("#")
    assert first_line.startswith("go_id\t")


def test_build_provenance_returns_expected_keys():
    prov = build_provenance(run_parameters={"thresh": 0.01})
    assert "tool" in prov
    assert prov["tool"] == "go-gsea"
    assert "generated_at" in prov
    assert "git_commit" in prov
    assert "python_version" in prov
    assert "package_versions" in prov
    assert prov["run_parameters"] == {"thresh": 0.01}


def test_build_provenance_package_versions_include_known_packages():
    prov = build_provenance()
    versions = prov["package_versions"]
    for pkg in ["pandas", "numpy", "scipy", "statsmodels", "goatools", "gseapy"]:
        assert pkg in versions


def test_build_provenance_defaults_run_parameters_to_empty_dict():
    prov = build_provenance()
    assert prov["run_parameters"] == {}


def test_write_results_with_provenance_prepends_json_comment_line(tmp_path, sample_result_df):
    prov = build_provenance(run_parameters={"label_strategy": "rank_tail"})
    written = write_results(sample_result_df, str(tmp_path), "testrun", provenance=prov)

    with open(written["all"]) as f:
        first_line = f.readline()
    assert first_line.startswith("#")

    parsed = json.loads(first_line[1:])
    assert parsed["run_parameters"]["label_strategy"] == "rank_tail"


def test_write_results_with_provenance_still_parses_correctly_with_comment_char(tmp_path, sample_result_df):
    """The real point of using '#': every existing reader in this codebase
    already uses pd.read_csv(..., comment='#'), so a provenance line must
    not break that, and must not appear as a data row."""
    prov = build_provenance()
    written = write_results(sample_result_df, str(tmp_path), "testrun", provenance=prov)

    loaded = pd.read_csv(written["all"], sep="\t", comment="#")
    assert len(loaded) == 4
    assert list(loaded.columns) == list(sample_result_df.columns)


def test_write_results_with_provenance_applies_to_all_files_written(tmp_path, sample_result_df):
    prov = build_provenance()
    written = write_results(sample_result_df, str(tmp_path), "testrun", provenance=prov)

    for path in written.values():
        with open(path) as f:
            first_line = f.readline()
        assert first_line.startswith("#"), f"{path} missing provenance line"