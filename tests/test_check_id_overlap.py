"""
tests/test_check_id_overlap.py
"""
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest


REPO_ROOT = Path(__file__).resolve().parent.parent
SCRIPT_PATH = REPO_ROOT / "scripts" / "check_id_overlap.py"

sys.path.insert(0, str(REPO_ROOT))
from scripts.check_id_overlap import (
    load_gene_table_ids, load_annotation_ids, check_id_overlap,
)


def run_cli(args):
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH)] + args,
        capture_output=True, text=True,
    )


def test_load_gene_table_ids_basic(tmp_path):
    path = tmp_path / "genes.tsv"
    pd.DataFrame({"gene_id": ["g1", "g2", "g3"], "val": [1, 2, 3]}).to_csv(path, sep="\t", index=False)
    ids, n_rows = load_gene_table_ids(path, "gene_id")
    assert ids == {"g1", "g2", "g3"}
    assert n_rows == 3


def test_load_gene_table_ids_applies_strip_suffix(tmp_path):
    path = tmp_path / "genes.tsv"
    pd.DataFrame({"gene_id": ["g1.v1.0", "g2.v1.0"]}).to_csv(path, sep="\t", index=False)
    ids, _ = load_gene_table_ids(path, "gene_id", strip_id_suffix=r"\.v\d+\.\d+$")
    assert ids == {"g1", "g2"}


def test_load_annotation_ids_with_header(tmp_path):
    path = tmp_path / "annot.tsv"
    path.write_text("locusName\tGO\ng1\tGO:001\ng2\t\ng3\tGO:002\n")
    ids = load_annotation_ids(path, "locusName", has_header=True)
    assert ids == {"g1", "g2", "g3"}


def test_load_annotation_ids_without_header_gaf_style(tmp_path):
    """Simulates a GAF-style file: no real header, '!' comment lines,
    ID in a specific column index."""
    path = tmp_path / "annot.gaf"
    path.write_text(
        "!gaf-version: 2.2\n"
        "!this is a comment line\n"
        "UniProtKB\tg1\tGENE1\t\tGO:0000001\n"
        "UniProtKB\tg2\tGENE2\t\tGO:0000002\n"
    )
    ids = load_annotation_ids(path, id_col=1, has_header=False, comment_char="!")
    assert ids == {"g1", "g2"}


def test_load_annotation_ids_with_comment_prefixed_header(tmp_path):
    """Real bug caught against the actual Chlamydomonas file: Phytozome's
    annotation_info.txt ships its header as '#pacId\tlocusName\t...',
    a real header that happens to start with '#'. It must be treated as
    the header (stripping the leading '#'), not skipped as a comment,
    which would silently shift the header to the first real data row."""
    path = tmp_path / "annot.tsv"
    path.write_text(
        "#pacId\tlocusName\tGO\n"
        "12345\tg1\tGO:001\n"
        "12346\tg2\tGO:002\n"
    )
    ids = load_annotation_ids(path, "locusName", has_header=True, comment_char="#")
    assert ids == {"g1", "g2"}


def test_load_annotation_ids_raises_on_missing_column(tmp_path):
    path = tmp_path / "annot.tsv"
    path.write_text("locusName\tGO\ng1\tGO:001\n")
    with pytest.raises(ValueError):
        load_annotation_ids(path, "nonexistent_col", has_header=True)


def test_check_id_overlap_full_match(tmp_path):
    gene_path = tmp_path / "genes.tsv"
    annot_path = tmp_path / "annot.tsv"
    pd.DataFrame({"gene_id": ["g1", "g2", "g3"]}).to_csv(gene_path, sep="\t", index=False)
    annot_path.write_text("locusName\tGO\ng1\tGO:001\ng2\tGO:002\ng3\tGO:003\n")

    report = check_id_overlap(
        gene_path, "gene_id", annot_path, "locusName",
    )
    assert report["match_rate"] == 1.0
    assert report["n_matched"] == 3
    assert report["unmatched_sample"] == []


def test_check_id_overlap_partial_match_and_sample(tmp_path):
    gene_path = tmp_path / "genes.tsv"
    annot_path = tmp_path / "annot.tsv"
    pd.DataFrame({"gene_id": ["g1", "g2", "g3", "g4"]}).to_csv(gene_path, sep="\t", index=False)
    annot_path.write_text("locusName\tGO\ng1\tGO:001\ng2\tGO:002\n")  # g3, g4 missing

    report = check_id_overlap(
        gene_path, "gene_id", annot_path, "locusName",
    )
    assert report["n_matched"] == 2
    assert report["match_rate"] == 0.5
    assert set(report["unmatched_sample"]) == {"g3", "g4"}


def test_check_id_overlap_zero_match_simulates_namespace_mismatch(tmp_path):
    """Simulates the real CC-503/CC-4532 case: completely different ID
    namespaces, near-zero overlap."""
    gene_path = tmp_path / "genes.tsv"
    annot_path = tmp_path / "annot.tsv"
    pd.DataFrame({"gene_id": ["Cre01.g000350_4532"]}).to_csv(gene_path, sep="\t", index=False)
    annot_path.write_text("locusName\tGO\nCre01.g000350\tGO:001\n")  # different suffix convention

    report = check_id_overlap(gene_path, "gene_id", annot_path, "locusName")
    assert report["match_rate"] == 0.0


def test_cli_end_to_end_high_match_rate(tmp_path):
    gene_path = tmp_path / "genes.tsv"
    annot_path = tmp_path / "annot.tsv"
    pd.DataFrame({"gene_id": [f"g{i}.v1.0" for i in range(10)]}).to_csv(gene_path, sep="\t", index=False)
    annot_lines = "locusName\tGO\n" + "".join(f"g{i}\tGO:{i:03d}\n" for i in range(10))
    annot_path.write_text(annot_lines)

    result = run_cli([
        "--gene-table", str(gene_path),
        "--gene-id-col", "gene_id",
        "--strip-id-suffix", r"\.v\d+\.\d+$",
        "--annotation-file", str(annot_path),
        "--annotation-id-col", "locusName",
    ])
    assert result.returncode == 0, result.stderr
    assert "Matched:   10 / 10 (100.0%)" in result.stdout
    assert "plausible real coverage" in result.stdout


def test_cli_low_match_rate_prints_warning(tmp_path):
    gene_path = tmp_path / "genes.tsv"
    annot_path = tmp_path / "annot.tsv"
    pd.DataFrame({"gene_id": [f"gene{i}" for i in range(20)]}).to_csv(gene_path, sep="\t", index=False)
    annot_path.write_text("locusName\tGO\nother_namespace_id\tGO:001\n")

    result = run_cli([
        "--gene-table", str(gene_path),
        "--gene-id-col", "gene_id",
        "--annotation-file", str(annot_path),
        "--annotation-id-col", "locusName",
    ])
    assert result.returncode == 0, result.stderr
    assert "WARNING" in result.stdout
    assert "namespace mismatch" in result.stdout


def test_cli_gaf_style_no_header(tmp_path):
    gene_path = tmp_path / "genes.tsv"
    annot_path = tmp_path / "annot.gaf"
    pd.DataFrame({"gene_id": ["g1", "g2"]}).to_csv(gene_path, sep="\t", index=False)
    annot_path.write_text(
        "!gaf-version: 2.2\n"
        "UniProtKB\tg1\tGENE1\t\tGO:0000001\n"
        "UniProtKB\tg2\tGENE2\t\tGO:0000002\n"
    )

    result = run_cli([
        "--gene-table", str(gene_path),
        "--gene-id-col", "gene_id",
        "--annotation-file", str(annot_path),
        "--annotation-id-col", "1",
        "--no-annotation-has-header",
        "--annotation-comment-char", "!",
    ])
    assert result.returncode == 0, result.stderr
    assert "Matched:   2 / 2 (100.0%)" in result.stdout