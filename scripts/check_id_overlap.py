#!/usr/bin/env python3
"""
scripts/check_id_overlap.py

Formalizes two real, high-stakes manual checks performed by hand for
Chlamydomonas, both bundled into one report here:

1. NAMESPACE ALIGNMENT: does the candidate annotation file's gene IDs
   actually overlap with your real gene table's IDs at all? This exact
   check caught a real strain/assembly mismatch (UniProt-GOA's
   CC-503/v5.6 file versus the actual CC-4532/v6.1 reference genome, see
   docs/examples/chlamydomonas.md) before it could silently corrupt every
   downstream GO enrichment result.

2. GO COVERAGE (optional, via --annotation-go-col): of the genes that DO
   match by ID, how many actually carry usable GO annotation? A file can
   score 100% on namespace alignment (nearly every gene in the genome is
   listed) while only a small fraction of those rows have real GO terms
   -- this is a DIFFERENT number from namespace alignment, and reporting
   only the first one risks false confidence. For Chlamydomonas, real
   ID alignment was ~100%, real GO coverage was ~18%; conflating the two
   would have been a genuine mistake.

Works against any tab-delimited annotation file, since format varies by
source: GAF files have no header and use '!' as a comment/header-line
prefix; Phytozome-style annotation_info.txt files have a header row that
is itself '#'-prefixed (a real header, not a comment to skip, confirmed
against the actual Chlamydomonas file, see load_annotation_ids's
docstring).

Example, matching the real Chlamydomonas Phytozome check, both numbers
at once:

    python scripts/check_id_overlap.py \\
        --gene-table data/raw/CR_3D.gene_data.tsv.gz \\
        --gene-id-col gene_id \\
        --strip-id-suffix '\\.v\\d+\\.\\d+$' \\
        --annotation-file data/go_reference/CR/CR.annotation_info.txt.gz \\
        --annotation-id-col locusName \\
        --annotation-go-col GO

Example against a GAF file (no header, ID typically in column 2, 0-indexed
column 1 -- confirm against the real GAF spec for the file in hand):

    python scripts/check_id_overlap.py \\
        --gene-table data/raw/AT_HS.gene_data.tsv.gz \\
        --gene-id-col gene_id \\
        --annotation-file some_species.gaf.gz \\
        --annotation-id-col 1 \\
        --no-annotation-has-header \\
        --annotation-comment-char '!'
"""
import sys
import os

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import argparse
import gzip
import re

import pandas as pd


def _opener(path):
    return gzip.open if str(path).endswith(".gz") else open


def load_gene_table_ids(path, id_col, strip_id_suffix=None):
    """
    Returns (set of gene IDs, total row count) from a gene-level table,
    matching the same loading convention scripts/run_pipeline.py uses.
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    if id_col not in df.columns:
        raise ValueError(f"'{id_col}' not found in {path}. Columns present: {list(df.columns)}")
    ids = df[id_col].astype(str)
    if strip_id_suffix:
        ids = ids.apply(lambda g: re.sub(strip_id_suffix, "", g))
    return set(ids), len(df)


def _resolve_header(header_line, comment_char):
    """
    The FIRST non-empty line of a has_header=True file is always treated
    as the header, even if it is itself comment-prefixed (a leading
    comment_char is stripped, not used to skip the line) -- Phytozome's
    own annotation_info.txt ships its header as "#pacId\\tlocusName\\t...",
    a real, load-bearing header line that happens to start with '#', not
    a comment to discard. Confirmed against the real Chlamydomonas file:
    treating it as a comment silently shifted the "header" to the first
    real data row instead.
    """
    if comment_char and header_line.startswith(comment_char):
        return header_line[1:].split("\t")
    return header_line.split("\t")


def load_annotation_ids(path, id_col, has_header=True, comment_char="#"):
    """
    Returns a set of unique IDs found in the given column of a tab-delimited
    annotation file, plain or gzipped. Namespace-alignment check only, use
    load_annotation_ids_with_go_status() for the coverage check too.

    id_col: column NAME if has_header=True, or a 0-indexed column NUMBER
    (as an int or a string convertible to int) if has_header=False.
    """
    ids = set()
    header = None
    idx = None if has_header else int(id_col)

    with _opener(path)(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            if has_header and header is None:
                header = _resolve_header(line, comment_char)
                if id_col not in header:
                    raise ValueError(
                        f"'{id_col}' not found in {path}'s header. Columns present: {header}"
                    )
                idx = header.index(id_col)
                continue

            if comment_char and line.startswith(comment_char):
                continue

            fields = line.split("\t")
            if idx is not None and idx < len(fields):
                ids.add(fields[idx])

    return ids


def load_annotation_ids_with_go_status(path, id_col, go_col, has_header=True, comment_char="#"):
    """
    Returns dict {id: has_go} for every ID found in id_col, where has_go
    is True if that ID has a non-empty value in go_col on ANY of its rows
    (a locus can appear on multiple rows, e.g. one row per transcript --
    if any row carries GO data, the gene counts as annotated). Requires
    has_header=True, since go_col is looked up by name.
    """
    if not has_header:
        raise ValueError("load_annotation_ids_with_go_status requires has_header=True, "
                          "go_col is looked up by name.")

    status = {}
    header = None
    id_idx = None
    go_idx = None

    with _opener(path)(path, "rt") as f:
        for line in f:
            line = line.rstrip("\n")
            if not line:
                continue

            if header is None:
                header = _resolve_header(line, comment_char)
                for needed in (id_col, go_col):
                    if needed not in header:
                        raise ValueError(
                            f"'{needed}' not found in {path}'s header. Columns present: {header}"
                        )
                id_idx = header.index(id_col)
                go_idx = header.index(go_col)
                continue

            if comment_char and line.startswith(comment_char):
                continue

            fields = line.split("\t")
            if id_idx >= len(fields):
                continue
            gene_id = fields[id_idx]
            has_go_here = go_idx < len(fields) and fields[go_idx].strip() != ""
            status[gene_id] = status.get(gene_id, False) or has_go_here

    return status


def check_id_overlap(gene_table_path, gene_id_col, annotation_file_path,
                      annotation_id_col, strip_id_suffix=None,
                      annotation_has_header=True, annotation_comment_char="#",
                      annotation_go_col=None, sample_size=10):
    """
    Runs the full check. Returns a dict report, does not print anything
    itself, so it can be used programmatically (e.g. from tests) as well
    as from the CLI below.

    If annotation_go_col is given, the report ALSO includes GO-coverage
    numbers (report['go_coverage_*']), answering a genuinely different
    question from namespace alignment alone, see this module's docstring.
    """
    gene_ids, n_total_rows = load_gene_table_ids(gene_table_path, gene_id_col, strip_id_suffix)

    if annotation_go_col:
        go_status = load_annotation_ids_with_go_status(
            annotation_file_path, annotation_id_col, annotation_go_col,
            has_header=annotation_has_header, comment_char=annotation_comment_char,
        )
        annotation_ids = set(go_status.keys())
    else:
        go_status = None
        annotation_ids = load_annotation_ids(
            annotation_file_path, annotation_id_col,
            has_header=annotation_has_header, comment_char=annotation_comment_char,
        )

    matched = gene_ids & annotation_ids
    unmatched = gene_ids - annotation_ids
    match_rate = len(matched) / len(gene_ids) if gene_ids else 0.0

    report = {
        "n_gene_table_rows": n_total_rows,
        "n_gene_table_ids": len(gene_ids),
        "n_annotation_ids": len(annotation_ids),
        "n_matched": len(matched),
        "n_unmatched": len(unmatched),
        "match_rate": match_rate,
        "unmatched_sample": sorted(unmatched)[:sample_size],
    }

    if go_status is not None:
        matched_with_go = {gid for gid in matched if go_status.get(gid, False)}
        go_coverage_rate = len(matched_with_go) / len(gene_ids) if gene_ids else 0.0
        report["go_coverage_n"] = len(matched_with_go)
        report["go_coverage_rate"] = go_coverage_rate

    return report


def print_report(report, gene_table_path, annotation_file_path):
    print(f"Gene table:      {gene_table_path}")
    print(f"  rows:          {report['n_gene_table_rows']}")
    print(f"  unique IDs:    {report['n_gene_table_ids']}")
    print(f"Annotation file: {annotation_file_path}")
    print(f"  unique IDs:    {report['n_annotation_ids']}")
    print()
    print("--- Question 1: does the ID namespace line up at all? ---")
    print(f"Namespace match: {report['n_matched']} / {report['n_gene_table_ids']} "
          f"({report['match_rate']:.1%})")
    print(f"Unmatched:       {report['n_unmatched']}")
    print()

    if report["match_rate"] < 0.05:
        print("WARNING: namespace match under 5%. This looks like a namespace mismatch, "
              "not sparse coverage -- check --strip-id-suffix, --gene-id-col, and "
              "--annotation-id-col before trusting this annotation source at all. "
              "This is the exact pattern a real strain/assembly mismatch produced "
              "for Chlamydomonas (see docs/examples/chlamydomonas.md).")
    else:
        print("Namespace match looks like real ID alignment, not a mismatch. "
              "This does NOT yet tell you how much usable GO coverage exists -- "
              "see Question 2 below if --annotation-go-col was given, or run this "
              "tool again with it if not.")

    if report["unmatched_sample"]:
        print()
        print(f"Sample of unmatched gene-table IDs (up to {len(report['unmatched_sample'])}):")
        for uid in report["unmatched_sample"]:
            print(f"  {uid}")

    if "go_coverage_rate" in report:
        print()
        print("--- Question 2: of the genes that match, how many actually carry ---")
        print("--- usable GO annotation? (a DIFFERENT number from namespace match) ---")
        print(f"GO coverage: {report['go_coverage_n']} / {report['n_gene_table_ids']} "
              f"({report['go_coverage_rate']:.1%})")
        print()
        print("This is the number that predicts how many genes reference/build_godb.py "
              "will actually annotate, and therefore how large your real enrichment "
              "population will be -- namespace match alone (Question 1) does not tell "
              "you this. A high namespace match with low GO coverage (as with the real "
              "Chlamydomonas/Phytozome case: ~100% namespace match, ~18% GO coverage) "
              "is a real, expected tradeoff for some annotation sources, not a bug.")
    elif report["match_rate"] >= 0.05:
        print()
        print("Note: GO coverage was not checked (no --annotation-go-col given). "
              "A high namespace match here does not by itself mean most genes will "
              "end up annotated once reference/build_godb.py runs -- rerun with "
              "--annotation-go-col to get that number too before committing to this "
              "source.")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Check gene-ID namespace alignment, and optionally GO coverage, "
                     "between a real gene table and a candidate GO annotation source.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--gene-table", required=True,
                         help="Path to a real gene-level TSV (plain or .gz).")
    parser.add_argument("--gene-id-col", default="gene_id")
    parser.add_argument("--strip-id-suffix", default=None,
                         help="Regex to strip from gene-table IDs before comparing "
                              "(e.g. a trailing version suffix).")

    parser.add_argument("--annotation-file", required=True,
                         help="Path to the candidate GO annotation file (GAF, "
                              "Phytozome-style, or any tab-delimited file).")
    parser.add_argument("--annotation-id-col", required=True,
                         help="Column NAME (if the file has a header) or 0-indexed "
                              "column NUMBER (if it does not) holding the gene/locus ID.")
    parser.add_argument("--annotation-go-col", default=None,
                         help="Optional column NAME holding GO term data. If given, "
                              "also reports what fraction of matched genes actually "
                              "carry non-empty GO data, a different question from "
                              "namespace alignment alone. Requires --annotation-has-header.")
    parser.add_argument("--annotation-has-header", action=argparse.BooleanOptionalAction,
                         default=True,
                         help="Whether --annotation-file has a real header row. GAF "
                              "files do not (use --no-annotation-has-header).")
    parser.add_argument("--annotation-comment-char", default="#",
                         help="Lines starting with this are skipped (after the header, "
                              "if any). Use '!' for GAF files.")

    parser.add_argument("--sample-size", type=int, default=10,
                         help="How many unmatched IDs to print.")

    return parser


def main():
    args = build_parser().parse_args()
    report = check_id_overlap(
        gene_table_path=args.gene_table,
        gene_id_col=args.gene_id_col,
        annotation_file_path=args.annotation_file,
        annotation_id_col=args.annotation_id_col,
        strip_id_suffix=args.strip_id_suffix,
        annotation_has_header=args.annotation_has_header,
        annotation_comment_char=args.annotation_comment_char,
        annotation_go_col=args.annotation_go_col,
        sample_size=args.sample_size,
    )
    print_report(report, args.gene_table, args.annotation_file)


if __name__ == "__main__":
    main()