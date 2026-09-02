#!/usr/bin/env python3
"""
scripts/check_id_overlap.py

Formalizes a real, high-stakes manual check performed by hand for
Chlamydomonas: before trusting a candidate GO annotation source for a new
species, confirm its gene ID namespace actually overlaps with the real
gene table's IDs. This exact check caught a real strain/assembly mismatch
(UniProt-GOA's CC-503/v5.6 file versus the actual CC-4532/v6.1 reference
genome, see docs/examples/chlamydomonas.md) before it could silently
corrupt every downstream GO enrichment result.

Works against any tab-delimited annotation file, since format varies by
source: GAF files have no header and use '!' as a comment/header-line
prefix; Phytozome-style annotation_info.txt files have a real header row.
You tell it which column holds the ID and whether there's a header,
rather than the tool guessing.

Example, matching the real Chlamydomonas Phytozome check:

    python scripts/check_id_overlap.py \\
        --gene-table data/raw/CR_3D.gene_data.tsv.gz \\
        --gene-id-col gene_id \\
        --strip-id-suffix '\\.v\\d+\\.\\d+$' \\
        --annotation-file data/go_reference/CR/CR.annotation_info.txt.gz \\
        --annotation-id-col locusName \\
        --annotation-has-header

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


def load_annotation_ids(path, id_col, has_header=True, comment_char="#"):
    """
    Returns a set of unique IDs found in the given column of a tab-delimited
    annotation file, plain or gzipped.

    id_col: column NAME if has_header=True, or a 0-indexed column NUMBER
    (as an int or a string convertible to int) if has_header=False.

    For GAF files (has_header=False), pass comment_char='!' -- the whole
    header block in a real GAF file is '!'-prefixed, there is no separate
    column-name header row to parse, comment lines are simply skipped.

    For has_header=True files, the FIRST non-empty line is always treated
    as the header, even if it is itself comment-prefixed (a leading
    comment_char is stripped from it, not used to skip it) -- Phytozome's
    own annotation_info.txt ships its header as "#pacId\tlocusName\t...",
    a real, load-bearing header line that happens to start with '#', not
    a comment to discard. Confirmed against the real Chlamydomonas file:
    treating it as a comment silently shifted the "header" to the first
    real data row instead. Any comment_char-prefixed line AFTER the
    header is still skipped normally.
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
                header_line = line[1:] if comment_char and line.startswith(comment_char) else line
                header = header_line.split("\t")
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


def check_id_overlap(gene_table_path, gene_id_col, annotation_file_path,
                      annotation_id_col, strip_id_suffix=None,
                      annotation_has_header=True, annotation_comment_char="#",
                      sample_size=10):
    """
    Runs the full check. Returns a dict report, does not print anything
    itself, so it can be used programmatically (e.g. from tests) as well
    as from the CLI below.
    """
    gene_ids, n_total_rows = load_gene_table_ids(gene_table_path, gene_id_col, strip_id_suffix)
    annotation_ids = load_annotation_ids(
        annotation_file_path, annotation_id_col,
        has_header=annotation_has_header, comment_char=annotation_comment_char,
    )

    matched = gene_ids & annotation_ids
    unmatched = gene_ids - annotation_ids
    match_rate = len(matched) / len(gene_ids) if gene_ids else 0.0

    return {
        "n_gene_table_rows": n_total_rows,
        "n_gene_table_ids": len(gene_ids),
        "n_annotation_ids": len(annotation_ids),
        "n_matched": len(matched),
        "n_unmatched": len(unmatched),
        "match_rate": match_rate,
        "unmatched_sample": sorted(unmatched)[:sample_size],
    }


def print_report(report, gene_table_path, annotation_file_path):
    print(f"Gene table:      {gene_table_path}")
    print(f"  rows:          {report['n_gene_table_rows']}")
    print(f"  unique IDs:    {report['n_gene_table_ids']}")
    print(f"Annotation file: {annotation_file_path}")
    print(f"  unique IDs:    {report['n_annotation_ids']}")
    print()
    print(f"Matched:   {report['n_matched']} / {report['n_gene_table_ids']} "
          f"({report['match_rate']:.1%})")
    print(f"Unmatched: {report['n_unmatched']}")
    print()

    if report["match_rate"] < 0.05:
        print("WARNING: match rate under 5%. This looks like a namespace mismatch, "
              "not sparse coverage -- check --strip-id-suffix, --gene-id-col, and "
              "--annotation-id-col before trusting this annotation source at all. "
              "This is the exact pattern a real strain/assembly mismatch produced "
              "for Chlamydomonas (see docs/examples/chlamydomonas.md).")
    elif report["match_rate"] < 0.5:
        print("Match rate is real but modest. Confirm this is plausible sparse "
              "annotation coverage for this source (check a few unmatched IDs below "
              "by hand) rather than assuming it, before committing to this source.")
    else:
        print("Match rate looks like plausible real coverage, not a namespace mismatch.")

    if report["unmatched_sample"]:
        print()
        print(f"Sample of unmatched gene-table IDs (up to {len(report['unmatched_sample'])}):")
        for uid in report["unmatched_sample"]:
            print(f"  {uid}")


def build_parser():
    parser = argparse.ArgumentParser(
        description="Check gene-ID overlap between a real gene table and a candidate "
                     "GO annotation source, before trusting the annotation source.",
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
    parser.add_argument("--annotation-has-header", action=argparse.BooleanOptionalAction,
                         default=True,
                         help="Whether --annotation-file has a real header row. GAF "
                              "files do not (use --no-annotation-has-header).")
    parser.add_argument("--annotation-comment-char", default="#",
                         help="Lines starting with this are skipped. Use '!' for GAF files.")

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
        sample_size=args.sample_size,
    )
    print_report(report, args.gene_table, args.annotation_file)


if __name__ == "__main__":
    main()