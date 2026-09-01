#!/usr/bin/env python3
"""
scripts/run_pipeline.py

CLI entry point tying reference/, filters/, labelers/, dataprep/, and
enrichment/ together. Wraps the same steps every real run so far has
done by hand in an interactive script.

Single-source example, matching the real Chlamydomonas PR_gene run:

    python scripts/run_pipeline.py \\
        --input-table /path/to/gene_data.tsv.gz \\
        --godb data/go_reference/CR/CR.godb.pkl \\
        --metric-col PR_gene \\
        --id-col gene_id \\
        --strip-id-suffix '\\.v\\d+\\.\\d+$' \\
        --exclude-id 'gene:Standard-R-luc' \\
        --label-strategy rank_tail --pct 10 \\
        --output-dir results/GO \\
        --dataset-name CR_3D.PR_gene \\
        --unknown-ratio-thresh 0.9 \\
        --thresh-type p --thresh 0.01

Multi-source (merged) example, for real multi-condition clustering, e.g.
PR_gene measured at two growth stages combined into one wide table
before clustering (see dataprep/merge_tables.py):

    python scripts/run_pipeline.py \\
        --merge-manifest merge_pr_gene.txt \\
        --godb data/go_reference/CR/CR.godb.pkl \\
        --metric-col cond_3D cond_6D \\
        --label-strategy cluster --n-clusters 4 \\
        --output-dir results/GO \\
        --dataset-name CR_3D-CR_6D.PR_gene_cluster

Exactly one of --input-table / --merge-manifest must be given.

--metric-col accepts one or more values. rank_tail, explicit_threshold,
and boolean_flag each require EXACTLY ONE (that is inherent to what
those strategies mean). cluster accepts ONE OR MORE, since
labelers.cluster() has supported multi-column input since it was first
written.

This script knows nothing about Chlamydomonas, PR, or any specific
species/condition -- every detail above is a command-line argument or a
manifest file, not a hardcoded assumption.
"""

import sys
import os

# Ensure the repo root (parent of this scripts/ folder) is importable,
# regardless of the current working directory this script is invoked from.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import argparse
import configparser
import pickle
import re

import pandas as pd

from labelers.labelers import rank_tail, explicit_threshold, boolean_flag, cluster
from enrichment.ora import run_ora
from enrichment.output import write_results
from dataprep.merge_tables import merge_gene_tables


LABEL_STRATEGIES = {
    "rank_tail": rank_tail,
    "explicit_threshold": explicit_threshold,
    "boolean_flag": boolean_flag,
    "cluster": cluster,
}

# Strategies that operate on exactly one column vs. strategies that can
# take one or more. Single-column strategies fail loudly if given more
# than one --metric-col value, rather than silently using only the first.
SINGLE_COLUMN_STRATEGIES = {"rank_tail", "explicit_threshold", "boolean_flag"}
MULTI_COLUMN_STRATEGIES = {"cluster"}


def build_parser():
    parser = argparse.ArgumentParser(
        description="Run GO over-representation analysis end to end.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--input-table",
                              help="Path to a single gene-level TSV (plain or .gz). "
                                   "Mutually exclusive with --merge-manifest.")
    input_group.add_argument("--merge-manifest",
                              help="Path to a plain-text, INI-style manifest listing "
                                   "multiple source files to merge into one wide table "
                                   "before labeling (see dataprep/merge_tables.py and this "
                                   "script's module docstring for the manifest format). "
                                   "Mutually exclusive with --input-table.")

    parser.add_argument("--godb", required=True,
                         help="Path to a cached .godb pickle from reference/build_godb.py.")
    parser.add_argument("--metric-col", required=True, nargs="+",
                         help="Column(s) to label genes by. With --input-table, these are "
                              "columns already present in that file. With --merge-manifest, "
                              "these should be the manifest's [section] names, which become "
                              "the merged table's column names. "
                              "rank_tail/explicit_threshold/boolean_flag require exactly one. "
                              "cluster accepts one or more.")
    parser.add_argument("--id-col", default="gene_id",
                         help="Gene ID column name in the input table(s).")
    parser.add_argument("--strip-id-suffix", default=None,
                         help="Regex to strip from gene IDs before matching against the godb "
                              "(e.g. a trailing version suffix like '.v6.1'). Applied BEFORE "
                              "--exclude-id, see that flag's help for what this means in practice. "
                              "Applied per-source when --merge-manifest is used, before merging.")
    parser.add_argument("--exclude-id", action="append", default=[],
                         help="Gene ID(s) to exclude before labeling (e.g. a spike-in control). "
                              "Applied AFTER --strip-id-suffix -- if both flags are used, give "
                              "the already-stripped form of the ID here, not the raw one, or the "
                              "exclusion will silently fail to match. Can be given multiple times. "
                              "Applied after merging when --merge-manifest is used.")

    parser.add_argument("--label-strategy", required=True, choices=list(LABEL_STRATEGIES.keys()))
    parser.add_argument("--pct", type=float, default=10,
                         help="Used by rank_tail: percent for each tail.")
    parser.add_argument("--high-thresh", type=float, default=None,
                         help="Used by explicit_threshold.")
    parser.add_argument("--low-thresh", type=float, default=None,
                         help="Used by explicit_threshold.")
    parser.add_argument("--n-clusters", type=int, default=None,
                         help="Used by cluster.")

    parser.add_argument("--output-dir", required=True)
    parser.add_argument("--dataset-name", required=True,
                         help="Used to name output files, e.g. summary.<dataset_name>.all.tsv")

    parser.add_argument("--unknown-ratio-thresh", type=float, default=0.9)
    parser.add_argument("--thresh-type", choices=["p", "q"], default="p")
    parser.add_argument("--thresh", type=float, default=0.01)
    parser.add_argument("--slim-godb", default=None,
                         help="Optional path to a cached slim .godb (from build_and_cache_slim_godb). "
                              "If given, also runs enrichment against it and writes to --slim-output-dir.")
    parser.add_argument("--slim-output-dir", default=None,
                         help="Required if --slim-godb is given.")

    return parser


def validate_metric_cols(metric_cols, label_strategy):
    """
    Raises ValueError if metric_cols' length doesn't match what
    label_strategy actually supports. Kept as its own function so this
    check is easy to unit test in isolation from argparse/CLI plumbing.
    """
    if label_strategy in SINGLE_COLUMN_STRATEGIES and len(metric_cols) != 1:
        raise ValueError(
            f"--label-strategy {label_strategy} requires exactly one --metric-col value, "
            f"got {len(metric_cols)}: {metric_cols}"
        )
    if label_strategy in MULTI_COLUMN_STRATEGIES and len(metric_cols) < 1:
        raise ValueError(
            f"--label-strategy {label_strategy} requires at least one --metric-col value"
        )


def read_merge_manifest(manifest_path):
    """
    Reads a plain-text, INI-style merge manifest. Each [section] becomes
    one entry in the sources list passed to dataprep.merge_gene_tables():
    (path, value_col, section_name) -- the section name IS the output
    column name, so --metric-col values should match the manifest's
    section names directly.

    Required keys per section: path, value_col

    Example:
        [cond_3D]
        path = data/raw/CR_3D.gene_data.tsv.gz
        value_col = PR_gene

        [cond_6D]
        path = data/raw/CR_6D.gene_data.tsv.gz
        value_col = PR_gene
    """
    parser = configparser.ConfigParser()
    parser.optionxform = str
    read_files = parser.read(manifest_path)
    if not read_files:
        raise FileNotFoundError(f"Merge manifest not found or empty: {manifest_path}")

    sources = []
    for section in parser.sections():
        section_dict = dict(parser[section])
        missing = [k for k in ("path", "value_col") if k not in section_dict]
        if missing:
            raise ValueError(f"Merge manifest section '{section}' missing key(s): {missing}")
        sources.append((section_dict["path"], section_dict["value_col"], section))

    return sources


def load_input_table(path, id_col, strip_id_suffix, exclude_ids, metric_cols):
    """
    Single-source load path. Order of operations: --strip-id-suffix is
    applied first, --exclude-id matching happens second, against the
    already-stripped IDs -- if a raw (unstripped) ID is passed to
    --exclude-id while --strip-id-suffix is also given, the exclusion
    will silently fail to match anything, since the comparison happens
    after stripping, not before.

    dropna is applied across ALL given metric_cols -- a row missing any
    one of them is dropped, matching merge_gene_tables()'s own inner-join
    philosophy (a gene needs a valid value in every column being used).
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    if strip_id_suffix:
        df[id_col] = df[id_col].apply(lambda g: re.sub(strip_id_suffix, "", g))
    if exclude_ids:
        df = df[~df[id_col].isin(exclude_ids)]
    df = df.dropna(subset=metric_cols)
    return df


def load_merged_table(manifest_path, id_col, strip_id_suffix, exclude_ids):
    """
    Multi-source load path. Reads the merge manifest, calls
    merge_gene_tables() to inner-join every listed source into one wide
    table, THEN applies --strip-id-suffix / --exclude-id to the merged
    result -- same order of operations as the single-source path, just
    after the merge instead of before, since the merge itself needs raw
    (unstripped) IDs to match correctly across source files that may all
    use the same ID convention consistently.
    """
    sources = read_merge_manifest(manifest_path)
    df = merge_gene_tables(sources, id_col=id_col, how="inner")
    if strip_id_suffix:
        df[id_col] = df[id_col].apply(lambda g: re.sub(strip_id_suffix, "", g))
    if exclude_ids:
        df = df[~df[id_col].isin(exclude_ids)]
    return df


def apply_label_strategy(df, args):
    fn = LABEL_STRATEGIES[args.label_strategy]
    if args.label_strategy == "rank_tail":
        return fn(df, col=args.metric_col[0], pct=args.pct)
    elif args.label_strategy == "explicit_threshold":
        return fn(df, col=args.metric_col[0], high_thresh=args.high_thresh, low_thresh=args.low_thresh)
    elif args.label_strategy == "boolean_flag":
        return fn(df, bool_col=args.metric_col[0])
    elif args.label_strategy == "cluster":
        return fn(df, cols=args.metric_col, n_clusters=args.n_clusters)
    raise ValueError(f"Unhandled label strategy: {args.label_strategy}")


def main():
    args = build_parser().parse_args()
    validate_metric_cols(args.metric_col, args.label_strategy)

    with open(args.godb, "rb") as f:
        godb = pickle.load(f)
    gene_go = godb["gene_go"]

    if args.merge_manifest:
        df = load_merged_table(args.merge_manifest, args.id_col,
                                args.strip_id_suffix, args.exclude_id)
        df = df.dropna(subset=args.metric_col)
        print(f"Loaded {len(df)} genes (merged from {args.merge_manifest})")
    else:
        df = load_input_table(args.input_table, args.id_col, args.strip_id_suffix,
                               args.exclude_id, args.metric_col)
        print(f"Loaded {len(df)} genes from {args.input_table}")

    labeled_df = apply_label_strategy(df, args)
    labeled_df = labeled_df.rename(columns={args.id_col: "gene_id"})
    for cls in labeled_df["class"].dropna().unique():
        n = (labeled_df["class"] == cls).sum()
        print(f"  class '{cls}': {n} genes")

    result = run_ora(
        gene_go, labeled_df,
        unknown_ratio_thresh=args.unknown_ratio_thresh,
        thresh_type=args.thresh_type, thresh=args.thresh,
    )
    print(f"Tested {len(result)} (class, GO term) rows, "
          f"{int(result['significance'].sum()) if len(result) else 0} significant")

    written = write_results(result, args.output_dir, args.dataset_name)

    if args.slim_godb:
        if not args.slim_output_dir:
            raise ValueError("--slim-output-dir is required when --slim-godb is given")
        with open(args.slim_godb, "rb") as f:
            slim_godb = pickle.load(f)
        slim_result = run_ora(
            slim_godb["gene_go"], labeled_df,
            unknown_ratio_thresh=args.unknown_ratio_thresh,
            thresh_type=args.thresh_type, thresh=args.thresh,
        )
        print(f"[GO-slim] Tested {len(slim_result)} rows, "
              f"{int(slim_result['significance'].sum()) if len(slim_result) else 0} significant")
        slim_written = write_results(slim_result, args.slim_output_dir, args.dataset_name)
        print("[GO-slim] Files written:")
        for key, path in slim_written.items():
            print(f"  {key}: {path}")

    print("Files written:")
    for key, path in written.items():
        print(f"  {key}: {path}")


if __name__ == "__main__":
    main()