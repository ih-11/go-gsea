#!/usr/bin/env python3
"""
scripts/run_pipeline.py

CLI entry point tying reference/, filters/, labelers/, and enrichment/
together. Wraps the same steps every real run so far has done by hand in
an interactive script.

Example, matching the real Chlamydomonas PR_gene run:

    python scripts/run_pipeline.py \\
        --input-table /path/to/gene_data.tsv.gz \\
        --godb data/go_reference/chlamy.godb.pkl \\
        --metric-col PR_gene \\
        --id-col gene_id \\
        --strip-id-suffix '\\.v\\d+\\.\\d+$' \\
        --exclude-id 'gene:Standard-R-luc' \\
        --label-strategy rank_tail --pct 10 \\
        --output-dir results/GO \\
        --dataset-name CR_3D.PR_gene \\
        --unknown-ratio-thresh 0.9 \\
        --thresh-type p --thresh 0.01

This script knows nothing about Chlamydomonas or PR specifically -- every
species/metric-specific detail above is a command-line argument, not a
hardcoded assumption.
"""

import sys
import os

# Ensure the repo root (parent of this scripts/ folder) is importable,
# regardless of the current working directory this script is invoked from.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import argparse
import pickle
import re

import pandas as pd

from labelers.labelers import rank_tail, explicit_threshold, boolean_flag, cluster
from enrichment.ora import run_ora
from enrichment.output import write_results


LABEL_STRATEGIES = {
    "rank_tail": rank_tail,
    "explicit_threshold": explicit_threshold,
    "boolean_flag": boolean_flag,
    "cluster": cluster,
}


def build_parser():
    parser = argparse.ArgumentParser(
        description="Run GO over-representation analysis end to end.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input-table", required=True,
                         help="Path to a gene-level TSV (plain or .gz). See README section 2.")
    parser.add_argument("--godb", required=True,
                         help="Path to a cached .godb pickle from reference/build_godb.py.")
    parser.add_argument("--metric-col", required=True,
                         help="Column in the input table to label genes by.")
    parser.add_argument("--id-col", default="gene_id",
                         help="Gene ID column name in the input table.")
    parser.add_argument("--strip-id-suffix", default=None,
                         help="Regex to strip from gene IDs before matching against the godb "
                              "(e.g. a trailing version suffix like '.v6.1'). Applied BEFORE "
                              "--exclude-id, see that flag's help for what this means in practice.")
    parser.add_argument("--exclude-id", action="append", default=[],
                         help="Gene ID(s) to exclude before labeling (e.g. a spike-in control). "
                              "Applied AFTER --strip-id-suffix -- if both flags are used, give "
                              "the already-stripped form of the ID here, not the raw one, or the "
                              "exclusion will silently fail to match. Can be given multiple times.")

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


def load_input_table(path, id_col, strip_id_suffix, exclude_ids, metric_col):
    """
    Order of operations matters here: --strip-id-suffix is applied first,
    --exclude-id matching happens second, against the already-stripped
    IDs. This is documented in --exclude-id's help text above -- if a raw
    (unstripped) ID is passed to --exclude-id while --strip-id-suffix is
    also given, the exclusion will silently fail to match anything, since
    the comparison happens after stripping, not before.
    """
    df = pd.read_csv(path, sep="\t", comment="#")
    if strip_id_suffix:
        df[id_col] = df[id_col].apply(lambda g: re.sub(strip_id_suffix, "", g))
    if exclude_ids:
        df = df[~df[id_col].isin(exclude_ids)]
    df = df.dropna(subset=[metric_col])
    return df


def apply_label_strategy(df, args):
    fn = LABEL_STRATEGIES[args.label_strategy]
    if args.label_strategy == "rank_tail":
        return fn(df, col=args.metric_col, pct=args.pct)
    elif args.label_strategy == "explicit_threshold":
        return fn(df, col=args.metric_col, high_thresh=args.high_thresh, low_thresh=args.low_thresh)
    elif args.label_strategy == "boolean_flag":
        return fn(df, bool_col=args.metric_col)
    elif args.label_strategy == "cluster":
        return fn(df, cols=[args.metric_col], n_clusters=args.n_clusters)
    raise ValueError(f"Unhandled label strategy: {args.label_strategy}")


def main():
    args = build_parser().parse_args()

    with open(args.godb, "rb") as f:
        godb = pickle.load(f)
    gene_go = godb["gene_go"]

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