#!/usr/bin/env python3
"""
scripts/run_pipeline.py

CLI entry point tying reference/, dataprep/, filters/, labelers/, and
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

Multi-source (merged) example, for real multi-condition clustering:

    python scripts/run_pipeline.py \\
        --merge-manifest merge_pr_gene.txt \\
        --godb data/go_reference/CR/CR.godb.pkl \\
        --metric-col cond_3D cond_6D \\
        --label-strategy cluster --n-clusters 4 \\
        --output-dir results/GO \\
        --dataset-name CR_3D-CR_6D.PR_gene_cluster

Exactly one of --input-table / --merge-manifest must be given.

--metric-col accepts one or more values. rank_tail, explicit_threshold,
and boolean_flag each require EXACTLY ONE. cluster accepts ONE OR MORE.

Optional Stage A population-eligibility filters (filters/population.py)
can be applied before labeling, in this fixed order: read-depth, then
usage, then minimum group size, then --drop-uninformative. All are off
by default -- omitting every --read-depth-col/--usage-col/--min-group-col/
--drop-uninformative flag reproduces the exact prior behavior of this
script unchanged:

    --read-depth-col T_n_reads --read-depth-thresh 10 \\
    --usage-col rTrans-usage.b --usage-thresh 0.05

--drop-uninformative is a special case: use it ONLY with
--label-strategy cluster, to avoid a real distance-clustering artifact
where genes tied at a shared floor value (e.g. all zero) mechanically
merge into a spuriously "perfect" cluster. It does not mean those genes
lack real signal -- pair it with a companion boolean_flag or
explicit_threshold run on the same population to test them on their own
terms, see filters/population.py's docstring for the full reasoning.

This script knows nothing about Chlamydomonas, PR, or any specific
species/condition -- every detail above is a command-line argument or a
manifest file, not a hardcoded assumption.
"""

import sys
import os

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
from filters.population import (
    chain_filters, read_depth_filter, usage_filter, min_group_size_filter,
    drop_uninformative_for_clustering,
)


LABEL_STRATEGIES = {
    "rank_tail": rank_tail,
    "explicit_threshold": explicit_threshold,
    "boolean_flag": boolean_flag,
    "cluster": cluster,
}

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
                                   "before labeling. Mutually exclusive with --input-table.")

    parser.add_argument("--godb", required=True,
                         help="Path to a cached .godb pickle from reference/build_godb.py.")
    parser.add_argument("--metric-col", required=True, nargs="+",
                         help="Column(s) to label genes by. rank_tail/explicit_threshold/"
                              "boolean_flag require exactly one. cluster accepts one or more.")
    parser.add_argument("--id-col", default="gene_id",
                         help="Gene ID column name in the input table(s).")
    parser.add_argument("--strip-id-suffix", default=None,
                         help="Regex to strip from gene IDs before matching against the godb. "
                              "Applied BEFORE --exclude-id.")
    parser.add_argument("--exclude-id", action="append", default=[],
                         help="Gene ID(s) to exclude before labeling. Applied AFTER "
                              "--strip-id-suffix -- give the already-stripped form if both "
                              "are used. Can be given multiple times.")

    # --- Stage A: population-eligibility filters, all optional, off by default ---
    parser.add_argument("--read-depth-col", default=None,
                         help="Column to apply filters.population.read_depth_filter on. "
                              "Requires --read-depth-thresh.")
    parser.add_argument("--read-depth-thresh", type=float, default=None,
                         help="Minimum value (inclusive) for --read-depth-col.")
    parser.add_argument("--usage-col", default=None,
                         help="Column to apply filters.population.usage_filter on. "
                              "Requires --usage-thresh.")
    parser.add_argument("--usage-thresh", type=float, default=None,
                         help="Minimum value (inclusive) for --usage-col.")
    parser.add_argument("--min-group-col", default=None,
                         help="Group-by column to apply filters.population.min_group_size_filter "
                              "on (e.g. a transcript ID, to require a minimum number of surviving "
                              "sites per transcript). Requires --min-group-n.")
    parser.add_argument("--min-group-n", type=int, default=None,
                         help="Minimum group size (inclusive) for --min-group-col.")
    parser.add_argument("--drop-uninformative", action="store_true",
                         help="ONLY use with --label-strategy cluster. Drop rows where every "
                              "--metric-col value is <= --drop-uninformative-thresh (default 0) "
                              "before clustering, to avoid a real distance-clustering artifact "
                              "where genes tied at the floor mechanically merge into one "
                              "spuriously 'perfect' cluster. Does NOT mean these genes are "
                              "biologically uninformative -- run a companion boolean_flag or "
                              "explicit_threshold pass on the same population to test them on "
                              "their own terms, see filters.population."
                              "drop_uninformative_for_clustering and the worked example.")
    parser.add_argument("--drop-uninformative-thresh", type=float, default=0.0,
                         help="Threshold used by --drop-uninformative.")

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
                         help="Optional path to a cached slim .godb. If given, also runs "
                              "enrichment against it and writes to --slim-output-dir.")
    parser.add_argument("--slim-output-dir", default=None,
                         help="Required if --slim-godb is given.")

    return parser


def validate_metric_cols(metric_cols, label_strategy):
    if label_strategy in SINGLE_COLUMN_STRATEGIES and len(metric_cols) != 1:
        raise ValueError(
            f"--label-strategy {label_strategy} requires exactly one --metric-col value, "
            f"got {len(metric_cols)}: {metric_cols}"
        )
    if label_strategy in MULTI_COLUMN_STRATEGIES and len(metric_cols) < 1:
        raise ValueError(
            f"--label-strategy {label_strategy} requires at least one --metric-col value"
        )


def build_population_filters(args):
    """
    Returns a list of (filter_fn, kwargs) tuples, ready to pass to
    filters.population.chain_filters(), in a fixed order: read-depth,
    then usage, then minimum group size. Each is only included if its
    matching *-col flag was given. Raises ValueError if a *-col flag is
    given without its required *-thresh/*-n counterpart, rather than
    silently skipping the filter.
    """
    population_filters = []

    if args.read_depth_col:
        if args.read_depth_thresh is None:
            raise ValueError("--read-depth-thresh is required when --read-depth-col is given")
        population_filters.append(
            (read_depth_filter, {"col": args.read_depth_col, "thresh": args.read_depth_thresh})
        )

    if args.usage_col:
        if args.usage_thresh is None:
            raise ValueError("--usage-thresh is required when --usage-col is given")
        population_filters.append(
            (usage_filter, {"col": args.usage_col, "thresh": args.usage_thresh})
        )

    if args.min_group_col:
        if args.min_group_n is None:
            raise ValueError("--min-group-n is required when --min-group-col is given")
        population_filters.append(
            (min_group_size_filter, {"group_col": args.min_group_col, "min_n": args.min_group_n})
        )

    return population_filters


def read_merge_manifest(manifest_path):
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
    df = pd.read_csv(path, sep="\t", comment="#")
    if strip_id_suffix:
        df[id_col] = df[id_col].apply(lambda g: re.sub(strip_id_suffix, "", g))
    if exclude_ids:
        df = df[~df[id_col].isin(exclude_ids)]
    df = df.dropna(subset=metric_cols)
    return df


def load_merged_table(manifest_path, id_col, strip_id_suffix, exclude_ids):
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
    population_filters = build_population_filters(args)  # validates *-thresh/*-n pairing early

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

    if population_filters:
        n_before = len(df)
        df = chain_filters(df, *population_filters)
        print(f"Population filter (Stage A): {n_before} -> {len(df)} genes")

    if args.drop_uninformative:
        n_before = len(df)
        df = drop_uninformative_for_clustering(
            df, args.metric_col, thresh=args.drop_uninformative_thresh
        )
        print(f"Dropped uninformative rows (all metric columns <= "
              f"{args.drop_uninformative_thresh}): {n_before} -> {len(df)} genes")

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