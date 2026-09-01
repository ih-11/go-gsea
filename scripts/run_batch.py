#!/usr/bin/env python3
"""
scripts/run_batch.py

Batch driver for scripts/run_pipeline.py: reads a manifest (CSV) listing
one row per (species, condition, metric) run and calls run_pipeline.py
once per row via subprocess, so any one row failing (bad annotation file,
missing input, etc.) does not stop the rest of the batch.

Handles every combinatorial case (n species x n condition x single/pair)
by just being a different number of manifest rows -- no new logic per
case. reference/, filters/, labelers/, and enrichment/ are all untouched
by this script.

Manifest CSV required columns:
    species, condition, input_table, godb, metric_col, label_strategy

Manifest CSV optional columns (blank/absent = use run_pipeline.py's own
default for that flag):
    slim_godb, id_col, strip_id_suffix, exclude_id (semicolon-separated
    for multiple), pct, high_thresh, low_thresh, n_clusters, output_dir,
    slim_output_dir, unknown_ratio_thresh, thresh_type, thresh,
    dataset_name

dataset_name, if not given, is auto-derived as
"<species>_<condition>.<metric_col>", matching the convention already
established for Chlamydomonas (e.g. CR_3D.PR_gene).

Example manifest row set covering several of the combinatorial cases:

    species,condition,input_table,godb,metric_col,label_strategy,pct
    CR,3D,data/raw/CR_3D.gene_data.tsv.gz,data/go_reference/CR/CR.godb.pkl,PR_gene,rank_tail,10
    CR,3D,data/raw/CR_3D.gene_data.tsv.gz,data/go_reference/CR/CR.godb.pkl,TPM,rank_tail,10
    AT,HS,data/raw/AT_HS.gene_data.tsv.gz,data/go_reference/AT/AT.godb.pkl,TPM,rank_tail,10
    OS,NC,data/raw/OS_NC.gene_data.tsv.gz,data/go_reference/OS/OS.godb.pkl,TPM,rank_tail,10
"""
import sys
import os
import subprocess
import argparse

import pandas as pd

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN_PIPELINE = os.path.join(REPO_ROOT, "scripts", "run_pipeline.py")

REQUIRED_COLS = ["species", "condition", "input_table", "godb", "metric_col", "label_strategy"]


def build_argv_from_row(row):
    """
    Builds the argv list (excluding the script path itself) to invoke
    run_pipeline.py for one manifest row. Pure function, no I/O, so it
    can be tested without actually running anything.

    row: a dict-like object (e.g. a pandas Series) with at least
    REQUIRED_COLS present and non-null.
    """
    missing = [c for c in REQUIRED_COLS if pd.isna(row.get(c, None))]
    if missing:
        raise ValueError(f"Manifest row missing required column(s): {missing}")

    species = str(row["species"])
    condition = str(row["condition"])
    metric_col = str(row["metric_col"])

    dataset_name = row.get("dataset_name")
    if pd.isna(dataset_name) or not dataset_name:
        dataset_name = f"{species}_{condition}.{metric_col}"

    argv = [
        "--input-table", str(row["input_table"]),
        "--godb", str(row["godb"]),
        "--metric-col", metric_col,
        "--label-strategy", str(row["label_strategy"]),
        "--output-dir", str(row.get("output_dir") or "results/GO"),
        "--dataset-name", dataset_name,
    ]

    if not pd.isna(row.get("slim_godb", None)) and row.get("slim_godb"):
        argv += ["--slim-godb", str(row["slim_godb"])]
        argv += ["--slim-output-dir", str(row.get("slim_output_dir") or "results/GO_slim")]

    if not pd.isna(row.get("id_col", None)) and row.get("id_col"):
        argv += ["--id-col", str(row["id_col"])]

    if not pd.isna(row.get("strip_id_suffix", None)) and row.get("strip_id_suffix"):
        argv += ["--strip-id-suffix", str(row["strip_id_suffix"])]

    if not pd.isna(row.get("exclude_id", None)) and row.get("exclude_id"):
        for ex_id in str(row["exclude_id"]).split(";"):
            ex_id = ex_id.strip()
            if ex_id:
                argv += ["--exclude-id", ex_id]

    if not pd.isna(row.get("pct", None)):
        argv += ["--pct", str(row["pct"])]
    if not pd.isna(row.get("high_thresh", None)):
        argv += ["--high-thresh", str(row["high_thresh"])]
    if not pd.isna(row.get("low_thresh", None)):
        argv += ["--low-thresh", str(row["low_thresh"])]
    if not pd.isna(row.get("n_clusters", None)):
        argv += ["--n-clusters", str(int(row["n_clusters"]))]

    if not pd.isna(row.get("unknown_ratio_thresh", None)):
        argv += ["--unknown-ratio-thresh", str(row["unknown_ratio_thresh"])]
    if not pd.isna(row.get("thresh_type", None)) and row.get("thresh_type"):
        argv += ["--thresh-type", str(row["thresh_type"])]
    if not pd.isna(row.get("thresh", None)):
        argv += ["--thresh", str(row["thresh"])]

    return argv


def run_batch(manifest_path, dry_run=False):
    """
    Reads manifest_path (CSV) and runs run_pipeline.py once per row via
    subprocess. One row failing does not stop the rest.

    Returns {"succeeded": [dataset_name, ...], "failed": [(dataset_name, stderr), ...]}.
    """
    manifest = pd.read_csv(manifest_path)
    succeeded, failed = [], []

    for i, row in manifest.iterrows():
        try:
            argv = build_argv_from_row(row)
        except ValueError as e:
            failed.append((f"row {i}", str(e)))
            continue

        dataset_name = argv[argv.index("--dataset-name") + 1]

        if dry_run:
            print(f"[dry run] {dataset_name}: python {RUN_PIPELINE} {' '.join(argv)}")
            succeeded.append(dataset_name)
            continue

        print(f"Running {dataset_name} ...")
        result = subprocess.run(
            [sys.executable, RUN_PIPELINE] + argv,
            capture_output=True, text=True,
        )
        if result.returncode == 0:
            print(result.stdout)
            succeeded.append(dataset_name)
        else:
            print(f"FAILED: {dataset_name}")
            print(result.stderr)
            failed.append((dataset_name, result.stderr))

    print()
    print(f"Batch complete: {len(succeeded)} succeeded, {len(failed)} failed")
    if failed:
        print("Failed runs:")
        for name, err in failed:
            last_line = err.strip().splitlines()[-1] if err.strip() else "(no stderr)"
            print(f"  {name}: {last_line}")

    return {"succeeded": succeeded, "failed": failed}


def main():
    parser = argparse.ArgumentParser(description="Batch-run run_pipeline.py from a manifest CSV.")
    parser.add_argument("manifest", help="Path to the manifest CSV.")
    parser.add_argument("--dry-run", action="store_true",
                         help="Print what would be run, without actually running anything.")
    args = parser.parse_args()
    run_batch(args.manifest, dry_run=args.dry_run)


if __name__ == "__main__":
    main()