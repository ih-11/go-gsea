#!/usr/bin/env python3
"""
scripts/run_batch.py

Batch driver for scripts/run_pipeline.py: reads a manifest (plain-text,
INI-style) listing one [section] per (species, condition, metric) run and
calls run_pipeline.py once per section via subprocess, so any one section
failing (bad annotation file, missing input, etc.) does not stop the rest
of the batch.

Manifest format chosen over CSV specifically because it is easier for a
beginner to hand-edit: no column alignment, no need to leave blank cells
for fields that do not apply to a given run, comments are supported (#
or ;), and each run's full set of settings is visually grouped together
instead of spread across one wide row.

Required keys per section:
    species, condition, input_table, godb, metric_col, label_strategy

Optional keys (omit the line entirely if not needed):
    slim_godb, id_col, strip_id_suffix, exclude_id (semicolon-separated
    for multiple), pct, high_thresh, low_thresh, n_clusters, output_dir,
    slim_output_dir, unknown_ratio_thresh, thresh_type, thresh,
    dataset_name

dataset_name, if not given, is auto-derived as
"<species>_<condition>.<metric_col>".

Example manifest.txt:

    [CR_3D_PR_gene]
    species = CR
    condition = 3D
    input_table = data/raw/CR_3D.gene_data.tsv.gz
    godb = data/go_reference/CR/CR.godb.pkl
    metric_col = PR_gene
    label_strategy = rank_tail
    pct = 10

    [CR_3D_TPM]
    species = CR
    condition = 3D
    input_table = data/raw/CR_3D.gene_data.tsv.gz
    godb = data/go_reference/CR/CR.godb.pkl
    metric_col = TPM
    label_strategy = rank_tail
    pct = 10
"""
import sys
import os
import subprocess
import argparse
import configparser

REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RUN_PIPELINE = os.path.join(REPO_ROOT, "scripts", "run_pipeline.py")

REQUIRED_KEYS = ["species", "condition", "input_table", "godb", "metric_col", "label_strategy"]


def read_manifest(manifest_path):
    """
    Reads a plain-text, INI-style manifest. Returns a list of
    (section_name, dict) pairs, one per run, in file order.
    """
    parser = configparser.ConfigParser()
    parser.optionxform = str  # keep keys as written, don't lowercase them
    read_files = parser.read(manifest_path)
    if not read_files:
        raise FileNotFoundError(f"Manifest file not found or empty: {manifest_path}")
    return [(section, dict(parser[section])) for section in parser.sections()]


def build_argv_from_row(row, row_label="run"):
    """
    Builds the argv list (excluding the script path itself) to invoke
    run_pipeline.py for one manifest run. Pure function, no I/O.

    row: a plain dict of string keys -> string values.
    """
    missing = [k for k in REQUIRED_KEYS if not row.get(k)]
    if missing:
        raise ValueError(f"Manifest section '{row_label}' missing required key(s): {missing}")

    species = row["species"]
    condition = row["condition"]
    metric_col = row["metric_col"]

    dataset_name = row.get("dataset_name") or f"{species}_{condition}.{metric_col}"

    argv = [
        "--input-table", row["input_table"],
        "--godb", row["godb"],
        "--metric-col", metric_col,
        "--label-strategy", row["label_strategy"],
        "--output-dir", row.get("output_dir") or "results/GO",
        "--dataset-name", dataset_name,
    ]

    if row.get("slim_godb"):
        argv += ["--slim-godb", row["slim_godb"]]
        argv += ["--slim-output-dir", row.get("slim_output_dir") or "results/GO_slim"]

    if row.get("id_col"):
        argv += ["--id-col", row["id_col"]]
    if row.get("strip_id_suffix"):
        argv += ["--strip-id-suffix", row["strip_id_suffix"]]
    if row.get("exclude_id"):
        for ex_id in row["exclude_id"].split(";"):
            ex_id = ex_id.strip()
            if ex_id:
                argv += ["--exclude-id", ex_id]

    if row.get("pct"):
        argv += ["--pct", row["pct"]]
    if row.get("high_thresh"):
        argv += ["--high-thresh", row["high_thresh"]]
    if row.get("low_thresh"):
        argv += ["--low-thresh", row["low_thresh"]]
    if row.get("n_clusters"):
        argv += ["--n-clusters", row["n_clusters"]]

    if row.get("unknown_ratio_thresh"):
        argv += ["--unknown-ratio-thresh", row["unknown_ratio_thresh"]]
    if row.get("thresh_type"):
        argv += ["--thresh-type", row["thresh_type"]]
    if row.get("thresh"):
        argv += ["--thresh", row["thresh"]]

    return argv


def run_batch(manifest_path, dry_run=False):
    """
    Reads manifest_path (plain-text, INI-style) and runs run_pipeline.py
    once per section via subprocess. One section failing does not stop
    the rest.
    """
    rows = read_manifest(manifest_path)
    succeeded, failed = [], []

    for section_name, row in rows:
        try:
            argv = build_argv_from_row(row, row_label=section_name)
        except ValueError as e:
            failed.append((section_name, str(e)))
            continue

        dataset_name = argv[argv.index("--dataset-name") + 1]

        if dry_run:
            print(f"[dry run] {section_name} -> {dataset_name}: "
                  f"python {RUN_PIPELINE} {' '.join(argv)}")
            succeeded.append(dataset_name)
            continue

        print(f"Running {section_name} -> {dataset_name} ...")
        result = subprocess.run(
            [sys.executable, RUN_PIPELINE] + argv,
            capture_output=True, text=True,
        )
        if result.returncode == 0:
            print(result.stdout)
            succeeded.append(dataset_name)
        else:
            print(f"FAILED: {section_name} -> {dataset_name}")
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
    parser = argparse.ArgumentParser(
        description="Batch-run run_pipeline.py from a plain-text, INI-style manifest."
    )
    parser.add_argument("manifest", help="Path to the manifest .txt file.")
    parser.add_argument("--dry-run", action="store_true",
                         help="Print what would be run, without actually running anything.")
    args = parser.parse_args()
    run_batch(args.manifest, dry_run=args.dry_run)


if __name__ == "__main__":
    main()