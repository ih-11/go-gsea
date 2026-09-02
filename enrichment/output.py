"""
enrichment/output.py

Writes run_ora() results to disk, matching the real precedent's output
shape: one "all" file (every term tested) plus, per class, an "over" file
(significant, positive fold_enrichment) and an "under" file (significant,
negative fold_enrichment), only written if non-empty.

Optionally embeds a provenance metadata footer (JSON, one line, '#'-
prefixed) at the top of every file written, matching the real precedent's
own output convention. This is fully backward compatible: every reader in
this codebase already opens result files with pd.read_csv(..., comment="#"),
so the provenance line is silently skipped by anything that doesn't know
about it, exactly like it already skips any other '#'-prefixed line.
"""
import os
import platform
import subprocess
import json
from datetime import datetime, timezone

import pandas as pd


def _get_git_commit_hash():
    """
    Returns the current go-gsea git commit hash, or None if unavailable
    (not a git repo, git not installed, etc.) -- never raises, provenance
    is a best-effort record, not something that should ever block a real
    run from completing.
    """
    try:
        repo_root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        result = subprocess.run(
            ["git", "rev-parse", "HEAD"],
            cwd=repo_root, capture_output=True, text=True, timeout=5,
        )
        if result.returncode == 0:
            return result.stdout.strip()
    except Exception:
        pass
    return None


def _get_package_versions():
    """
    Returns {package_name: version_string_or_'not installed'} for the
    key packages this pipeline's correctness actually depends on. Never
    raises on an individual package failing to import.
    """
    packages = ["pandas", "numpy", "scipy", "statsmodels", "goatools", "gseapy"]
    versions = {}
    for pkg in packages:
        try:
            mod = __import__(pkg)
            versions[pkg] = getattr(mod, "__version__", "unknown")
        except Exception:
            versions[pkg] = "not installed"
    return versions


def build_provenance(run_parameters=None):
    """
    Builds a provenance metadata dict: when this ran, what code produced
    it (git commit, if available), what package versions were active, and
    whatever run-specific parameters the caller supplies (e.g. thresh_type,
    thresh, label_strategy, metric_col, dataset_name).

    run_parameters: any JSON-serializable dict, passed through as-is under
    the "run_parameters" key. Not validated or interpreted by this
    function, just recorded.
    """
    return {
        "tool": "go-gsea",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "git_commit": _get_git_commit_hash(),
        "python_version": platform.python_version(),
        "package_versions": _get_package_versions(),
        "run_parameters": run_parameters or {},
    }


def _write_tsv_with_provenance(df, path, provenance):
    """
    Writes df to path as a TSV. If provenance is given, prepends one
    '#'-prefixed line containing it as JSON, before the header row.
    """
    if provenance is None:
        df.to_csv(path, sep="\t", index=False)
        return

    with open(path, "w") as f:
        f.write("#" + json.dumps(provenance) + "\n")
        df.to_csv(f, sep="\t", index=False)


def write_results(result_df, output_dir, dataset_name, provenance=None):
    """
    Writes result_df (as produced by enrichment.ora.run_ora) to output_dir.

    Layout:
        output_dir/summary.<dataset_name>.all.tsv      -- every row tested
        output_dir/summary.<dataset_name>.<class>.over.tsv    -- per class, if any
        output_dir/summary.<dataset_name>.<class>.under.tsv   -- per class, if any

    provenance: optional dict (e.g. from build_provenance()). If given,
    written as a '#'-prefixed JSON line at the top of every file. If
    omitted (the default), output is byte-identical to every prior
    version of this function -- fully backward compatible.

    Returns a dict of the paths actually written (all.tsv always included,
    even if result_df is empty -- over/under only present if non-empty).
    """
    os.makedirs(output_dir, exist_ok=True)
    written = {}

    all_path = os.path.join(output_dir, f"summary.{dataset_name}.all.tsv")
    _write_tsv_with_provenance(result_df, all_path, provenance)
    written["all"] = all_path

    if len(result_df) == 0:
        return written

    for cls in sorted(result_df["class"].unique()):
        cls_df = result_df[result_df["class"] == cls]

        over_df = cls_df[cls_df["significance"] & (cls_df["fold_enrichment"] > 0)]
        if len(over_df) > 0:
            over_path = os.path.join(output_dir, f"summary.{dataset_name}.{cls}.over.tsv")
            _write_tsv_with_provenance(over_df, over_path, provenance)
            written[f"{cls}.over"] = over_path

        under_df = cls_df[cls_df["significance"] & (cls_df["fold_enrichment"] < 0)]
        if len(under_df) > 0:
            under_path = os.path.join(output_dir, f"summary.{dataset_name}.{cls}.under.tsv")
            _write_tsv_with_provenance(under_df, under_path, provenance)
            written[f"{cls}.under"] = under_path

    return written