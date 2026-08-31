"""
enrichment/output.py

Writes run_ora() results to disk, matching the real precedent's output
shape: one "all" file (every term tested) plus, per class, an "over" file
(significant, positive fold_enrichment) and an "under" file (significant,
negative fold_enrichment), only written if non-empty.
"""
import os
import pandas as pd


def write_results(result_df, output_dir, dataset_name):
    """
    Writes result_df (as produced by enrichment.ora.run_ora) to output_dir.

    Layout:
        output_dir/summary.<dataset_name>.all.tsv      -- every row tested
        output_dir/summary.<dataset_name>.<class>.over.tsv    -- per class, if any
        output_dir/summary.<dataset_name>.<class>.under.tsv   -- per class, if any

    Returns a dict of the paths actually written (all.tsv always included,
    even if result_df is empty -- over/under only present if non-empty).
    """
    os.makedirs(output_dir, exist_ok=True)
    written = {}

    all_path = os.path.join(output_dir, f"summary.{dataset_name}.all.tsv")
    result_df.to_csv(all_path, sep="\t", index=False)
    written["all"] = all_path

    if len(result_df) == 0:
        return written

    for cls in sorted(result_df["class"].unique()):
        cls_df = result_df[result_df["class"] == cls]

        over_df = cls_df[cls_df["significance"] & (cls_df["fold_enrichment"] > 0)]
        if len(over_df) > 0:
            over_path = os.path.join(output_dir, f"summary.{dataset_name}.{cls}.over.tsv")
            over_df.to_csv(over_path, sep="\t", index=False)
            written[f"{cls}.over"] = over_path

        under_df = cls_df[cls_df["significance"] & (cls_df["fold_enrichment"] < 0)]
        if len(under_df) > 0:
            under_path = os.path.join(output_dir, f"summary.{dataset_name}.{cls}.under.tsv")
            under_df.to_csv(under_path, sep="\t", index=False)
            written[f"{cls}.under"] = under_path

    return written