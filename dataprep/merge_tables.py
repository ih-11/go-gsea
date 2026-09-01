"""
dataprep/merge_tables.py

Generic N-source table merge: joins any number of named columns from any
number of gene-level files into one wide table, keyed on gene_id.

This is the one new primitive needed for real multi-condition (or
multi-technology) clustering -- e.g. PR_gene measured at CR_3D and CR_6D,
combined into a single wide table with one column per condition, before
handing it to labelers.cluster(). Knows nothing about species, condition,
or technology: it only ever sees paths and column names, the caller
decides what those represent.

Default join is INNER, matching the real precedent's own stated rule
(hierarchical clustering methods, section 2.3): "only genes with
[data] available across all datasets used in the analysis were selected."
"""
import pandas as pd


def merge_gene_tables(sources, id_col="gene_id", how="inner"):
    """
    sources: list of (path, value_col, output_col_name) triples.
        path: path to a gene-level TSV (plain or .gz)
        value_col: the column in that file to keep
        output_col_name: what to call it in the merged output (e.g. a
            condition or dataset label, so the wide table's columns are
            distinguishable even if value_col is the same name in every
            source file, e.g. "PR_gene" everywhere)
    id_col: gene ID column name, assumed the same across all sources.
    how: join type passed to pandas.merge. "inner" (default) matches the
        real precedent: only genes present in every source survive.
        "outer" is available if a caller genuinely wants the union
        instead, but that is a deliberate deviation from the precedent,
        not the default.

    Returns one wide DataFrame: id_col, <output_col_name_1>, <output_col_name_2>, ...

    Raises ValueError if sources is empty, or if any two entries share
    the same output_col_name (which would silently collide in the merge).
    """
    if not sources:
        raise ValueError("sources must contain at least one (path, value_col, output_col_name) triple")

    output_names = [s[2] for s in sources]
    if len(output_names) != len(set(output_names)):
        dupes = {n for n in output_names if output_names.count(n) > 1}
        raise ValueError(f"Duplicate output_col_name(s) in sources: {dupes}")

    merged = None
    for path, value_col, output_col_name in sources:
        opener_kwargs = {"compression": "infer"}
        df = pd.read_csv(path, sep="\t", comment="#", **opener_kwargs)
        if id_col not in df.columns:
            raise ValueError(f"'{id_col}' not found in {path}. Columns present: {list(df.columns)}")
        if value_col not in df.columns:
            raise ValueError(f"'{value_col}' not found in {path}. Columns present: {list(df.columns)}")

        sub = df[[id_col, value_col]].rename(columns={value_col: output_col_name})

        if merged is None:
            merged = sub
        else:
            merged = merged.merge(sub, on=id_col, how=how)

    return merged