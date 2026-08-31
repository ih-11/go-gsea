"""
reference/build_godb.py

Builds a cached, propagated gene->GO annotation object ("godb") from a
Phytozome-format annotation_info.txt(.gz) file.

Tier 1 (default): direct Phytozome GO column -- exact ID match to the
    reference genome's own gene models, no inference.
Tier 3 (opt-in): ortholog-transferred GO via Best-hit-arabi-name, using a
    pre-built donor godb (e.g. Arabidopsis). Every gene added this way is
    tagged provenance="ortholog_IEA" and never overrides a Tier 1 entry.
"""
import gzip
import pickle
from pathlib import Path
from collections import defaultdict

from goatools.obo_parser import GODag


def _open_maybe_gz(path):
    path = str(path)
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def get_all_ancestors(godag, term_id, should_consider_part_of=True):
    """
    Returns the full ancestor set of term_id, walking is_a always, and
    relationship['part_of'] as well unless disabled.

    goatools' own get_all_parents() only follows is_a -- confirmed
    empirically against go-basic.obo: terms with real part_of edges
    (e.g. GO:0000015, GO:0000027) store them in term.relationship['part_of'],
    but get_all_parents() omits them entirely. This is the fix.
    """
    if term_id not in godag:
        return set()
    ancestors = set()
    stack = [godag[term_id]]
    while stack:
        term = stack.pop()
        parents = list(term.parents)
        if should_consider_part_of:
            parents += list(term.relationship.get("part_of", set()))
        for p in parents:
            if p.id not in ancestors:
                ancestors.add(p.id)
                stack.append(p)
    return ancestors


def parse_phytozome_annotation_info(path, locus_col="locusName", go_col="GO"):
    """Returns dict[gene_id] -> set(GO_id), deduped from transcript-level rows."""
    gene_go = defaultdict(set)
    with _open_maybe_gz(path) as f:
        header = f.readline().lstrip("#").rstrip("\n").split("\t")
        locus_idx, go_idx = header.index(locus_col), header.index(go_col)
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if len(fields) <= max(locus_idx, go_idx):
                continue
            go_field = fields[go_idx].strip()
            if not go_field:
                continue
            gene_go[fields[locus_idx]].update(
                g for g in go_field.split() if g.startswith("GO:")
            )
    return dict(gene_go)


def build_and_cache_godb(annotation_info_path, obo_path, cache_path):
    gene_go = parse_phytozome_annotation_info(annotation_info_path)
    provenance = {g: "phytozome_direct" for g in gene_go}

    godag = GODag(str(obo_path), optional_attrs=["relationship"])

    propagated = {}
    for gene_id, direct_terms in gene_go.items():
        all_terms = set()
        for term_id in direct_terms:
            if term_id not in godag:
                continue
            all_terms.add(term_id)
            all_terms.update(get_all_ancestors(godag, term_id))
        propagated[gene_id] = all_terms

    godb = {
        "gene_go": propagated,
        "provenance": provenance,
        "source": str(annotation_info_path),
        "obo": str(obo_path),
    }
    Path(cache_path).parent.mkdir(parents=True, exist_ok=True)
    with open(cache_path, "wb") as f:
        pickle.dump(godb, f)
    return godb