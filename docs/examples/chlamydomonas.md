# Worked example: Chlamydomonas reinhardtii, polysome ratio (PR)

The first real, end-to-end application of `go-gsea`. This document records
the specific choices made for this species and dataset, and why, following
the general principles in the main `README.md` (section 4), but making
decisions a different organism, condition, or metric might reasonably make
differently.

---

## Research question

Translational efficiency (Polysome Ratio, PR, from `Project1_polysomelongread_PR`)
compared between the top and bottom 10% of genes by `PR_gene`: which GO terms
are over or under represented among the best and worst translated genes?

## Input data used

A gene-level table, `CR_3D.PR.gene_data.tsv.gz`, produced by the upstream
long-read pipeline described in the main README (section 2). Relevant
columns for this run: `gene_id`, `PR_gene`. One row is excluded before
labeling: the luciferase spike-in control (`gene:Standard-R-luc`), since it
is a synthetic construct with no meaningful GO annotation.

## GO annotation source: Phytozome, not UniProt-GOA

The obvious first choice for a UniProt-GOA-format GAF file for this species,
`32602.C_reinhardtii_CC-503.goa`, from EBI's GOA proteomes index, turned out
to be annotated against strain CC-503, on the older v5.6 genome assembly.
This project's actual reference genome (and every gene ID in
`gene_data.tsv.gz`) is strain CC-4532, assembly v6.1, a separately
reassembled genome, confirmed from JGI's own documentation to share "close
ancestry" with CC-503 but not be a renumbering of it. Using the CC-503 file
risked a large, silent gene-ID mismatch.

Phytozome's own `CreinhardtiiCC_4532_707_v6.1.P14.annotation_info.txt.gz`
(a per-gene functional annotation table, Pfam, Panther, KOG, KO, GO,
best-hit orthologs, bundled with the same genome release the reference
genome itself came from) uses the exact same gene ID namespace by
construction. Confirmed directly: real gene IDs from `gene_data.tsv.gz`
match this file's `locusName` column after stripping only a trailing
`.vX.Y` version suffix (for example `Cre01.g000350_4532.v6.1` becomes
`Cre01.g000350_4532`).

**Tradeoff accepted:** Phytozome's own GO column has substantially sparser
coverage than UniProt-GOA (confirmed: about 19% of genes genome-wide, 3,203
of about 17,700 loci). Phytozome's GO column is the direct output of their
own annotation pipeline, not a GO-specialist resource layering in orthology
transfer and electronic annotation the way UniProt-GOA does. Given the
project's own methodology explicitly prioritizes population and ID
correctness over raw coverage (main README section 4), this tradeoff was
accepted deliberately rather than reaching for a larger but ID-mismatched
source.

**A same-species crosswalk was checked and ruled out.** Phytozome also
ships `CreinhardtiiCC_4532_707_v6.1.synonym.txt`, which might plausibly
have bridged CC-503-era IDs to CC-4532 IDs, making the UniProt-GOA file
usable after all as a supplementary tier. Checked directly: every ID in
this file is already `_4532`-suffixed, it is a gene-symbol and alias table
(mapping v6.1 transcript IDs to gene symbols and pre-Phytozome internal
gene-caller IDs like `g4.t1`), not a cross-assembly ID map. This path was
dropped rather than left as an unverified assumption.

**Designed but not built:** an ortholog-transfer supplement using
Phytozome's `Best-hit-arabi-name` column (available in the same
`annotation_info.txt` file) plus a donor `.godb` (for example Arabidopsis).
For genes with zero direct Phytozome GO annotation, the idea is to
transfer GO terms from the best-hit Arabidopsis ortholog, tagged with a
distinct provenance (`ortholog_IEA`) so it is never conflated with direct
annotation. Not yet implemented as reusable `go-gsea` code.

## `is_a` and `part_of` propagation

Verified empirically against `go-basic.obo`, not assumed from
documentation: `goatools.GODag.get_all_parents()` only follows `is_a`,
even with `optional_attrs=["relationship"]` loaded. Confirmed via two real
terms known to have no `part_of` edges (`GO:0006914` autophagy,
`GO:0071456` cellular response to hypoxia, both returned an empty
`part_of` set correctly) versus terms confirmed via raw grep of the OBO
file to have real `part_of` edges (`GO:0000015`, `GO:0000027`, and others,
both returned real parent IDs via `term.relationship["part_of"]`).
`reference/build_godb.py`'s custom `get_all_ancestors()` walker was built
and verified against this distinction.

Sanity check on one real gene: `Cre01.g000650_4532` has 4 direct GO terms
(copper ion binding, primary methylamine oxidase activity, amine metabolic
process, quinone binding) leading to 19 terms after propagation, including
both namespace roots and no unrelated terms, manually checked against
`go-basic.obo`'s term names, not just trusted by count.

## GO-slim

Built with `goslim_generic.obo` (206 terms), from
`reference/build_godb.py`'s `build_and_cache_slim_godb()`, which intersects
the already-propagated full-GO term set per gene against the slim
vocabulary rather than re-propagating from scratch. 2,100 of the 3,203
full-GO-annotated genes have at least one slim term (average 2.8 slim
terms per gene, versus 16.3 full terms per gene). Verified with a direct
subset check before trusting the result: every one of the 2,100 slim
genes' term sets was confirmed to be a strict subset of that same gene's
full term set, zero violations.

## Real-data integration numbers

- `build_godb.py` on the real Phytozome file: 3,203 annotated genes,
  matching an independently predicted count from the raw file before any
  code existed to build it.
- `gene_data.tsv.gz`: 3,594 total genes (after dropping the luciferase
  spike-in row and rows with missing `PR_gene`: 3,593).
- ID overlap after `.vX.Y`-suffix normalization: 650 of 3,594 matched
  (18.1%), consistent with the about-19% genome-wide Phytozome coverage
  figure above, confirming the ID-matching logic is correct rather than
  coincidentally non-crashing.
- Unknown-gene ratio for this run: about 82% (2,943 of 3,593 unannotated
  genes), a direct consequence of the coverage tradeoff above. `run_ora()`'s
  `unknown_ratio_thresh` was set to `0.9` for this run (main README section
  4 explains why this threshold is a tunable parameter, not a fixed
  constant: two copies of the original precedent script disagreed on this
  exact value, `0.2` versus `0.9`, and which was current was never
  definitively resolved).

## Labeling strategy used

`labelers.rank_tail(df, col="PR_gene", pct=10)`: top and bottom 10% by
`PR_gene`, matching the confirmed real-data convention from prior lab work
on Arabidopsis and rice (`AT-T87-HS-PR.tsv` and similar files were
confirmed to use exactly 10%/10% splits, not a median split). 359 genes in
each of `High`/`Low`.

## Result, and how to read it

At `thresh_type="p", thresh=0.01`, full-GO:

- **Low-PR (worst translated) genes:** enriched for cellular component and
  organelle organization and assembly terms (`GO:0016043`, `GO:0022607`,
  `GO:0071840`, `GO:0044085`, all p<0.01). Notably, `GO:0006412`
  (translation) had zero observed genes in this group against an expected
  about 5.05 (fold_enrichment = -2.60, p = 0.0103, just outside the p<0.01
  cutoff). Core translation-machinery genes essentially never appearing
  among the worst translated genes is the expected biological direction,
  and a reassuring sign this is not a pipeline artifact.
- **High-PR (best translated) genes:** `GO:0046982`/`GO:0046983` (protein
  dimerization activity) over represented; `GO:0008152` (metabolic
  process) under represented.

**GO-slim, same labeled genes, same threshold:** 172 (class, term) rows
tested (versus 2,334 for full-GO), 0 rows reaching p<0.01. Checked
directly rather than assumed correct: the actual p-value spread in
`summary.CR_3D.PR_gene.all.tsv` ranges from about 0.016 up through 1.0,
not a suspicious wall of identical values, so this reflects a real
absence of a slim-vocabulary-detectable signal at this threshold, not a
broken computation. The closest term, `GO:0005694` (chromosome), sits at
p = 0.021, fold_enrichment = 1.40. Plausible explanation: the 206-term
slim vocabulary is coarse enough that the specific translation-machinery
signal visible in full-GO (`GO:0006412`, itself only just outside the
cutoff) does not survive being folded into a broader category.

Per the main README's language-constraint principle, these are described
as genes that tended to include these terms, not as statistically
significant findings, since this threshold applies no multiple-testing
correction.

**Caveat to carry into any write-up using this result:** every full-GO
number above is built on the about-19%-coverage Phytozome annotation
layer, not a comprehensive one. A `GO:0006412` population of 49 genes is
out of about 650 GO-annotated genes in this run, not the full 3,593-gene
analyzed set.

## Reproducing this run

Through the CLI (`scripts/run_pipeline.py`), both full-GO and GO-slim in
one call:

```
python scripts/run_pipeline.py \
    --input-table /path/to/CR_3D.PR.gene_data.tsv.gz \
    --godb data/go_reference/chlamy.godb.pkl \
    --slim-godb data/go_reference/chlamy.slim.godb.pkl \
    --metric-col PR_gene \
    --id-col gene_id \
    --strip-id-suffix '\.v\d+\.\d+$' \
    --exclude-id 'gene:Standard-R-luc' \
    --label-strategy rank_tail --pct 10 \
    --output-dir results/GO \
    --slim-output-dir results/GO_slim \
    --dataset-name CR_3D.PR_gene \
    --unknown-ratio-thresh 0.9 \
    --thresh-type p --thresh 0.01
```

Confirmed to reproduce the exact numbers above (359 High, 359 Low, 2,334
full-GO rows tested, 7 significant, 172 slim rows tested, 0 significant)
against an earlier, independent hand-written script that called the same
library functions directly, before `scripts/run_pipeline.py` existed.

## Open items specific to this example

- Only `PR_gene` (gene-level PR) has been run. `TPM` and `PTPM` (which
  fraction bare `TPM` represents in `gene_data.tsv.gz` was never fully
  confirmed) and transcript- or variant-level PR are unexercised.
- Ortholog-transfer supplementation and provenance-stamped output files
  are still open across the whole tool (main README section 8), nothing
  here is Chlamydomonas-specific about those gaps.