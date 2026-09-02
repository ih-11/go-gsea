# Worked example: Chlamydomonas reinhardtii

Seven real, end-to-end applications of `go-gsea` on real data: two
growth-stage conditions (3 and 6 days), three metrics, two labeling
strategies with real lab precedent, and both single-file and genuine
merged-multi-condition scopes. This document records the specific choices
made for this species and dataset, and why, following the general
principles in the main `README.md` (section 4), but making decisions a
different organism, condition, or metric might reasonably make
differently.

---

## Demo gallery

| # | Strategy | Metric(s) | Scope | Genes (per class) | Full-GO sig. | Slim sig. |
|---|---|---|---|---|---|---|
| A | `rank_tail` | `PR_gene` | CR_3D alone | 359 / 359 | 7 | 0 |
| B | `rank_tail` | `TPM` | CR_3D alone | 359 / 359 | 25 | 5 |
| C | `cluster`, k=4 | `PR_gene` | CR_3D alone | 470/431/1612/1080 | 13 | 4 |
| D | `rank_tail` | `PTPM` | CR_6D alone | 305 / 305 | 24 | 4 |
| E | `cluster`, k=3 | `PTPM`, `TPM` jointly | CR_6D alone | 1156/689/1212 | 21 | 4 |
| F | `cluster`, k=2 | `PR_gene` | CR_3D + CR_6D merged | 332 / 2326 | 28 | 6 |
| G | `cluster`, k=2 | `TPM` | CR_3D + CR_6D merged | 1374 / 1284 | 18 | 3 |

**Every one of the seven demos independently surfaces the same core
biological signal**: ribosome and translation-machinery genes
(`GO:0003735` structural constituent of ribosome, `GO:0005840` ribosome,
`GO:0006412` translation, `GO:0160307` protein biosynthetic process),
over-represented in high-expression/high-PR/consistently-well-translated
groups and under-represented in the opposite groups, across two
conditions, three metrics, both labeling strategies, and both scopes.
That convergence, not any single demo, is the strongest evidence in this
project that the pipeline is finding something real rather than an
artifact of one dataset or one method.

---

## Shared setup, applies to all seven demos

### Input data

Gene-level tables `CR_3D.PR.gene_data.tsv.gz` and `CR_6D.PR.gene_data.tsv.gz`,
produced by the upstream long-read pipeline (main README section 2).
Relevant columns: `gene_id`, `PR_gene`, `PTPM`, `TPM`. The luciferase
spike-in control (`gene:Standard-R-luc`) is excluded in every demo, it is
a synthetic construct with no meaningful GO annotation.

**`PTPM`/`TPM`/`PR_gene` are not three independent metrics.** Confirmed
from the upstream pipeline's own source code: `PTPM` is Polysome-fraction
TPM, `TPM` is Total-fraction TPM, and:

```
PR_gene = (PTPM / luciferase_P_TPM) / (TPM / luciferase_T_TPM)
```

Demos comparing across these three are useful checks that the pipeline
generalizes across different label distributions, not independent
validations of unrelated questions (main README section 4).

### GO annotation source: Phytozome, not UniProt-GOA

The obvious UniProt-GOA file for this species (`32602.C_reinhardtii_CC-503.goa`)
is annotated against strain CC-503, assembly v5.6, a different, older
assembly than the actual reference genome used here (strain CC-4532,
assembly v6.1). Using it risked a large, silent gene-ID mismatch.
Phytozome's own `CreinhardtiiCC_4532_707_v6.1.P14.annotation_info.txt.gz`,
bundled with the same genome release the reference genome itself came
from, uses the exact same gene ID namespace by construction, confirmed
directly (gene IDs match after stripping only a trailing `.vX.Y` version
suffix).

**Tradeoff accepted**: Phytozome's own GO column has substantially sparser
coverage than UniProt-GOA (confirmed: about 19% of genes genome-wide,
3,203 of about 17,700 loci). Given the project's methodology explicitly
prioritizes population/ID correctness over raw coverage, this was accepted
deliberately.

**A same-species crosswalk was checked and ruled out.** Phytozome's
`CreinhardtiiCC_4532_707_v6.1.synonym.txt` might plausibly have bridged
CC-503-era IDs to CC-4532 IDs. Checked directly: every ID in it is already
`_4532`-suffixed, it is a gene-symbol/alias table, not a cross-assembly ID
map.

**Designed but not built**: an ortholog-transfer supplement using
Phytozome's `Best-hit-arabi-name` column plus a donor `.godb` (e.g.
Arabidopsis), for genes with zero direct Phytozome annotation.

### `is_a` and `part_of` propagation

`goatools.GODag.get_all_parents()` was confirmed to only follow `is_a`,
even with `optional_attrs=["relationship"]` loaded. Confirmed via terms
known to have no `part_of` edges (returning correctly empty) versus terms
confirmed via raw grep of the OBO file to have real `part_of` edges
(returning real parent IDs via `term.relationship["part_of"]`).
`reference/build_godb.py`'s custom `get_all_ancestors()` walker combines
both, verified against this distinction, and against one real gene
(4 direct terms -> 19 propagated terms, manually checked as a coherent,
correct generalization).

### GO-slim

Built with `goslim_generic.obo` (206 terms), intersecting the
already-propagated full-GO term set per gene against the slim vocabulary.
2,100 of 3,203 full-GO-annotated genes have at least one slim term.
Verified with a direct subset check: every slim gene's term set is a
strict subset of that gene's full term set, zero violations across 2,100
genes.

### Population filter check (real data)

`filters/population.py`'s `read_depth_filter` and `usage_filter`, checked
against real `T_n_reads`/`rTrans-usage.b` columns:

| Filter | Genes passing (of 3,594) |
|---|---|
| `read_depth_filter(col="T_n_reads", thresh=10)` alone | 2,232 |
| `usage_filter(col="rTrans-usage.b", thresh=0.05)` alone | 3,346 |
| Both, chained | 1,984 |

Both filters independently removed a real, distinct subset (neither a
no-op), chained result stricter than either alone, consistent with
correct AND composition.

### Determinism check

Ward's-method clustering (demo C) was rerun back-to-back with identical
parameters; the two output files were confirmed byte-identical via `diff`
(exit code 0). Clustering is deterministic within a session on this
codebase and library versions. **This does not guarantee a result
documented at one point in time will still reproduce exactly later** — see
demo C below for a real case where it did not.

---

## Demo A: `rank_tail` on `PR_gene`, CR_3D alone

Top/bottom 10% by `PR_gene` (359 genes each), matching the confirmed
real-data convention from prior lab work on Arabidopsis/rice (exact
10%/10% splits, not a median split).

**Low-PR (worst translated):** enriched for cellular component/organelle
organization and assembly terms (`GO:0016043`, `GO:0022607`, `GO:0071840`,
`GO:0044085`, all p<0.01). `GO:0006412` (translation) had zero observed
genes against an expected ~5.05 (fold_enrichment=-2.60, p=0.0103, just
outside the cutoff).

**High-PR (best translated):** `GO:0046982`/`GO:0046983` (protein
dimerization activity) over-represented; `GO:0008152` (metabolic process)
under-represented.

**GO-slim: 0 significant** (172 rows tested). Checked directly: p-values
range from ~0.016 to 1.0, not a suspicious uniform wall, a real absence of
slim-detectable signal at this threshold, not a broken run.

**q-value check**: smallest observed p-value (0.0025) carries
q=1.0 with 1,167 terms tested per class. A q<0.01 cutoff would erase all 7
real findings.

---

## Demo B: `rank_tail` on `TPM`, CR_3D alone

Same split convention, on `TPM` instead.

**High-TPM: a textbook ribosome/translation signature.**
`GO:0003735` (p=9.7e-15), `GO:0005840` (p=3.9e-14), `GO:0006412`/`GO:0160307`
(p=6.1e-14), `GO:0005198` (p=7.6e-12), plus broader organelle/gene-expression
terms. 24 of 1,167 terms tested reached significance for this class alone,
21 over-represented, 3 under-represented.

**Low-TPM: the mirror image.** `GO:0003723` (RNA binding) significant at
p=0.0064 (under-represented). `GO:0006412`/`GO:0160307` both near-miss
under-represented (p=0.0105).

**GO-slim: 5 significant**, contrasting directly with demo A's 0. Real
pileup near zero in the High class's p-value distribution, not only near
1.0.

---

## Demo C: `cluster`, k=4, `PR_gene`, CR_3D alone

Ward's-method hierarchical clustering (Yeo-Johnson -> Z-score -> Euclidean
-> Ward's), single column (`PR_gene`), single condition. Verifies the
clustering mechanism itself, not the original multi-condition use case
its lab precedent was designed for, see demos F/G for the genuine
multi-condition case.

**k selection**: elbow curve showed a sharp bend around k=3-5. The first
k tried, k=13 (borrowed from an unrelated precedent), was rejected after
its per-sample silhouette showed cluster 1 (470 genes) with a suspicious
*perfect* silhouette (mean=min=1.0) — confirmed to be exactly the 470
genes with `PR_gene == 0`, a tied-value artifact, not real structure.
k=4 was chosen instead: real elbow bend, cluster 1 (470 genes, the
zero-floor group) cleanly isolated on its own rather than fragmenting the
rest of the distribution into artifacts, cluster 3 (1,612 genes)
strongly separated (mean silhouette 0.65, no negative minimum).

Cluster sizes: 1=470, 2=431, 3=1612, 4=1080.

**Result, current and reproducible: 13 significant full-GO, 4 significant
slim.** An earlier version of this document reported 16/3 for this exact
command; that number could not be reproduced in a later session (cluster
sizes matched exactly, significant-term counts did not) and, after
confirming this session's own result is internally deterministic (byte-
identical on rerun, see "Determinism check" above), the earlier 16/3 is
treated as superseded rather than carried forward. The likely cause is a
library-version or environment difference between sessions, not a bug
in this codebase as currently verified — but the exact cause was not
tracked down, and the correct practice going forward is: **do not assume
a previously-documented clustering result still reproduces without
re-verifying it in the current environment.**

---

## Demo D: `rank_tail` on `PTPM`, CR_6D alone

First real run of any kind on CR_6D. Same convention as demos A/B, on
`PTPM` this time.

**High-PTPM**: the same ribosome/translation signature as demo B,
independently reproduced on a different condition (6 days vs 3) with a
different metric (polysome-fraction, not total-fraction expression).
`GO:0043229`/`GO:0043226` (organelle terms, p~5-8e-15), `GO:0003735`
(p=9.8e-14), `GO:0005198` (p=1.3e-13), `GO:0005840`/`GO:0160307`/`GO:0006412`
(p~3-4e-13). 21 over-represented, 3 under-represented of 24 significant.

**Low-PTPM**: mirror pattern. `GO:0005840`/`GO:0003735` both at zero
observed genes (p~0.10, near-miss, not significant at this threshold but
the same directional pattern as every other Low-expression group tested).

**GO-slim: 4 significant.**

---

## Demo E: `cluster`, k=3, `[PTPM, TPM]` jointly, CR_6D alone

First multi-column clustering demo, no merge needed since both columns
already live in one file.

**k selection**: elbow bend at k=3 (6,100 -> 2,350 -> 900), silhouette
peaked cleanly at k=3 (0.541), dropping at k=4 (0.445) with no competing
secondary peak at higher k. Both diagnostics agreed.

Cluster sizes: 1=1156, 2=689, 3=1212. Silhouette means: 0.63 / 0.67 / 0.39.

**Result, corrected an initial hypothesis.** Expected "two well-expressed
populations + one diffuse mass"; actual finding was one dominant
ribosome/translation axis split into three positions: cluster 1
under-represents it, cluster 2 over-represents it (and is the standout:
best-separated, min silhouette 0.15 is the only positive minimum of the
three, and the strongest, cleanest GO signal), cluster 3 under-represents
it from a third angle (`translation`, `protein biosynthetic process`
specifically). 21 significant full-GO, 4 slim, 0 full-GO/slim
discrepancies across all 3 classes.

---

## Demo F: `cluster`, k=2, `PR_gene`, CR_3D + CR_6D merged

**The real multi-condition use case** demo C explicitly could not test.
Built via `dataprep.merge_gene_tables()`, inner join on `gene_id`: 2,658
genes present in both conditions (of 3,593 in CR_3D, 3,057 in CR_6D).

**k selection was genuinely ambiguous, and the ambiguity was resolved by
the artifact check, not by picking whichever number looked better.**
Elbow showed no sharp single bend; silhouette had a global max at k=2
(0.503) and a smaller secondary peak at k=4-5 (~0.43). Checked whether
k=2's peak was a repeat of the zero-floor artifact: only 42 of the
resulting 332-gene cluster (12.6%) were near-zero-in-both-conditions
genes, ruling out a pure artifact. Per-sample silhouette: both clusters
positive mean (0.31, 0.53), moderate but real separation, weaker than
demo E's but not degenerate.

Cluster sizes: 1=332 (small, the ribosome/translation-enriched group),
2=2326.

**Result: 28 significant full-GO, 6 significant slim — the strongest
result in the whole worked example.** Cluster 1 over-represents the same
signature every other demo has found (`GO:0003735` p=5.2e-11,
`GO:0005840` p=1.4e-10, `GO:0006412`/`GO:0160307` p=3.9e-10), cluster 2
under-represents the identical term list, near-perfect mirror symmetry.
This is the first result requiring a gene's translational status to be
*consistent across two independent growth-stage conditions simultaneously*
to land in the enriched cluster, a stronger biological claim than any
single-condition result: not "looks this way once," but "stable across
conditions." 0 full-GO/slim discrepancies; the one near-threshold mismatch
(`GO:0003824`, catalytic activity, slim-significant/full-GO-not) is
explainable as the coarser slim vocabulary aggregating related full-GO
terms.

---

## Demo G: `cluster`, k=2, `TPM`, CR_3D + CR_6D merged

Same merge mechanism, `TPM` instead of `PR_gene`. 2,658 genes (same inner
join, confirmed identical merge behavior).

**k selection**: cleaner than demo F. Elbow bend sharply at k=2->3 (62%
then 37% drop); silhouette single clean peak at k=2 (0.511), no competing
secondary peak, monotonic decline after. Artifact check: 615 of 1,374
genes in cluster 1 (44.8%) are near-zero-in-both, a real but partial
fraction, not a pure artifact. Both clusters positive mean silhouette
(0.61, 0.41), the strongest separation of the two merged demos.

Cluster sizes: 1=1374 (low-expression), 2=1284 (high-expression).

**Result: 18 significant full-GO, 3 significant slim.** Cluster 1
strongly under-represents ribosome/translation (`GO:0160307`/`GO:0006412`
p=5.7e-5, `GO:0005840` p=1.7e-4, `GO:0003735` p=1.75e-4, the strongest
under-representation p-values for these terms across all seven demos),
cluster 2 over-represents the identical terms (`GO:0006412`/`GO:0160307`
p=0.0084). 0 full-GO/slim discrepancies.

---

## Reproducing these runs

Through the CLI (`scripts/run_pipeline.py`, full flag reference in the
main README section 6). Single-file demos (A-E) only differ in
`--metric-col`, `--label-strategy`/`--n-clusters`/`--pct`, and
`--dataset-name`; merged demos (F-G) additionally use `--merge-manifest`
in place of `--input-table`.

Example, demo F, using a merge manifest (`merge_pr_gene.txt`):

```ini
[cond_3D]
path = data/raw/CR_3D.gene_data.tsv.gz
value_col = PR_gene

[cond_6D]
path = data/raw/CR_6D.gene_data.tsv.gz
value_col = PR_gene
```

```bash
python scripts/run_pipeline.py \
    --merge-manifest merge_pr_gene.txt \
    --godb data/go_reference/CR/CR.godb.pkl \
    --slim-godb data/go_reference/CR/CR.slim.godb.pkl \
    --metric-col cond_3D cond_6D \
    --strip-id-suffix '\.v\d+\.\d+$' \
    --exclude-id 'gene:Standard-R-luc' \
    --label-strategy cluster --n-clusters 2 \
    --output-dir results/GO \
    --slim-output-dir results/GO_slim \
    --dataset-name CR_3D-CR_6D.PR_gene_cluster \
    --unknown-ratio-thresh 0.9 \
    --thresh-type p --thresh 0.01
```

All seven demos confirmed to reproduce the exact numbers in the gallery
table above.

## Open items specific to this example

- Only gene-level data has been used, both species-wide and per-condition.
  Transcript- or variant-level PR, and a second species, are unexercised.
- `rank_tail` and `cluster` have both been run extensively on real data
  (single-file, multi-column, and genuinely merged multi-condition cases).
  `explicit_threshold` and `boolean_flag` remain unverified outside
  synthetic tests, no natural real use case for either has presented
  itself on this dataset.
- `k` selection for the merged demos (F, G) required a short standalone
  script rather than `notebooks/select_cluster_k.ipynb`, which currently
  only supports single-file input (main README section 9).
- Ortholog-transfer supplementation and provenance-stamped output files
  are still open across the whole tool (main README section 9).