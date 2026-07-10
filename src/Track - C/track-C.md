# Track C — Alterations (GeneTraceAI)

**Owner:** Fiona (Chaithali) — Track C: Mutations, Fusions, Signatures, miRNA
**Repo:** UoB-AstraZeneca-25-26 · Branch: `Chai`
**Source notebooks:** `01_mutations_cleaning.ipynb`, `02_fusions_cleaning.ipynb`, `03_signatures_cleaning.ipynb`, `04_mirna_cleaning.ipynb`
**Status:** All four datasets cleaned and written. Track C complete.

This document is verified against the final, executed notebooks — numbers and logic below reflect what was actually run, not the exploratory profiling that preceded each one.

---

## 1. Overview

Track C covers alteration-type evidence for the GeneTraceAI composite scorer: somatic mutations, gene fusions, genome-wide instability signatures, and miRNA expression. All four were cleaned independently at native grain and joined to the shared reference tables owned by Track A (Thiruvel) — `gene_lookup.parquet` and `cell_line_lookup.parquet` — following the project's hub-and-spoke architecture. No physical merging across tracks; datasets stay at native grain and join on demand via `ensg_id` and `model_id`.

Each notebook follows the same six-step recipe: **profile → decide → resolve identifiers → attach keys → log unmapped → write and test.**

### Summary table

| Dataset | Grain | Final rows | Genes | Models | Gene axis? | Status |
|---|---|---|---|---|---|---|
| Mutations | `(ensg_id, model_id)` | 636,568 | 18,122 | 1,744 | Yes | ✅ Complete |
| Fusions | `(ensg_id, model_id)` | 145,842 | 16,198 | 1,699 | Yes | ✅ Complete |
| Signatures | `model_id` only | 1,955 | — | 1,955 | **No** | ✅ Complete |
| miRNA | `(mirna_id, model_id)` | 695,832 | — | 948 | **No** | ✅ Complete |

Two of four Track C datasets (signatures, miRNA) have **no gene axis** and cannot join to `gene_lookup` — this is a structural finding, not a cleaning gap, and directly affects the scoring architecture (§6.1).

---

## 2. Mutations

**Notebook:** `01_mutations_cleaning.ipynb`
**Output:** `mutations_collapsed.parquet`, `mutations_variant_detail.parquet`, `mutations_resolution_log.csv`, `mutations_out_of_universe.csv`
**Grain:** `(ensg_id, model_id)` — 636,568 rows · 18,122 genes · 1,744 models

### Pipeline

1. **Load + column pruning.** Raw `mutations.csv` (1,066,869 rows × 70 columns) had 22 low-value annotation columns dropped up front (GTEx, DIDA, PharmGKB, BRCA1 functional score, GWAS disease/PMID, CIViC description/ID/score, Hess signature, dbSNP rsID, rescue reason, VEP clinsig, intron, molecular consequence, NMD, PS, and the LOF-transcript detail fields). These are third-party annotation-database fields, not primary evidence — kept the file lean before doing any join work.
2. **Bridge resolution.** `ensemblgeneid` uppercased directly (source ships clean Ensembl IDs, no symbol fallback needed here — unlike fusions). `profileid` merged against `depmap_profiles.csv` on `profileid → modelid`, then `model_id` uppercased.
3. **Variant dedup.** A `variant_key` (`chrom:pos:ref:alt`) was built, then rows deduped on `(ensg_id, model_id, variant_key)` — 913,984 → 715,123 rows (198,861 removed), collapsing repeated annotation rows for the same physical variant (e.g. multi-transcript VEP output) before any aggregation.
4. **Per-variant signal engineering.** `vep_rank` (ordinal: modifier < low < moderate < high), `path_max_variant` (max of AlphaMissense and REVEL scores), `is_high_mod` (VEP rank ≥ moderate).
5. **Collapse to `(ensg_id, model_id)`.** Engineered fields: `max_vep_rank`, `max_pathogenicity`, `variant_count`, `multi_hit_high_impact` (≥2 high/moderate-impact variants), `oncogene_hit` / `tsg_hit` (direction-preserving, from `oncogenehighimpact` / `tumorsuppressorhighimpact`), `max_vaf`, `variant_ids` (list), `any_driver` (any of 5 driver-flag columns true), `variant_burden` (`log1p(variant_count)`). Cell line name attached from `cell_line_lookup` for readability.

### Decisions locked

| # | Decision | Resolution | Why | Effect downstream |
|---|---|---|---|---|
| 1 | Drop 22 annotation columns up front | Removed before any join work | These are external-database cross-references, not primary evidence signals for scoring | Smaller, faster pipeline; no loss of scorer-relevant information |
| 2 | Bridge via `profileid → modelid`, not direct `model_id` in source | Left-merge against `depmap_profiles.csv`, unmatched rows logged rather than dropped silently | Raw mutation rows are keyed by sequencing profile, not model, and DepMap models can have profiles that don't resolve (deprecated/superseded profiles) | 152,885 rows (14.3% of raw) lost `model_id` and were excluded from scoring — logged in `mutations_resolution_log.csv` as `unmatched_profile` entries. This is a real, sizeable gap worth stating explicitly in the report rather than letting it hide in a row-count delta. |
| 3 | Variant-key dedup before collapse | Drop duplicates on `(ensg_id, model_id, variant_key)` | Same physical variant can appear multiple times per transcript/annotation; without this, `variant_count` and `variant_burden` would be inflated by annotation multiplicity rather than true biological burden | Ensures `variant_burden` reflects genuinely distinct mutations, not annotation artefacts |
| 4 | `variant_burden = log1p(variant_count)`, driver flags retained (not dropped for near-zero variance) | Confirmed prior architectural principle | Log-transform controls for skew in variant count; class-weighting at training time is the correct fix for imbalanced driver flags, not feature removal | Keeps a defensible, low-variance-but-real signal available to the scorer rather than pre-emptively discarding it |
| 5 | No symbol→Ensembl fallback (unlike fusions) | `ensemblgeneid` used directly, uppercased | Source mutation calls already ship clean Ensembl gene IDs — no `.`-placeholder problem like fusions had | Simpler pipeline than fusions' identifier resolution step |
| 6 | **Universe filter added (patched after initial review)** | `ensg_id` validated against `gene_lookup`; out-of-universe rows logged to `mutations_out_of_universe.csv` and dropped from both `mutations_collapsed.parquet` and `mutations_variant_detail.parquet` | Brings mutations in line with the fusions pipeline; without it, 130 genes' worth of drift accumulates silently and risks join failures against `gene_lookup`/Neo4j | 642,491 → 636,568 rows (0.92% dropped), 19,576 → 18,122 genes — consistent, defensible universe scoping across both alteration datasets |

### ✅ Resolved: universe filter added, matching the fusions pattern

The gap flagged in the initial notebook review (19,576 genes vs. the 19,446-gene universe) has been fixed. A universe-validation step, mirroring fusions' decision #5, was added between the collapse and write-output stages: `ensg_id` checked against `gene_lookup`, out-of-universe rows logged to `mutations_out_of_universe.csv` and dropped from both `mutations_collapsed.parquet` and `mutations_variant_detail.parquet` (the detail table shares the same join key and needed the same filter or it would carry genes the collapsed table no longer has).

**Result:** 642,491 → 636,568 rows (5,923 dropped, 0.92%), 19,576 → 18,122 unique genes (1,454 genes removed). The row-level impact is small relative to the gene-level impact — the dropped genes were thin, low-mutation-count entries rather than systematically important ones, consistent with non-canonical/non-protein-coding loci being swept up by the caller the same way DepMap's fusion caller was.

Mutations now lands at 18,122 genes, **1,324 short of** the 19,446 ceiling (unlike fusions, which landed close to it) — expected, since many protein-coding genes are simply never mutated across the 1,744 profiled models, not a data quality signal.

---

## 3. Fusions

**Notebook:** `02_fusions_cleaning.ipynb`
**Output:** `fusions_gene_level.parquet`
**Grain:** `(ensg_id, model_id)` — 145,842 rows · 16,198 genes · 1,699 models

### Pipeline (three stages, matching the notebook's own structure)

1. **Event-level dedup + identifier resolution** — filter to `isdefaultentryformodel == 'yes'`, uppercase `model_id` and Ensembl IDs, symbol fallback via `gene_lookup['hgnc_symbol']` for missing IDs, dedupe repeated breakpoint calls down to one row per `(model_id, canonicalfusionname)` event (keeping the highest-`ffpm`, then highest-confidence call).
2. **Explode + collapse to gene grain** — each fusion event split into its two gene participants (`side1`/`side2`), concatenated, collapsed via `groupby(['ensg_id','modelid']).apply(collapse_group, include_groups=False)` into `fusion_count`, `max_confidence`, `best_ffpm`, `any_in_frame`, `top_partner_ensg_id`, `top_partner_gene`.
3. **Universe validation** — only `ensg_id` (the join key) is checked against `gene_lookup`; `top_partner_ensg_id` is deliberately **not** filtered, since it's explanation-layer metadata, not a join key, and a partner gene outside the universe is still legitimate biological signal.

### Decisions locked

| # | Decision | Resolution | Why | Effect downstream |
|---|---|---|---|---|
| 1 | Filter to `isdefaultentryformodel == 'yes'` before any collapse | Applied first | Confirmed clean: 178,041 rows, exactly 1 `sequencingid` per model post-filter | Prevents double-counting fusion burden for models with multiple sequencing profiles |
| 2 | Column inclusion | Kept `ffpm`, `confidence`, `reading_frame`; dropped all breakpoint/coverage/strand/tag columns | `confidence` is monotonic with `ffpm` (verified: high mean 0.76 → medium 0.39 → low 0.08) — it's a validated signal, not noise; the dropped columns are its technical inputs | Leaner evidence table; nothing lost that isn't already summarized in `confidence`/`ffpm` |
| 3 | Two-stage collapse (event dedup → gene explosion) | Breakpoint multiplicity per gene-pair-per-model deduped to one row per event before exploding to gene level | Median 1 breakpoint call per gene pair, max 156 — re-detections of the same event, not independent evidence | Without stage 1, `fusion_count` would be wildly inflated for genes with many re-detected (not re-occurring) fusion calls |
| 4 | Symbol fallback + resolution log | `gene1_ens_id`/`gene2_ens_id` case-fixed; missing IDs (`.` placeholder) resolved via `gene_lookup['hgnc_symbol']`; still-unmapped logged | Fallback cut unmapped from 7.6%/16.1% to 4.57%/10.87%; remainder confirmed as lncRNAs, pseudogenes, small RNAs, and viral integration coordinates — correctly outside the protein-coding universe | Clean, defensible unmapped tail; nothing silently dropped without a paper trail |
| 5 | **Universe filter on `ensg_id` only, not on `top_partner_ensg_id`** | Explicitly scoped in the notebook's own markdown | The join key must be universe-valid or downstream joins silently fail; the partner is metadata for the SHAP explanation layer, where an out-of-universe partner is real signal worth keeping, not a defect | 197,943 → 145,842 rows (26.32% dropped, 9,424 genes removed), landing at 16,198 — safely under the 19,446 ceiling. This is the pattern mutations should also apply (see §2's open item). |

### Defensibility note
26.32% row loss is large in isolation, but it's the same underlying story as `gene_lookup`'s own construction gap, at higher volume: DepMap's fusion caller (and viral integration detection specifically) is genome-wide, not protein-coding-scoped.

### Resolved during this pass
`clean_fusions.py:132` `FutureWarning` — the notebook already uses `include_groups=False` in the `groupby().apply()` call. Confirmed patched.

---

## 4. Signatures

**Notebook:** `03_signatures_cleaning.ipynb`
**Output:** `signatures_model_level.parquet`
**Grain:** `model_id` only — **no gene axis** — 1,955 rows · 1,955 models

### Fields
`model_id`, `msiscore`, `lohfraction`, `wgd` (nullable boolean), `cin`, `ploidy`, `aneuploidy` — genome-wide instability phenotypes, not SBS-style mutational signature decomposition.

### Decisions locked

| # | Decision | Resolution | Why | Effect downstream |
|---|---|---|---|---|
| 1 | Filter to `isdefaultentryformodel == 'yes'` | Applied first | Confirmed exact 1:1 `model_id` grain post-filter (1,955 rows = 1,955 models), no collapse math needed | Simplest pipeline of the four — filter, rename, cast, write |
| 2 | Column inclusion | Kept `model_id` + 6 metrics only | Everything else was sequencing/profile-technical metadata | Compact context table |
| 3 | **Confirmed no `ensg_id` axis exists** | No gene column anywhere in source | Signatures are computed per cell line as a whole, not per gene | **Cannot join to `gene_lookup` or feed the per-gene scorer directly** — see §6.1 for the open architectural decision this creates |
| 4 | Missing values left as NaN | 333 rows (17.0%) null across all 5 WGS-derived columns as a clean all-or-nothing block; `msiscore` has zero nulls | Structural gap (insufficient sequencing depth/assay type for WGS-derived metrics), not missing-at-random | Imputation would fabricate biological signal; NaN correctly propagates "not measured" downstream |
| 5 | `wgd` cast to nullable boolean (`'boolean'` dtype) | Explicit cast, not left as float | Preserves the binary semantic while keeping NaN representable | Cleaner type for any downstream filtering (`wgd == True`) without float comparison quirks |

### Open architectural question (unchanged from prior review)
Because signatures has no gene axis, its role in the composite score is still an open team decision: **(a)** post-hoc cell-line-level modifier, **(b)** pure Neo4j metadata, or **(c)** excluded from v1 scoring, reserved for the explanation layer.

---

## 5. miRNA

**Notebook:** `04_mirna_cleaning.ipynb`
**Output:** `mirna_model_level.parquet`
**Grain:** `(mirna_id, model_id)` — **no gene axis** (Option B confirmed) — 695,832 rows · 734 miRNAs · 948 models

### Pipeline

1. **Column → `model_id` resolution.** Source is wide (734 miRNA rows × 954 cell-line columns, headers like `mcf7_breast`). Built a normalized lookup (`re.sub(r'[^a-z0-9]', '', ...)`, lowercased) from `cell_line_lookup['cell_line_name']` and `['stripped_cell_line_name']`, **deduped by `model_id` into a `set`** before counting matches. Column headers matched via iterative tissue-suffix stripping (trying 0, 1, 2… trailing underscore segments) since tissue suffixes are themselves multi-word (`haematopoietic_and_lymphoid_tissue`).
2. **Manual override + drop holdouts.** `d341med_central_nervous_system` mapped by hand to `ACH-000095` (confirmed via `stripped_cell_line_name` match to `D341`). `ncih1339_lung` and `colo699_lung` confirmed absent from the lookup after checking **all three** available name fields — `cell_line_name`, `stripped_cell_line_name`, and `synonyms`. Four genuinely ambiguous columns (`tt_oesophagus`, `tt_thyroid`, `u251mg_central_nervous_system`, `kmh2_haematopoietic_and_lymphoid_tissue`) dropped rather than guessed, since `cell_line_lookup` has no tissue/lineage field to break the tie.
3. **Melt to long format.** `(mirna_id, mirna_symbol, model_id, expression_value)`, verified zero duplicate `(mirna_id, model_id)` pairs.

### Decisions locked

| # | Decision | Resolution | Why | Effect downstream |
|---|---|---|---|---|
| 1 | **Grain / architecture confirmed as Option B** | No gene mapping columns anywhere in source | miRNA-to-target-gene mapping is a separate concern (e.g. via miRTarBase/TargetScan), not present in this expression matrix | miRNA is its own entity type in the knowledge graph, keyed on `(mirna_id, model_id)`, not exploded onto protein-coding genes — same structural category as signatures |
| 2 | Deduped-by-`model_id` matching (not raw match count) | Lookup built as `dict[key] = set(model_ids)`, not a list | An earlier version of this logic iterated two near-identical name fields per row without deduping, which double-counted matches and misclassified 928/954 columns as ambiguous | Caught before shipping by comparing against a simpler prefix-only match that returned 951/954 clean resolutions — a concrete example of validating derived output against a second method rather than trusting the first result |
| 3 | Iterative tissue-suffix stripping (0 to 5 segments) | Try progressively longer strips until a lookup hit | Naive last-underscore splitting fails on multi-word tissue names (`large_intestine`, `central_nervous_system`) | Without this, dozens of legitimately resolvable columns would have been wrongly dropped as unmatched |
| 4 | 6 of 954 columns dropped (0.63%) | 1 manual override, 2 genuinely absent from lookup, 4 real name collisions (no tissue field to disambiguate) | Checked all 3 name fields (including `synonyms`) before concluding absence; ambiguous ties left unresolved rather than guessed, since a wrong guess here silently corrupts expression data for the wrong cell line | 948 of 954 cell lines retained (99.4%) — small, fully logged, fully explained gap |
| 5 | Melt to long, verify zero duplicate `(mirna_id, model_id)` pairs | Explicit post-melt check | Confirms the column-to-model mapping was injective — a silent many-to-one mapping here would have meant two different source columns overwriting each other's expression values for the same model | Confidence that every row in the final output is a genuinely distinct, correctly-attributed observation |

### Note on `cell_line_lookup.synonyms`
This field exists in the lookup table and **was** used — for the manual verification of the 2 genuinely-absent holdouts — but is **not** wired into the automated matching pipeline itself (only `cell_line_name` and `stripped_cell_line_name` feed `lookup_norm`). Worth considering for Track A as a future enhancement to the automated matching path, since it may resolve some of the current ambiguous/no-match cases without manual intervention.

---

## 6. Cross-cutting findings

### 6.1 Not all Track C datasets have a gene axis
Mutations and fusions produce `(ensg_id, model_id)` grain and feed the per-gene composite score directly. **Signatures (`model_id` only) and miRNA (`mirna_id, model_id`) do not**, and cannot join to `gene_lookup`. Both need a team-level decision on scoring role — post-hoc modifier, Neo4j-only metadata, or excluded from v1/reserved for the explanation layer. Flagged for Week 5 team discussion alongside the flat-merge vs. hub-and-spoke resolution and the P1/P2/P3 scoring proposition decision.

### 6.2 Universe-filtering is now applied consistently across both gene-axis datasets
Fusions validates `ensg_id` against the 19,446-gene `gene_lookup` universe (dropping 26.32% of rows) and documents why the partner-gene column is deliberately exempt. Mutations initially skipped this check (19,576 genes, 130 over the universe) — patched after review to bring it in line, dropping 0.92% of rows and 1,454 genes (18,122 final). The underlying cause is the same in both datasets: DepMap callers operate genome-wide, not scoped to protein-coding genes, so any dataset built from a DepMap caller should be assumed to need this check unless proven otherwise.

### 6.3 Identifier resolution difficulty scaled with dataset structure
- **Signatures**: no resolution needed beyond `model_id` casing (no gene axis at all)
- **Mutations**: `profileid → model_id` bridge merge, with a real 14.3% unmatched-row gap logged rather than hidden, and no symbol fallback needed (Ensembl IDs shipped clean)
- **Fusions**: symbol→Ensembl fallback for a `.`-placeholder gap, plus universe validation
- **miRNA**: hardest — free-text column headers required multi-segment tissue-suffix stripping, and a real matching-logic bug (duplicate-field double-counting) was caught and fixed before trusting the output

This escalation is a strong narrative thread for the reflective account, and the miRNA matching-bug catch in particular is a concrete, citable example of validating tool/pipeline output against a second method rather than accepting first-pass results.

### 6.4 Structural nulls handled consistently
Signatures' 17.0% null block (WGS-derived metrics absent for models without adequate sequencing depth) and mutations' variant-level predictor nulls are both coverage/type-conditional, not missing-at-random. Both left as NaN rather than imputed, consistent with the project-wide principle that imputation would fabricate signal.

---

## 7. Open items carried forward from Track C

| Item | Owner | Status |
|---|---|---|
| ~~Mutations gene universe not validated against `gene_lookup`~~ | Fiona | **Closed** — universe filter added; 18,122 genes final, 636,568 rows |
| Signatures scoring role (post-hoc modifier / Neo4j metadata / excluded) | Team | Open — Week 5 |
| Mutations bridge: 152,885 rows (14.3%) lost `model_id` via unmatched `profileid` | Fiona / Thiruvel | Logged, worth a one-line note in the report on scope |
| `cell_line_lookup` missing lineage/tissue field (blocked 4 miRNA columns from disambiguation) | Thiruvel | Flagged |
| `cell_line_lookup.synonyms` not wired into automated miRNA matching | Thiruvel / Fiona | Noted, low priority |
| Flat-merge vs. hub-and-spoke resolution | Team | Open — Week 5 |
| Scoring proposition (P1/P2/P3) | Team | Open |

---

## 8. Outputs summary

```
Track - C/outputs/
├── mutations/
│   ├── mutations_collapsed.parquet
│   ├── mutations_variant_detail.parquet
│   ├── mutations_resolution_log.csv
│   └── mutations_out_of_universe.csv
├── fusions/
│   ├── fusions_gene_level.parquet
│   ├── fusions_unmapped_gene1.csv
│   ├── fusions_unmapped_gene2.csv
│   └── fusions_out_of_universe.csv
├── signatures/
│   └── signatures_model_level.parquet
└── mirna/
    ├── mirna_model_level.parquet
    └── mirna_unresolved_columns.csv
```

**Track C status: complete.** All four datasets are cleaned, resolved to canonical keys where a gene axis exists, and consistently validated against the protein-coding gene universe where a gene axis applies. Documented for handoff into the composite scoring layer once the remaining open architectural items (§7) — all team-level, not Track C-internal — are resolved.

[TRACK-C][MUTATIONS] Fix: strip versioned ENSG IDs at source before universe filter

Inserted strip cell at position [4] — before variant dedup, collapse, and
detail output. Removes ENSG version suffixes (e.g. ENSG00000157764.14 →
ENSG00000157764) so the universe filter join operates on clean base IDs.

Result: 64 non-protein-coding genes correctly excluded (18,122 → 18,058);
out-of-universe row count rises by 3,649 as expected (~57 variants/gene × 64
genes now correctly filtered). mutations_collapsed.parquet and
mutations_variant_detail.parquet both rebuilt from clean IDs.

[TRACK-C][FUSIONS] Fix: strip versioned ENSG IDs and redo collapse in single cell

Strip cell at position [6] strips exploded['ensg_id'] then immediately redoes
fusions_clean via groupby — initial fix had groupby before strip (silent no-op);
caught from log output and corrected before commit.

Result: 181 non-protein-coding genes correctly excluded (16,198 → 16,017);
32 versioned ID pairs merged into base form before filter. fusions_gene_level.parquet
rebuilt from clean IDs.