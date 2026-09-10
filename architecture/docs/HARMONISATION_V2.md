# Harmonisation v2 — the two-wave architecture

What Stage 0 became, why, and what it means for everything downstream.

Companion documents: [DUCKDB_LAYER.md](DUCKDB_LAYER.md) (the query layer),
[GLOSSARY.md](GLOSSARY.md) (what the quantities mean),
[REBUILD_RUNBOOK.md](REBUILD_RUNBOOK.md) (how to re-run it).

**Every number below was measured on this repository on 2026-08-03**, by running
the two notebooks end to end. None is quoted from a doc.

---

## 1. What changed

Stage 0 used to be one pandas notebook (`00_harmonisation.ipynb`, now archived at
`discarded/code/00_harmonisation_v1_pandas.ipynb`) that resolved 15 cleaned tables
onto `model_id` and wrote a single artefact, `outputs/harmonised.parquet`.

It is now **two notebooks and a database**:

```
        15 cleaned tables                    5 newer sources
        (DepMap/CCLE, GEO, HPA,              (GDSC, ProCan, COSMIC CNA,
         Cellosaurus)                         COSMIC CGC, HGNC)
                │                                     │
                ▼                                     ▼
   00_harmonisation.ipynb  ───────────►  00b_enriched_harmonisation.ipynb
   identity hub + gene dim              extends the hub with the Sanger axis,
   + coverage matrix                    resolves the newer sources onto it
                │                                     │
                └──────────────┬──────────────────────┘
                               ▼
               outputs/celllineselector.db   (DuckDB, 1.5 GB, 28 tables)
                               │
                               ├── harmonised.parquet          → Stage 1
                               ├── gene.parquet                → gene-axis consumers
                               ├── coverage_matrix.parquet     → Stage 2+
                               └── *_enriched.parquet          → the second wave
```

Three structural differences from v1:

| | v1 | v2 |
|---|---|---|
| storage | one parquet | DuckDB warehouse + parquet exports |
| gene axis | none — genes were never a first-class entity | `gene` dimension table, one row per ENSG |
| ambiguity | `model_id_set` columns, handled per-table | uniform flag → explode → unwrap, `is_ambiguous` on every table |

The parquet exports exist so **nothing downstream had to change**. Stage 1 reads
`harmonised.parquet` exactly as before; its column set is a strict superset of v1's.

---

## 2. Wave 1 — `00_harmonisation.ipynb`

Authored by Triveni; adapted here only for repo-relative paths and the parquet
exports in §9. The method is unchanged, and is worth stating because wave 2 reuses
it verbatim:

- **Flag, don't drop** — ambiguous identities are retained and marked
- **Credit both** — a row resolving to 2 candidate ACHs is duplicated, one row each
- **Grain before joining** — event tables are aggregated before any join
- **Identity hubs hold lists** — one row per entity, alternative ids as list columns,
  so `model_id` stays unique and no join inflates a count

### Measured output

| | v1 | v2 |
|---|---:|---:|
| ACHs in the identity hub | 1,840 | **2,127** |
| genes in the gene dimension | — | **20,163** |
| lines with all 9 modalities | — | **237** |

The hub grew because v2 bases it on the **union of every `model_id` seen anywhere**
rather than on `sample_info`. `sample_info` is authoritative but not complete —
`fusions` carries lines whose ACH ids post-date its release. `in_roster` records the
distinction rather than hiding it.

### The protein-coding filter

v2 restricts `depmap_expr` from 53,961 gene columns to 19,897, and `geo_expr` from
19,914 to 16,063, using a whitelist **inferred** by unioning `mutations`
(`vepbiotype == 'protein_coding'`) with `hpa_rna`. The notebook is explicit that this
is a proxy for a real biotype table.

Wave 2 checks that inference against HGNC — see §4.3. It holds up.

---

## 3. Wave 2 — `00b_enriched_harmonisation.ipynb`

Five sources that arrived after v1 was written and had never been harmonised:

| source | rows | what it adds |
|---|---:|---|
| `gdsc_models` | 2,266 | the SIDM ↔ ACH ↔ CVCL ↔ CCLE ↔ COSMIC bridge, plus tissue / cancer type / ploidy / MSI |
| `procan_proteomics` | 952 × 8,453 | DIA-MS proteomics |
| `cosmic_cna` | 168,411 | gene-level total copy number |
| `hgnc` | 19,273 | the authoritative protein-coding reference |
| `cosmic_cgc` | 768 | Cancer Gene Census role and tier |

### 3.1 The Sanger axis — why wave 2 needed its own hub extension

v1's hub holds CVCL, profile, GEO, CCLE and name axes. It has no **SIDM** (Sanger
model id) axis, because none of the 15 original tables carried one in usable form.

ProCan is keyed on SIDM and nothing else. Without a SIDM axis it can only be reached
by name — the weakest possible identifier. So wave 2 builds the axis first, from two
independent routes, and compares them:

```
route A  gdsc_models.broad_ach          1,770 pairs
route B  sample_info.sanger_model_id    1,153 pairs

  agreed by both routes   1,138
  route A only              632
  route B only               15
  ACHs where they disagree    0
```

Zero disagreements between Sanger's curation and DepMap's is a meaningful result, not
a trivial one — the machinery to record a conflict is in place (`sanger_id_conflict`)
and simply never fired.

The hub widens from 2,127 to **2,145** ACHs: 18 lines that only the newer sources know
about, added with `in_roster = False` and `source_wave = 'enriched'`.

### 3.2 The multi-path resolver — the one genuinely new mechanism

v1 resolved each table through **one** identifier axis. Each newer source carries
several at once, so wave 2 adds a priority-ordered chain — strongest identifier first,
falling back down the rungs, recording which rung actually fired in `matched_via`:

```
ach  →  cvcl  →  ccle_name  →  sidm  →  plain name
```

A match that came from a name is not as trustworthy as one from an accession, and
`matched_via` keeps that difference inspectable instead of flattening it into an
undifferentiated "matched". Measured:

| table | resolved | by what |
|---|---|---|
| `gdsc_models` | 1,781 / 2,266 (78.6%) | ach 1,770 · cvcl 11 |
| `procan_proteomics` | 948 / 948 (**100%**) | sidm 948 |
| `cosmic_cna` (per sample) | 996 / 1,013 (98.3%) | gdsc_name 974 · name 22 |

The 485 unresolved `gdsc_models` rows are models with no ACH at all — 256 of them are
organoids, which are **flagged (`is_organoid`), not dropped**.

---

## 4. What the second wave actually bought

### 4.1 Proteomics reach — the headline

| | lines | % of hub |
|---|---:|---:|
| CCLE/Gygi only | 84 | 3.9% |
| ProCan only | 661 | 30.8% |
| both platforms | 291 | 13.6% |
| **either** | **1,036** | **48.3%** |

Proteomic coverage goes from 375 to 1,036 lines — a **176% increase** — with 291 lines
measured on both platforms available for cross-platform comparison.

The two are stored as **separate layers and never concatenated**. That constraint is
not stylistic: `test_run_proteomics_platform_overlap.py` measured the batch effect
between the TMT and DIA-MS platforms before the merge was permitted, precisely because
Stage 2 percentile-ranks each protein across cell lines, and a naive merge would let a
protein look "high" merely because it was measured on the more sensitive instrument.

**Note on Stage 2.** `02_core_score.ipynb` Cell 7d already merges ProCan ad-hoc,
resolving 941 of 948 lines from the CSV directly. Wave 2 resolves 948 of 948 through
the hub. Stage 2 can now read `procan_proteomics` from the warehouse instead of
re-deriving the bridge, which would remove the 7-line discrepancy and one duplicated
resolver. That change is **not** made here — it belongs to a Stage 2 re-run.

### 4.2 Modality coverage

`coverage_matrix_enriched` is a strict superset of `coverage_matrix`: the nine
original modality columns carried through unchanged, plus `procan_proteomics` and
`cosmic_cna`, so `n_modalities` is now out of 11.

| | of 9 (core) | of 11 (enriched) |
|---|---:|---:|
| lines at full coverage | 237 | 195 |
| lines at ≥ 8 modalities | 577 | 548 |

### 4.3 Was the v1 protein-coding inference sound?

This is the check v1 explicitly deferred. HGNC is the rigorous source it stood in for:

```
Stage 0 inferred whitelist : 20,163 genes
HGNC protein-coding (ENSG) : 19,236 genes
  agreed by both           : 19,215
  inferred only            :    948   (HGNC says not protein-coding, or no ENSG)
  HGNC only                :     21   (real coding genes the inference missed)

precision vs HGNC : 95.3%
recall    vs HGNC : 99.9%
```

**The inference was sound.** It missed 21 genes out of 19,236 and admitted 948 extras.
`gene_enriched.is_protein_coding` now carries HGNC's verdict as its own column beside
the inference rather than replacing it, so Stage 2 can choose which definition to
filter on and the choice stays visible in the data.

### 4.4 Gene roles

581 of 20,163 genes carry a Cancer Gene Census role (257 oncogene, 254 TSG, 70 both).

Symbol matching goes through an alias-aware table built from HGNC's `symbol`,
`alias_symbol` and `prev_symbol` columns — 64,382 symbol→ENSG pairs. COSMIC still uses
names HGNC has retired, and without the alias table those genes silently fail to join;
10 census genes matched on a previous symbol and 1 on an alias.

Alias expansion cuts both ways, so the match is taken at the **strongest symbol tier
that fires** for each census gene, and anything still ambiguous at that tier is
dropped rather than picked from: 760 of 768 census symbols resolve to exactly one
ENSG, 8 have no HGNC route, 0 remain ambiguous.

The same discipline applies to `cosmic_cna`, which is symbol-keyed: of 17,605 distinct
symbols, 15,828 resolve to one ENSG, 474 are ambiguous and **left unresolved rather
than guessed**, 1,303 have no HGNC symbol at all. Guessing would put copy number on
the wrong gene, which is worse than leaving it unjoined.

---

## 5. Warehouse contents

`outputs/celllineselector.db` — 1.5 GB, 28 tables, regenerated by re-running the two
notebooks. Gitignored; the parquet exports are the portable interface.

| wave | tables |
|---|---|
| core (18) | the 15 source tables + `cell_line_connection`, `gene`, `coverage_matrix` |
| enriched (10) | `gdsc_models`, `procan_proteomics`, `procan_protein_map`, `cosmic_cna`, `cosmic_cgc`, `hgnc`, `hgnc_symbols`, `cell_line_connection_enriched`, `coverage_matrix_enriched`, `gene_enriched` |

No wave-1 table is ever overwritten by wave 2. Every wave-2 output is either a new
source table or a `*_enriched` variant sitting beside the original, so the two waves
stay separable and wave 1 can be re-run without wave 2 having to be re-run first.

---

## 5.1 Querying it — `build_warehouse_views.py`

```bash
python src/pipeline/build_warehouse_views.py --verify
```

Exposes the Stage 0-7 parquet outputs, `reference/`, `cleaned_track_data/` and
`data/parquet/data_clean/` as **views** inside the warehouse, so one database
answers every question. Measured: 28 materialised tables + 76 views across 4 schemas.

Views, not tables, because those files are already parquet on disk and
`core_score.parquet` alone would add ~1 GB if copied in. A view is re-bound on every
query, so a rebuilt parquet is picked up with no re-import.

Two convenience views sit on top:

- `v_core_score` — `core_score` with `gene_id` aliased for joining to the gene dimension
- `v_score_annotated` — score + gene identity + CGC role + modality coverage + GDSC
  tissue, in one object

```sql
SELECT model_id, hugo_symbol, gene_role, core_score, n_modalities, tissue
FROM v_score_annotated
WHERE ensg_id = 'ensg00000157764' AND n_layers = 2
ORDER BY core_score DESC LIMIT 5;
```

0.92 s over 29.8 M rows with three joins.

---

## 5.2 Key casing — one convention, lowercase

**Canonical key case across the whole project is lowercase**, matching the warehouse:
`ach-000219`, `ensg00000157764`.

It was not always so. Stage 0/0b/1 wrote lowercase while Stage 2-7 wrote UPPERCASE,
because `02_core_score.ipynb` carried an engineering rule reading *".str.upper()
model_id before every join"*. That was self-consistent within each half and broke at
the boundary — `explain_pair.py` had a function holding `pred` (upper) and `harm`
(lower) at once, with a comment explaining the discrepancy rather than fixing it, and
a cross-schema SQL join written the obvious way returned **zero rows, looking like
missing data rather than a key mismatch**.

Converted 2026-08-03 across 26 files: a 174-edit mechanical sweep flipping every key
normalisation from `.upper()` to `.lower()` and lowering the `ENSG…` / `ACH-…` literals
compared against them, plus 16 hand-written fixes where a read had **no**
normalisation at all and had been relying on every source happening to be uppercase.

**Why this was safe to do without a flag day.** Every one of those calls sits
immediately after a read, so `.lower()` is *idempotent*: it is a no-op on files
already written lowercase, and it silently upgrades the legacy UPPERCASE artefacts
that cannot currently be rebuilt. The Stage 7 CLI reads a stale UPPERCASE
`predictions_with_confidence.parquet` and works unchanged.

Verified behaviour-preserving: Stage 2 and Stage 4 both re-ran to **identical shapes**
(29,781,274 × 5 and 885,844 × 13) before and after.

Two latent bugs surfaced during the conversion, both pre-existing:

- **`ensg_id` was never normalised at all** on the Track-C reads in Stage 4
  (`mutations_collapsed`, `fusions_gene_level`, `mutations_scores`, `fusions_scores`).
  It only ever worked because every table happened to be uppercase — an accident, not
  a rule. Now explicit.
- **`cna_flags.parquet` is merged with `how="outer"`.** A key mismatch there does not
  join to nothing and fail visibly; it *appends* every CNA row as an unjoinable orphan
  and inflates the table. Row growth therefore cannot detect the fault, so the guard
  asserts on matched pairs instead. Measured: 5,216 of 117,540 CNA pairs (4.4%) match
  an existing flag row; the other 112,324 are CNA-only pairs the outer join is
  designed to keep.

---

## 5.3 Stale artefacts in `outputs/` — resolved 2026-08-04

Until the `gene_regime` producer existed, five files predated the v2 chain and
**could not be rebuilt**, because every route to them ran through the Stage 5
blocker. All five have now been regenerated; the table is kept as a record of what
"stale" looked like and how to spot it again.

| file | was | now |
|---|---|---|
| `predictions_with_confidence.parquet` | 28,403,110 rows (old `core_score` count) | 29,781,274, matches `core_score` |
| `evidence_ledger.parquet` | derived from the above | rebuilt |
| `chronos_validation.parquet` | 16,866 rows, old scores | rebuilt, 16,866 |
| `cna_flags.parquet` | needed `gene_regime` | rebuilt, 117,540 |
| `stage4_eval_consistent_denom.parquet` | metrics on the old score | see §6 |

**The tell to remember:** a row count of 28,403,110 anywhere means the *old*
`core_score`, and UPPERCASE keys mean an artefact written before 2026-08-03. Either
one identifies a stale file at a glance. A join between a stale artefact and the
current `core_score` does not error — it silently drops ~1.38 M pairs.

`outputs/variants/` (18 parquets, 3.5 GB, 2026-07-20, pre-ProCan, zero readers) was
**deleted** on 2026-08-04. It had been registered as 18 `out.variants_*` views by
`build_warehouse_views.py`, which made a stale scoring experiment look as legitimate
as current output; the script drops views whose backing file is gone, so re-running it
clears them.

---

## 5.4 `00c_harmonisation_eda.ipynb`

Read-only profiling of the warehouse along four axes — completeness, scale, shape
and agreement. Writes nothing; it is a companion to the two builders, not a stage.

Adapted from the handed-over version in two places:

- `DB_PATH` resolved to `outputs/celllineselector.db` rather than the CWD
- **`source_agreement` had a grain bug.** It aligned `depmap_expr` and `geo_expr` on
  `model_id` as though that were a unique index. It is not: `depmap_expr` is
  profile-grain (one model can have RNA, WES and WGS profiles) and both tables are
  exploded during harmonisation so an ambiguous line is credited to every candidate
  ACH. Aligning on a duplicated index fanned the join out and the per-gene masks
  stopped matching, raising `IndexError: boolean index did not match`. Both sides are
  now collapsed to one value per line (mean over replicate profiles) before
  correlating — *grain before joining*, applied where it had been skipped.

Headline result: **median Pearson r = 0.552** between DepMap and GEO per gene, over
the lines they share. Same axis, independent sources, moderate agreement — which is
the evidence for treating them as two columns rather than averaging them into one.

---

## 6. Downstream impact

### Stage 1 — verified working, no change needed

`01_lookup_metadata_join.ipynb` runs unmodified against the new
`harmonised.parquet` and writes `harmonised_enriched.parquet` at (2,127 × 26).

One number moved, for a legitimate reason:

| | v1 | v2 |
|---|---:|---:|
| `canonical_name` / `gap_class` match rate | ~98.8% of 1,840 | 85.5% of 2,127 |

The **absolute** match count is unchanged at 1,818. The rate fell because the hub now
carries 2,127 lines instead of 1,840, and `union_lookup_entity.parquet` — built when
the roster was DepMap-only — has no row for the 309 lines outside it. This is the
denominator changing, not a join breaking. The notebook's `<-- INVESTIGATE KEY
MISMATCH` warning fires on the rate and can be safely read past here.

`biotype` / `hgnc_status` remain at 0%, unchanged from v1: `harmonised.parquet` is
line-grain and has no `ensg_id`, so the gene-level join was always skipped. Wave 2 now
supplies exactly these fields properly on the gene axis, in `gene_enriched`.

### Stage 2 — re-run, passed

`02_core_score.ipynb` was re-run end to end against the v2 hub. A new **Cell 7e**
sources ProCan from the warehouse instead of `validation/prepared/id_bridge.parquet`
(a file with no builder in the repo), following the notebook's supersede-don't-modify
convention — Cell 7d is left intact and runs first, so both are comparable in one
session.

| | Cell 7d (id_bridge) | Cell 7e (warehouse) |
|---|---:|---:|
| ProCan lines resolved | 941 / 948 | **948 / 948** |
| protein-layer lines (union with CCLE) | 1,025 | **1,036** |

`core_score.parquet` rewritten: **29,781,274 rows** (was 28,403,110), keys lowercase. Two-layer cells
5,277,064; the correlation-penalised override moved that stratum's mean from 0.7371 to
0.5170 and cells above 0.95 from 22.0% to 2.0%, consistent with the pre-existing
Cell 8b correction.

### Stage 4 — re-run, passed after two fixes

`flags_with_driver.parquet` rebuilt at (885,844 × 13), keys lowercase. The second fix
was the `ensg_id` normalisation and CNA outer-merge guard described in §5.2.

Cell 2 loaded `layer3_continuous.parquet`, archived to `discarded/outputs/` on
2026-08-03, and raised `FileNotFoundError` before anything else ran. The notebook's
own markdown already established that this file "has no producer anywhere in this
repository" and added Cell 2b to reconstruct `flags_cont` from Track-C's reproducible
scored outputs — Cell 2b **reassigns `flags_cont` outright**, so the loaded value was
always discarded. Only the read itself was load-bearing, and only in the sense that
its absence stopped the notebook.

The read is now conditional rather than deleted, so the original data path stays
visible and the file is still used if it is ever restored.

### Stage 5 onwards — unblocked 2026-08-04

`gene_regime.parquet` (141 genes) had no producer anywhere in the repo. It is now
built by `src/pipeline/build_gene_regime.py` — see that module's docstring for what
is *certain* about the derivation (the gene set and `n_sensitive_gdsc` reproduce
exactly, 141/141) and what is *reconstructed* (the threshold).

The full chain now runs end to end:

| stage | output | rows |
|---|---|---:|
| 02 `core_score` | `core_score.parquet` | 29,781,274 |
| `build_gene_regime` | `gene_regime.parquet` | 141 |
| `build_cna_layer` | `cna_flags.parquet` | 117,540 |
| 04 driver routing | `flags_with_driver.parquet` | 885,844 |
| 05 confidence tiers | `predictions_with_confidence.parquet` (curated) | 223,678 |
| `build_full_predictions` ×2 | `predictions_with_confidence.parquet` (all genes) | 29,781,274 |
| `build_chronos_validation` | `chronos_validation.parquet` | 16,866 |
| `build_evidence_ledger` | `evidence_ledger.parquet` | 29,781,274 |
| 06 / `stage4_eval_save` | `stage4_eval_*` | 28 |
| `compute_stats` + 4 figure scripts | `_stats.json`, 17 figures | — |

**Ordering is load-bearing**: `build_cna_layer` reads `gene_regime`, and Stage 4 reads
`cna_flags`. Running Stage 4 first silently uses the previous CNA layer.

Held-out evaluation (28 test genes, 26 activation-driven): driver gating lifts hit@20
from 0.01402 to 0.01773 (**+26.5%**) and cuts the per-gene spread by 26.7%.
`abundance_tracking` shows delta exactly 0.000, which is correct by construction —
those genes always rank on `core_score` alone, so the flag cannot move them.

### A second orphan, same shape as the first

`outputs/coverage_tables/` — `coverage_per_cell_line.parquet` (1,485 rows) and
`coverage_per_gene.parquet` (19,176) — is read in four places by
`src/figures/compute_stats.py` and **written by nothing**. Both files date from
2026-07-26 and are derived from the old `core_score`.

The staleness is visible in the published statistics: `scored_cell_lines` reports
**1,485** while `core_score` now covers **1,746**. Any figure quoting the number of
scored cell lines is understating it by 261 lines. This is the same failure mode
`gene_regime` had and it has not been fixed here — writing a producer means
re-deriving the schema's intent, which is a decision for the author.
