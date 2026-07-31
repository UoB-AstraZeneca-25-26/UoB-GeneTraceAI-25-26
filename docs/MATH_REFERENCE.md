# GeneTraceAI — Mathematics Reference

Catalog of every formula, statistical test, and decision rule used in the GeneTraceAI
pipeline: what it is, where it lives, why it's there, and how its output is consumed
downstream. Organized by pipeline stage, in execution order. Line numbers are
approximate (notebooks renumber cells on edit) but function/file names are exact.

Paths are relative to the repo root (`C:\Disertation\UoB-GeneTraceAI-25-26`).

---

## Stage 0 — Harmonisation

`src/pipeline/00_harmonisation.ipynb`, `src/new-arch/00_data_harmonisation_pipeline.ipynb`,
`src/Track - B/TrackB_final_resolution_pipeline.ipynb`

No scoring math at this stage — it's identity resolution (joining cell-line/gene
identifiers across sources) plus data-quality auditing. The only quantitative
operations are diagnostic counts and percentages, listed here because they gate
decisions made later (which normalisation scheme to use, which columns to drop).

**Missing-data audit**
```
missing_pct = df[col].isna().mean() * 100
```
- Why: surfaces which columns/tables have severe gaps before choosing a per-modality
  normalisation strategy. Columns at 100% missing are dropped.
- How: computed per (table, column) across all 15 pipeline tables, sorted descending.

**Numeric scale summary** — min/max/mean/median/std per numeric column.
- Why: "useful for spotting which columns are raw counts, log-ratios, z-scored,
  percentages, etc. before deciding on a per-modality normalisation strategy"
  (informs the percentile-vs-zscore choice in Stage 2).

**Modality coverage count**
```
n_modalities = sum(coverage_bool_columns, axis=1)
complete_cycle = (n_modalities == 8)
```
- Why: quantifies how many of the 8 modalities (expression, mutation, fusion,
  protein, HPA-RNA, metabolomics, miRNA, signatures) each cell line has. Diagnostic
  only — not fed into any score.

**Track B identity-bridge QA**: set intersection/union counts and match-rate
percentages (`mapped / before * 100`) at each join step of the RRID↔CVCL↔ACH↔profile
crosswalk. Gates whether a row is accepted, logged as "disputed," or "orphaned" —
no numeric scoring.

**Name matching**: no fuzzy/edit-distance similarity score anywhere in this
pipeline. Matching is either exact-match after deterministic string normalisation
(`normalise_base()` in Track B — lowercase, strip brackets/hyphens/underscores via
regex) or substring containment (`fuzzy_model_id_for()` in `02_core_score.ipynb`,
using `.str.contains()`). "Fuzzy" here means *permissive substring match*, not a
Levenshtein/Jaccard score.

---

## Stage 1 — Lookup / Metadata Join

`src/pipeline/01_lookup_metadata_join.ipynb`

Pure additive left-join (no scoring). One QA gate:
```
pct = 100 * matched / len(enriched)
flag "<-- INVESTIGATE KEY MISMATCH" if pct < 90
```
- Why: "Low match rates = key mismatch (case, version suffix, wrong grain), not
  genuinely missing metadata." Below 90% signals a join-key bug, not a data gap.

---

## Stage 2 — Core Score

`src/pipeline/02_core_score.ipynb`, `src/pipeline/scoring_variants.py`

This is the mathematical center of the pipeline: turning raw expression and
proteomics measurements into one comparable `core_score` per (gene, cell line).

### 2.1 Percentile-rank normalisation

```
E = expr.rank(axis=0, pct=True, na_option="keep")   # per gene, across cell lines
P = prot.rank(axis=0, pct=True, na_option="keep")
```
*Where:* `scoring_variants.normalise(method="percentile")` ([scoring_variants.py:30-31](src/pipeline/scoring_variants.py)); `percentile_rank()` in `02_core_score.ipynb`.

*Why:* expression (log2 TPM) and proteomics (log-ratio) live on different, non-comparable
native scales. Percentile rank maps both onto [0,1] using only within-gene ordering,
which is robust to outliers and distributional differences between the two assays.
An internal audit found min-max scaling ("2.5 Alternative variants" below) to be the
weakest normalisation by comparison.

*How:* applied independently to the deduplicated expression matrix and the
proteomics matrix (both `model_id × ensg_id` wide form) before combination.

### 2.2 Duplicate-profile deduplication

```
expr_dedup = expr.groupby(level=0).mean()   # collapse duplicate model_id rows
```
*Where:* `dedup_mean()`, `02_core_score.ipynb`.

*Why:* 16 ACH cell lines have two RNA profiles each. Mean is chosen because it's
"the least lossy option (keeps magnitude; max would inflate, drop would lose data)."

### 2.3 Noisy-OR combination (superseded — kept for context)

```
core_score = 1 - (1 - E)(1 - P)      when both E and P present
core_score = E                        when only E present
core_score = P                        when only P present
core_score = NaN                      when neither present
```
*Where:* `scoring_variants.combine(method="noisy_or")` ([scoring_variants.py:59-60](src/pipeline/scoring_variants.py)); original `02_core_score.ipynb` cell.

*Why:* treats expression and proteomics as independent evidence sources for
"protein present/abundant." Noisy-OR is the standard way to combine independent
probabilistic evidence sources into one combined probability of at least one
being true.

*Problem:* mRNA and protein abundance are **not** independent (ρ ≈ 0.46–0.58 per
Nusinow et al. 2020), so Noisy-OR systematically inflates `core_score` wherever
both layers agree — it double-counts correlated evidence. Measured impact on
production data: **22.5% of two-layer rows scored > 0.95** under Noisy-OR, and
**379 rows were pinned at exactly 1.0**. This motivated 2.4.

### 2.4 Correlation-penalised weighted sum (production formula)

The fix: weight each layer inversely to how correlated it is with the others, so a
redundant (highly-correlated) layer is down-weighted rather than allowed to inflate
the combination.

**Two-layer form:**
```
w_E = 1 / (1 + |ρ_EP|)
w_P = 1 / (1 + |ρ_EP|)          (symmetric in the 2-layer case)
w_E', w_P' = normalise(w_E, w_P) so w_E' + w_P' = 1

core_score = w_E' · E + w_P' · P
```

**Three-layer form** (`combine_three`, with expression E, proteomics P, Chronos C):
```
mean_ρ_E = (|ρ_EP| + |ρ_EC|) / 2
mean_ρ_P = (|ρ_EP| + |ρ_PC|) / 2
mean_ρ_C = (|ρ_EC| + |ρ_PC|) / 2

w_E_raw = 1 / (1 + mean_ρ_E)
w_P_raw = 1 / (1 + mean_ρ_P)
w_C_raw = 1 / (1 + mean_ρ_C)

w_E, w_P, w_C = normalise(w_E_raw, w_P_raw, w_C_raw)     # sum to 1

core_score = w_E · E + w_P · P + w_C · C     (all three layers present)
```
2-layer and 1-layer subsets of the 3-layer function re-normalise weights over
whichever layers are present (e.g. `w_EP_E = w_E_raw / (w_E_raw + w_P_raw)`).

Default empirical correlations used as constants: `ρ_EP = 0.46` (expression–protein,
prior literature/validation), `ρ_EC = ρ_PC = 0.005` (expression/protein–Chronos,
measured directly in this project, see Stage 2.7).

*Where:* `scoring_variants.combine_three()` ([scoring_variants.py:75-136](src/pipeline/scoring_variants.py)) — canonical implementation; applied as a post-hoc override in `02_core_score.ipynb`.

*Why:* directly fixes the Noisy-OR inflation problem. Verified on production data:
the same 22.5%-scored->0.95 stratum drops to **2.3%** under this formula, and the
379 rows previously pinned at 1.0 score **< 0.5**.

*How:* `override_mask = both_layers_present & weight_available`; overrides
`core_score` for the two-layer stratum only (1-layer fallback cells are untouched,
since there's nothing to weight). This is the value written to
`outputs/core_score.parquet` and consumed by every stage downstream: driver-gated
routing, confidence tiers, held-out evaluation, and the explain/rank/CLI layer.

### 2.5 Alternative normalisation/combination variants (validation harness only)

`scoring_variants.py` implements a 4×3 grid (4 normalisations × 3 combinations = 12
variants) used *only* to empirically compare against GDSC ground truth — none of
these except 2.1+2.4 are in the production path.

| Normalisation | Formula |
|---|---|
| `minmax` | `(x - min) / (max - min)`, span 0 → NaN. Flagged as the weakest by audit. |
| `percentile` | `rank(pct=True)` — production choice, see 2.1 |
| `zscore` | `z = (x - μ) / σ`, then `sigmoid(z) = 1 / (1 + e^-z)` to map into [0,1] |
| `hill` | percentile rank, then `pct^k / (pct^k + p0^k)` with placeholder `k=2, p0=0.5` — explicitly marked "uncalibrated" |

| Combination | Formula |
|---|---|
| `noisy_or` | `1 - (1-E)(1-P)` — see 2.3 |
| `weighted_sum` | `0.5·E + 0.5·P` (equal weight, no correlation penalty) |
| `product_floor` | `max(E·P, 0.5·min(E,P))` — "prevents one weak layer from annihilating a strong one" |

### 2.6 Stratum-aware rank

```
stratum_rank = core_score.groupby(["ensg_id", "n_layers"]).rank(pct=True, method="average")
```
*Where:* `scoring_variants.score()` ([scoring_variants.py:172-175](src/pipeline/scoring_variants.py)); `02_core_score.ipynb`.

*Why:* a cell line with `n_layers=2, core_score=0.7` isn't comparable to one with
`n_layers=1, core_score=0.7` — they rest on different amounts of evidence.
Cross-stratum calibration is explicitly deferred pending a held-out regression
against CRISPR/GDSC anchors.

*How:* percentile rank of `core_score` computed **within** the (gene, n_layers)
group. This is the metric surfaced to end users as `stratum_rank_percentile` in
`explain_pair.py` / `rank_cell_lines.py` — it's the safe way to compare cell lines
that may have different evidence coverage.

### 2.7 Chronos three-layer blend (tested, rejected)

```
chronos_pct = 1 - rank(chronos_raw, pct=True)   # invert: more negative raw = more essential
core_score_3 = w_EP · core_score_2layer + w_C · chronos_pct
```
built only on a 50/50 train split of cell lines (to avoid circularity), and gated
per gene by `ρ(chronos_pct, ground_truth) > 0.05` (Stage "Confound tests" below).

*Where:* `src/pipeline/build_chronos_layer.py:113-137`; re-derived in `validate_chronos.py:54-77`.

*Why:* tests whether genome-wide CRISPR essentiality (Chronos) improves ranking as
a third evidence layer.

*Outcome:* **rejected.** "BCL2's hit@20 collapsed to 0 when Chronos was blended
into core_score as a 3rd layer" — not used in production.

---

## Stage 3 — Cell-line Name Resolution Audit

`02_core_score.ipynb` (audit cell)

Compares the deterministic RRID→CVCL→ACH resolver (`model_id_for()`) against naive
substring matching (`fuzzy_model_id_for()`) across 30 test names, classifying each
into `OK / RRID_RESCUED / WRONG / FUZZY_ONLY / AMBIGUOUS_BUT_AGREE / BOTH_MISSING`.
Result: 18 rescued, 9 OK, 1 wrong, 1 both-missing, 1 fuzzy-only. Diagnostic only —
not consumed by scoring, no numeric similarity metric involved.

---

## Stage 4 — Driver-Gated Routing

`src/pipeline/04_driver_gated_routing.ipynb`, `src/pipeline/build_cna_layer.py`

### 4.1 Driver flag (Boolean logic, no numeric formula)
```
mut_driver = any_driver | oncogene_hit | tsg_hit
fusion_driver = (max_confidence == "high")
has_driver_alteration = mut_driver | fusion_driver
```
*Why:* replaces the earlier "any alteration present" gate (`p_mutation > 0.5`),
which was found to add variance without improving hit rate. Restricting the routing
tie-break to *driver* events only reduces noise (quantified in 4.3/6.4).

*How:* used as the primary sort key ahead of `core_score` for `activation_driven`
and `unknown` gene classes:
```
sort_values(["has_driver_alteration", "core_score"], ascending=[False, False])
```

### 4.2 Copy-number amplification/deletion thresholds
```
AMP_THRESHOLD = 2.5    # TOTAL_CN > 2.5 -> amplification
DEL_THRESHOLD = 1.5    # TOTAL_CN < 1.5 -> deletion
```
*Where:* `build_cna_layer.py:30-31, 88-89`.

*Why:* standard copy-number cutoffs around diploid (CN=2), gated by gene role —
deletion is only meaningful evidence for a tumour-suppressor gene (TSG), and
amplification only for an oncogene:
```
has_cna_alteration = (
    (is_amplification & role in {oncogene, both}) |
    (is_deletion       & role in {tsg, both}) |
    (is_amplification & regime_class in {activation_driven, abundance_tracking}) |
    (is_deletion       & regime_class == loss_of_function)
)
```
*How:* collapsed to one row per (model_id, ensg_id); feeds `flags_with_driver.parquet`
and confidence-tier upgrades in Stage 5.

### 4.3 / 6.4 Variance-reduction metric

```
reduction_pct = (1 - std(Δ_driver) / std(Δ_any)) * 100
```
where `Δ_driver = hit_rate_driver_gated - hit_rate_flat` and `Δ_any = hit_rate_any_alteration - hit_rate_flat`.

*Where:* `04_driver_gated_routing.ipynb` / `06_held_out_eval.ipynb`; `stage4_eval_save.py:76-79`.

*Why/How:* quantifies how much switching from "any alteration" gating to
"driver-only" gating reduces the *spread* (noise) of hit-rate lift across genes in
the `activation_driven` class. Reported as a headline number, not thresholded into
an automatic accept/reject rule — the routing change was accepted based on this and
the held-out hit@20 results together (Stage 6).

---

## Stage 5 — Confidence Tiers

`src/pipeline/05_confidence_tiers.ipynb`, `src/pipeline/build_full_predictions.py`

Deliberately **not** a numeric 0–1 confidence score — categorical, deterministic
rules with an explicit precedence order:
```
tier = "high"     if regime_source == "measured" & class == "abundance_tracking" & n_layers == 2
tier = "moderate" (default fallback)
tier = "low"      if class == "activation_driven" & rank_basis == "score_only"
tier = "unknown"  if regime_source != "measured"
# then, applied last (highest precedence):
tier = CNA-based upgrade                          if a CNA mechanistic hit is present
tier = "moderate" for unknown-class genes         if chronos_check == "validated" (Stage "Confound tests")
```
*Where:* `05_confidence_tiers.ipynb`; production version `build_full_predictions.py:91-134`.

*Why:* "Confidence is categorical and deterministic — not a fabricated 0–1 score."
Encodes an explicit priority order between four evidence types: empirical GDSC
validation, COSMIC gene-role priors, CNA mechanistic confirmation, and the
Chronos-essentiality cross-check for otherwise-unclassified genes.

*How:* surfaced as the user-facing `confidence` field (CLI, `explain_pair.py`,
`rank_cell_lines.py`) and feeds `signal_quality` classification
(`np.select` over confidence/regime_source/chronos_check, [build_full_predictions.py:189-194](src/pipeline/build_full_predictions.py)).

### 5.1 Tumour-suppressor score inversion
```
core_score_TSG    = 1.0 - core_score
stratum_rank_TSG  = 1.0 - stratum_rank
```
*Where:* `build_full_predictions.py:65-67`.

*Why:* for a tumour-suppressor gene (TSG), *low* abundance is the biologically
relevant signal (loss of the suppressor), the opposite direction from oncogenes or
abundance-tracking genes where high abundance is relevant. Inverting at build time
means every downstream sort can consistently sort descending.

---

## Stage 6 — Held-Out Evaluation

`src/pipeline/06_held_out_eval.ipynb`, `src/pipeline/stage4_eval_save.py`

### 6.1 Train/test gene split
```
test_genes = rng.choice(sorted(curated_genes), size=len(curated_genes)//5, replace=False, seed=42)
```
*Why:* reproducible 20% held-out gene set for evaluating ranking strategies without
leakage. Genes are sorted before sampling specifically "to ensure the same test
genes regardless of Python's set ordering" (set iteration order is not guaranteed
stable across runs/versions).

### 6.2 Vectorised routing sort key
```
sort_key = has_driver_flag.astype(float) * 1e6 + core_score
rank = sort_key.rank(ascending=False, method="first")
```
*Why:* the `1e6` multiplier guarantees the boolean flag always dominates
`core_score` in the ordering — this collapses a two-column sort
(`[flag, score]`) into one sortable numeric key, which vectorises cleanly across
all genes at once instead of looping.

### 6.3 Hit-rate@20 (hit@k)
```
hit_rate = |top_20_predicted ∩ GDSC_sensitive| / |GDSC_sensitive ∩ scored_in_core_score|
```
*Where:* `06_held_out_eval.ipynb`; `stage4_eval_save.py:54-64`; also `validate_chronos.py hit_at_k()` and `test_run_metabolomics_mirna.py hit_at_k()`.

*Why:* "the only honest denominator — we can only rank what we have a score for."
Standard top-k recall metric: of the cell lines known (from GDSC) to be sensitive to
a gene's inhibition, what fraction land in the model's top 20 predicted lines for
that gene?

*How:* computed per gene, averaged by class (`abundance_tracking` vs
`activation_driven`);
```
driver_lift_pct = (mean(hit_rate_driver) - mean(hit_rate_flat)) / mean(hit_rate_flat) * 100
```
used as the primary accept/reject evidence for the Stage 4 routing change.

---

## Stage 7 — Explain / Rank / CLI Query Layer

`src/pipeline/explain_pair.py`, `src/pipeline/rank_cell_lines.py`, `src/pipeline/cli.py`

No new math — pure retrieval and sorting over already-computed columns
(`core_score`, `n_layers`, `stratum_rank`, `confidence`).

- **Sort logic** (`sort_gene_rows()`): `abundance_tracking` genes sort by
  `core_score` descending; other classes sort by
  `[has_driver_alteration, core_score]` descending (TSG inversion from 5.1 means
  descending is always the correct direction, no per-class branching needed at
  query time).
- **Consistency self-check**: `rank_cell_lines.py` asserts the known BRAF/A375
  stratum-rank percentile improved after the 2.4 scoring-formula change
  (0.767 → 0.814) and stays within the top 5%:
  ```
  top5pct = max(1, round(total_n * 0.05))
  ```
  This is a regression-test threshold, not a scoring formula.

---

## Track C — Mutations & Fusions Scoring

`src/Track - C/03_mutations_scoring.ipynb`, `src/Track - C/04_fusions_scoring.ipynb`

Track C produces **two distinct kinds of output, with different fates**. This
distinction matters and was previously stated too loosely here:

| Track C output | Columns | Consumed by | Affects ranking? |
|---|---|---|---|
| **Cleaning** — `mutations_collapsed.parquet` | `any_driver`, `oncogene_hit`, `tsg_hit` | `mut_driver` → `has_driver_alteration` (Stage 4) | **Yes — this is the production gate** |
| **Cleaning** — `fusions_gene_level.parquet` | `max_confidence` | `fusion_driver` → `has_driver_alteration` | **Yes** |
| **Scoring** — `mutations_scores.parquet` | `p_mutation` | `explain_pair.py:137-153` (display); `has_alteration` (rejected arm) | **No** |
| **Scoring** — `fusions_scores.parquet` | `p_fusion` | `explain_pair.py:157-162` (display); `has_alteration` | **No** |

The Boolean driver columns from Track C's *cleaning* stage drive every
driver-gated ranking decision in the pipeline. The continuous Hill/Noisy-OR
*scores* (`p_mutation`, `p_fusion`) never enter the sort order: the production
gate is `has_driver_alteration = mut_driver | fusion_driver`, which is purely
Boolean and never reads them. Their only live consumer is `explain_pair.py`,
which surfaces them to the user tagged in-code as
`"CONFIDENCE MODIFIER — does not affect stratum_rank"`.

They also feed `has_alteration = (p_mutation > 0.5) | (p_fusion > 0.5)` — the
any-alteration gate that was **rejected** at Stage 4 (§4.1) and now survives only
as the comparison arm `hr_any_alt` in the held-out evaluation.

*Consequence for calibration:* the uncalibrated Hill parameters below matter less
than they appear to, because the scores they produce are display-only. Calibrating
them is only worthwhile if this layer is ever promoted into ranking. This is also
why the signature-discount test returned null (see Confound tests) — it tuned a
component whose output does not reach the metric being measured.

### Hill function (shared primitive)
```
Hill(x, p0, k) = x^k / (x^k + p0^k)
```
An S-shaped map from a raw or percentile feature onto (0,1), with midpoint at
`x = p0`. **Explicitly flagged as a placeholder**: "TEMP SCORING — p0 and k values
are placeholders pending calibration against CRISPR/GDSC anchors." Not yet fit to
any ground truth.

### Mutations scoring
```
p_vep    = Hill(max_vep_rank,      p0=3.0, k=2.0)   # 0-3 ordinal VEP impact
p_path   = Hill(max_pathogenicity, p0=0.5, k=2.0)   # max(AlphaMissense, REVEL)
p_burden = Hill(variant_burden,    p0=3.0, k=1.5)   # variant_burden = log1p(variant_count)

p_base   = 1 - (1-p_vep)(1-p_path)(1-p_burden)          # Noisy-OR
p_mutation = min(p_base + 0.15 * driver_flag, 1.0)       # DRIVER_BOOST=0.15 (locked)

p_mutation_quality = 1 - (1-p_vep)(1-p_path)             # burden excluded
p_mutation_burden  = p_burden                             # kept separate
```
*Why:* each raw feature is bounded into a (0,1) "evidence strength" via Hill, then
combined by Noisy-OR (same rationale as Stage 2.3 — treats VEP impact,
pathogenicity, and variant burden as independent evidence of a damaging mutation).
A flat +0.15 additive boost rewards a known oncogene/TSG driver hit specifically.
Quality (`p_vep`+`p_path`) is kept separate from burden because a downstream
signature-discount step scales burden but not quality
(`p_mutation_adj = NoisyOR(p_mutation_quality, m_mut · p_mutation_burden) + driver_boost`).

*How:* max-aggregated to one row per (ensg_id, model_id) → `mutations_scores.parquet`
→ consumed by Stage 4 (`p_mutation`, driver flags).

### Fusions scoring
```
conf_ord   = {low: 1, medium: 2, high: 3}[max_confidence]
ffpm_pctl  = rank(best_ffpm,     pct=True) within ensg_id
recur_pctl = rank(fusion_count,  pct=True) within ensg_id

p_conf  = Hill(conf_ord,   p0=2.0, k=2.0)    # midpoint at "medium"
p_ffpm  = Hill(ffpm_pctl,  p0=0.5, k=1.5)    # right-skewed distribution, gentle k
p_recur = Hill(recur_pctl, p0=0.5, k=1.5)    # 87% mass at percentile 1, gentle k

p_base   = 1 - (1-p_conf)(1-p_ffpm)(1-p_recur)
p_fusion = min(p_base + 0.15 * any_in_frame, 1.0)   # INFRAME_BOOST=0.15 (locked)
```
*Why:* same Hill→Noisy-OR→driver-boost pattern as mutations. In-frame fusions are
more likely to produce a functional (translated) chimeric protein, hence the boost.
Fusions has no quality/burden split — the whole score is scaled by the signature
discount (`s_fus' = m_fus · s_fus`).

*How:* max-aggregated to (ensg_id, model_id) → `fusions_scores.parquet` → consumed
by Stage 4.

---

## Confound / Diagnostic Tests

`test_run_metabolomics_lineage_confound.py`, `test_run_metabolomics_mirna.py`,
`test_run_mirna_lineage_confound.py`, `validate_chronos.py`, `build_chronos_validation.py`

These are validation experiments, not production code — but they determine which
candidate layers get accepted or rejected, and are the direct source of the
narrative in [DECISIONS.md](../DECISIONS.md).

### Spearman rank correlation (ρ) — the primary validation metric throughout

```
r1, r2 = rank(x), rank(y)                       # per gene, average-rank ties
ρ = Σ(r1-r̄1)(r2-r̄2) / sqrt(Σ(r1-r̄1)² · Σ(r2-r̄2)²)     # Pearson correlation on ranks
t  = ρ · sqrt((n-2) / (1-ρ²))
p  = 2 · (1 - T_cdf(|t|, df=n-2))                 # two-sided Student-t test
```
*Where:* vectorised manual implementation in `build_chronos_validation.py:65-87`
(avoids a per-gene Python loop across thousands of genes); `scipy.stats.spearmanr`
used directly everywhere else.

*Why:* measures whether `core_score` (or a candidate variant of it) agrees in
ranking with an independent ground truth — Chronos CRISPR essentiality, or GDSC
drug-sensitivity. This is the load-bearing metric for every "does adding layer X
help" decision in the project.

### Validated / Inverted / None classification
```
RHO_THRESHOLD = 0.1 ;  P_THRESHOLD = 0.05 ;  MIN_N = 30
"validated"  if ρ >  0.1 and p < 0.05
"inverted"   if ρ < -0.1 and p < 0.05
"none"       otherwise
```
*Where:* `build_chronos_validation.py:41-43, 89-92`; identical thresholds in
`test_run_metabolomics_mirna.py classify():43-45, 100-105`.

*Why:* per "unknown"-class gene (no curated regime, no COSMIC role), classifies
whether abundance-based `core_score` ranking tracks real functional dependency.
"Validated" genes get promoted from `confidence="low/unknown"` to `"moderate"`
(Stage 5); "inverted" genes are flagged as likely running backwards.

*How:* written to `chronos_validation.parquet` (`ensg_id, n, rho, pval, chronos_check`),
consumed by `build_full_predictions.py` (lines 82-89, 127-133, 158-168).

### Chronos correlation gate
```
use_chronos = ( ρ(chronos_pct, ground_truth) > 0.05 )     # per gene
```
*Where:* `validate_chronos.py chronos_correlation_gate():121-128`.

*Why:* "only add Chronos if it aligns with drug sensitivity" — a per-gene safety
gate so a disagreeing Chronos signal isn't blended in even experimentally.

### Bootstrap confidence interval on ρ
```
for b in 1..1000:
    sample = resample_with_replacement(cell_lines)
    ρ_b = spearman(x[sample], y[sample])
CI_95 = [percentile(ρ_b, 2.5), percentile(ρ_b, 97.5)]
```
*Where:* `validate_chronos.py bootstrap_rho():101-110`.

*Why:* quantifies uncertainty around a point-estimate ρ (e.g. 3-layer vs 2-layer
score against GDSC ground truth for BRAF/BCL2/PTEN) — a single ρ number without a
CI can't distinguish a real effect from sampling noise at this cell-line count.

### Paired bootstrap on Δρ (confound tests)
```
for b in 1..2000:
    sample = resample_with_replacement(rows)     # same rows for both scores
    ρ_baseline = spearman(baseline[sample],  ground_truth[sample])
    ρ_candidate = spearman(candidate[sample], ground_truth[sample])
    Δ_b = ρ_candidate - ρ_baseline
CI_95 = [percentile(Δ, 2.5), percentile(Δ, 97.5)]
significant  if 0 ∉ CI_95
```
*Where:* `test_run_metabolomics_lineage_confound.py:100-114`;
`test_run_mirna_lineage_confound.py:95-113`.

*Why:* pairing the resample (same rows feed both the baseline and candidate score)
removes shared sampling noise, giving a more powerful test than two independent
bootstraps for "did adding this candidate layer change correlation with ground
truth?"

*How:* `candidate = 0.5·baseline_core_score + 0.5·global_pct` (see below).
**Result in both tests (metabolomics, miRNA): CI included 0 — not significant, and
the point estimate moved in the wrong direction** (lower correlation). Both were
rejected from production — see [DECISIONS.md](../DECISIONS.md).

### Constant-term global-scalar blend (rejected candidate layers)
```
global_pct = rank(mean_z_across_all_metabolites_or_mirnas, pct=True)   # per cell line
candidate  = 0.5 · baseline_core_score + 0.5 · global_pct
```
*Why:* tests whether a single gene-agnostic "global metabolic activity" or "global
miRNA burden" scalar, added as a constant term, changes correlation with GDSC
ground truth for BCL2/MET — specifically probing for a **lineage-confound trap**:
a scalar that looks predictive only because it's secretly proxying tissue-of-origin,
which is itself strongly linked to certain drug responses (e.g. BCL2i in
blood/lymphoid lineages).

*Outcome:* rejected for both modalities — see z-score/ANOVA below for why.

### Z-score standardisation (per feature, before averaging into a global scalar)
```
z = (x - mean(x)) / std(x)             # per metabolite or per miRNA
global_z = mean(z across all features, axis=1)   # per cell line
```
*Where:* `test_run_metabolomics_lineage_confound.py:39-40`;
`test_run_mirna_lineage_confound.py:35-38`.

*Why:* standardises each metabolite/miRNA to comparable scale before collapsing
hundreds of them into one average — otherwise features with larger raw variance
would dominate the mean.

### One-way ANOVA + eta-squared (lineage effect size)
```
F, p = f_oneway(group_1, group_2, ..., group_k)     # groups = per-lineage arrays of global_z
R² (eta-squared) = SS_between / SS_total
SS_between = Σ_lineage n_lineage · (mean_lineage - grand_mean)²
SS_total   = Σ_i (x_i - grand_mean)²
```
*Where:* `test_run_metabolomics_lineage_confound.py:55-63`;
`test_run_mirna_lineage_confound.py:52-59`.

*Why:* directly tests whether a cell line's lineage (tissue of origin) explains
variance in the global scalar. A large R² here would mean the scalar is mostly a
lineage proxy, which would make the constant-term blend above dangerous to adopt
(inflating a gene's apparent correlation for the wrong reason — via lineage, not
biology).

*Results (recorded in [DECISIONS.md](../DECISIONS.md)):*
- Metabolomics: R² = 0.042 (4.2%), F=1.713, p=0.0198, n=928, 24 lineage groups
- miRNA: R² = 0.072 (7.2%), F=2.975, p=2.5e-06, n=948, 25 lineage groups

Both statistically significant (large n) but small in effect size — lineage does
*not* dominate either global scalar. Combined with the null/negative Δρ result
above, the conclusion in both cases was that the global scalar is closer to noise
than to a usable signal of any kind (lineage or gene), most likely from averaging
across features spanning unrelated biological pathways — **not** folded into
`core_score` either way.

### Metabolomics enzyme-flux proxy (tested, null result)
```
proxy = product_metabolite_percentile - substrate_metabolite_percentile
```
*Where:* `test_run_metabolomics_mirna.py:244-250` (`METABOLITE_MAP`).

*Why:* "higher proxy = more flux through that enzyme" — literature-curated
substrate/product metabolite pairs per gene (e.g. NAMPT: niacinamide → NAD+), tested
against Chronos essentiality via the same ρ + classify() rule as above.

*Outcome:* null — "no gene cleared ρ>0.1, p<0.05." Not used in production.

### miRNA-as-inhibitor dampening (tested, mixed/weak result)
```
burden      = mean(expression of curated inhibitory miRNAs)     # per cell line
burden_pct  = rank(burden, pct=True)
dampener    = 1.0 - ALPHA · fillna(burden_pct, 0.0)              # ALPHA = 1.0
E_effective = E_raw · dampener
```
*Where:* `test_run_metabolomics_mirna.py:174-187`.

*Why:* models curated, experimentally-known miRNA regulators (e.g. BCL2 ←
miR-34a/15a/16) as multiplicative inhibitors of expression, for the two genes
(BCL2, MET) that have both a curated miRNA map and GDSC ground truth.

*Outcome:* mixed/weak — "MET improved slightly, BCL2's already-odd baseline
direction shrank toward zero." Not adopted into production.

---

## `new-arch` — Continuous vs Discrete Diagnostic

`src/new-arch/01_long_form_gene_axis_diagnostic.ipynb`

```
score_continuous = (expression - min) / (max - min)     # per gene, min-max
score_discrete   = (expression > 1.0).astype(float)      # HPA_THRESHOLD = 1.0 log2-TPM
disagreement     = |score_continuous - score_discrete|
```
*Why:* compares the project's own continuous percentile-style scoring paradigm
against the Human Protein Atlas's binary present/absent convention (`>1.0 log2-TPM`
is the standard HPA cutoff), to see where the two paradigms diverge most (e.g. top
disagreement rows for BRAF). Explicitly diagnostic — "no fix is proposed here, that
is a separate, later decision." Not part of the production formula.

---

## Cross-Pipeline Formula Reuse

| Formula | Canonical source | Reused in |
|---|---|---|
| Median | **0.0426** | 0.0413 |
| Mean | 0.0509 | 0.0490 |
| SD | 0.0867 | 0.0850 |
| IQR | [−0.001, +0.099] | [−0.002, +0.097] |
| Fraction > 0.2 | **4.9%** | 4.4% |
| Fraction > 0.3 | 0.7% | 0.5% |

Mean shared lines per pair: 1,042 (median 978, range 678–1,493).

L2-only lines (two layers, tighter evidence): median r = 0.009 — more independent
than L1-only (median r = 0.050). Two-layer genes show less cross-gene correlation,
the opposite of what correlation-driven inflation would predict.

Rank displacement (top-20, selectivity formula vs correlation-adjusted): mean 3.51,
median 3.0 positions out of 20. Driven by the ~5% of pairs with r > 0.2.

Pre-declared: median corr > 0.2 → **WRONG**. Measured median is 0.043.

**Conclusion.** Core scores of randomly selected gene pairs are near-independent
across cell lines. The selectivity formula `score_A × (1 − score_B)` is
empirically justified — the independence assumption holds on average. Tail pairs
(r > 0.2, ~5%) contribute modest rank displacement (3–4 positions in top-20);
no systematic correction is warranted.

*Source: `Ranking/cli.py`, selectivity computation; `diagnostics/X3_selectivity_corr_fast.py` (2026-08-14).*

---

## 11.10 Stage 5 co-selection — tie inflation from min() (X4)

The Stage 5 joint co-selection score for a gene pair is:

```
joint_score = min(score_A, score_B, …)
```

with floor `JOINT_FLOOR = 0.50` (`cli.py:172`). The concern is that `min()` is
idempotent only when both scores are equal — when they differ it projects to the
lower value, potentially tying lines that have different profiles.

**X4 measured (100 pairs, seed=99, floor=0.50):**

| Quantity | Measured |
|---|---|
| Joint min() tie rate | median **17.7%** (IQR 14.2–21.9%) |
| Single-gene tie rate | median **16.3%** |
| Excess (joint − single) | median **+0.7 pp** |
| Top-20: lines at min value | mean=1.0, median=1.0 |
| Non-min score spread | median=0.065 |

Pre-declared: excess ≥ 10 pp → **WRONG**. The actual excess is 0.7 pp — the min()
operator creates negligible additional tie inflation relative to single-gene
scoring at this floor. Lines above the floor have broadly similar score values,
making the tie rates comparable.

**Conclusion.** No tie-inflation correction is warranted.

---

## 11.11 Fusion channel dominance (X5)

The Noisy-OR fusion score combines three channels:

```
p_base   = 1 − (1−p_conf)(1−p_ffpm)(1−p_recur)
p_fusion = min(p_base + 0.15·any_in_frame, 1.0)
```

**X5 measured (144,221 fusion records, 64.8% pass gate p_fusion ≥ 0.5):**

| Category | Count | Share of passing |
|---|---|---|
| 0 channels ≥ 0.4 (inframe-boost only) | 9,296 | 10.0% |
| Exactly 1 channel ≥ 0.4 | 51,825 | **55.5%** |
| 2 channels ≥ 0.4 | 25,073 | 26.8% |
| 3 channels ≥ 0.4 | 7,223 | 7.7% |

Among single-channel passes: p_conf dominates — **93.0%** of single-channel passes
occur via p_conf alone. The p_conf channel reflects `max_confidence ∈ {low,
medium, high}` mapped through Hill(conf_ord; p₀=2, k=1.5), with median
p_conf = 0.500 across all fusions (55.8% have p_conf ≥ 0.5).

**Genuine Noisy-OR** (no single channel ≥ 0.5): 10.3% of passes → CONFIRMED < 20%.
**Pre-declared** (>60% single-channel) → **WRONG** at 55.5%.

**Implication.** The fusion channel is effectively a confidence-tier classifier
with p_ffpm and p_recur as secondary modifiers. The ordinal-to-continuous mapping
via Hill(conf_ord; ·) is the structural issue noted in §C.3 — the choice of
encoding {low=1, medium=2, high=3} determines the score for the majority of
passing fusions.

---

## 11.12 CNA attrition waterfall (X6)

Starting from all rows in the `main.cosmic_cna` DuckDB table:

| Step | Rows | Lost |
|---|---|---|
| All rows (Step 1) | 143,052 | — |
| After sample join to pipeline lines (Step 2) | 142,572 | 480 |
| After is_ambiguous = FALSE (Step 3) | 142,412 | 160 |
| After gene_role join (Step 5) | 142,412 | 0 |
| After direction gate (Step 7) | 3,070 | 139,342 |
| Final cna_flags has_cna_alt | 2,925 | 145 |

**Dominant loss: unknown gene role.** Of 142,412 post-ambiguity rows, 138,185
(97%) carry `gene_role = unknown` — they belong to genes not in the curated
580-gene set (`oncogene=256, tsg=254, both=70`). These are structurally
unevaluable; the direction gate applies only to role-annotated genes.

Among role-annotated rows (4,227): 1,157 (27.4%) fail the direction gate
(amplification for a TSG, deletion for an oncogene, or neutral).

**Pipeline coverage.** 992 of 1,840 pipeline lines appear in COSMIC CNA (54%).
The remaining 46% of lines have no CNA data and receive `has_cna_alteration = FALSE`.

**Ploidy artefact (X6, ploidy correlation):**

```
Spearman r(per_line_mean_CN, per_line_alt_rate) = 0.2488   p = 1.97 × 10⁻¹⁵
```

Lines with higher mean copy number have higher CNA alteration rates. Pre-declared
|corr| < 0.15 → **WRONG**.

**Stratified analysis (diagnostics/D1_ploidy_stratified.py, 2026-08-15):**

The pipeline uses COSMIC's pre-computed `cna_call` (ploidy-aware), not raw absolute
CN thresholds. The stratification reveals whether the residual correlation is artefact
or biology:

| cna_call | n lines | r(mean_CN, correct-dir rate) | r(mean_CN, correct-dir n) | r(mean_CN, total calls) |
|---|---|---|---|---|
| amplification | 795 | +0.154 (p=1.2×10⁻⁵) | +0.396 (p=2.7×10⁻³¹) | +0.721 (p=1.1×10⁻¹²⁸) |
| deletion | 993/997¹ | +0.090 (p=4.6×10⁻³) | −0.063 (95% CI [−0.124, −0.001], p=0.048) | −0.452 (p=4.6×10⁻⁵¹) |

¹ D1 n=993 lines; Y9 re-measurement n=997 (slightly different filtering; values consistent).

`correct-dir rate` = fraction of amplifications that are oncogene-direction (or
deletions that are TSG-direction) — **per-gene** fraction. `correct-dir n` = absolute
count of direction-gated events **per line**. `total calls` = raw amp or del call count.

**Interpretation:** The r=+0.72 for total amplification count is driven by
chromosomal instability — highly aneuploid lines have more focal amp events genome-wide,
including non-oncogene loci. For ONCOGENE amplifications specifically (correct-direction
rate), r=+0.154. For TSG deletions: the correct-direction rate (per gene) is r=+0.090,
while the correct-direction count (per line) is r=−0.063 (95% CI [−0.124, −0.001],
p=0.048 — barely excludes zero; marked **UNEXPLAINED**: the CI sign is negative but
the mechanism producing TSG-direction deletion depletion in high-CN lines, independently
of the total deletion depletion at r=−0.45, has not been demonstrated). The D1 value
of −0.056 (p=0.08) was within sampling uncertainty of the Y9 value and is now
superseded.

**Conclusion:** The ploidy correlation in the pipeline is primarily chromosomal-instability
biology, not a threshold artefact. COSMIC's ploidy-aware calls already prevent the naive
absolute-threshold problem (CN>2.5/CN<1.5). The residual r≈0.15 is biologically expected.
No code change required. §9 gap 21 is **resolved by the existing implementation**.

---

## 11.13 p_burden confound (X7)

The mutation score includes `p_burden = Hill(variant_burden; p₀=3.0, k=1.5)` where
`variant_burden = variant_count.clip(upper=5)`. The concern is that
hypermutated lines receive elevated p_mutation regardless of biology.

**X7 measured (632,919 gene×line records, 1,744 lines):**

| Quantity | Measured | Pre-declared | Verdict |
|---|---|---|---|
| Spearman r(per_line_TMB, p_mutation) | **0.159** | > 0.20 | WRONG |
| Hypermutator rate / rest rate | **1.06×** | ≥ 1.50× | WRONG |

TMB distribution: mean=404, median=205, p95=1,392 variants.
Hypermutators (top 5%, TMB ≥ 1,392): 88 of 1,744 lines.
Lineage: colorectal (25), uterus (19), blood (11), ovary (8), lung (6).

**Why pre-declarations are wrong.** `mut_driver` rate is 88.9% in non-hypermutators
and 94.1% in hypermutators — compressed near the ceiling. With 89% baseline nearly
everything is already classified as a driver, leaving little room for hypermutators
to differ (ratio 1.06× vs pre-declared 1.5×). The ceiling is a consequence of the
`any_driver ∨ oncogene_hit ∨ tsg_hit` disjunction driving most records to true.

**Ablation (remove p_burden entirely):** 170,334 / 632,919 records (26.91%) change
from `mut_driver = TRUE` to `FALSE`. All are losses (none gained). This shows
p_burden is a material contributor to driver calls, despite the low TMB correlation.
The mechanism: even with modest TMB (median=205 variants), some lines have
`variant_burden = 5` (clipped ceiling) for particular (gene, line) pairs, which
Hill-saturates p_burden. The saturation is not proportional to TMB in the way the
pre-declaration assumed.

**Conclusion.** p_burden does not produce a genome-wide TMB confound at the
line level. Its effect is gene-specific and ablation-measurable but not
lineage-predictable from TMB alone.

---

## 11.14 Chronos re-run delta (X8)

The Chronos validation was re-run using real Chronos scores
(`chronos_validation.parquet`) in place of negated Project Score
(`chronos_validation_OLD.parquet` = Project Score × −1).

**Distribution comparison:**

| | OLD (Project Score negated) | NEW (real Chronos) |
|---|---|---|
| Genes | 16,866 | 17,226 |
| Mean ρ | −0.0052 | **+0.0103** |
| SD ρ | 0.0686 | 0.0725 |
| p5 | −0.118 | −0.105 |
| p95 | +0.096 | +0.122 |

The old data introduced a systematic negative bias (mean ρ = −0.005 vs +0.010
with real Chronos). This is consistent with using −1 × Project Score: any
positive dependency in Project Score became negative.

**Label changes (16,802 genes in both):**

| Tier change | Count |
|---|---|
| none → weak_positive | 944 |
| weak_negative → none | 633 |
| weak_positive → none | 221 |
| weak_negative → weak_negative (stay) | 547 |
| Total changed | **2,253 (13.4%)** |
| Direction flips (validated↔inverted) | **0** |

No gene moved between validated and inverted — the label changes are in the
weak/none region where the effect size (ρ ≈ 0.04–0.10) crosses the floor in one
direction. The validated gene count rose from 31 to 37 (+6), inverted from 6 to 8
(+2).

*Source: `diagnostics/T12_delta.py` against `src/pipeline/outputs/chronos_validation_OLD.parquet`
and `chronos_validation.parquet`.*

---

## §9 additions (from COWORK v4.0, X1–X8)

The following items extend the open-gaps list of §9:

**Gap 20. W_RNA is calibrated to ρ_implied = 0.25, but measured ρ_avg = 0.548 — RESOLVED 2026-08-15.**
`W_RNA = sqrt(1.431)` now applied in `core_score.py:39`. A/B test confirmed flat hit@20
improves (oncogene +10%, tsg +101%) with no loss in driver hit@20. See §11.3.

**Gap 21. Ploidy artefact in CNA direction gate — PARTIALLY RESOLVED.**
X6 measured r(per_line_mean_CN, per_line_alt_rate) = 0.248 (p = 2×10⁻¹⁵). Initial
interpretation: absolute threshold artefact. Stratified analysis (D1_ploidy_stratified.py,
2026-08-15) shows the pipeline already uses COSMIC's ploidy-aware `cna_call` — not raw
absolute thresholds. For oncogene-direction amplifications, correct-direction rate r=+0.154
per gene (modest). For TSG-direction deletions: correct-direction rate r=+0.090 per gene
(modest, biologically plausible); correct-direction count per line r=−0.063 (95% CI
[−0.124, −0.001], p=0.048, Y9 2026-08-15) — **UNEXPLAINED**: the CI barely excludes
zero on the negative side; the mechanism producing this residual effect after direction
gating has not been demonstrated. The total deletion count correlates at r=−0.45 (ploidy
suppression); the 6× reduction to r=−0.063 after direction gating is consistent with
direction-gating filtering to TSG-specific events. No code change required. See §11.12.

**Gap 22. MNAR shrink-direction defect (protein arm, 83.4% of pairs — updated).**
For low-abundance proteins absent from one source (MNAR), shrink moves z toward
zero when the true absent z would be strongly negative. The direction of shrink is
wrong for this specific missingness mechanism. Affected pairs: 83.4% (5,041,853 of
6,043,430 protein-arm pairs have n_sources=1; Y3 measurement 2026-08-15). Earlier
documentation (~23%) was wrong. Three resolution options escalated to Fiona (COWORK
v5.0 Y3): remove shrink, carry detected as feature, or document as limitation.
Architecture decision needed from Fiona (§11.5).

---

## 11.15 p_mutation distribution and ceiling (Decision 5)

`p_mutation` is computed per (gene, cell-line) pair in `mutations_scoring.py` via:
```
p_base   = 1 − (1 − p_vep)(1 − p_path)(1 − p_burden)     [Noisy-OR]
p_mutation = min(p_base + 0.15 × any_driver, 1.0)
mut_driver = (p_mutation ≥ 0.50)
```

**Measured on `mutations_scores.parquet` (632,919 pairs, 2026-08-15):**

| Statistic | Value |
|---|---|
| mean | 0.698 |
| std | 0.155 |
| min | 0.161 |
| 25th pctile | 0.531 |
| median | 0.710 |
| 75th pctile | 0.822 |
| At ceiling (=1.0) | 1.67% |
| ≥ 0.85 | 21.4% |
| mut_driver (≥ 0.50) | **90.6%** |

For role-annotated genes (oncogene/tsg/both): mut_driver rate = **92.3%**,
ceiling fraction 8.25%.

**Channel decomposition (among mut_driver=TRUE records):**

| Channel | Fraction single-channel |
|---|---|
| Only p_vep ≥ 0.50 | 1.4% |
| Only p_path ≥ 0.50 | 31.1% |
| Only p_burden ≥ 0.50 | 0.3% |
| ≥2 channels (multi-channel) | 2.2% |
| None ≥ 0.50 individually (Noisy-OR combination) | 65.0% |
| Boost-only (p_base < 0.50, passes only via DRIVER_BOOST) | **0.006%** (32 of 573,600) |

**Y2 finding — DRIVER_BOOST is effectively a dead code path.** Of 573,600 driver
calls, only 32 (0.006%) are boost-only passes that would not survive without the
`any_driver` annotation. Pre-declared >5%: WRONG by a factor of ~800. All 32 pairs
have p_base ∈ [0.35, 0.50) with any_driver=True; the near-gate region (p_base ∈
[0.35, 0.50)) contains 29,781 pairs of which only 32 are driver-annotated. The gate
passes on evidence, not annotation, for >99.99% of driver calls.

**Why 65% pass without any single channel ≥ 0.50:** Noisy-OR sums complementary
evidence. With p_vep=0.40, p_path=0.40, p_burden=0.40 (each below threshold),
p_base = 1 − 0.6³ = 0.784 >> 0.50. The gate is not binary-OR — it rewards
moderate multi-channel evidence. This is the intended architecture.

**Why 90.6% overall is high:** The input (`mutations_collapsed.parquet`) is
pre-filtered to non-synonymous, functionally annotated variants. The population
already excludes synonymous and intronic mutations. The 90.6% rate reflects that
most remaining mutations are pathogenic at some level in this filtered set, not
a broken gate.

**The 9.4% that fail the gate (p_mutation < 0.50)** have p_vep=1 (synonymous
VEP rank ~1–2), low CADD/REVEL score, and fewer than 3 variants — records that
passed the initial variant-type filter but have very weak functional evidence.

**Conclusion.** The gate is not architecturally broken. p_mutation operates as
a soft filter on pre-filtered mutations; the effective discrimination is on the
continuous score (0.16–1.0) used in ranking, not on the binary 90.6% pass rate.
No code change is applied for this decision (see §9 gap 23 for the open question
on whether to raise the gate to 0.65 for selectivity).

---

## §9 gap 23 (from D2 analysis, 2026-08-15)

**Gap 23. p_mutation gate selectivity: 90.6% pass rate on pre-filtered mutations.**
The binary `mut_driver` flag captures 90.6% of all (gene, line) pairs with any
observed non-synonymous mutation. The continuous p_mutation score (mean=0.70, std=0.16)
still varies and feeds the ranking. Open question: would raising the gate to 0.65
or using the continuous score directly as a ranking key (rather than binary promotion)
improve hit@20? Requires a Stage 6 A/B test against the current driver-promotion
architecture.

---

## 11.16 Code ownership boundary (Y1 audit, COWORK v5.0, 2026-08-15)

The production pipeline contains two categories of code:

**Fiona's independent implementations** — protein_scorer.py, platform_tier.py,
core_score.py, all of Stage 3 (altercations), Stage 4b (driver routing), Stage 4c
(confidence tiers), Stage 5 (ranking/CLI), Stage 6 (eval). These are not in
`02_transcriptonomics.ipynb` in any form.

**Adapted from the shared transcriptomics notebook** (`src/pipeline/02_transcriptonomics.ipynb`) —
`final_pipeline/utils/common.py::robust_z_matrix` and
`final_pipeline/01_Transcriptomics/rna_scorer.py`. Key simplification: the
notebook's `robust_z` returns `(z, reason)` with three guards (too-few-peers,
gene-silent-in-lineage, no-spread/MAD-floor); production `robust_z_matrix` retains
two guards (too-few-peers, MAD-floor) and drops the gene-silent-in-lineage guard
and the reason return entirely.

**Features absent from production:**

| Feature | Notebook function | Production status |
|---|---|---|
| `calibrate_constants()` | Per-gene MAD_FLOOR/MIN_PEERS derivation | ABSENT — frozen constants |
| Silent-lineage guard | `SILENT_FRAC=0.20` threshold | **RESTORED 2026-08-15** — added to `common.py:score_source_lineage()` after `robust_z_matrix`; sets Z column to NaN where < 20% of lineage lines express above 1.0 log2 TPM+1 |
| `z_reason` return | `(z, reason)` tuple | ABSENT — matrix only |
| Cochran's Q / I² | `heterogeneity()`, `I2_CUT=50` | ABSENT |
| FDR control | `false_discovery_control()` | ABSENT |
| `source_floor_check` | Per-run detection from data | ABSENT |
| `exclusion_gate` | 3-valued pass/fail/unknown | ABSENT |
| `shrink(z, k) = z×k/(k+1)` | EB shrink with posterior weight k | PROTEIN SHRINK DIFFERS: √(n_obs/N) |
| `collapse_isoforms` | Isoform deduplication | ABSENT |

**z_reason gap (T10/D26) — CLOSED 2026-08-15:** The silent-lineage guard fires
when fewer than 20% of a lineage expresses a gene above 1.0 log2 TPM+1 (EXPRESSED_MIN).
Previously omitted, this meant a line with near-zero expression received
`z = (value − median) / MAD_FLOOR` rather than NaN — an absence-as-signal error.
The guard is now present in `common.py:score_source_lineage()`, restoring the
notebook's third guard. Affected (gene, lineage) pairs produce NaN z-scores, which
propagate to NaN core_score; these orphan pairs are excluded from the hit@20 ranking
by the `core_score.notna()` filter in `eval.py`.

**Effect of guard restoration:** ~87K additional (gene, model) pairs moved from
"scored with uninformative z-score" to "unscored NaN" at the 2026-08-15 run.
Total unscored pairs at the final run: 364,470 (of which 87K are silence-guard-masked;
the remainder are alteration-only orphans with no expression measurement at all).

*Consequence for dissertation:* The adaptation decisions (which guards to retain)
are Fiona's; the notebook is the teammate's. No production file should be described
as containing the notebook's methodology without noting the simplification. The
silence guard restoration is documented in the §C.2 change and in ALTERATION_DEFENCE §13.
