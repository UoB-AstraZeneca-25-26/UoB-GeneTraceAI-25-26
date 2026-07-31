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
| Percentile rank `rank(pct=True)` | `scoring_variants.normalise()` | core_score (E, P), fusions (ffpm/recur), Chronos pct, metabolomics flux proxy, global scalars |
| Noisy-OR `1-(1-a)(1-b)...` | `scoring_variants.combine()` | core_score (superseded), mutations `p_mutation`, fusions `p_fusion` |
| Correlation-penalised weighted sum | `scoring_variants.combine_three()` | production core_score (2-layer), 3-layer Chronos blend (rejected) |
| Hill function `xᵏ/(xᵏ+p0ᵏ)` | Track C `hill()` | mutations (vep/path/burden), fusions (conf/ffpm/recur) — all placeholders pending calibration |
| Spearman ρ + p-value | `scipy.stats.spearmanr` / manual rank impl | Chronos validation, GDSC validation, all confound tests |
| hit@k | `validate_chronos.hit_at_k()` | Chronos validation, held-out eval (Stage 4/6), confound tests |
| Bootstrap CI (percentile method) | `validate_chronos.bootstrap_rho()` | 3-layer ρ CI, paired Δρ CI in confound tests |
| Stratum-aware percentile rank | `scoring_variants.score()` | production ranking metric, `explain_pair.py`/`rank_cell_lines.py` |

## What's Calibrated vs. Placeholder

For dissertation purposes, worth flagging explicitly which formulas rest on fitted/
measured parameters vs. which are acknowledged placeholders:

- **Measured/fitted:** `ρ_EP=0.46` (literature + internal validation), `ρ_EC=ρ_PC=0.005`
  (measured directly), `AMP_THRESHOLD/DEL_THRESHOLD` (standard diploid-based cutoffs),
  `RHO_THRESHOLD=0.1`/`P_THRESHOLD=0.05` classification rule (chosen, not fit, but
  consistently applied throughout).
- **Explicit placeholders, not yet calibrated:** every Hill function `(p0, k)` pair
  in Track C mutations/fusions scoring, and the `hill` normalisation variant in
  `scoring_variants.py` (`k=2, p0=0.5`) — all flagged in-source as pending
  calibration against CRISPR/GDSC anchors.
- **Tested and explicitly rejected:** Chronos as a 3rd core_score layer (Stage 2.7),
  metabolomics global scalar, miRNA global scalar, metabolomics enzyme-flux proxy,
  miRNA-inhibitor dampening (all documented with reasoning in
  [DECISIONS.md](../DECISIONS.md)).
