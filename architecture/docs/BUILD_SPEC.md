# BUILD_SPEC — target architecture C1–C6

> ## ⚠ RETRACTION — the gate-region pAUC in this document is withdrawn
>
> **Issued 2026-08-10.** This is an EMPIRICAL retraction: the formula in
> §C.15 is mathematically correct and is still used in `gate_audit/`. What
> is withdrawn is the *claim* the formula was used to support. No arithmetic
> error is involved.
>
> This document cites **gate-region pAUC 0.596–0.608**, p = 2.8e-76, as verified
> and replicated. **That result does not stand.** It was measured before
> `gate_audit/` ran.
>
> | | |
> |---|---|
> | decile-1 ratio, as published | 0.8322 |
> | decile-1 ratio, screened lines only | **0.9875** |
> | placebo with dependency labels randomised | 0.8350 — **98% of the published depletion** |
>
> The effect was an artefact of counting never-screened cell lines as confirmed
> non-dependencies (`gate_audit/03_denominator.py`). It reproduces at DepMap 24Q4,
> so it was never a small-sample effect (`gate_audit/08`). What survives is
> CRISPR-specific: flat on GDSC2 drug response, and **opposite** on RNAi at
> matched prevalence (`gate_audit/09`).
>
> **Do not cite this number.** The claim it supports is withdrawn. Statements
> below about the gate being "the only place this architecture has replicated
> above-chance evidence" are superseded — the validated result in this project is
> the cell-line similarity graph (+0.186 [+0.149, +0.227] drug-response
> concordance over a same-tissue baseline, `cell_similarity/08`) and the Kleene
> filter (80.9% on 53 published cell-line facts,
> `docs/MULTIGENE_GROUND_TRUTH.md`).
>
> Retained rather than deleted so the record shows what was believed, when, and
> what changed it.

For each requirement: **what can be built from current assets**, **what is blocked and on what**, **what needs a decision**, and the **recommended design with its evidential basis**.

Evidence tags: `PROVEN` (published, cited) / `VERIFIED` (measured in `MEASUREMENTS.md`) / `CONTESTED` / `UNVALIDATED`. Designs with no published precedent are labelled **NOVEL FORMALISATION**.

Standing rule 2 applies throughout: **every free parameter below is specified with a sweep and a held-out selection procedure, never a picked value.**

---

## C1. All expression data used

### Buildable now

| Component | Status |
|---|---|
| DepMap layer (1,479 lines × 19,173 genes, `log2(TPM+1)`) | ready |
| HPA layer (1,104 × 19,209, `nTPM`, raw `TPM` also present) | ready |
| GEO layer (588 × 15,962) | ready |
| `model_id` mapping for all three | ready — `harmonised_enriched.parquet` |

### Can DepMap and HPA be placed on a common scale?

**Yes, on rank; no, on level.** VERIFIED.

- DepMap ships only `log2(TPM+1)`; the pre-log values are not on disk (recoverable as `2^x − 1`, but the original file is absent — a provenance gap, see C3).
- HPA ships `TPM`, `pTPM` and `nTPM`. `nTPM` is HPA's own cross-sample normalisation and is **not** the same transform as DepMap's.
- The two are therefore not on a common absolute scale and cannot be made so without re-deriving both from counts, which the repo cannot do.
- They **are** directly comparable after a within-gene rank or quantile transform, which is what the pipeline already does (`scoring_variants.normalise`, `percentile`).

**Recommendation:** keep the within-gene percentile as the common scale. Do not attempt a level-based harmonisation; there is no shared reference to anchor it to.

### Is an independence-assuming combination valid? **No.**

Measured (B3.1): median per-gene Spearman **ρ(DepMap, HPA) = 0.817**, IQR 0.700–0.903, 99.5% of genes positive. VERIFIED.

An uncorrected Stouffer combination of two sources at ρ = 0.817 inflates the combined *z* by a factor of

```
sqrt(2) / sqrt(2 + 2·0.817)  =  1.414 / 1.907  =  0.741
```

i.e. the naive combined z is **1/0.741 = 1.35× too large**, and the corresponding p-value is overstated accordingly. Noisy-OR fails worse, because it assumes the two sources can independently "fire".

**Required correction — Strube (1985) / Hartung (1999) generalised Stouffer** (`PROVEN`, published):

```
Z_combined = Σ wᵢ zᵢ  /  sqrt( Σᵢ Σⱼ wᵢ wⱼ Rᵢⱼ )
```

where `R` is the source-by-source correlation matrix. With equal weights and the measured R:

```
        DepMap   HPA     GEO
DepMap   1.000   0.817   0.467
HPA      0.817   1.000   0.523
GEO      0.467   0.523   1.000
```

**Estimate R per gene, not globally.** The IQR on ρ(DepMap,HPA) is 0.700–0.903 — a spread wide enough that a single global 0.817 will be wrong for the tails. Where a gene has too few shared lines for a stable per-gene estimate, shrink toward the global median (James–Stein / empirical-Bayes shrinkage).

**Free parameters and their sweep:**

| Parameter | Sweep | Selection |
|---|---|---|
| shrinkage weight λ on per-gene ρ toward global | λ ∈ {0, 0.1, 0.25, 0.5, 0.75, 1.0} | held-out 28 genes, gate-region pAUC, paired bootstrap over genes |
| minimum shared lines for a per-gene ρ | {20, 30, 50, 100} | same |

### Why GEO is held out

The brief cites τ = 0.666 (GEO) vs 0.748 (DepMap) / 0.741 (HPA) from the teammate's notebook. **That notebook is not in this repo** — `UNVERIFIED — could not locate`. The claim is not reproducible here.

What *is* verified here supports the same conclusion for partly different reasons:
- GEO's observed minimum value is **1.74**, not 0. It floors at background, consistent with RMA microarray. VERIFIED.
- GEO correlates with DepMap at ρ = 0.467 and with HPA at ρ = 0.523 — far below the 0.817 between the two RNA-seq sources. VERIFIED.
- GEO adds only **37 lines** (2.3% of the union) not already covered. VERIFIED.

So GEO is a genuinely different measurement, which is exactly what makes it useful as *corroboration* and dangerous as a *score component*: it is different partly because it is worse.

**Recommendation: DepMap + HPA in the score with the Strube/Hartung correction; GEO as categorical corroboration only.** Evidential basis: VERIFIED ρ matrix, VERIFIED coverage. Note the awkward fact that GEO agrees better with HPA (0.523) than with DepMap (0.467); this should be explained before GEO is trusted as an arbiter.

### GEO corroboration representation

A three-valued categorical, **never a score contribution**:

| Value | Condition |
|---|---|
| `corroborated` | GEO covers this (gene, line) **and** its within-gene percentile agrees with the DepMap+HPA percentile within a tolerance τ |
| `contradicted` | GEO covers it and disagrees by more than τ |
| `not_covered` | GEO has no measurement for this (gene, line) |

`not_covered` must be its own value and must never be rendered as "no corroboration found". Tolerance τ swept over {0.1, 0.2, 0.3} percentile units, selected on held-out genes; report the corroborated/contradicted/not-covered split at each τ so the choice is visible.

---

## C2. Final output = ranked cell lines per gene

### The conflict, restated with this audit's numbers

The brief cites ranker-region pAUC 0.500–0.518. Confirmed from `test_run_abundance_gate_chronos_results.json`: **ranker_region_Z = 0.5002, ranker_region_RNA = 0.5015**, against gate_region 0.5962 / 0.6076. VERIFIED.

This audit adds a second, independent confirmation that is stronger than the pAUC result because it is on held-out genes with an explicit chance baseline (B1.4):

> On the 28 held-out genes, **every** metric at **every** k has a 95% bootstrap CI on excess-over-chance that includes zero. The only above-chance cell in the whole table is `recall@20` on the 113 **tuning** genes.

A score-sorted list is not defensible. Confirmed, and now with an interval.

### Can each proposed level be computed from current data?

| Level | Sort key | Computable? | Evidential basis |
|---|---|---|---|
| 1 | Gate tier | **Yes** | gate pAUC 0.5962 (Z) / 0.6076 (RNA), p = 2.77e-76, `gate_frac_above_chance` = 0.722. Replicated on a second screen (`test_run_screen_replication`, `test_run_geo_gate_replication`). **supported**, VERIFIED |
| 2 | Count of independent confirmations | **Yes, but the count is not what it appears** | see below. **defensible with a correction** |
| 3 | Absolute protein copies/cell | **Partly** — see C6. Available for at most 952 of 1,746 scored lines, ploidy-correctable for 70% of the panel | **pending feasibility** |
| 4 | Alteration demotion | **Yes** — touches 1.6% of pairs (B3.6) | **pending decision**, see C4 |
| 5 | RNA percentile | **Yes** | ranker pAUC ≈ 0.500. **explicitly unvalidated** |

### Level 2 needs a correction before it is defensible

"Count of independent confirmations" presumes independence. Measured: DepMap and HPA agree at ρ = 0.817. **Counting DepMap and HPA as two confirmations double-counts one measurement.** VERIFIED.

**Recommendation.** Define the count over *sources that are actually independent*, not over source files:

| Confirmation | Counts as |
|---|---|
| RNA consensus (DepMap + HPA, combined under C1) | **1** |
| Protein confirmed (ProCan and/or TMT) | 1 |
| GEO corroborated | 1 |
| Essentiality confirmed | 1 |

Maximum 4, not 5. **NOVEL FORMALISATION** — I am aware of no published precedent for an evidence-count tie-break of this form; it is proposed here as a defensible construct, not as an established method.

### Output schema

```
gene_id, model_id,
-- level 1
gate_tier                  ENUM('survived','excluded_by_1','excluded_by_2plus')  NOT NULL
gate_exclusions            LIST<STRING>          -- which gates fired, empty if none
-- level 2
n_independent_confirmations  TINYINT NOT NULL     -- 0..4, per the table above
confirmations              STRUCT(rna, protein, geo, essentiality : ENUM('positive','negative','not_measured'))
-- level 3
protein_copies_per_cell    DOUBLE                -- NULL where not estimable
protein_copies_state       ENUM('estimated','not_measured','not_ploidy_correctable')  NOT NULL
-- level 4
alteration_state           ENUM('altered','not_altered','not_measured')  NOT NULL
alteration_detail          STRUCT(mutation, fusion, cna : ENUM('positive','negative','not_measured'))
-- level 5
rna_percentile             DOUBLE
rna_percentile_validated   BOOLEAN NOT NULL DEFAULT FALSE   -- always FALSE; forces the label
-- provenance (C3)
matched_via                STRING  NOT NULL
source_versions            MAP<STRING,STRING>  NOT NULL
n_used_in_denominator      MAP<STRING,INT>     NOT NULL
-- ordering
final_rank                 INT NOT NULL
ordering_validated_to_level TINYINT NOT NULL     -- 2 today; 3 if C6 lands; never 5
```

The renderer must refuse to display a ranking without `ordering_validated_to_level`, and must print a visible notice that ordering below that level is not validated.

### How this ranking must be evaluated

Non-negotiable, per standing rules 3 and 4:

1. **Paired bootstrap over genes**, ≥10,000 resamples, resampling *genes* not pairs, reporting the CI on `(observed − chance)`, not on `observed`.
2. **Hypergeometric chance baseline** alongside every @k figure: `E[recall@k] = k/N`; `P[any-hit@k] = 1 − C(N−K,k)/C(N,k)`. Without it, an @k number is uninterpretable — as B1.3 demonstrates, the headline 93.6% sits below its own 90.8% baseline.
3. **Lineage-stratified splits.** The current split (`06_held_out_eval.ipynb` cell 3) is an unstratified `rng.choice` over genes. With 382 haematopoietic lines against 1 adrenal, unstratified resampling understates variance (gap 7/17).
4. **Cluster bootstrap over lineage**, since cell lines within a lineage are not independent.
5. **Report the label-reliability ceiling** (`label_fpr = 0.00928`) beside every figure. A measured effect below the label noise floor is not an effect.
6. **28 genes is not enough for a point estimate.** Quote intervals only. If an interval is too wide to be useful, that is the finding.

**Recommended primary metric:** gate-tier **recall@k with its hypergeometric baseline**, because the gate is the only level with replicated above-chance evidence. Do not lead with a metric that measures the ordering.

---

## C3. Complete traceability

### Which tables currently have `matched_via`

| Has it | Lacks it |
|---|---|
| `cosmic_cna.parquet` | `core_score.parquet` |
| `gdsc_models.parquet` | `predictions_with_confidence.parquet` |
| | `evidence_ledger.parquet` |
| | `flags_with_driver.parquet` |
| | `cna_flags.parquet` |
| | `chronos_validation.parquet` |
| | `gene_dispersion.parquet` |
| | `gene_regime.parquet` |
| | `harmonised.parquet` / `harmonised_enriched.parquet` |
| | `gene.parquet` / `gene_enriched.parquet` |
| | `gene_ambiguity_flags.parquet` |
| | `stage4_eval_consistent_denom.parquet` |

**2 of 14 audited tables carry `matched_via`.** VERIFIED. Neither of the two is on the scoring path.

### Every place "not measured" is currently treated as a value — exhaustive

This is the most important list in this document. Ordered by how much damage each does.

| # | Site | What is collapsed | Scale of the damage |
|---|---|---|---|
| **1** | `build_full_predictions.py:122` — `c["has_cna_alteration"].fillna(False)` | *CNA not measured* → *no CNA* | **29,689,395 of 29,781,274 pairs (99.69%)**. CNA context exists for only 0.31% of pairs. VERIFIED |
| **2** | `build_evidence_ledger.py:68` — `c["has_cna_context"].fillna(False)` | same, in the ledger | same scale |
| **3** | `build_full_predictions.py:121` — `has_driver_alteration.fillna(False)` | *no mutation/fusion call* → *no driver* | pairs outside `flags_with_driver`'s 885,844 rows, i.e. ~28.9M of 29.78M |
| **4** | `build_evidence_ledger.py:66-67` — `mut_driver`/`fusion_driver` `.fillna(False)` | same | same |
| **5** | `stage4_eval_save.py:35-36`, `93` and `06_held_out_eval.ipynb` cell 4 — `has_driver_alteration`/`has_alteration` `.fillna(False)` | same, **inside the evaluation harness** | the held-out metrics in `stage4_eval_summary.json` are computed on this |
| **6** | `stage4_eval_save.py:41`, `06_held_out_eval.ipynb` cell 4 — `is_sensitive.fillna(False)` | *not tested in GDSC* → *not sensitive* | every (gene, line) not in GDSC2 is scored as a true negative. This is the denominator of every @k figure the project quotes |
| **7** | `build_gene_regime.py:150` — `is_sensitive…fillna(False)` | same, in regime assignment | decides `abundance_tracking` vs `activation_driven` for all 141 curated genes |
| **8** | `03_mutations_scoring.ipynb` / `test_run_signature_discount.py:85-87` — `hill(...).fillna(0.0)` on `p_vep`, `p_path`, `p_burden` | *VEP/pathogenicity/burden unavailable* → *probability 0* | feeds `p_mutation` for every variant with a missing annotation |
| **9** | `03_mutations_scoring.ipynb` — `oncogene_hit`/`tsg_hit` `.fillna(False)` (`test_run_signature_discount.py:90`) | *not in COSMIC CGC* → *not a driver* | COSMIC CGC has 768 genes; the other ~18,400 are silently "not drivers" |
| **10** | `04_fusions_scoring.ipynb` — `hill(df['conf_ord'].fillna(0), …)`, same for `ffpm_pctl`, `recur_pctl` | *fusion confidence unknown* → *lowest confidence* | all fusions lacking a confidence tier |
| **11** | `build_gene_regime.py:164-167` — `n_known`/`hits_at_k`/`n_scored` `.fillna(0)`, then `hit_at_20_train…fillna(0.0)` | *gene never evaluable* → *hit rate 0.0* | genes with no GDSC overlap get a real-looking 0.0 and are classified from it |
| **12** | `build_paralog_flags.py:66,79-80` — paralog counts `.fillna(0)` | *gene absent from the BioMart export* → *has no paralogs* | genes not in the 3.55M-row table |
| **13** | `test_run_tsg_two_hit.py:86` — `max_vaf.fillna(0) >= VAF_LOH_THRESHOLD` | *VAF not called* → *VAF 0* | biases the two-hit test toward "no second hit" |
| **14** | `test_run_metabolomics_mirna.py:184` — `1.0 - ALPHA * burden_pct.fillna(0.0)` | *burden unknown* → *no dampening* | |
| **15** | `test_run_stouffer_regime_split.py:356`, `test_run_stratum_centring.py:287` — `d.pE.fillna(0) + d.pP.fillna(0)` | *layer absent* → *contributes 0* | directly contradicts `scoring_variants.py`'s stated invariant that "missing layer contributes no term (not zero)" |
| **16** | `track_c_eda.py:188-189` — `mut_lines`/`fus_lines` `.reindex(...).fillna(0)` | *gene not in Track C* → *0 altered lines* | drives the 596 "genes with 0% of lines altered" count in B3.6 |
| **17** | `scoring_variants.py:169` — `out.dropna(subset=["core_score"])` | not-measured pairs are **dropped**, so downstream cannot distinguish "not scored" from "not in the panel" | 24.5M of 29.78M rows are one-layer; the absent layer leaves no trace |
| **18** | `harmonised_enriched.parquet` — `biotype` and `hgnc_status` typed `null`, 100% missing | a column that exists but carries nothing reads as "checked, no value" | 2,127 rows |

Items **1, 2, 6 and 15** are the load-bearing ones. Item 6 is the most consequential of all: **the evaluation harness treats every untested (gene, line) pair as a confirmed negative**, which is why an @k figure can look high while carrying no information — precisely the effect B1.3 measured.

Item 15 is a direct contradiction between a stated invariant and the code that implements it, in two separate scripts.

### Assertion set

Build-time barriers. Each fails the build loudly; none is a warning.

```python
# --- identity -------------------------------------------------------------
A1  every table has a non-null `matched_via` on every row
A2  no identifier is composite:
      no ';' or '/' in any gene symbol column          # 15 ProCan + 1 CCLE today
      no ';' in any model_id column                    # 400 dry-run rows today
A3  every ensg_id matches ^ENSG\d{11}$ (case-normalised)
A4  every model_id resolves to exactly one row of harmonised_enriched

# --- joins ----------------------------------------------------------------
A5  every merge asserts the post-merge row count equals the expected count;
    fan-out is declared explicitly or it is a failure
A6  no left join silently introduces NaN into a column declared NOT NULL

# --- three-state ----------------------------------------------------------
A7  no `.fillna(False)` or `.fillna(0)` on any column whose absence means
    "not measured"  -- enforced by a lint rule over src/, allowlisted per site
    with a written justification
A8  every measured column has a companion `<col>_state` in
    {measured_positive, measured_negative, not_measured}, and
    `<col>_state != 'not_measured'` wherever `<col>` is non-null
A9  aggregate reports state coverage: for every score component, the count of
    each of the three states, printed at build time

# --- provenance -----------------------------------------------------------
A10 every score component carries source_dataset, dataset_version,
    transform_applied, n_used_in_denominator
A11 dataset_version is non-null -- today DepMap's release is unrecorded, so
    this assertion fails on day one, by design
A12 no all-null column ships (catches harmonised_enriched.biotype)

# --- evaluation -----------------------------------------------------------
A13 no @k metric may be emitted without its chance baseline in the same record
A14 no comparison may be emitted without an interval
```

**A11 will fail immediately.** That is the intended behaviour: the DepMap release quarter is genuinely unknown (A1.1) and should block a traceability claim until it is recorded.

---

## C4. Alterations demote — **DECISION REQUIRED, NOT IMPLEMENTED**

### Every site where alteration status currently affects a score or ordering

| Site | Effect | Direction |
|---|---|---|
| `03_mutations_scoring.ipynb:62` | `p_mutation = (p_base + DRIVER_BOOST·driver_flag).clip(upper=1.0)`, `DRIVER_BOOST = 0.15` | **PROMOTE** |
| `04_fusions_scoring.ipynb:68` | `p_fusion = (p_base + INFRAME_BOOST·inframe_flag).clip(upper=1.0)`, `INFRAME_BOOST = 0.15` | **PROMOTE** |
| `stage4_eval_save.py:45-47` | `sort_driver = where(class=='abundance_tracking', core_score, has_driver_alteration·1e6 + core_score)` | **PROMOTE** (absolute precedence — 1e6 dominates any score) |
| `stage4_eval_save.py:51-52` | same with `has_alteration` | **PROMOTE** |
| `06_held_out_eval.ipynb` cell 5 | identical to the two above | **PROMOTE** |
| `stage4_eval_save.py:97` | `sort_values(["has_driver_alteration","core_score"], ascending=[False,False])` | **PROMOTE** |
| `explain_pair.py:109-112` | `driver_is_primary_sort_key(class)` → `sort_keys = ["has_driver_alteration","core_score"]` | **PROMOTE** (primary key) |
| `explain_pair.py:114-116` | else → `sort_keys = ["core_score","has_driver_alteration"]` | **PROMOTE** (tie-break only) |
| `build_full_predictions.py:140-145` | `rank_basis` ∈ {`core_score`, `lof_flag+inverted`, `lof_inverted`, `score+driver_tiebreak`, `driver_flag+score`, `score_only`} | selects which of the above applies |
| `build_full_predictions.py:155` | `activation_driven & rank_basis=='score_only'` → confidence `"low"` | **PROMOTE** (absence of alteration lowers confidence) |
| `build_cna_layer.py:30-31` | `AMP_THRESHOLD=2.5`, `DEL_THRESHOLD=1.5` → `has_cna_alteration` | role-gated, feeds the above |

**Every site promotes. There is no demotion anywhere.** VERIFIED. The requirement is a full reversal of a documented, load-bearing design decision that appears at 11 sites.

### Fraction demoted under a global rule

**1.597% of scored pairs** (475,695 of 29,781,274) carry any alteration. Driver alterations: 0.482%. Per-gene median 1.16% of lines; 87 genes exceed 10%; 596 genes have zero. VERIFIED (B3.6).

A global demotion is therefore a near-no-op for the median gene and a large re-ordering for ~87 genes.

### Is the existing role classification sufficient to condition on?

**Two independent classifications already exist, and they are not the same thing:**

| Classification | Values | Source | Coverage |
|---|---|---|---|
| `gene_role` | oncogene / TSG / both / neither | COSMIC CGC via `build_gene_roles.py` | 768 CGC genes; everything else is `neither` **by fillna, not by evidence** (C3 item 9) |
| `class` (regime) | `abundance_tracking` / `activation_driven` / `loss_of_function` / `unknown` | `build_gene_regime.py`, cut at `LIFT_DEFAULT = 2.0` | **141 curated genes only**; all others are `unknown` |

Neither is sufficient alone:
- `gene_role` covers only 768 genes with evidence, and its "neither" is indistinguishable from "not in COSMIC".
- `class` covers 141 genes, was fitted on the training half, and is decided by an unswept threshold.

**Answer: no. A use-case switch at query time is required.** The role classification can *inform* the default, but it cannot carry the decision for 19,177 genes when it has evidence for at most 768.

### The max-aggregation hazard

`03_mutations_scoring.ipynb` and `04_fusions_scoring.ipynb` both aggregate with `.groupby(["ensg_id","model_id"]).max()`. For m independent draws from a uniform, `E[max] = m/(m+1)`, so a gene with more variants in a line scores higher **mechanically**, independent of biology. §C.4 records 8.7% exposure. VERIFIED that the pattern is live.

**Any continuous alteration score inherits this, in either direction.** Under a demotion rule the bias flips sign but does not disappear: multi-variant pairs would be demoted *more*, again mechanically.

**Recommended form that does not inherit it — a boolean gate on a pre-specified variant class**, e.g.

```
altered := (∃ variant with vep_impact == 'HIGH')  OR  (∃ in-frame fusion with confidence ≥ tier)
```

A boolean is invariant to m. If a continuous score is required, use a count-calibrated statistic — the max's p-value under a null of m independent draws (`1 − F(x)^m`), or a fixed-rank order statistic — not the raw max.

### Recommendation: **a gate, not a continuous demotion**

Evidential basis (`PROVEN` by precedent within this project, VERIFIED here):
- CNA deletion is already a measured, independent exclusion signal: `test_run_cna_gate_results.json` gives deleted-rate 0.0397 vs base 0.1044, **ratio 0.380**, usable on the **69 genes with ≥20 deleted lines**. (The brief's OR 0.373 / ≥10 genes corresponds to an earlier run.)
- The gate region is the only place this architecture has replicated above-chance evidence (pAUC 0.596–0.608 vs ranker 0.500).
- 1.6% prevalence is gate-shaped, not tie-break-shaped.

**Free parameters and their sweep** (standing rule 2):

| Parameter | Sweep | Selection |
|---|---|---|
| variant class admitted to the gate | {HIGH only, HIGH+MODERATE, HIGH+hotspot} | held-out 28 genes, gate-region pAUC, paired bootstrap |
| fusion confidence floor | {high, medium, any} | same |
| minimum altered lines per gene for the gate to apply | {5, 10, 20, 30} | same; report coverage at each |
| role conditioning | {none, `gene_role` only, `class` only, both} | same |

**Do not implement until the user answers `DECISIONS_REQUIRED.md` Q1.** Both directions are defensible and the choice is not the auditor's.

---

## C5. Tissue lineage

### The decisive experiment is blocked

The brief asks for the teammate's guarded scorer to be run through the §C.15 harness. **`02_transcriptonomics.ipynb` is not in this repository.** A repo-wide search for the file, and for `silent_lineage`, `MIN_PEERS_FOR_LINEAGE`, `mad_floor` and `MAD_FLOOR`, returns nothing. `UNVERIFIED — could not locate`.

**Blocked on:** the notebook itself, or a specification of the three guards precise enough to reimplement (silent-lineage test, MAD floor value, minimum-peers rule).

The EGFR-in-blood observation the brief calls decisive (n=104, median 0.12, 38% exactly zero, MAD 0.178) is likewise not reproducible here. It is a **plausible and well-argued** claim — a MAD floor genuinely cannot catch a distribution that is degenerate but dispersed — but it remains `UNVALIDATED` in this repo.

### What the repo does contain

`test_run_improvement_scan.py` implements lineage-conditioned percentiles with `MIN_LINEAGE = 15`, and has been run. Result (`test_run_improvement_scan_chronos_results.json`, n = 1,792 genes):

| | pAUC |
|---|---|
| panel-wide percentile | 0.51756 |
| lineage-conditioned percentile | 0.51673 |
| **Δ** | **−0.0018**, p = 0.540, better in **47.8%** of genes |
| verdict | `lineage_conditioning_neutral` |

VERIFIED. This reproduces §C.15 almost exactly and confirms the brief's characterisation. **But it has no guards** — so it tests what §C.15 tested, not what the teammate proposes. The brief's central point stands: the guards are the novel component and remain untested.

**No paired bootstrap interval is available** for this comparison; the stored result reports a p-value only. Standing rule 4 means this Δ is not quotable as-is. Re-running with a bootstrap is cheap and should be done regardless of the guard question.

### Is the silent-lineage guard applicable to the pooled percentile as a standalone improvement?

**Yes, and this is the most valuable salvageable piece.** The guard asks: *is this gene's distribution within this stratum degenerate enough that a percentile is meaningless?* Nothing in that question requires the stratum to be a lineage. Applied to the pooled panel it becomes a per-gene validity check — and the repo already has the raw material in `gene_dispersion.parquet` (`dynamic_range`, `iqr_log2tpm`, `signal_spread`, `top_vs_median_fold`) and in `evidence_state.py`'s `FLAT_FOLD_CEILING = 3.0`.

**Recommendation:** implement the silent-lineage guard as a **pooled-percentile validity gate**, independent of any lineage decision. It converts a meaningless ordering into an honest `UNINFORMATIVE` abstention, which is a claim type this architecture already supports (`docs/QUERY_LAYER_STATES.md`). This is testable now, on held-out genes, without the missing notebook.

### Panel coverage at ≥15 lines

**1,587 of 1,746 scored models = 90.9%** sit in one of 21 lineages with ≥15 lines. On the full 2,145-model panel it is 1,725 = 80.4%, and **364 models (17.0%) carry no tissue label at all**. VERIFIED (B3.5).

### Recommendation

**Lineage as a stratification and annotation variable, not a scoring variable** — unless the guarded test comes back positive, which cannot currently be run.

Lineage is nonetheless **required** for:
- lineage-stratified CV splits (gap 17) — the current split is unstratified,
- cluster bootstrap over lineage (gap 7) — currently absent, and the reason held-out intervals are probably too narrow,
- output display,
- the 17% unlabelled models, which need an explicit `lineage_not_recorded` state rather than a default.

---

## C6. Absolute protein quantification — feasibility assessment

### Which dataset governs the answer

**Both are on disk, and they give opposite answers.** VERIFIED (A1.2).

| | CCLE/Gygi **TMT** | ProCan **DIA** |
|---|---|---|
| Values | log-ratios: median −0.042, **52.6% negative** | log2 intensities: median 3.44, **1.4% negative** |
| Verdict | **Absolute quantification impossible** — ratios to a bridge channel; no absolute scale to recover, and ratio compression distorts even relative values | **Possible in principle** |

The TMT verdict is exactly as the brief predicted. Everything below concerns ProCan.

### The three feasibility questions

**Q1 — Are histones present and quantified? YES.** VERIFIED.

17 histone proteins in the ProCan matrix, including several near-universally quantified across 952 lines:

| Accession | Symbol | Lines quantified | Median log2 |
|---|---|---|---|
| `P62805` | H4C16 | **952 / 952 (100.0%)** | 13.491 |
| `O75367` | MACROH2A1 | **952 / 952 (100.0%)** | 7.450 |
| `Q5TEC6` | H3-7 | 951 (99.9%) | 8.321 |
| `Q8IUE6` | H2AC21 | 947 (99.5%) | 8.182 |
| `Q92522` | H1-10 | 947 (99.5%) | 5.537 |
| `P16104` | H2AX | 907 (95.3%) | 4.737 |
| `P10412` | H1-4 | 894 (93.9%) | 2.433 |
| `P16401` | H1-5 | 878 (92.2%) | 10.925 |
| `P07305` | H1-0 | 847 (89.0%) | 7.707 |

Plus MACROH2A2 (87.6%), H1-1 (70.7%), H2BC26 (51.3%), H2BC1 (48.3%) and four sparse entries.

**Caveat:** `P62805` is an ambiguous accession mapping to 14 H4C genes (and `P68431` to 10 H3C genes in the TMT matrix). For the ruler this is *correct behaviour* — the ruler wants total histone mass, so collapsed paralogues are what you want. The same accessions must never be used for per-gene histone statements.

**Q2 — Has cross-sample normalisation destroyed the scale? PARTLY, AND THE REST IS UNDOCUMENTED.** VERIFIED + `UNVERIFIED`.

- Per-cell-line median log2 intensity: median 3.436, **sd 0.212**, range [2.667, 4.811], p1–p99 spread **0.98 log2 units ≈ 2× in linear intensity**.
- A hard median-scaled matrix would have sd ≈ 0 here. It does not. **Between-sample intensity differences survive.**
- Compare the TMT matrix, which *is* centred: per-line median sd 0.088.

**But:** the file is named `Protein_matrix_averaged_20250211.tsv`. An averaging step over replicates has been applied, and **what normalisation ProCan applied upstream is not recorded anywhere in this repository** — no README, no provenance JSON, no manifest. `UNVERIFIED — could not locate`. A 2× residual spread is consistent with *no* median scaling, but also with a partial or per-batch normalisation. **This must be resolved against the ProCan publication/portal before any absolute claim is made.** Less-processed intensities, if obtainable, would remove the ambiguity.

**Q3 — Is coverage deep enough? CANNOT BE ASSESSED ON THE PUBLISHED CRITERION.** VERIFIED.

- ProCan depth: median **5,218 proteins per line** (range 0–6,226; one all-missing row).
- The published guidance (Wiśniewski et al. 2014, *MCP* 13:3497) is that the histone-to-total ratio stabilises from **~12,000 peptides**. This matrix is **protein-level only** — peptide counts are not in the file. The criterion cannot be applied directly. `UNVERIFIED`.
- What *can* be measured: the histone/total linear-intensity ratio has median **0.0606**, IQR 0.0454–0.0816, **CV = 0.411**.
- Depth dependence is weak — Spearman(ratio, proteins quantified) = **−0.066, p = 0.042** — and the by-quintile profile is non-monotone (0.0687 / 0.0604 / 0.0561 / 0.0579 / 0.0625). So the variation is **not** primarily a depth artefact.

**A 41% CV on the ruler's internal standard is the binding constraint.** It propagates directly and multiplicatively into every copy-number estimate. An estimate carrying ±41% is still useful for an order-of-magnitude gate ("below ~1,000 copies") but is not usable as a tie-break, which is exactly how the brief proposes to use it.

### Ploidy correction — mandatory, and available for 70% of the panel

Median panel ploidy **2.653**; **55.8%** of lines exceed 2.5; only ~26% are near-diploid. VERIFIED (B3.4).

- Assuming diploidy would misestimate DNA per cell by a median factor of **1.33×**, up to 2.7× at the extreme, **and the error correlates with lineage** — exactly the failure mode the brief warns about.
- Ploidy is available for **1,502 of 2,145** models (70.0%). The remaining **30% cannot be ploidy-corrected** and must carry `protein_copies_state = 'not_ploidy_correctable'`, never a diploid default.
- **One correction closes gap 8 as well**: `AMP_THRESHOLD = 2.5` sits *below* the median line's average copy number of 2.65, so CNA calls should move to ploidy-relative thresholds in the same change.

### Overall verdict on C6

**Feasible in principle, blocked in practice on one item, and constrained by two.**

| | Status |
|---|---|
| Histones present | ✅ VERIFIED — not a blocker |
| Ploidy available | ✅ VERIFIED for 70%; explicit state for the other 30% |
| Scale survives normalisation | ⚠️ **BLOCKED** on ProCan's upstream normalisation being documented |
| Depth adequate | ⚠️ Cannot be assessed on the published peptide criterion; 41% CV on the standard |
| Coverage representative | ❌ No — Cliff's δ = +0.251 selection bias (B3.3) |

**Recommendation.** Pursue the ruler, but as a **second gate with an absolute threshold**, not as a tie-break — which is also what the brief argues, and the CV finding strengthens the case. Frame the output as an order-of-magnitude band with an explicit uncertainty, e.g. `<10³ copies (95% CI spans one order of magnitude)`, never a point estimate.

**Sequence the work:**
1. Resolve Q2 — document ProCan's upstream normalisation, or obtain less-processed intensities. Everything else is contingent on this.
2. Obtain peptide-level depth, or state explicitly that the published stabilisation criterion could not be checked.
3. Implement ploidy correction with the 30% not-estimable state.
4. Validate the gate on held-out genes with a paired bootstrap and the hypergeometric baseline.

**Fallback if the ruler fails.** `estimated_copies(gene, line) ≈ anchor_copies(gene) × relative_ratio(line)`, anchored to PaxDb (ppm, gene-level), Geiger et al. 2012 *MCP* 11:M111.014050 (11 lines), Bekker-Jensen et al. 2017 *Cell Systems* 4:587 (deep HeLa), Nagaraj et al. 2011 *MSB* 7:548 (HeLa) — all `PROVEN`. This inherits all anchor error and should be labelled as a rescaling of a relative profile, not a measurement.

### Two caveats that must travel with any C6 output

1. **Absolute quantification does not fix MNAR.** It covers at most 952 of 1,746 scored lines, and those lines are selected with Cliff's δ = +0.251 toward higher expression. Precision on a biased subset is not representativeness. VERIFIED.
2. **It reverses §2.1 deliberately.** §2.1 trades magnitude away for cross-assay comparability — that is what makes log₂-TPM and proteomic log-ratios combinable at all. Going absolute forfeits it, and **a new argument is required for how RNA and protein sit on comparable scales**. No such argument exists yet. This is a design debt, not a detail.
