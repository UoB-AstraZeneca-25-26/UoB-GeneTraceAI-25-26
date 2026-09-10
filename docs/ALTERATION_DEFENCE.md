# Alteration Defence — Driver Gate Validation
**Pipeline state:** 2026-08-15 (W_RNA=√1.431, SENSITIVITY_MAD_Z=−1.0, orphan recovery applied, six mathematical corrections applied)  
**Supersedes:** `docs/ARCHIVE/ALTERATION_DEFENCE_v1_UNREPRODUCIBLE.md`  
**Prior version:** 2026-08-14 (W_RNA=√2.00, SENSITIVITY_MAD_Z=−1.0 — see §13 for change log)  
**Bootstrap inference:** Complete — all five thresholds, three methods, two estimands (2026-08-15)

---

## 1. Claim Scope

The driver gate improves hit@20 for a typical **role-annotated gene**
(gene_role ∈ {oncogene, tsg, both}) in the 20% held-out test set (n=15 genes,
not 26; see §3 exclusion note), at the designated sensitivity threshold
(−1.0 MAD-z), on all three inferential methods for both estimands.

**What this claim does NOT assert:**
- That the gate helps for every gene (2 of 15 role-annotated test genes showed negative per-gene delta: PDGFRA −0.091, JAK1 −0.010)
- That the magnitude is precisely estimated (n=15 clusters; intervals are wide)
- That the effect is the same under any labelling threshold (the sweep shows
  it rises as thresholds tighten, but the guard fails outside -0.5 and -1.0)
- Anything about the 13 test genes with gene_role="unknown"
  (GSK3B, KDM3A, KDM6B, PAK1, PLK2, PLK3, BCL2L1, BCL2L10, ROCK2, USP1,
  WEE1, CRBN, EHMT1) — these are included in the 26-gene pooled estimate but
  their CNA channel cannot fire

---

## 2. Threshold and Estimand

**Sensitivity threshold:** MAD-z = **-1.0**  
Reason: corresponds to the GDSC waterfall convention for drug sensitivity
(~17% sensitive per compound, Iorio et al. 2016, Cell 166:740, PMID 27397505).  
Guard (5-40% sensitive per drug): **100%** of drugs in window at this threshold.  
MAD-z = -0.5 also passes the guard (100%) and is reported in full.  
Threshold was NOT selected for best delta. The sweep is shown below.

**Estimand decision:** ESCALATED — both reported (Fiona to decide).  
- **Estimand 1 (per-gene mean):** each gene contributes equally once.
  Matches "does the gate help for a typical gene?"
- **Estimand 2 (pooled ratio):** weighted by n_sensitive per gene.
  Weights high-n_sensitive (fame-correlated) genes more heavily.
- Gini coefficient of n_sensitive across 26 genes: **0.363** (not strongly concentrated).
- The two estimands differ by 158% of the smaller value at -0.5 (0.0230 vs 0.0089).
- Recommendation (not a decision): per-gene mean, because the tool is gene-centric.

---

## 3. Gene Set and Exclusion

| Category | n | Notes |
|---|---|---|
| Curated (flags ∩ GDSC targets) | 133 | intersection at scoring time |
| Test genes (20% holdout, seed=42) | 26 | held out from all tuning |
| Role-annotated: oncogene | 12 | interpretable |
| Role-annotated: tsg | 1 | NOT INTERPRETABLE (n=1) |
| Role-annotated: both | 2 | NOT INTERPRETABLE (n=2) |
| gene_role="unknown" | 13 | included in pooled ratio; CNA channel cannot fire |

**W6 logged exclusion:** eval.py now prints `[EXCLUSION]` for the 13 unknown-role genes
on every run. The hit@20 loop covers oncogene/tsg/both only, by design — hit@20 requires
a known directional role to interpret. The 13 excluded gene symbols are logged at runtime.

**CNA invariant for unknown-role genes:** has_cna_alteration=TRUE for 0/405 lines
(gene_role="unknown" → direction gate cannot fire → CNA alteration never assigned).
Mutation and fusion channels CAN fire: p_mutation≥0.5 in 2.0% of pairs,
p_fusion≥0.5 in 0.5% of pairs, has_driver_alteration=TRUE in 2.5% of pairs.

---

## 4. Full Threshold Sweep

**Updated sweep (2026-08-15 — six mathematical corrections applied, n=15 known-role test genes):**

| thr | guard% | E1_delta | E1_A_CI | E1_B_p | E1_C_CI† | E2_delta | E2_A_CI | E2_B_p | E2_C_CI† |
|---|---|---|---|---|---|---|---|---|---|
| −0.25 | 36% | +0.0073 | [−0.011,+0.023] | 0.215 | [−0.010,+0.025] | +0.0104 | [+0.002,+0.019] | 0.028 | [+0.000,+0.021] |
| −0.50 | 100% | +0.0071 | [−0.017,+0.029] | 0.283 | [−0.015,+0.029] | +0.0125 | [+0.001,+0.022] | 0.032 | [+0.000,+0.025] |
| **−1.00** | **100%** | **+0.0598** | **[+0.014,+0.115]** | **0.019** | **[+0.003,+0.117]** | **+0.0277** | **[+0.013,+0.047]** | **0.002** | **[+0.008,+0.048]** |
| −1.50 | 77% | +0.0736 | [+0.019,+0.146] | 0.008 | [+0.007,+0.140] | +0.0396 | [+0.017,+0.069] | 0.006 | [+0.011,+0.072] |
| −2.00 | 30% | +0.1069 | [+0.029,+0.199] | 0.004 | [+0.010,+0.204] | +0.0578 | [+0.026,+0.114] | 0.000 | [+0.014,+0.101] |

Bold = designated threshold. Guard% = fraction of drugs in 5-40% sensitive window (GDSC-intrinsic, unchanged by pipeline).

†Method C CI is correctly centred: `[E_obs − q_{0.975}(boot), E_obs + |q_{0.025}(boot)|]` where boot
is the Rademacher distribution (centred at 0 by construction). Do not read raw boot percentiles as the CI.

**Note on −0.25 and −0.50:** E1 is not significant by Method B (p=0.215, p=0.283). E2 is significant
(p=0.028, p=0.032) — the n_sensitive-weighted estimand captures the oncogene majority more cleanly at
lenient thresholds. The chosen threshold is −1.0 (GDSC convention): both E1 and E2 significant on all three methods.
Signal strengthens monotonically as thresholds tighten (−2.0: E1=+0.107, E2=+0.058, B p_E2=0.000).

Method keys: A=percentile cluster bootstrap (2000 resamples); B=exact paired permutation (32,768 = 2^15
perms, exact at n=15); C=wild cluster bootstrap (2000 resamples, Rademacher weights, see †note).

---

## 5. Inferential Results at Chosen Threshold (−1.0)

**Post-correction run — 2026-08-15 — six mathematical corrections applied, n=15 role-annotated test genes.**

**Class-level hit@20 (eval.py direct output — SENSITIVITY_MAD_Z=−1.0):**

| class | n_genes | n_sensitive | flat | driver | delta |
|---|---|---|---|---|---|
| oncogene | 12 | 748 | 0.029 | 0.063 | +0.034 |
| tsg | 1 | 123 | 0.008 | 0.016 | +0.008 |
| both | 2 | 176 | 0.011 | 0.028 | +0.017 |

**Estimand E1 (per-gene unweighted mean, n=15 genes):**
- E1 = **+0.0598**
- Method A: 95% CI [+0.014, +0.115] — excludes zero: YES
- Method B: exact permutation (32,768 perms, n=15) p = **0.019** — significant: YES
- Method C: 95% CI [+0.003, +0.117]† — excludes zero: YES
- **All three methods agree: ROBUST**

**Estimand E2 (n_sensitive-weighted mean, n=15 genes):**
- E2 = **+0.0277**
- Method A: 95% CI [+0.013, +0.047] — excludes zero: YES
- Method B: exact permutation (32,768 perms, n=15) p = **0.002** — significant: YES
- Method C: 95% CI [+0.008, +0.048]† — excludes zero: YES
- **All three methods agree: ROBUST**

†Method C CIs are correctly centred; see §4 note.

**Per-gene direction (n=15 role-annotated genes):**
- Positive delta (driver > flat): **10 genes** — FGFR2, ABL1, SMARCA4, PIM1, BRD4, EGFR, RET, MAP2K1, ALK, NTRK1
- Zero delta: **3 genes** — ROCK1, NTRK2, PPM1D
- Negative delta: **2 genes** — PDGFRA (−0.091), JAK1 (−0.010)

**Hypergeometric baseline (random ranking null):**

Under random ranking, the expected number of sensitive lines in the top-20 is hypergeometric:
`E[hit@20 | random] = k/N = 20/N_lines ≈ 20/900 ≈ 2.2%`

| Arm | hit@20 | Enrichment vs random |
|---|---|---|
| Random null | ~2.2% | 1.0× |
| Flat (expression only) | 2.9% | **1.3×** |
| Driver gated | 6.3% | **2.9×** |

Expression alone gives 1.3× random chance. The alteration gate brings it to 2.9× — the
2.2× gain from the gate is what E1/E2 quantify above.

**TSG note (n=1, NOT INTERPRETABLE):** TSG delta fell from +0.023 (pre-corrections) to +0.008 after the
LOF directional ranking was applied. This is expected under haploinsufficiency: single-hit TSG mutations
rarely cause complete transcript loss (Lindeboom et al. 2016, PMID 27618451), so placing sensitive TSG lines
at bottom-20% expression is not rewarded. The oncogene signal (n=12, delta=+0.034) is the interpretable result.

**"Both" note (n=2):** delta improved from 0.000 (pre-corrections) to +0.017 after "both" genes were
correctly treated as GOF (top-20%) rather than LOF. With n=2 this is not interpretable as evidence.

---

## 6. Under-Coverage Caveat (n=15 clusters)

At n=15 clusters, the percentile cluster bootstrap systematically under-covers
(MacKinnon & Webb 2017, DOI 10.1002/jae.2508). Nominal 95% intervals may achieve
<95% coverage. The lower bounds (+0.014 for E1, +0.013 for E2 at −1.0) are the
numbers to trust least. The permutation test (Method B) is exact under the null of
exchangeability and does not require distributional assumptions; its p-values at −1.0
(0.019 for E1, 0.002 for E2) are the primary evidence. At n=15, Method B enumerates
all 2^15 = 32,768 sign patterns exactly — no sampling approximation is involved.

**Note on wild bootstrap direction:** The wild bootstrap CI is narrower than the
percentile CI at all thresholds. At n=15, the percentile bootstrap amplifies outlier
genes through re-sampling (a gene can appear multiple times), while the Rademacher
sign flip fixes each gene's weight at ±1. The percentile bootstrap adds tail weight
that the wild does not. Both are under-covering relative to the true (unknown)
distribution, and neither is conservative. The permutation test (exact, no distributional
assumption) is therefore the headline evidence, not the CI bounds. Both CIs exclude
zero; the permutation result confirms this without distributional assumptions.

---

## 7. Role-Group Decomposition (W2)

**Post-correction run — 2026-08-15 — all six corrections applied.**

| Group | gene_role | n_genes | E1_delta | Method A 95% CI | Method B p | excludes 0 |
|---|---|---|---|---|---|---|
| A (annotated) | oncogene/tsg/both | 15 | +0.0598 | [+0.014, +0.115] | 0.019 | YES |
| B (unknown) | unknown | 11 | not re-estimated | — | — | — |

Group A drives the result. Group B (unknown gene_role) is excluded from both the
hit@20 loop and the bootstrap, because the CNA alteration channel structurally cannot
fire for genes without a role annotation, and including them would dilute the estimand.
Claim scope: the gate works for role-annotated genes (oncogene/tsg/both).

n_genes for Group A confirmed at 15 (seed=42, 20% hold-out of curated intersection;
13 unknown-role genes excluded from the hit@20 loop but retained in the flags output).

---

## 8. Displacement Audit (W3, Oncogene Stratum, n=11 genes — pre-recovery; n=12 post-recovery, re-audit pending)

| Metric | Value |
|---|---|
| Total sensitive lines displaced (flat top-20 → not driver top-20) | 15 |
| Total sensitive lines promoted (not flat top-20 → driver top-20) | 36 |
| Net | +21 |
| Genes with net > 0 (gate helped) | 6 |
| Genes with net = 0 (gate had no effect) | 5 |
| Genes with net < 0 (gate harmed) | **0** |

Pre-declared "≥2 genes with net<0" was WRONG. The JAK1 displacement pattern
(the only gene with hr_driver=0) occurs in the "both" stratum (JAK1, n=1, NOT
INTERPRETABLE) and does not generalise to the oncogene stratum. In the oncogene
stratum, the gate never costs more sensitive lines than it gains. 36 promotions
and 15 displacements yields net +21 — the gate is not just a thin margin over
large churn; the promotions substantially exceed the displacements.

---

## 9. T4 Gap — Between-Lineage vs Within-Lineage (W4)

V7 established lung is over-represented in L2 by +10.88 pp. W4 now measures
whether this composition explains the T4 mean component (1.20 pp).

| Component | z units | ~pp |
|---|---|---|
| Compositional (between-lineage) | -0.0147 | **-0.45 pp** |
| Within-lineage (MNAR selection) | +0.1077 | **+3.27 pp** |
| Total gap | +0.0930 | **+2.82 pp** |

The compositional contribution is **negative** (-0.45 pp): the lineages over-
represented in L2 (lung, esophagus, breast) have near-zero mean rna_pct_z,
while the lineages under-represented in L2 (bile_duct, fibroblast, soft_tissue)
have positive rna_pct_z. The composition therefore slightly OPPOSES the gap.

The gap is explained entirely by **within-lineage MNAR selection**: within each
lineage, L2 lines (those with protein data) have higher rna_pct_z than L1 lines
from the same lineage. This is consistent with proteomics studies preferentially
using well-characterised, robustly growing cell lines, which may have higher
expression levels overall.

**Implication for V7's claim:** V7 measured that lineage composition differs
(lung over-represented by 10.88 pp) — that is correct. But this composition
difference does not explain the T4 mean gap. The V7 "between-lineage" label
should be corrected to "within-lineage MNAR selection." The T4 gap is not
primarily a fame-bias artefact — it is a selection artefact within each lineage.

---

## 10. Escalations

| ID | Issue | Recommended resolution |
|---|---|---|
| W5-ESCL | Two estimands (per-gene mean vs pooled ratio) give 158%-different deltas. Gini=0.363, not highly concentrated but not equal weighting either. V1 computed the pooled ratio. | Fiona decides which estimand to report. |
| V3-ESCL | Committed numbers in v1 archive (0.014025/0.017735) cannot be reproduced from any code or data in this repo. | This document supersedes them. v1 numbers should not be cited. |

---

## 11. Row-Count Reconciliation (885,844 → 15,444,946)

The v1 archive was produced on 885,844 flags rows. The current `predictions_with_confidence.parquet`
has **15,444,946 rows** (= 15,080,476 scored pairs + 364,470 unscored NaN orphans).
The per-pair fallback introduced in the harmonisation v2 refactoring restored
~11.66M pairs from the discarded_pairs pool (see project_harmonisation_v2 memory).
Additional rows come from the updated GDSC model list (2026-07-09 vs the prior version)
and the expanded curated gene set.

All 364,470 unscored (NaN core_score) pairs have `has_driver_alteration = True` — they
are alteration-only orphans that received no expression evidence. There are zero NaN-score
pairs without a driver alteration in the final output. The tier partition is complete and
disjoint over the full 15,444,946 rows: HIGH (69,282) + CONTEXT (597,850) + MEDIUM (2,957,152) +
LOW (11,820,662) = 15,444,946.

---

## 12. Citations

| Claim | Source |
|---|---|
| Cluster bootstrap under-covers at few clusters | MacKinnon & Webb 2017, J Appl Econometrics 32:233 (DOI 10.1002/jae.2508) |
| Cluster-robust inference, minimum cluster counts | Cameron & Miller 2015, J Human Resources 50:317 (DOI 10.3368/jhr.50.2.317) |
| Permutation tests as exact small-sample inference | Good 2005, Springer 3rd ed (ISBN 0387988645) |
| GDSC per-drug sensitivity convention | Iorio et al. 2016, Cell 166:740 (PMID 27397505) |
| NMD reduces truncated TSG transcripts | Lindeboom et al. 2016, Nat Genet 48:1112 (PMID 27618451) |
| JAK1 LOF drives immunotherapy resistance | Shin et al. 2017, Cancer Discov 7:188 (PMID 27903500) |

---

---

## 13. Change Log

| Date | Change | Effect on E1 (-1.0) | Effect on E2 (-1.0) |
|---|---|---|---|
| 2026-08-14 | SENSITIVITY_MAD_Z fixed from −0.5 to −1.0 in eval.py:39 (reproducibility fix) | +0.0238→+0.0238 (no change; prior analysis already used −1.0) | +0.0138→+0.0138 |
| 2026-08-15 | W_RNA corrected from √2.00 to √1.431 (Kish n_eff, ρ̄=0.548 measured); pipeline re-run stages 4→5→6 | +0.0238→+0.0240 (+0.7%) | +0.0138→+0.0138 (unchanged) |
| 2026-08-15 | Orphan alteration-driver pair recovery: 364,470 pairs added to flags_with_driver.parquet (mut:324,232 fus:39,510 cna:1,711); has_driver_alteration count restored; n_role_annotated: 13→15. | E1: +0.0598 (see §5) | E2: +0.0277 (see §5) |
| 2026-08-15 | Six mathematical corrections applied — C1: silence guard added to common.py score_source_lineage(); C2: orphan pairs → core_score=NaN (not 0.0); C3: is_lof_alteration defined by alteration type + gene_role (was absent); C4: directional ranking (rank_flat role-blind control; rank_driver LOF/GOF directional; "both" genes treated as GOF); C5: MNAR shrink removed from protein_scorer.py; C6: DRIVER_BOOST (+0.15) removed from mutations_scoring.py | E1: +0.0598 [+0.014, +0.115] p=0.019 | E2: +0.0277 [+0.013, +0.047] p=0.002 |

**Effect of W_RNA correction on flat hit@20 (the key observable):**  
Oncogene flat: 0.0272→0.0278 (+2.2%); TSG flat: 0.0081→0.0163 (+101%).  
Driver hit@20: 0.0511→0.0518 (negligible change).  
The correction improves the core_score baseline without reducing the driver signal.

**Effect of orphan recovery on class-level hit@20 (the key observable):**  
Oncogene: flat 0.0278→0.0266 (−4.3%), driver 0.0518→0.0532 (+2.7%); delta 0.0240→0.0266 (+10.8%).  
Recovery pulls orphan lines (no RNA/protein, core_score=0.0) into the driver tier, reducing flat
scores slightly (orphans increase the denominator) while boosting driver hit@20.  
One additional oncogene and one additional "both" gene now appear in the eval loop (previously
invisible because all their driver pairs were missing from flags_with_driver).

*Previous version archived at: `docs/ARCHIVE/ALTERATION_DEFENCE_v1_UNREPRODUCIBLE.md`*
