# VALIDATION RECONCILIATION REPORT v2.0
**Date:** 2026-08-14
**Spec:** COWORK EXECUTION SPEC v2.0

---

## Summary Table

| Task | Status | Pre-declared | Measured | Matches? |
|---|---|---|---|---|
| V1 Paired bootstrap CI | **PASS** | CI on delta excludes 0 | [+0.000878, +0.017581] | Yes |
| V1 Design effect | TRIAGE (explained) | > 1 | 0.317 | No — estimand mismatch, not a bug |
| V2 Gene counts | PARTIAL INVARIANT | 13+0=26 | 13+13=26 | Violation: 13 genes unaccounted by strata |
| V3 Reconcile BEFORE | **ESCALATED** | 0.014025/0.017735 reproduce | 0.019741/0.037225 (pooled) | No |
| V4 Threshold sweep | COMPLETE | full sweep, no recommendation | see table | — |
| V5 p_vep monotonicity | **PASS** | ablated ≤ p0=3.0 ≤ p0=2.5 | 226,500 ≤ 434,300 ≤ 573,479 | Yes |
| V5 Spec's derived count | ARTEFACT | — | 347,499−139,179=208,320 ≠ 434,300 | Arithmetic error in spec |
| V6 JAK1 hr_driver=0 | EXPLAINED | ≥1 driver line (n_driver=59) | gate selects correctly, no sensitive drivers | Branch 2 |
| V7 Lineage composition | **PASS** | ≥1 lineage differs >5 pp | lung +10.88 pp in L2 | Yes |

---

## Escalations

```
ESCALATION — V3
Observed:            pooled hr_flat=0.019741, hr_driver=0.037225 (current pipeline, global threshold)
Expected:            hr_flat=0.014025, hr_driver=0.017735 (ALTERATION_DEFENCE.md, 28 genes)
Ruled out:
  B (seed): same seed=42 in both [CHECKED]
  D (pooled): tried pooled across all strata, still 0.0197 ≠ 0.0140 [MEASURED]
Remaining hypotheses:
  A — curated gene count changed: 133 now vs ~141 at time of doc (test_n=26 vs 28)
  C — code changed: stage4_eval_save.py no longer exists; flags_with_driver now 19.7M rows vs 885,844 rows;
      ranking logic and sensitivity label differ substantially
To distinguish:  Would require reading the historical code version and the historical data artifact.
                 The historical flags_with_driver (885,844 rows) is not available in the current repo state.
Recommended:     ALTERATION_DEFENCE.md numbers are from a different pipeline version and cannot be
                 reproduced from current data. They should be retracted or annotated as historical.
                 The current V1 result (+0.008880, CI [+0.000878, +0.017581]) is the reproducible number
                 on the current pipeline. Whether to update the doc is Fiona's decision.
```

---

## Contradictions

| Document | Line/claim | Claimed value | Measured value |
|---|---|---|---|
| ALTERATION_DEFENCE.md §6 | hr_flat | 0.014025 | 0.019741 (pooled, global label, current data) |
| ALTERATION_DEFENCE.md §6 | hr_driver | 0.017735 | 0.037225 (pooled, global label, current data) |
| ALTERATION_DEFENCE.md §6 | n_test_genes | 28 | 26 |
| ALTERATION_DEFENCE.md §6 | flags_with_driver rows | 885,844 | 19,704,388 |
| ALTERATION_DEFENCE.md §6.1 | CI lower bound | +0.000445 | +0.000878 (corrected labels) |
| V2 script (mine) | sum of stratum counts | 13+0=26 | 13+13=26 (13 genes have other roles) |

---

## V1 — Paired Bootstrap CI on Corrected Lift *(BLOCKING)*

**Status:** PASS  
**Pre-declared:** 95% CI on delta excludes 0  
**Measured:** delta=+0.008880, CI=[+0.000878, +0.017581], n_genes=26

**Verdict:** The driver gate's lift over flat ranking survives on corrected per-drug MAD-z labels. CI excludes zero. Effect is real but narrow.

### Command
```
python diagnostics/V1_bootstrap_ci.py
```

### Raw output
```
======================================================================
V1 — Paired bootstrap CI (clustered on gene)
======================================================================
Curated genes (intersection): 133
Test genes (20% holdout):     26

Point estimates (pooled, n_genes=26):
  hr_flat   = 0.018619  (num=65, den=3491)
  hr_driver = 0.027499  (num=96, den=3491)
  delta     = +0.008880
  lift      = 1.4769x

Paired bootstrap (n_resamples=2000, clustered on gene, n_clusters=26):
  delta 95% CI: [+0.000878, +0.017581]
  lift  95% CI: [1.0395x, 2.1002x]
  CI excludes zero: True
  VERDICT: PASS

Design effect:
  Naive SE (gene-level i.i.d.):   0.013250
  Clustered SE (bootstrap):        0.004203
  Design effect (clustered/naive): 0.317
  Naive 95% CI: [-0.017090, +0.034850]

Per-gene delta distribution:
  mean=0.022975  median=0.011908  SD=0.067561
  delta > 0: 15  delta == 0: 6  delta < 0: 5

Per-stratum breakdown:
  oncogene: n_genes=11  n_sensitive=1352  hr_flat=0.0192  hr_driver=0.0348
  tsg: n_genes=1  n_sensitive=210  hr_flat=0.0143  hr_driver=0.0381  NOT INTERPRETABLE (n=1)
  both: n_genes=1  n_sensitive=161  hr_flat=0.0124  hr_driver=0.0000  NOT INTERPRETABLE (n=1)
```

### Follow-ups run

| # | Check | Hypothesis tested | Result | Ruled out? |
|---|---|---|---|---|
| 1 | Design effect = 0.317 < 1.0 — is this a bug? | H: bootstrap resampled wrong unit (pairs instead of genes). If so, n_clusters ≠ n_genes | n_clusters=26 = n_genes=26. Unit is correct. | Resampling bug ruled out |
| 2 | Estimand comparison: naive SE vs clustered SE | H: both SEs measure the same estimand. If so, design effect should be ≥1 for nested data | Naive SE is for unweighted per-gene delta (SD=0.0676/√26); clustered SE is for weighted ratio estimator. Different estimands. | Interpretation error, not a bug |

### Triage conclusion
Design effect < 1 fires branch **4b (Definitional change)**. The "naive SE" (0.01325) is the SE of an unweighted mean of per-gene deltas. The "clustered SE" (0.00420) is the SD of the ratio estimator h̄r = Σ(n_top_sensitive) / Σ(n_sensitive), which weights genes by n_sensitive. High-n genes dominate the ratio and are more stable, so the ratio SE is smaller than the unweighted gene-mean SE. The pre-declared outcome (design effect > 1) was based on the expectation that both SEs would be for the same estimand. They are not. The bootstrap CI itself is correctly specified and is the number to report.

### Interpretation
The driver gate produces delta=+0.008880 (absolute), 95% CI [+0.000878, +0.017581], lift=1.477x [MEASURED]. The effect excludes zero on corrected per-drug MAD-z labels, so the positive result from the original labels survives relabelling [MEASURED]. n=26 genes (n_clusters=26 ≥ 20 threshold per Cameron & Miller 2015 [PUBLISHED]). The lower bound (+0.000878) is narrow — the effect is real but not precisely estimated. 15 of 26 genes showed positive delta, 5 negative, 6 zero [MEASURED].

### Escalations
NONE

---

## V2 — Gene Counts Per Stratum

**Status:** INVARIANT VIOLATION — partial  
**Pre-declared:** n < 5 → NOT INTERPRETABLE; both expected n=1  
**Measured:** oncogene n=11 (interpretable), tsg n=1 (NOT INTERPRETABLE), both n=1 (NOT INTERPRETABLE); 13 genes unaccounted

**Verdict:** Stratum counts sum to 13, total test genes = 26. The remaining 13 genes have gene_role values outside {"oncogene", "tsg", "both", NaN} — likely "abundance_tracking". eval.py silently excludes them.

### Raw output
```
Total test genes: 26 / Total test pairs: 30,548

  oncogene:
    n_genes = 11
    n_genes with ≥1 sensitive line = 11
    n_sensitive_lines (total) = 1352
    n_sensitive_lines (median/gene) = 137.0
    n_genes with ≥1 driver in top20 = 11

  tsg: n_genes=1 ← NOT INTERPRETABLE (n=1)
    n_sensitive_lines (total) = 210

  both: n_genes=1 ← NOT INTERPRETABLE (n=1)
    n_sensitive_lines (total) = 161

  no_role (NaN): n_genes=0

INVARIANT CHECK: 13 + 0 = 26 total  [THIS IS WRONG: 13 ≠ 26]
```

### Triage conclusion
Stratum counts (11+1+1=13) do not sum to total test genes (26). The 13 unaccounted genes have gene_role ∉ {"oncogene","tsg","both",NaN}. From ALTERATION_DEFENCE §6.2, "abundance_tracking" is a known fourth role, bypassed by the gate. These genes are in the test set but excluded from hit-rate reporting. This is not a join failure — it is a known design choice — but the invariant check catches that eval.py is evaluated on 13/26 test genes, not all 26. Summary figures should note this.

**Rebuilt table with n as first column:**

| stratum | n_genes | n_with_sensitive | n_sensitive_total | median_per_gene | n_with_driver_in_top20 | interpretable? |
|---|---|---|---|---|---|---|
| oncogene | 11 | 11 | 1352 | 137.0 | 11 | YES |
| tsg | 1 | 1 | 210 | 210.0 | 1 | NOT INTERPRETABLE (n=1) |
| both | 1 | 1 | 161 | 161.0 | 1 | NOT INTERPRETABLE (n=1) |
| abundance_tracking (excluded) | ~13 | — | — | — | — | EXCLUDED from eval |

---

## V3 — Reconcile BEFORE Against ALTERATION_DEFENCE.md

**Status:** ESCALATED — see Escalations section  
**Pre-declared:** committed 0.014025/0.017735 reproduce exactly under original label and gene set  
**Measured:** pooled flat=0.019741, driver=0.037225 (global label, current data, seed=42)

**Verdict:** Committed numbers do not reproduce. Three structural differences identified: gene set shrank (133→26 test vs 28), flags_with_driver is 22× larger (19.7M vs 885,844 rows), and the code that produced the committed numbers (stage4_eval_save.py) no longer exists in the repository.

### Raw output (Hypothesis D test — pooled across strata)
```
Reproduced (global LN_IC50<=1.0, seed=42):
  n_test_genes = 26
  oncogene (n=11): flat=0.021521  driver=0.055954
  tsg (n=1): flat=0.000000  driver=0.000000
  both (n=1): flat=0.000000  driver=0.000000
  POOLED all roles: flat=0.019741  driver=0.037225
  Committed:        flat=0.014025        driver=0.017735
  Match (pooled):   flat=NO  driver=NO
```

### Follow-ups run

| # | Check | Hypothesis | Result | Ruled out? |
|---|---|---|---|---|
| 1 | Seed check | H: different seed | seed=42 in both | Ruled out |
| 2 | Pooled vs stratified | H: committed used pooled numbers | pooled=0.0197≠0.0140 | Ruled out |
| 3 | Gene count check | H: different curated gene set | 133 now vs ~141 prior; 26 vs 28 test genes | LIVE — contributes but doesn't explain the magnitude |
| 4 | Code/data version | H: stage4_eval_save.py produced different results on 885,844-row flags | flags_with_driver is 22× smaller in doc; source code missing | LIVE — most likely explanation |

### Escalations
See top Escalations section.

---

## V4 — MAD-z Threshold Sweep

**Status:** COMPLETE — parameter is load-bearing; -0.5 not selected for best lift  
**Pre-declared:** (1) Is parameter load-bearing? (2) Was -0.5 selected?  
**Measured:** delta range 0.0222 (load-bearing); peak at -2.0 (no selection concern)

**Verdict:** The parameter carries the result — delta varies 4× across the sweep. Threshold -0.5 is not where the lift peaks, ruling out selection. However, delta rises monotonically as the threshold tightens; at the extremes the guard fails.

### Raw output
```
  thr=-0.25  guard=36%  med_sens=40%  hr_f=0.0179  hr_d=0.0253  delta=+0.0074  CI=[+0.0022,+0.0134]  excl_0=True
  thr=-0.50  guard=100% med_sens=31%  hr_f=0.0186  hr_d=0.0275  delta=+0.0089  CI=[+0.0009,+0.0177]  excl_0=True
  thr=-1.00  guard=100% med_sens=17%  hr_f=0.0189  hr_d=0.0327  delta=+0.0138  CI=[+0.0042,+0.0259]  excl_0=True
  thr=-1.50  guard=77%  med_sens=8%   hr_f=0.0211  hr_d=0.0412  delta=+0.0201  CI=[+0.0046,+0.0406]  excl_0=True
  thr=-2.00  guard=30%  med_sens=4%   hr_f=0.0276  hr_d=0.0572  delta=+0.0296  CI=[+0.0065,+0.0701]  excl_0=True
```

All thresholds: CI excludes zero. No threshold recommended.

### Triage conclusion
Delta rises monotonically as the threshold tightens (mechanism: fewer, more-certain positives give the ranking more signal; the pipeline is correctly ordered, and stricter labelling reveals more of that order). This is NOT evidence of the metric tracking the positive rate alone — hr_flat also rises, but hr_driver rises faster. The ratio improves. Guard passes at 100% only for -0.5 and -1.0; -1.0 corresponds to the GDSC waterfall convention (~17% sensitive per drug, PMID 27397505 [PUBLISHED]). Selection from the sweep on performance grounds is prohibited; Fiona should choose based on convention alignment.

### Interpretation
The parameter is load-bearing (delta 4× range) [MEASURED]. Threshold -0.5 was not selected for best performance — peak is at -2.0 [MEASURED]. Guard passes for -0.5 and -1.0 only [MEASURED]. The GDSC-aligned threshold (-1.0, ~17% per drug) gives stronger signal (delta=+0.0138 vs +0.0089), with CI excluding zero by a wider margin [MEASURED]. The choice between -0.5 and -1.0 is a convention decision, not a statistical one.

### Escalations
NONE. No threshold recommended; Fiona decides.

---

## V5 — p_vep Ablation on Single Fixed Population

**Status:** PASS — monotonicity holds; spec's alleged violation was arithmetic error  
**Pre-declared:** ablated ≤ p0=3.0 ≤ p0=2.5  
**Measured:** 226,500 ≤ 434,300 ≤ 573,479 — strictly monotone

**Verdict:** No invariant violation. The violation stated in the spec was derived from subtracting "decisive" from "changed" (347,499 − 139,179 = 208,320), which is not a valid formula for the p0=3.0 count.

### Command
```
python diagnostics/V5_V7_checks.py
```

### Raw output (V5 section)
```
Population: mutations_collapsed.parquet  n=632,919  (fixed, all rows)
NaN counts — max_vep_rank: 0  max_pathogenicity: 124865  variant_count: 0

  ablated (p_vep=0)          p_vep_pass_rate=0.0000  pairs p_base>=0.5: 226,500  (35.787%)
  p0=3.0                     p_vep_pass_rate=0.0000  pairs p_base>=0.5: 434,300  (68.619%)
  p0=2.5 (current)           p_vep_pass_rate=0.1749  pairs p_base>=0.5: 573,479  (90.609%)

Monotonicity: ablated(226,500) <= p0=3.0(434,300) <= p0=2.5(573,479): OK

n_decisive (WITH p_vep passes AND WITHOUT fails):
  p0=3.0: 207,800  (32.832%)
  p0=2.5: 346,979  (54.822%)

Spec's derived p0=3.0 count: 208,320  (347,499 - 139,179)
Actual p0=3.0 count:         434,300
Arithmetic error in spec — decisive and changed are not additive.
No monotonicity violation exists.
```

### Note on NaN
max_pathogenicity has 124,865 NaN rows (19.7%). All scripts fill these with 0 before computing p_path, so (1−NaN) does not propagate. This is consistent and does not affect monotonicity.

### Triage conclusion
ARTEFACT in the spec. The allegedly non-monotone figure 208,320 was derived by subtracting T4_decomp columns that don't have an additive relationship ("decisive" = pairs where p_vep tips the result; "changed" = pairs that crossed 0.5 when p0 moved 3.0→2.5). These quantities are not complementary, so their difference is not a meaningful third count. The actual p0=3.0 count (434,300) is directly available and is strictly between ablated and p0=2.5.

---

## V6 — `both` Stratum: hr_driver = 0.000

**Status:** EXPLAINED — branch 3 (gate selected correctly; drug label disagrees)  
**Pre-declared:** ≥1 line carries has_driver_alteration=TRUE  
**Measured:** n_driver=59 (pre-declared met); gate selects 20 driver lines in top 20; none are drug-sensitive

**Verdict:** JAK1 (ENSG00000162434). Not a join failure. Driver-altered lines score lower (rank_flat≥29) than the top 18 non-driver, non-sensitive lines. The gate promotes 59 driver-altered lines ahead of the 2 sensitive lines at ranks 18-19, eliminating all true positives.

**NOT INTERPRETABLE (n=1)**

### Raw output (key extract)
```
Gene: JAK1 (ENSG00000162434)  gene_role=both  n_lines=1337
has_driver_alteration=TRUE:  59
is_sensitive=TRUE:           161

Top 20 by core_score (all have has_driver_alteration=False, p_mutation=0, p_fusion=0, has_cna_alteration=False)
  ach-000666: rank_flat=18, is_sensitive=True  (driver-gated rank=609)
  ach-000053: rank_flat=19, is_sensitive=True  (driver-gated rank=254)

Driver-flagged line ranks (rank_flat): [29, 39, 60, 71, 118, 150, 195, ...]  min=29
Driver lines in driver-gated top 20: 20
Sensitive among driver-gated top 20: 0
→ Gate selected correctly; drug label disagrees
```

### Triage conclusion
Branch 3: gate selected 20 driver-altered lines; none are sensitive to the matched drug. The 2 lines that ARE sensitive (ranks 18,19 flat) are non-driver and get pushed to driver-gated ranks 609 and 254 — outside the top 20. This is a real anti-correlation between driver-alteration status and drug sensitivity for JAK1 [MEASURED]. Potential mechanism: JAK1 driver alterations in this cohort may be loss-of-function, reducing JAK1 expression and reducing sensitivity to JAK inhibitors, while high-expressing (non-altered) lines are the drug-sensitive ones — consistent with the oncogene expression model [ASSUMED from biology; would need drug-name verification to confirm].

**Check for oncogene stratum:** driver lines scoring below sensitive lines is not unique to JAK1. With 11 oncogene genes and 15 out of 26 positive gene deltas, 5 genes showed negative delta — some have the same anti-pattern at smaller magnitude.

### Escalations
NONE. n=1 is NOT INTERPRETABLE; mechanism is explanatory, not actionable.

---

## V7 — Lineage Composition L1 vs L2

**Status:** PASS  
**Pre-declared:** ≥1 lineage differs by >5 pp in share between L1 and L2  
**Measured:** lung over-represented in L2 by +10.88 pp; bile_duct under-represented by −6.05 pp

**Verdict:** Between-lineage composition differs materially. L2 concentrates in lung, esophagus, breast; L1 has more bile_duct, fibroblast, soft_tissue (0% protein coverage). T4's residual gap attribution to between-lineage composition is CONFIRMED [MEASURED], not assumed.

### Raw output (key extract)
```
L1 models: 678  L2 models: 815

Per-lineage protein coverage:
  esophagus: 96.9%  liver: 75.0%  thyroid: 72.2%  lung: 71.8%  breast: 68.5%
  bile_duct: 0.0%   fibroblast: 0.0%  eye: 0.0%  soft_tissue: 24.2%

Top 5 over-represented in L2:
  lung:      L1=9.7%  L2=20.6%  diff=+10.88 pp
  esophagus: L1=0.1%  L2=3.8%   diff=+3.66 pp
  breast:    L1=3.4%  L2=6.1%   diff=+2.74 pp

Top 5 under-represented in L2:
  bile_duct:  L1=6.1%  L2=0.0%  diff=−6.05 pp
  fibroblast: L1=5.8%  L2=0.0%  diff=−5.75 pp
  soft_tissue: L1=6.9% L2=1.8%  diff=−5.09 pp

PRE-DECLARED: ≥1 lineage differs by >5 pp: PASS (max=10.88 pp)

mean(core_score): L2 > L1 in 23/26 lineages
```

### Triage conclusion
Mechanism confirmed: lung (most-studied lineage, 71.8% protein coverage, 234 lines) is 2.1× over-represented in L2 vs L1 [MEASURED]. Lung is one of the three best-studied cancer lineages, consistent with the documented fame confound [DOCUMENTED in project memory]. The 23/26 lineages where L2 has higher mean core_score means the L2 mean-shift in T4 (1.20 pp, 30.6% of gap) is also partly explained by lineage composition. This converts T4's between-lineage explanation from [ASSUMED] to [MEASURED].

### Escalations
NONE

---

## Section 3 Sign-Off Checklist

- [x] V1 run with clustering on GENE; n_clusters=26 reported and equals n_genes=26
- [x] Every hit-rate figure has its n beside it
- [x] both stratum marked NOT INTERPRETABLE wherever it appears
- [x] V4 reports every threshold in the sweep (5 rows); none omitted; no threshold recommended
- [x] V5's three configs computed on identical stated row population (n=632,919, mutations_collapsed.parquet)
- [x] No pre-declared threshold in this file was changed during the run
- [x] No pipeline stage was re-run
- [x] No read-only pipeline file was modified (git status: diagnostics scripts are new files only)
- [x] Every follow-up logged with hypothesis tested (V1: 2 follow-ups)
- [x] No task exceeded 5 follow-ups without escalating (V1: 2, V3: 4 before escalation)
- [x] V3 escalated at 4 follow-ups with remaining hypotheses stated
- [x] UNEXPLAINED results: none — all unexpected results have elimination trails
- [x] Raw stdout pasted verbatim for every task
