# RISK_REGISTER

Every place the target architecture (C1–C6) would make a claim the current evidence does not support. Nothing here should be built on without the assumption being declared in writing first.

**Severity:** **CRITICAL** = the claim is measurably false today · **HIGH** = unsupported and load-bearing · **MEDIUM** = unsupported but contained · **LOW** = known, bounded.

---

## CRITICAL

### R1 — "The tool ranks cell lines well" · C2
**Claim implied:** an ordered list carries information about which lines are best.
**Evidence against:** ranker-region pAUC **0.5002 / 0.5015** (chance). On the 28 held-out genes **every** metric at **every** k has a 95% CI on excess-over-chance including zero. The only above-chance cell in the table is `recall@20` on the 113 **tuning** genes. VERIFIED.
**Risk if built anyway:** the product's headline function is not supported by any held-out measurement.
**Mitigation:** ship the lexicographic ranking of `BUILD_SPEC.md` C2 with `ordering_validated_to_level = 2` and a visible notice. Never sort by score alone without that label.

### R2 — "top-20 hit rate 93.6%" · C2, reporting
**Claim:** the `.docx` headline performance figure.
**Evidence against:** the metric is *any-hit@20*, whose **hypergeometric chance baseline is 90.8%**. Recomputed on the current build the observed value is **89.4% — below chance**. VERIFIED.
**Risk if built anyway:** a published figure that a reviewer can falsify in one calculation.
**Mitigation:** withdraw the figure. Standing rule 3: no @k number without its baseline. See `DECISIONS_REQUIRED.md` Q2.

### R3 — "not measured" is reported as "measured negative" · C3
**Claim:** absence of a flag means absence of the thing.
**Evidence against:** **99.69%** of scored pairs (29,689,395 of 29,781,274) have no CNA measurement at all, yet `has_cna_alteration = False` via `.fillna(False)` at `build_full_predictions.py:122` and `build_evidence_ledger.py:68`. Eighteen such sites are enumerated in `BUILD_SPEC.md` C3. VERIFIED.
**Risk if built anyway:** every "this line has no alteration in this gene" statement is, for the overwhelming majority of pairs, a statement that nothing was checked. This is §9 gap 13 in its most concrete form.
**Mitigation:** three-state representation everywhere plus assertions A7/A8. Non-negotiable prerequisite for C3.

### R4 — the evaluation harness scores untested pairs as true negatives · C2
**Claim:** the @k figures measure ranking quality.
**Evidence against:** `is_sensitive.fillna(False)` at `stage4_eval_save.py:41`, `06_held_out_eval.ipynb` cell 4 and `build_gene_regime.py:150`. Every (gene, line) pair not tested in GDSC2 enters the denominator as a confirmed negative. VERIFIED.
**Risk if built anyway:** the metric is partly measuring GDSC2's testing pattern, not the ranker. This is *why* an @k number can look high while carrying no information.
**Mitigation:** restrict denominators to tested pairs, or model the missingness. Report the tested fraction beside every figure.

---

## HIGH

### R5 — combining DepMap and HPA as independent evidence · C1, C2 level 2
**Claim:** two RNA sources agreeing is two confirmations.
**Evidence against:** median per-gene Spearman **ρ = 0.817** (IQR 0.700–0.903, 99.5% of genes positive). An uncorrected Stouffer combination inflates the combined *z* by **1.35×**. VERIFIED.
**Risk if built anyway:** the same error class as §2.3, at a larger magnitude than the ρ_EP = 0.46 already acknowledged. Level 2 of the proposed ranking would systematically over-count.
**Mitigation:** Strube/Hartung generalised Stouffer with per-gene R and shrinkage; count RNA consensus as **one** confirmation, capping the count at 4 not 5.

### R6 — no multiple-testing correction anywhere · B2.3, gap 3
**Claim:** the per-gene labels in `chronos_validation.parquet` are significant findings.
**Evidence against:** **16,866** simultaneous tests at uncorrected p < 0.05. ~843 genes expected to pass by chance; 4,332 pass. Zero correction code exists repo-wide. VERIFIED.
**Risk if built anyway:** `weak_positive` / `weak_negative` (1,918 genes) are assigned on significance alone with no effect-size floor and no correction. The 31 `validated` / 6 `inverted` survive an effect floor and are safer.
**Mitigation:** BH-FDR across genes; or state plainly that only the effect-size-gated labels are used, and stop emitting the weak bands.

### R7 — abundance→dependency runs weakly *negative*, contradicting §C.15 · B1.1
**Claim:** §C.15's gate finding implies a positive abundance→dependency relationship.
**Evidence against:** of 4,332 significant genes, **57.0% negative** (2,468 vs 1,864), binomial **p = 4.46e-20**. On a pan-essential control set the median ρ is **−0.0655**, 79.3% negative. The screen-quality confound is ruled out (Spearman = −0.019, p = 0.58). VERIFIED.
**Risk if built anyway:** an unexplained sign contradiction sits under the gate result the whole architecture leans on.
**Note:** this is **not** a sign bug. The convention was checked against DepMap's own reference sets (essentials median −5.198 vs non-essentials +3.492, p = 3.29e-229) and is correct and consistently applied at all three sites. The contradiction is real and unexplained, which is worse than a bug.
**Mitigation:** declare as an open limitation. Do not describe the abundance→dependency relationship as positive anywhere.

### R8 — `n_layers == 2 → "high"` confidence · B2.1, gap 14
**Claim:** two-layer coverage means higher confidence.
**Evidence against / for:** §C.18 measures the protein layer as a net negative **on global AUROC** (−0.0086 GDSC / −0.0042 Chronos) but a net *positive* on **pAUC**, the metric that matches a top-N recommender (+0.0234 GDSC / +0.0171 Chronos — `PROTEOMICS_DEFENCE.md` §6.2). The rule is live at `build_full_predictions.py:151` and `05_confidence_tiers.ipynb` cell 5. VERIFIED, recharacterised — the direction the claim was checked against was the wrong metric for this product.
**Risk if built anyway:** the confidence tier is directionally consistent with the pAUC result but has never been calibrated — no check exists that "high"-tier rows are actually more often correct than "moderate"-tier ones, across 5,277,064 two-layer rows. `n_layers == 2` is also measured as 1.37 effective layers (Kish), not 2, which the tier rule does not account for.
**Mitigation:** calibrate the tier rule against top-N hit-rate outcomes by tier before relying on it; do not drop `n_layers` from the rule on the basis of the AUROC number alone.

### R9 — protein-derived claims generalise to the panel · C6, gap 13
**Claim:** protein evidence describes the cell-line panel.
**Evidence against:** lines with proteomics have higher mean RNA — **Cliff's δ = +0.251**, Cohen's d = +0.445, MWU p = 6.68e-17. Both platforms show it at nearly the same magnitude, so it is a property of *line selection for MS*, not of the assay. VERIFIED.
**Risk if built anyway:** absolute quantification (C6) would deliver precision on a biased subset and invite it to be read as representativeness.
**Mitigation:** every protein-derived output carries a coverage statement and a `not_measured` state. Never impute.

### R10 — the ruler's internal standard is imprecise · C6
**Claim:** copies-per-cell estimates are precise enough to order lines.
**Evidence against:** histone/total intensity ratio **CV = 0.411** across 948 lines (median 0.0606, IQR 0.0454–0.0816). Not a depth artefact — Spearman(ratio, depth) = −0.066, non-monotone by quintile. VERIFIED.
**Risk if built anyway:** ±41% propagates multiplicatively into every estimate. Usable for an order-of-magnitude gate; **not** usable as level 3 of a tie-break, which is how the brief proposes to use it.
**Mitigation:** frame C6 as a second gate with an absolute threshold. Report bands, never point estimates.

---

## MEDIUM

### R11 — ProCan's upstream normalisation is undocumented · C6
Per-line medians spread **0.98 log2 units ≈ 2× linear**, so scale partly survives — but the file is `..._averaged_...` and no README, provenance JSON or manifest records what ProCan applied. `UNVERIFIED — could not locate`. Blocks any absolute claim. See `DECISIONS_REQUIRED.md` Q3.

### R12 — depth adequacy cannot be checked on the published criterion · C6
Wiśniewski et al. give a **~12,000 peptide** stabilisation threshold. The matrix is protein-level only (median 5,218 proteins/line); peptide counts are absent. The criterion cannot be applied. Do not claim it is met. `UNVERIFIED`.

### R13 — 30% of the panel cannot be ploidy-corrected · C6, gap 8
Ploidy available for **1,502 of 2,145** models (70.0%). Median ploidy **2.653**; **55.8%** above 2.5; only ~26% near-diploid. Assuming diploidy misestimates DNA per cell by a median **1.33×** and up to 2.7×, **correlated with lineage**. VERIFIED. The other 30% must carry `not_ploidy_correctable`, never a diploid default. Note `AMP_THRESHOLD = 2.5` currently sits *below* the median line's average copy number.

### R14 — the guarded lineage scorer is unavailable and untested · C5
`02_transcriptonomics.ipynb` is **not in this repo**; nor is any `silent_lineage` / `MAD_FLOOR` / `MIN_PEERS_FOR_LINEAGE` implementation. The unguarded result is neutral (Δ = **−0.0018**, p = 0.540, better in 47.8% of genes), reproducing §C.15. The guards may well work — the EGFR-in-blood argument is sound — but the claim is `UNVALIDATED` here. Do not adopt lineage conditioning for scoring on the strength of an untested guard.

### R15 — the §C.15 lineage comparison has no interval · C5
The stored result reports p = 0.540 and a Δ, with **no paired bootstrap CI**. Standing rule 4 makes it unquotable as-is. Cheap to fix; should be done regardless of R14.

### R16 — 17% of panel models have no lineage label · C5
**364 of 2,145** panel models (6.4% of the 1,746 scored) carry no tissue label. Two unreconciled lineage vocabularies coexist (`gdsc_models.tissue`, `sample_info.lineage`). VERIFIED. Any lineage-stratified split, cluster bootstrap or display must handle `lineage_not_recorded` explicitly.

### R17 — composite identifiers are still present · B2.4
**15** ProCan accessions map to semicolon-joined composite symbols (`a8mya2 → CXorf49;CXorf49B`, `p12532 → CKMT1A;CKMT1B`, `p23610 → F8A1;F8A2;F8A3`, …), plus 1 in CCLE. **1,633** ProCan and **420** CCLE proteins are lost at the join. The detecting audit is read-only; no build-time assertion blocks them. VERIFIED. C3's assertion A2 fails today.

### R18 — the protein-coding universe is four different numbers · B2.5
19,213 / 19,187 / 20,163 / 19,177 depending on which script is speaking, with **no version column on either source**. VERIFIED. Any "N protein-coding genes" claim is ambiguous by up to 986. See `DECISIONS_REQUIRED.md` Q5.

### R19 — DepMap release version is unrecorded · C3
No release/quarter string in the filename, parquet metadata or any manifest. `UNVERIFIED — could not locate`. C3's assertion A11 fails on day one — deliberately. A traceability claim cannot be made until this is recorded.

### R20 — "Chronos" in this codebase means "Project Score" · A1.6
`chronos_long.parquet` has median +3.05, range −55.1…+29.1 — a Bayes-factor scale, not Chronos gene effect. The real Chronos HDF5s (155 MB) are on disk and **unused**. VERIFIED. The sign convention is nonetheless correct. Every mention of "Chronos" in code and docs is mislabelled, and an independent second screen is being left on the table. See `DECISIONS_REQUIRED.md` Q6.

### R21 — `chronos_pct` cannot support the claims read from it · B1.1
`chronos_pct` is ranked **within gene**, so mean `chronos_pct` is **0.4995 for reference essentials and 0.4995 for reference non-essentials** — it carries **zero** between-gene information by construction. VERIFIED. Any downstream reading of `chronos_check` as evidence about *whether a gene is essential* is reading a quantity that cannot carry that information.

---

## LOW

### R22 — max-aggregation extreme-value bias · C4, gap 4
`.groupby(["ensg_id","model_id"]).max()` live in both Track C scoring notebooks. `E[max] = m/(m+1)`, 8.7% exposure. Any continuous alteration score inherits it in either direction. Mitigation: boolean gate, or a count-calibrated statistic.

### R23 — unswept constants decide the headline results · A3
`DEP_THRESHOLD = −0.5` defines the positive class for **every** gate and ranker figure the project quotes, at **8 independent sites**, unswept. `LIFT_DEFAULT = 2.0` decides the routing class of all 141 curated genes, unswept. `MIN_LINES` and `MIN_POS/MIN_NEG` are **inconsistent across scripts** (100 vs 150; 10/10 vs 8/8). The superseded `RHO_THRESHOLD = 0.1` survives in two live scripts after being deliberately replaced by 0.30 elsewhere. `FLOOR_RATE = 0.0139` / `FLOOR_BASE = 0.0169` are **hardcoded results of a previous run used as constants**. VERIFIED. Standing rule 2 is currently not met anywhere.

### R24 — the held-out split is unstratified · gaps 7, 17
`rng.choice` over genes with no lineage stratification and no cluster bootstrap, against a panel ranging from 382 haematopoietic lines to 1 adrenal. VERIFIED. Held-out intervals are probably narrower than the truth — which matters, because several already include zero.

### R25 — dead code carries a latent asymmetry · B1.2
`scoring_variants.combine_three()` has **zero callers**; `core_score.parquet` is a hardcoded `0.5·E + 0.5·P` (5,277,064 two-layer rows; 24,504,210 one-layer; **no** three-layer stratum). Had it been called, the E–C and P–C pairs would **not** be (½,½) — with the stored ρ values they work out to ≈(0.448, 0.552). The asymmetry would appear silently the moment a third layer is wired in. VERIFIED.

### R26 — a stated invariant is violated in two scripts · C3
`scoring_variants.py` documents "missing layer contributes no term (not zero)", but `test_run_stouffer_regime_split.py:356` and `test_run_stratum_centring.py:287` both compute `d.pE.fillna(0) + d.pP.fillna(0)`. VERIFIED.

### R27 — `04_driver_gated_routing.ipynb` cannot run as committed · A2.2
It reads `layer3_continuous.parquet`, whose only copy on disk is in `discarded/outputs/`. The Stage 3 chain is broken. VERIFIED.

### R28 — regression testing rests on one anchor · gap 12
`rank_cell_lines.py:122-131` asserts on A375 only, at a top-5% tolerance. A single anchor cannot distinguish recalibration from regression.

### R29 — CNA residual exceeds the label noise floor · gap 19
`test_run_cna_gate_results.json`: residual 0.0304 against label FPR 0.00928, verdict `residual_above_noise__check_calls`. The CNA gate's effect is real but partly unexplained; the calls themselves need checking before the gate is leaned on.

### R30 — noisy-OR probabilities are user-visible · B2.2, gap 2
`explain_pair.py:178-206` emits raw `p_mutation` / `p_fusion` floats into user output via `cli.py:89`, tagged "does not affect stratum_rank". They are display-only for *ordering*, but they are user-facing probabilities produced by an untested independence assumption. VERIFIED.

---

## Assumptions that would be undeclared unless written here

1. That within-gene percentiles are comparable across sources with different transforms (`log2(TPM+1)`, `nTPM`, RMA). Rank-based combination makes this defensible; **level-based** combination would not be.
2. That GDSC2 drug sensitivity is a valid proxy for "this cell line is a good model for this gene". Never tested in this repo.
3. That `DEP_THRESHOLD = −0.5` marks a biologically meaningful dependency boundary. Unswept, at 8 sites.
4. That the 141 curated genes are representative of the 19,177 scored genes. Never tested; they are drug-target genes, which is a strongly selected set.
5. That COSMIC CGC membership is a complete driver annotation. 768 genes have evidence; ~18,400 are "not drivers" by `fillna`.
6. That the 28 held-out genes support any point estimate. They do not — standing rule 4.
