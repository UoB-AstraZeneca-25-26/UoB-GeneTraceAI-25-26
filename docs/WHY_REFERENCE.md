# The "Why" Reference — Every Clip, Threshold, and Design Choice

This document answers, for every clip bound, shrinkage constant, threshold, and
combination-weight in `transcriptomics.py`, `06_protein_score_M2.py` /
`platform_tier.py`, and `core_score.py`: **what it is, why this value, what
happens if you change it, and whether it has been empirically tested.**

It supersedes the parameter-by-parameter content of `docs/SCORING_METHODS_FULL.md`
§1–§3 (that document's narrative walkthrough and citation table remain the
reference for *how the pipeline works*; this document is the reference for
*why each number is what it is*). It is organized file-by-file, stage-by-stage,
to mirror the structure of the companion mathematical-audit report — cross-reference
by stage name.

**Status tags**, consistent with this project's existing discipline
(`docs/RISK_REGISTER.md`, `docs/MATH_REFERENCE.md`):
- **PROVEN** — derived from this project's own data via a stated, reproducible procedure, or backed by a cited external result that applies directly.
- **CONTESTED** — has a stated rationale, but the audit or this task found a reason to doubt it fully applies here (wrong reference population, self-acknowledged as unresolved, or a plausible alternative wasn't ruled out).
- **UNVALIDATED** — a first-pass convention with no calibration or citation behind the specific value; stated plainly as such, not dressed up.

**Update (this task).** Two corrections are folded in below that change the record from `SCORING_METHODS_FULL.md`'s existing text:
1. A claim that RNA's inter-source agreement uses Pearson while protein's inter-platform agreement uses Spearman is **not accurate** — both use Spearman (§3.5 below). `SCORING_METHODS_FULL.md` did not itself contain this claim, but nothing in it corrects the record either; §3.5 here is the authoritative statement.
2. Finding **C7** (the `n_layers=1`/`n_layers=2` scale mismatch) has a verified, ready-to-promote fix (`core_score_v7_fallback_fix.py`, not yet promoted) — §3.6 reflects the corrected formula and reports the fix's status as **candidate, not live**.

**Update (promotion task, 2026-08-29).** Both C7 and C3 were promoted into live `core_score.py` in this task, each with the full backup/replace/run/regenerate-downstream-chain/spot-check-via-real-HTTP-call procedure. §3.2 and §3.6 below, and the Part C summary table, are updated to reflect final, live status rather than candidate status. Key result: **C3 did NOT become a fifth instance of the "upstream improvement fails to compose downstream" pattern** — the full downstream spike check (the step explicitly left undone in the prior task) came back flat-to-slightly-improved, not worse. This is the pattern's first *clean* composition among the cases checked in this project, worth noting precisely because the pattern had held in every prior case.

**Update (promotion task, 2026-08-31).** Two further findings promoted into live `core_score.py` and, for the first time, `confidence_tiers.py` — the latter extends this document's scope beyond its original file list (§4 is new). §3.7 (new) covers a `.stack()` NaN-row-survival bug in `core_score.py`, discovered independently of the mathematical audit that produced C3/C7 (found while auditing this project's own test/validation suite, not the scoring code directly). §4 (new) covers a companion default-tier fix in `confidence_tiers.py`, made necessary in part by §3.7's finding and in part by an entirely separate, legitimate NaN source (`driver_routing.py`'s orphan-recovery step) that §3.7's fix does not and should not touch.

---

## 1. `transcriptomics.py`

### 1.1 Detection floor
| | |
|---|---|
| **What** | `TPM_FLOOR = 1.0` (HPA/DepMap axes: absolute floor `log2(1+1)=1.0` in log2 space). GEO axes: `ax.floor = quantile(vv, floor_q)`, a per-method, data-driven floor. |
| **Why** | HPA's own convention: TPM ≥ 1 as the detection cutoff (Uhlén et al. 2015). GEO methods have no comparable convention (they aren't TPM-scaled), so each gets its own empirical floor instead. |
| **Sensitivity** | Not tested in isolation here; the HPA/DepMap floor is a literature constant, not a tuned one. |
| **Status** | **PROVEN** for HPA/DepMap (cited). **CONTESTED** for the mechanism generally — see `floor_q` below, which sets *where in the GEO distribution* that floor sits. |

### 1.2 `floor_q`, `n0`, `upper_cut`, `lower_cut` (calibrated `Params`)
| | |
|---|---|
| **What** | Live values (`chosen_params_v2_widened.json`): `floor_q=0.3, n0=25, upper_cut=0.5, lower_cut=0.3`. |
| **Why** | Chosen by `calibrate_params()` — grid search scored by held-out cross-source Spearman ρ (refit on non-held-out sources, compared to a fixed-reference held-out score, ambiguous detection band only), 1-SE rule tie-broken toward the conservative (more rank-based) end. |
| **Sensitivity** | Per this project's own earlier overtuning-check work: `upper_cut` was found **well-identified** by the calibration surface (the held-out concordance score is genuinely sensitive to it); `n0`, `floor_q`, and `lower_cut` were found **weakly- or edge-identified** (the scoring surface is comparatively flat across a range of these, so the specific chosen value carries less evidential weight than `upper_cut`'s does). |
| **Status** | **CONTESTED**, upgraded from the audit's plain PASS on the calibration *mechanism*: the live values (`floor_q=0.3`, `lower_cut=0.3`) fall **outside** `PARAM_GRID`'s own searchable range (`floor_q∈{0.02,…,0.20}`, `lower_cut∈{0.02,…,0.20}`, `transcriptomics.py:61-66`). Either the grid was widened for the run that actually produced these values and never synced back into the constant, or the constant is stale — either way, `PARAM_GRID`/`calibrate_params()` cannot currently reproduce or refine the parameters actually in production, and `05_transcriptomics_stats_layer.py:145` confirms the live pipeline never calls `calibrate_params()` at all (it loads the JSON directly). The calibration *procedure* is sound; the *infrastructure* is inconsistent with what's live. |

### 1.3 `MIN_PEERS=5`, `MAD_FLOOR=0.384`
| | |
|---|---|
| **What** | `MIN_PEERS`: minimum non-NaN lines in a (source, lineage) group before any z-score is computed for it (`transcriptomics.py`'s `method_z`, imported live from `utils/common.py`). `MAD_FLOOR`: absolute floor on the blended scale to prevent division-by-near-zero. |
| **Why** | `common.py`'s own comment: `MAD_FLOOR` is "calibrated against the shared transcriptomics notebook." No equivalent statement exists for why `MIN_PEERS=5` specifically (vs. 3, 10, etc.). |
| **Sensitivity** | Not swept in this task. |
| **Status** | **CONTESTED**, not because the values are unreasonable but because of their provenance: `utils/common.py`'s own computational surface (`robust_z_matrix`, `score_source_lineage`, `fetch_depmap/hpa/geo`) is dead code relative to the live pipeline — it serves only `protein_scorer.py`, which `run_all.py:44-47` confirms was superseded by `run_protein_m2_stage.py`. `MAD_FLOOR`'s calibration notebook is external and unverifiable from code; the module that housed its derivation context is no longer exercised by anything live. |

### 1.4 Scale-detection heuristic thresholds (100, 50, 25)
| | |
|---|---|
| **What** | `detect_scale_from_values`: `linear` if `hi>100` or `hi/med>50`; else `log` if `hi<25`, else `linear`. |
| **Why** | No comment, citation, or calibration anywhere in the file. |
| **Sensitivity** | Untested. |
| **Status** | **UNVALIDATED.** Plausible domain reasoning (log2(TPM) rarely exceeds ~25) but presented with no more support than that. |

### 1.5 Robust-z shrinkage (`w = nb/(n0+n0)`, `n0=25`)
| | |
|---|---|
| **What** | Lineage MAD blended toward the method-global MAD: `scale = sqrt(w·local² + (1−w)·global²)`, `w = nb/(nb+n0)`. |
| **Why** | James-Stein-style partial pooling: small lineages should defer to the global scale, large ones should trust their own. |
| **Sensitivity** | Verified explicitly in the audit: `w→0` as `nb→0` (fully global), `w→1` as `nb→∞` (fully local) — the formula behaves exactly as claimed at both limits. |
| **Status** | **PROVEN** as a shrinkage *mechanism*. **CONTESTED** in two specific respects found by the audit: (a) the weight uses `nb` (total lineage cell-line count) rather than `nfin` (the gene's own per-lineage finite-value count, computed separately in the same function for other purposes) — a gene with heavy missingness within an otherwise-large lineage is shrunk as if fully observed; (b) the *center* (median) uses a binary `nb>=MIN_PEERS` switch rather than the same continuous blend applied to scale — an unexplained asymmetry between how the two parameters of the same distribution are protected against small samples. |

### 1.6 Rank inverse-normal transform
| | |
|---|---|
| **What** | `Zb = norm.ppf((R-0.5)/nfin)`, commented "(Blom / Van der Waerden form)". |
| **Why** | Standard technique for converting ranks to normal quantiles when the parametric (robust-z) branch isn't trusted. |
| **Sensitivity** | Not separately tested; converges to the same limit as Blom/van der Waerden for large n. |
| **Status** | **PROVEN** as a technique, **CONTESTED** as documented: the formula `(r−0.5)/n` is the *Hazen* plotting position, not literally Blom (`(r−3/8)/(n+1/4)`) or van der Waerden (`r/(n+1)`) as the comment claims. Practically immaterial (all three agree closely for reasonably large n); a documentation-accuracy defect, not a correctness one. |

### 1.7 Stouffer combination (two-stage: GEO methods → sources → final)
| | |
|---|---|
| **What** | `z_combined = Σ(w·z)/√Σw²`; GEO stage weights by `√n_samples`; final 3-source stage uses equal weight (1.0 each). |
| **Why** | Coverage-based weighting for combining differently-sized studies (Mosteller & Bush 1954, per `SCORING_METHODS_FULL.md`'s citation table); equal weight across the 3 final sources is a deliberate "one source, one vote" choice, stated rationale: more robust to any single source's artefacts than a flat pooled Stouffer would be. |
| **Sensitivity** | Verified algebraically: this normalization produces exactly unit variance when inputs are independent standard normals — the general property the whole pipeline leans on repeatedly. |
| **Status** | **PROVEN** (both the math and the coverage-weighting citation). The equal-weight choice for the final 3-source stage is a **reasoned design decision** with a stated rationale, not independently benchmarked against a coverage-weighted alternative at that specific stage — flagged as untested-alternative, not wrong. |

---

## 2. `06_protein_score_M2.py` / `platform_tier.py`

### 2.1 Isoform collapse (max vs. mean; concordance threshold)
| | |
|---|---|
| **What** | Min pairwise Spearman ρ across all isoform pairs; if it clears a data-derived threshold, gene value = **max** abundance across isoforms; else **mean**. |
| **Why** | Weakest-link concordance (every pair must agree, not just the best pair) — consistent with this project's MIN-based combination philosophy elsewhere. Threshold: 95th percentile of a permutation null (random cross-gene isoform-pair correlations), a stated statistical convention. |
| **Sensitivity** | Not swept here. |
| **Status** | Concordance-detection mechanism: **PROVEN** (sound null construction, `n_pairs=1500` adequately sampled, power-gated minimum overlap). The **max-over-sum** choice for the concordant case: **UNVALIDATED** — no stated reason to prefer max over sum/mean when concordant isoforms may represent genuinely distinct co-expressed protein products. |

### 2.2 `robust_z_matrix` scale fallback chain (MAD → std → `1.0`)
| | |
|---|---|
| **What** | `scale = 1.4826·MAD`, falling back to `std` if MAD≈0, falling back to a bare `1.0` if both are degenerate. |
| **Why** | No comment given for the final `1.0`. |
| **Sensitivity** | Fires only on fully degenerate (near-constant) columns — rare, but not quantified in this task. |
| **Status** | **UNVALIDATED.** Serves the identical purpose as RNA's `MAD_FLOOR=0.384`, but with an invented, uncalibrated number instead of a data-derived one. |

### 2.3 Platform-agreement thresholds (`agree_low=0.266`, `agree_high=0.425`) and `MIN_SHARED`
| | |
|---|---|
| **What** | Confidence ramp bounds from `derive_thresholds()`: null (mismatched-gene) vs. matched (same-gene) ProCAN/CCLE correlation distributions, percentile-based. Separately, `platform_tier.py` hardcodes its own `MIN_SHARED=44` for the tiering pass, distinct from `06_protein_score_M2.py`'s own power-derived `min_shared=55` (`m2_thresholds.json`) for the *residualization*-feeding agreement threshold. |
| **Why** | `agree_low`/`agree_high`: data-derived, matched-vs-null percentile design (documented, sound — see audit B4). `MIN_SHARED=44`: **no comment or derivation given anywhere in `platform_tier.py`.** |
| **Sensitivity** | Not swept here. |
| **Status** | `agree_low`/`agree_high`/`min_shared=55`: **PROVEN** (derived, with a stated procedure). `platform_tier.py`'s `MIN_SHARED=44`: **UNVALIDATED** — a same-concept, different-value, uncited sibling of the `min_shared=55` used one file over for a similar purpose. Not previously flagged in the audit; found while assembling this document. |

### 2.4 Tier cutoffs (`consistent≥0.50`, `cautious 0.30–0.50`, `conflicting<0.30`)
| | |
|---|---|
| **What** | `platform_tier.py`'s `_assign_tier`. |
| **Why** | Round-number conventions, not derived from the null/matched distributions the way `agree_low`/`agree_high` are. |
| **Sensitivity** | `MATH_REFERENCE.md` §11.4 already measured this directly: at the minimum shared-line count, the 95% CI on ρ=0.30 spans **[−0.068, +0.596]** — wider than the entire 0.20-unit tier bin. ~2,000 proteins sit within ±0.10 of a tier boundary — tier assignment there is dominated by sampling noise, not a real difference in agreement. |
| **Status** | **CONTESTED** — reasonable round-number conventions, but the project's own sensitivity analysis already shows the bins are frequently mis-assignable at the point of minimum eligible data. |

### 2.5 Coverage-weighted combination + confidence scoping
| | |
|---|---|
| **What** | `z_comb = Σ(√src_n·z)/√Σsrc_n`; confidence multiplier applied only when `n_sources>=2`. |
| **Why** | Same Stouffer-normalization principle as RNA's combination; single-source rows have no second platform to disagree with, so confidence down-weighting doesn't apply to them. |
| **Sensitivity** | Verified algebraically correct (proper `Σwz/√Σw²` form). |
| **Status** | **PROVEN.** `n_sources>=2` scoping reconfirmed correctly implemented (audit B7) — single-source rows pass `z_raw` through unmodified. |

---

## 3. `core_score.py`

### 3.1 RNA per-stratum standardization (`RNA_STRATUM_MIN_N=5`, `RNA_STRATUM_N0=25`, `RNA_WINSOR_BOUND=5.0`)
| | |
|---|---|
| **What** | Per-(gene, n_sources) median/SD, thin strata (`n<5`) blended toward the gene-global median/SD with `w=n/(n+25)`, clipped at ±5. |
| **Why (N0=25, MIN_N=5)** | Explicitly acknowledged in-file as first-pass, not independently derived for this context — reused from `transcriptomics.py`'s own calibrated `n0` "by convention." |
| **Why (WINSOR_BOUND=5.0)** | **This one has real cited evidence**, unlike most bounds in this pipeline: the file's own docstring states it was "observed to bind on <0.2% of rows post-standardization, so this rarely activates" — an empirical winsorization-rate check, not a guess. |
| **Sensitivity** | This task tested the scale choice specifically (below, §C1 test): switching stratum SD from `std(ddof=1)` to a MAD-based scale gives `corr=0.995` against the live version, median `|diff|=0.044`, but a real tail (`p99=0.53`, `1.13%` of rows differ by `>0.5`, `0.19%` by `>1.0`); MAD-based standardization also hits the ±5 winsor bound **3× more often** (0.14% vs 0.04%). |
| **Status** | Winsor bound: **PROVEN** (cited empirical rate). `N0=25`/`MIN_N=5`: **UNVALIDATED**, self-acknowledged. Non-robust `std(ddof=1)` scale choice: **CONTESTED** — the audit flagged the inconsistency with RNA's own upstream MAD-based philosophy; this task's isolated test confirms the practical difference is small-but-real (not negligible enough to dismiss on measurement alone, but not large enough that the current choice is clearly wrong either). See Part C recommendation below. |

### 3.2 Protein residualization (`RHO_PRIOR=0.353`, `N_MIN_RHO=30`, `K_RHO=30`, `beta_g` formula)
| | |
|---|---|
| **What** | `rho_raw = pearsonr(prot_z, rna_z)` per gene; `rho_g = lam·RHO_PRIOR + (1−lam)·rho_raw`, `lam=K_RHO/(n+K_RHO)`; `beta_g = rho_g·(sd_prot/sd_rna)`. |
| **Why (RHO_PRIOR=0.353)** | Stated as "median per-gene RNA-protein correlation across the panel" — a genuine, data-derived prior, confirmed via `SCORING_METHODS_FULL.md`: **this uses Pearson correlation**, correctly described. |
| **Why (N_MIN_RHO=K_RHO=30)** | No specific derivation beyond "minimum shared lines to estimate a per-gene rho" — a round, plausible, first-pass minimum, not derived from a power calculation the way protein's own `seed_overlap`/`min_shared` are. |
| **Sensitivity** | EB blend limits verified explicitly: `n→0` gives `rho_g→RHO_PRIOR` (fully prior), `n→∞` gives `rho_g→rho_raw` (fully data) — correct James-Stein-style behavior. `beta_g` formula verified as the standard, correct OLS-slope-from-correlation identity (`β=ρ·σ_Y/σ_X`), not an invented substitute. |
| **Status** | `RHO_PRIOR`: **PROVEN.** `N_MIN_RHO`/`K_RHO=30`: **UNVALIDATED**, plausible convention. `beta_g` formula: **PROVEN** algebraically. The unshrunk-SD issue ("C3") is **FIXED AND PROMOTED (2026-08-29)** — `sd_prot`/`sd_rna` are now shrunk toward their panel-wide values with the same `w=n/(n+n0)` shape used elsewhere in this file (`SD_SHRINK_N0=25`, reused convention, itself **UNVALIDATED** for this specific application — shrinking an SD that feeds a ratio inside a slope, not a direct standardization — flagged the same way every other `n0` in this file is). See §3.6-adjacent promotion note below and the Part C summary table. |

### 3.3 Protein-residual shrinkage (`PROT_RESID_N0=25`, `PROT_RESID_WINSOR_BOUND=5.0`)
| | |
|---|---|
| **What** | Gene's own residual SD blended toward the panel-wide residual SD, `w=n/(n+25)`, continuous (not gated behind a hard cutoff the way RNA's is). |
| **Why (N0=25)** | Explicitly "SAME first-pass value as RNA_STRATUM_N0, reused by convention, NOT independently re-derived." |
| **Why (WINSOR_BOUND=5.0)** | Matches RNA's bound because `prot_resid_z` is standardized to the same SD=1.00 target — "a 5-sigma bound is an equally rare event on either arm under that shared calibration," per the file's own docstring, and confirmed (not blindly assumed) by the printed post-standardization SD check. |
| **Sensitivity** | Not independently swept in this task beyond the C3 test above. |
| **Status** | Winsor bound: **PROVEN** (checked, not assumed). `N0=25`: **UNVALIDATED**, self-acknowledged. This fix (continuous shrinkage) is itself an already-evidenced, settled correction — audit finding C4, PASS. |

### 3.4 Regime-3 detection and substitution (`REGIME3_SIG_THRESH=1.0`)
| | |
|---|---|
| **What** | `\|rna_std_z\|>1 ∧ \|prot_std_z\|>1 ∧ opposite sign` → substitute the direct-average (Variant C) formula for that row only. |
| **Why** | Matches the original regime-discovery task's own "signaled" cutoff exactly, reused for consistency between discovery and promotion rather than re-tuned. `|z|>1` is a natural, round threshold, not independently optimized via ROC/power analysis for this specific decision. |
| **Sensitivity** | Verified before promotion (per `core_score.py`'s own docstring, re-confirmed in the discovery-rerun task): population 75,117/6,762,902 rows reproduced exactly; spike improved 3.0389%→3.0313% (n_layers=2), 2.8473%→2.8444% (full panel). |
| **Status** | **PROVEN** as a promoted, verified fix. The specific threshold value (`1.0`) is a **reasonable convention**, not independently re-optimized — worth naming, not worth re-litigating (already settled per this task's own scope boundary). |

### 3.5 `W_RNA=√1.431`, `W_PROT=√1.45` combination weights — Pearson/Spearman correction
| | |
|---|---|
| **What** | Kish effective-sample-size weights: `n_eff = N/(1+(N−1)·ρ̄)`. RNA: `N=3, ρ̄=0.548 → n_eff=1.431`. Protein: `N=2, ρ̄=0.373 → n_eff=1.45`. |
| **Correlation estimator used (verified this task)** | **Both use Spearman, not a Pearson/Spearman split.** RNA's `ρ̄=0.548` comes from `diagnostics/X2_rna_rho.py` — a **one-time, non-live** diagnostic (not part of `run_all.py`), which computes `scipy.stats.spearmanr` pairwise between DepMap/HPA/GEO on 2,000 sampled genes (median of the three pairwise ρ's: DepMap×HPA 0.793, DepMap×GEO 0.417, HPA×GEO 0.435 → average 0.548). Protein's `ρ=0.373` comes from `platform_tier.py:97` — **live**, run every `run_all.py` execution as part of Stage 2 — which also computes `scipy.stats.spearmanr(val_procan, val_ccle)` per protein, median across ~5,957 proteins. **Neither value is read back into `core_score.py` at runtime** — both `W_RNA`/`W_PROT` are frozen historical literals baked into the source, not live-recomputed. |
| **Where Pearson actually appears** | `core_score.py`'s own `rho_raw = stats.pearsonr(prot_z, rna_z)` (§3.2 above) — this is a **different quantity entirely**: the RNA-vs-protein correlation feeding `RHO_PRIOR`/`beta_g`, unrelated to `W_RNA`/`W_PROT`. A prior explainer conflated these three distinct correlations (RNA-internal-source agreement, protein-internal-platform agreement, and RNA-vs-protein agreement) into an incorrect "RNA=Pearson, protein=Spearman inconsistency" framing. **That framing is retracted here.** The recommendation built on it ("switch to Spearman for consistency") is moot, since there is no Pearson-vs-Spearman split between the two arms to reconcile — both already use the same estimator. (This task searched `docs/SCORING_METHODS_FULL.md`, `docs/MATH_REFERENCE.md`, and `docs/Testing_and_Validation_Report.docx` for the specific document making this claim and did not find it verbatim in any of them; if it exists in a document outside this repo, this section is the corrected text to carry into it.) |
| **Sensitivity** | Verified algebraically: `(W_RNA²+W_PROT²)/W_NORM² = 1` exactly, under the assumption that `rna_std_z` and `prot_resid_z` are independent unit-variance normals. |
| **Status** | Kish arithmetic itself: **PROVEN** (correctly computed, correctly verified here). **CONTESTED** on two separate grounds: (a) both weights measure *within-arm* source/platform agreement, then get repurposed as the weight for a *different* combination problem (arm-vs-arm) with no independent derivation for that specific question — `core_score.py`'s own docstring already flags this as unresolved, not a decision this task revisits; (b) the "independent unit-variance normals" assumption underlying the unit-variance proof is only approximate in practice, because `prot_resid_z` is a residual built from an EB-*shrunk* `beta_g`, not the true OLS slope, so exact orthogonality to `rna_std_z` isn't guaranteed at the sample level. |

### 3.6 `n_layers=1` fallback scaling — finding C7, fix candidate (not yet promoted)
| | |
|---|---|
| **What (live, unchanged)** | `core_score = Φ(rna_std_z)` for protein-free lines / RNA-only genes — no `W_RNA/W_NORM` scaling. |
| **What (fix candidate, `core_score_v7_fallback_fix.py`)** | `core_score = Φ((W_RNA/W_NORM)·rna_std_z)`, treating "no protein arm" as "protein residual = 0" rather than "skip arm weighting". |
| **Why the live version is wrong** | `n_layers=2` rows downweight `rna_std_z` by `W_RNA/W_NORM≈0.7048` before Φ; `n_layers=1` rows don't — the same RNA evidence scores systematically more extreme via the fallback path, and both row types are ranked against each other in the same `(ensg_id, lineage)` stratum. |
| **Fix verification (this task)** | `n_layers=2` rows bit-identical (max\|diff\|=0.0) between live and fix. `n_layers=1` scores are pulled toward 0.5 symmetrically (not uniformly "lower" — high old scores move down, low old scores move up equally, since the correction is a shrinkage toward the center, not a one-directional cut). Full-panel spike (`<0.01` or `>0.99`): **2.84% → 1.48%**. `n_layers=1`-only spike: **2.73% → 0.49%** — more than a 5× reduction. In a 150-gene stratum-rank sample, 5,738 `n_layers=1` rows previously outranked the best `n_layers=2` row in their own stratum; 24.5% no longer do after the fix (the remainder are winsor-capped, genuinely very strong RNA signals that legitimately stay near the top even after correct weighting — not a residual bug). Concrete example: ENSG00000214309/breast/ach-001396 — old 0.886 (ranked above the stratum's best n_layers=2 row at 0.803), new 0.802 (now ranked just below it). |
| **Status** | **CONFIRMED BUG, FIXED, PROMOTED 2026-08-29.** Final live-verified numbers (re-confirmed against unchanged upstream inputs immediately before promotion, then reconfirmed bit-identical after promotion): `n_layers=2` bit-identical (max diff 0.0); full-panel spike **2.8444% → 1.4790%**; `n_layers=1`-only spike **2.7252% → 0.4894%**. Worked example ENSG00000214309/breast/ach-001396: **0.8858 → 0.8021**, confirmed end-to-end via a real `/gene/detail` HTTP call against the promoted, server-loaded data (`score: 0.802077`, gene symbol MBLAC1). Downstream chain (`driver_routing.py` → `confidence_tiers.py`) regenerated and structurally sound (row counts, null counts, tier-partition assertions all unchanged/passing). Pre-promotion state backed up as `core_score_v4_pre_fallback_fix_20260829_041807.py.bak` + `final_pipeline/outputs/pre_fallback_fix_backup_20260829_041807/`. Lambda redeploy is a separate, deliberate decision not made in this task. |

### 3.7 `.stack()` NaN-row-survival bug — finding, fix, promoted 2026-08-31
| | |
|---|---|
| **What (bug)** | `shared_models`/`shared_genes`, and the two `n_layers=1` fallback blocks' own row/column choices, are each built from the FULL row/column universe of a wide-pivoted table, not the true per-(model, gene) overlap of real measurements. This is harmless as long as a cell lacking real data produces NaN and that NaN drops out before the row is written — which `.stack()` used to do by default. In the pandas version this project now runs (3.0.1), `.stack()`'s default implementation no longer drops NaN rows. `core_score` itself was never fabricated — traced end to end for a concrete example (`ach-001333` × CD86), every intermediate value (`rna_std_z`, `prot_resid_z`, `core_z`, `core_pct`) was correctly `NaN`; only row *survival* was the bug, not score *computation*. |
| **Discovery path** | Not found by the mathematical audit that produced C3/C7. Found while auditing this project's own `Testing/run_tests.py` / `Validation/eval.py` / `run_gene_tests.py` test suite for unexamined assumptions (a separate task) — specifically, while checking whether `run_gene_tests.py`'s stale "CD86: no protein layer" / "CD3E: sparse data" expectations were still accurate against current data. |
| **Exact scope (measured before the fix)** | 9,188,147 / 26,558,496 rows (34.60%) of the live output. Correlation between a gene's real-row coverage and its phantom fraction: Pearson r = **−1.000** (mechanically deterministic, not stochastic) — coverage-decile phantom rate ranged 97.4% (least-covered decile) to 0.29% (best-covered). Confirmed **NOT** reaching real users: `ranking.py` and `Validation/eval.py` already filter on `core_score.notna()`/`.dropna()`. Confirmed reaching `run_gene_tests.py`'s printed tier breakdowns and `predictions_with_confidence.parquet` for any ungated `confidence_tier`/`n_layers` grouping (e.g. CD86 showed 1,504 "LOW" rows, of which 1,397 were phantom). |
| **Why (`future_stack` behavior change)** | Confirmed empirically, not assumed: `.stack(dropna=True)` under the current default (`future_stack=True`) raises `ValueError: dropna must be unspecified as the new implementation does not introduce rows of NA values` — that message refers only to alignment-introduced NaNs, not the genuine pre-existing NaN cells this grid has. `.stack(future_stack=False, dropna=True)` (the old implementation) does restore the drop, but pandas itself flags that code path for removal in a future version, so it was rejected as trading one version-fragile assumption for another. |
| **Fix** | `.dropna()` immediately on the stacked Series, before `.reset_index()`, at all three sites (main combination; protein-free-line fallback; RNA-only-gene fallback). Works identically regardless of which `.stack()` implementation pandas defaults to next, since it operates on the Series `.stack()` already produced. |
| **Verification** | Row count **17,370,349** (exact match to the independently-derived "real-row" count from the diagnostic task). All 9,188,147 previously-phantom rows confirmed absent (not reset to another value). Every real row bit-identical across `core_score`, `n_layers`, `combination_variant`, `lineage`, `stratum_rank` — outer-merge indicator 100% "both", 0 left/right-only, max abs `core_score` diff 0.0. Gene/model universe unchanged (16,992 genes, 1,563 models, 0 genes lost entirely). Spot-checked via a real `/gene/detail` HTTP call: `ach-001333`×CD86 now correctly returns HTTP 422 ("No prediction for 'CD86' + 'ACH-001333'") instead of a stale/absent-data ambiguity; a real CD86 row (`ach-000702`) correctly returns `"total":134` — CD86's true population size, not the old phantom-inflated count. |
| **Status** | **CONFIRMED BUG, FIXED, PROMOTED 2026-08-31.** Pre-promotion state backed up as `core_score_v6_pre_nan_fix_20260831_112601.py.bak` + matching `.parquet.bak` snapshots of `core_score.parquet`, `gene_dispersion.parquet`, `flags_with_driver.parquet`, `predictions_with_confidence.parquet`. Downstream chain (`driver_routing.py` → `confidence_tiers.py`) regenerated; see §4 below for the companion fix this promotion also required. Lambda redeploy is a separate, deliberate decision not made in this task. |

---

## 4. `confidence_tiers.py`

This file was outside this document's original scope (stated in the opening paragraph above: `transcriptomics.py`, `06_protein_score_M2.py`/`platform_tier.py`, `core_score.py`); added here because its default-tier defect was found and fixed in the same promotion as §3.7, and is closely coupled to it.

### 4.1 `NO_EVIDENCE` default-tier fix — finding, fix, promoted 2026-08-31
| | |
|---|---|
| **What (bug)** | `confidence_tier = "LOW"` was an unconditional default before any score-based override. For a `core_score.isna()` row, `score_rank_pct` is also NaN, so `top` (`score_rank_pct >= 0.80`) evaluates `False` — comparisons against NaN are always False, never an error — so a NaN row was indistinguishable, downstream, from a row genuinely measured and found NOT in the top 20%: "no evidence" silently collapsed into "measured and unremarkable." |
| **Two independent NaN sources feeding this file** | (1) §3.7's `.stack()` bug — fixed at the source, no longer produces NaN rows as of that fix. (2) `driver_routing.py`'s own, separate, **intentional** recovery step: it adds back alteration-driver pairs with no expression score at all, deliberately setting `core_score=NaN, n_layers=0`. **Confirmed directly, not assumed**: re-running `driver_routing.py` against the §3.7-fixed `core_score.parquet` still produces 333,417 such rows, 100% correctly tagged `n_layers=0` (0 residual-bug rows) — this fix is necessary regardless of §3.7, not made redundant by it. |
| **Why the naive early-default fix would have failed** | Because `has_driver_alteration` is `True` by construction for every recovery-step row, an early `NO_EVIDENCE` default would be silently **overwritten** by the later `CONTEXT` assignment — `~directional_top & drv` evaluates `True` for every driver-flagged NaN row, since `directional_top` is `False` for NaN inputs by the same NaN-comparison-is-False mechanism as `top`. The fix applies the `NO_EVIDENCE` override **last**, after the existing HIGH/MEDIUM/CONTEXT/LOW logic runs completely unchanged. |
| **Verification** | Verified in isolation (against the pre-§3.7 live `flags_with_driver.parquet`): 9,326,483 `NO_EVIDENCE` rows = 9,188,147 (§3.7's phantom count) + 138,336 (orphans already recovered even against the buggy file) — exact arithmetic match. Verified combined (both fixes together): `NO_EVIDENCE` = 333,417 (exactly the legitimate orphan-recovery count). **LOW/MEDIUM/CONTEXT/HIGH counts for real rows are byte-identical between isolation and combined** (13,630,183 / 3,406,451 / 255,900 / 77,815 both times) — confirms the fix isolates NaN rows correctly regardless of source. Partition assertions extended (not loosened): HIGH+CONTEXT and MEDIUM+LOW invariants now scoped to scored rows only; new assertions confirm `NO_EVIDENCE` is exactly (nothing more, nothing less) the `core_score.isna()` population. |
| **Status** | **CONFIRMED BUG, FIXED, PROMOTED 2026-08-31.** Live re-run (unmodified `run_gene_tests.py`, real promoted pipeline) confirms CD86's tier breakdown went from `{CONTEXT:36, HIGH:3, LOW:1504, MEDIUM:24}` (1,397 of the 1,504 "LOW" rows phantom) to `{HIGH:3, LOW:107, MEDIUM:24, NO_EVIDENCE:36}` — every real row's tier unchanged, every phantom row either correctly relabeled (if a legitimate driver-evidence orphan) or simply absent (if not). CD3E: `{CONTEXT:22, LOW:1507, MEDIUM:34}` → `{LOW:133, MEDIUM:34, NO_EVIDENCE:22}`. GAPDH (spot-checked as an unaffected control): mean/std unchanged to 3 decimal places (0.479/0.253 both before and after); only 2 of its rows moved from CONTEXT to NO_EVIDENCE. Pre-promotion state backed up as `confidence_tiers_v1_pre_no_evidence_fix_20260831_112601.py.bak` + matching `predictions_with_confidence...parquet.bak`. Lambda redeploy is a separate, deliberate decision not made in this task.

---

## Part C summary — promotion outcomes

| Finding | Measured impact | Outcome |
|---|---|---|
| **C1** — non-robust `std(ddof=1)` scale in RNA stratum standardization | corr=0.995 vs. MAD-based; median \|diff\|=0.044; but 1.13% of rows differ by >0.5 and MAD-based hits the winsor bound 3× more often | **NOT PROMOTED — recommendation-support only, unchanged.** Modest measured impact; lower priority than C3 was. Still open. |
| **C3** — shrunk `rho_g` × unshrunk `sd_prot`/`sd_rna` in `beta_g` | Isolated: near `N_MIN_RHO=30` (30–60 shared lines), `beta_g` shifted by median 19.9%, mean 72.1%, p99 1074%. **Reproduced against the real implementation** in `core_score_v8_beta_shrinkage.py`. Full downstream check (the step explicitly left undone previously): full-panel spike **1.4790% → 1.4769%**, `n_layers=2` spike **3.0313% → 3.0257%** — both flat-to-slightly-improved, not worse. Full-panel `core_score` mean/quantiles unchanged to 3 decimal places (aggregate effect is small because only ~339/11,460 qualifying genes sit near the threshold), while individual near-threshold genes shifted meaningfully per-line (e.g. `ENSG00000137463`/ach-000939: 0.8609→0.8244, confirmed via real `/gene/detail` call returning `score: 0.824384`, gene symbol MGARP). | **PROMOTED 2026-08-29 — downstream spike improved from 1.4790% to 1.4769% (full panel) and 3.0313% to 3.0257% (n_layers=2).** This is the *fourth confirmed instance* of the "one half of a paired estimate gets shrinkage, the other doesn't" pattern in this project, and notably **not** a fifth instance of the separate "upstream improvement fails to compose downstream" pattern (four prior instances existed before this task; C3 composed cleanly). Pre-promotion state backed up as `core_score_v5_pre_beta_shrinkage_fix_20260829_050249.py.bak` + `final_pipeline/outputs/pre_beta_shrinkage_fix_backup_20260829_050249/`. |
| **§3.7** — `.stack()` no longer drops NaN by default (pandas 3.0.1) | 9,188,147 / 26,558,496 rows (34.60%) were phantom (NaN `core_score`, real-looking `n_layers`); coverage-vs-phantom-fraction correlation r=−1.000 (deterministic). Fix verified bit-identical on all 17,370,349 real rows (max diff 0.0); all phantom rows confirmed absent, not reset. | **PROMOTED 2026-08-31** — see §3.7 above for full detail. Backed up as `core_score_v6_pre_nan_fix_20260831_112601.py.bak` + 4 matching `.parquet.bak` outputs. |
| **§4.1** — `confidence_tiers.py` defaulted `core_score.isna()` rows to `LOW` | 9,326,483 rows mislabeled in isolation (pre-§3.7); 333,417 legitimate `NO_EVIDENCE` rows remain even combined with §3.7 (driver_routing.py's own orphan-recovery is a separate, valid NaN source). Real-row tier counts byte-identical isolation vs. combined. | **PROMOTED 2026-08-31** — see §4.1 above. Backed up as `confidence_tiers_v1_pre_no_evidence_fix_20260831_112601.py.bak` + matching `predictions_with_confidence...parquet.bak`. |

---

*Companion documents: `docs/SCORING_METHODS_FULL.md` (mechanism walkthrough), `docs/MATH_REFERENCE.md` (broader project-wide statistical reference), the mathematical-audit report (this task's predecessor, file-by-file pass/fail/flagged table).*
