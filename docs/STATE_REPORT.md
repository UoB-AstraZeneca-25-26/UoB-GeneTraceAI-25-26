# STATE_REPORT — GeneTraceAI inventory and document-vs-code verification

**Scope:** Parts A and B of the Cowork brief. Repo root `C:\Disertation\UoB-GeneTraceAI-25-26`, branch `chai`, audited 2026-08-09.

**Evidence tags** — `PROVEN` (published, cited) / `VERIFIED` (computed here) / `CONTESTED` (sources disagree) / `UNVALIDATED` (asserted, no measurement) / `UNVERIFIED — could not locate`.

**Reproduction:** every measured number below comes from a script in `MEASUREMENTS.md`. No pipeline file was modified.

---

## 0. Headline findings

| # | Finding | Tag |
|---|---|---|
| 1 | The `.docx` headline "top-20 93.6%" uses an *any-hit@k* metric whose hypergeometric chance baseline is **90.8%**. The figure is at chance. | VERIFIED |
| 2 | On the 28 held-out genes, **every** metric at every k has a 95% CI on excess-over-chance that includes zero. | VERIFIED |
| 3 | The Chronos sign convention is **correct** and applied identically at all three sites. The brief's hypothesis of an inverted sign is not supported. | VERIFIED |
| 4 | `scoring_variants.combine_three()` — the `w = 1/(1+ρ)` code path — has **zero callers**. It is dead code; `core_score.parquet` is built by a hardcoded `0.5·E + 0.5·P`. | VERIFIED |
| 5 | DepMap and HPA RNA correlate at **median per-gene Spearman ρ = 0.817**. Any independence-assuming combination of them is invalid. | VERIFIED |
| 6 | ProCan **DIA** proteomics is on disk, carries **log2 intensities** (not ratios), and contains **17 histone proteins**, several quantified in 100% of 952 lines. The proteomic ruler is not blocked at step 1. | VERIFIED |
| 7 | Per-line **ploidy exists** (`gdsc_models.parquet:ploidy_wes/ploidy_wgs`, 1,502 of 2,145 models). §C.17's premise that it must be sourced is stale. | VERIFIED |
| 8 | **Paralog annotation exists** (3,552,265-row Ensembl BioMart table + `build_paralog_flags.py`). §C.17's claim that it is absent is stale. | VERIFIED |
| 9 | No multiple-testing correction exists anywhere in `src/`. Gap 3 still open. | VERIFIED |
| 10 | The teammate's `02_transcriptonomics.ipynb` is **not in this repo**. C5's decisive experiment cannot be run. | UNVERIFIED — could not locate |

---

# PART A — Inventory

## A1. Data assets

### A1.1 Expression sources

| Source | Path | Shape | Scale | Raw available? |
|---|---|---|---|---|
| DepMap | `data/parquet/2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.parquet` | 1,495 × 53,962 | `log2(TPM+1)`, observed range 0.00–17.78 | **No** — only the transform is on disk. Untransformed TPM is recoverable arithmetically (`2^x − 1`), but the pre-log source file is not present. |
| HPA | `data/parquet/1_4_hpa_rna_celline.parquet` | 24,315,372 × 6 (long) | Three columns shipped: `TPM`, `pTPM`, `nTPM`. Observed `nTPM` range 0.00–421,173.10 | **Yes** — raw TPM is present alongside the normalised nTPM. |
| GEO | `data/parquet/3_GEOexpression.parquet` | 19,914 genes × 3,267 samples | Observed range 1.74–43,460.70. Non-zero floor at 1.74 is consistent with RMA/background-floored microarray. | **Mixed** — GEO is an archive of heterogeneous submissions; per-series processing state is not recorded in a single column. |

- DepMap release: `UNVERIFIED — could not locate` a release/version string in the filename, in the parquet metadata, or in a manifest. The file is `OmicsExpressionAllGenesTPMLogp1Profile`, a DepMap Omics naming convention, but the quarter is not recorded anywhere in the repo. **This is a traceability defect (C3).**
- Coverage after mapping to `model_id` (VERIFIED): DepMap **1,479** lines × 19,173 protein-coding genes; HPA **1,104** × 19,209; GEO **588** × 15,962.
- DepMap headers are `SYMBOL (ENSG…)`; 53,961 of 53,962 columns carry an ENSG, 1 does not.

### A1.2 Proteomics — **both datasets are on disk** (high-priority question)

| | Nusinow / CCLE Gygi **TMT** | ProCan **DIA** |
|---|---|---|
| Path | `data/parquet/4_Harmonized_MS_CCLE_Gygi_subsetted.parquet` | `data/proteomics_procan/Protein_matrix_averaged_20250211.tsv` (pipeline copy `src/pipeline/outputs/procan_proteomics.parquet`) |
| Shape | 375 lines × 12,558 proteins | 950 lines × 8,453 proteins (pipeline copy 952 × 8,459) |
| **Ratios or intensities** | **Log-ratios.** median −0.042, **52.6% of values negative**, sd 0.96, range −23.05…12.78 | **Log2 intensities.** median 3.44, only **1.40% negative**, range −5.45…15.65 |
| Cross-sample normalisation | **Yes, median-centred.** Per-line median = −0.041, sd **0.088**, range [−0.636, 0.109] | **Not hard-scaled.** Per-line median = 3.436, sd 0.212, p1–p99 spread **0.98 log2 units ≈ 2× linear** |
| Missing | 27.8% | 38.8% |
| Depth | median **9,072** proteins/line (8,096–10,425) | median **5,218** proteins/line (0–6,226; one all-missing row) |
| Histones present | **23** histone proteins | **17** histone proteins |
| Absolute quantification | **Not possible** — ratios to a bridge channel, no absolute scale to recover | **Possible in principle** — see C6 in `BUILD_SPEC.md` |

All tags VERIFIED.

Histone detail (ProCan, quantified fraction of 952 lines): `P62805`→H4C16 **100.0%**, `O75367`→MACROH2A1 **100.0%**, `Q5TEC6`→H3-7 99.9%, `Q8IUE6`→H2AC21 99.5%, `Q92522`→H1-10 99.5%, `P16104`→H2AX 95.3%, `P10412`→H1-4 93.9%, `P16401`→H1-5 92.2%, `P07305`→H1-0 89.0%, `Q9P0M6`→MACROH2A2 87.6%, `Q02539`→H1-1 70.7%, plus six sparser entries.

**Caveat (VERIFIED, load-bearing):** several histone accessions are *ambiguous* — `P62805` maps to 14 H4C genes, `P68431` to 10 H3C genes, `P04908` to H2AC4/H2AC8 (`test_run_identifier_ambiguity_audit_results.json`, `T2_Q1`). For the ruler this is tolerable (the ruler wants *total* histone signal, so collapsing paralogues is correct), but the same accessions must never be used for per-gene histone claims.

### A1.3 Ploidy — **present**

`src/pipeline/outputs/gdsc_models.parquet`, columns `ploidy_wes` (n=1,248) and `ploidy_wgs` (n=254), sourced from `data/GDSC/model_list_20260709.csv`. Any estimate available for **1,502 of 2,145** panel models (70.0%). Median 2.653. VERIFIED. Distribution in `MEASUREMENTS.md` §B3.4.

### A1.4 Mutations, fusions, CNA

| Table | Path | Rows | Notes |
|---|---|---|---|
| Mutations (raw) | `data/parquet/6_OmicsSomaticMutationsProfile.parquet` | 1,066,869 × 70 | carries VEP consequence, pathogenicity, hotspot, VAF |
| Mutations (collapsed) | `cleaned_track_data/mutations_collapsed.parquet` | 632,919 × 13 | `max_vep_rank`, `max_pathogenicity`, `variant_burden`, `oncogene_hit`, `tsg_hit`, `max_vaf` |
| Mutations (variant detail) | `cleaned_track_data/mutations_variant_detail.parquet` | 704,713 × 18 | `proteinchange`, `hotspot`, `vepimpact` |
| Fusions | `data/parquet/5_OmicsFusionFilteredSupplementary.parquet` | 184,237 × 30 | `ffpm`, `confidence`, `totalreadssupportingfusion` |
| Fusions (gene-level) | `cleaned_track_data/fusions_gene_level.parquet` | 144,221 × 9 | |
| CNA (COSMIC) | `src/pipeline/outputs/cosmic_cna.parquet` | 168,411 × 14 | `total_cn`, `minor_allele`, `cna_call`, **`matched_via`** |
| CNA flags | `src/pipeline/outputs/cna_flags.parquet` | 117,540 × 7 | `max_cn`, `min_cn`, `cna_type` |
| Combined flags | `src/pipeline/outputs/flags_with_driver.parquet` | 885,844 × 13 | `p_mutation`, `p_fusion`, `mut_driver`, `fusion_driver` |

**REVEL and AlphaMissense: `UNVERIFIED — could not locate`.** No column, file, or reference to either scorer exists in `src/` or `data/`. The pathogenicity signal in use is `max_pathogenicity` from the DepMap mutation table, not a named external predictor.

Cell-line coverage (from `coverage_matrix_enriched.parquet`, VERIFIED): mutations **1,744** models, fusions **1,698**, COSMIC CNA **997**, of 2,145 panel models.

### A1.5 Lineage / tissue annotation

Source: `gdsc_models.parquet:tissue`, derived from `data/GDSC/model_list_20260709.csv`. **30 distinct labels.** Full distribution in `MEASUREMENTS.md` §B3.5.

- Of the 2,145 panel models, **1,781 (83.0%)** carry a tissue label; **364 (17.0%) have none**.
- Of the 1,746 models actually scored in `core_score.parquet`, **1,634 (93.6%)** are labelled; 112 (6.4%) are not.
- A second lineage column (`sample_info.lineage`, DepMap-derived) exists in the warehouse and is surfaced by `metadata.py:30`. **Two lineage vocabularies coexist and no code reconciles them** — a traceability defect.

### A1.6 Labels / ground truth

| Label set | Path | Coverage | Consumed by |
|---|---|---|---|
| Chronos / Project Score | `validation/prepared/chronos_long.parquet` | 16,409,850 rows; 17,645 genes × ~930 lines | `build_chronos_validation.py:86`, `test_run_abundance_gate.py`, `test_run_gate_mechanism.py`, `test_run_improvement_scan.py`, `test_run_stratum_centring.py`, `test_run_metabolomics_mirna.py` |
| Real DepMap Chronos | `data/DepMap_Chronos/GeneFitnessEffect_Chronos_Achilles.hdf5` (110 MB), `…_Score.hdf5` (44 MB) | not profiled | **No consumer** — orphan input |
| DepMap reference essentials | `data/DepMap_Chronos/ReferenceEssentials.csv` (1,079) / `ReferenceNonEssentials.csv` (576) | 1,043 / 558 map to ENSG | `test_run_abundance_gate.py`, `test_run_cna_gate.py`, `test_run_detection_vs_percentile.py`, `test_run_geo_gate_replication.py`, `test_run_improvement_scan.py`, `test_run_screen_replication.py` |
| GDSC2 | `data/GDSC/GDSC2_fitted_dose_response_27Oct23.xlsx` → `validation/prepared/gdsc_long.parquet` (229,934) | | `06_held_out_eval.ipynb`, `stage4_eval_save.py`, `build_gene_regime.py` |
| Project Score matrices | `data/GDSC/Project_score_combined_*.tsv` (4 files, 508 MB) | | `test_run_project_score_replication.py`, `test_run_screen_replication.py` |
| 141-gene curated set | derived: `sorted(set(gene_regime.ensg_id) & set(gdsc.target_ensg))` | 141 genes → 113 train / 28 test, seed 42 | `06_held_out_eval.ipynb:cell 3`, `stage4_eval_save.py` |

**Important provenance note (VERIFIED).** `chronos_long.parquet` is **not** DepMap Chronos. Its `essentiality` values have median +3.05 and range −55.1…+29.1 — a Bayes-factor-like scale from Project Score, not Chronos gene effect (which centres on 0 with essentials near −1). The genuine Chronos HDF5 files sit unused in `data/DepMap_Chronos/`. Every file and variable naming this column "chronos" is mislabelled. This does **not** break the sign convention (see B1.1) but it does mean *"Chronos"* in the codebase means *"Project Score"* throughout.

### A1.7 Paralog annotation — **present, contradicting §C.17**

`data/ensembl_paralogy/hsapiens_paralogy_biomart_20260806.parquet` — 3,552,265 rows × 4 (`gene`, `paralog`, `perc_id_target`, `perc_id_query`), with `PROVENANCE.json`, dated 2026-08-06. Consumed by `build_paralog_flags.py`, which writes `max_paralog_identity`, `n_high_identity_paralogs`, `n_identical_paralogs`, `paralog_identical_flag`, `paralog_ambiguous_flag` into `src/pipeline/outputs/gene_ambiguity_flags.parquet` (20,163 × 12). Thresholds `IDENTICAL = 100.0`, `AMBIGUOUS = 95.0` (`build_paralog_flags.py:44-45`). VERIFIED. **§C.17 and Appendix gap 16 are stale — paralog data is no longer blocked on external sources.**

---

## A2. Code inventory and dependency graph

### A2.1 Execution-order graph (statically extracted from `read_parquet`/`to_parquet` calls)

```
Stage 0  00_harmonisation.ipynb → 00b_enriched_harmonisation.ipynb → 00c_harmonisation_eda.ipynb
             └→ harmonised.parquet, harmonised_enriched.parquet, gene.parquet, gene_enriched.parquet,
                coverage_matrix{,_enriched}.parquet, procan_proteomics.parquet, cosmic_cna.parquet,
                gdsc_models.parquet, celllineselector.db
Stage 1  01_lookup_metadata_join.ipynb   reads gene_lookup, union_lookup_entity → harmonised.parquet
Stage 2  02_core_score.ipynb             reads depmap_expr_clean, proteomics, harmonised_enriched
                                          → core_score.parquet, gene_dispersion.parquet
Stage 3  build_cna_layer.py              → cna_flags.parquet
         build_gene_roles.py             → gene_lookup.parquet (in-place)
         build_gene_regime.py            → gene_regime.parquet
         04_driver_gated_routing.ipynb   → flags_with_driver.parquet
Stage 4  05_confidence_tiers.ipynb  /  build_full_predictions.py
                                          → predictions_with_confidence.parquet
Stage 5  build_chronos_validation.py     → chronos_validation.parquet
         build_evidence_ledger.py        → evidence_ledger.parquet
Stage 6  06_held_out_eval.ipynb / stage4_eval_save.py
                                          → stage4_eval_consistent_denom.parquet, stage4_eval_summary.json
Query    cli.py → rank_cell_lines.py, explain_pair.py, multi_gene.py, export_web.py
```

### A2.2 Orphans and broken edges (VERIFIED)

| Issue | Detail |
|---|---|
| **Broken input** | `04_driver_gated_routing.ipynb` reads `layer3_continuous.parquet`. The only copy on disk is `discarded/outputs/layer3_continuous.parquet`. The notebook cannot run from `src/pipeline/outputs/`. |
| **Duplicated artefact** | `gene_regime.parquet` exists in both `src/pipeline/outputs/` (141 × 10) and `discarded/outputs/`. Consumers resolve to the former by relative path only. |
| **Orphan input** | `data/DepMap_Chronos/*.hdf5` (155 MB, the real Chronos) — no consumer. |
| **Orphan input** | `data/DepMap_Chronos/CCLEGeneCopyNumber20Q2.hdf5` (202 MB) — no consumer; the CNA layer uses COSMIC instead. |
| **Dead code** | `scoring_variants.combine_three()` — zero callers (see B1.2). |
| **Circular-looking self-reads** | Several `test_run_*` scripts read their own `*_per_gene.parquet` output; these are resume/cache paths, not cycles. |
| **All-null columns** | `harmonised_enriched.parquet` ships `biotype` and `hgnc_status` typed `null` — 100% missing, 2,127 rows. |
| **Three lineages of the same table** | `data/parquet/raw_data/`, `data/parquet/data_clean/` and `cleaned_track_data/` hold near-duplicate copies of the same 14 sources with differing row counts (e.g. mutations 1,066,869 / 1,066,869 / 632,919). |

### A2.3 End-to-end execution

**Not attempted, by decision.** Every `build_*.py` and Stage 0–6 notebook writes into `src/pipeline/outputs/`, overwriting the exact artefacts this audit measures. Running them would destroy the evidence base mid-audit and silently "fix" findings that Standing Rule 7 says must be written up rather than patched. What was verified instead:

- All read-only `test_run_*` result JSONs on disk are internally consistent with the parquet artefacts they describe (spot-checked on `stage4_eval_summary.json`, `test_run_cna_gate_results.json`, `test_run_abundance_gate_chronos_results.json`).
- The Stage 6 split was **re-executed independently** and reproduces the shipped 28-gene test set exactly (`split reproduces shipped stage4 test set: True`). VERIFIED.
- `04_driver_gated_routing.ipynb` **cannot** run end-to-end as committed — its `layer3_continuous.parquet` input is only in `discarded/`. VERIFIED.

To run the chain destructively, snapshot `src/pipeline/outputs/` first; `docs/REBUILD_RUNBOOK.md` documents the intended order.

---

## A3. Constants register

`Swept?` = a sensitivity analysis over this value exists in the repo. `Justified?` = a written rationale exists in code comments or `docs/`.

### A3.1 Constants named in the brief

| Constant | Value | File:line | Swept? | Justified? |
|---|---|---|---|---|
| `AMP_THRESHOLD` | 2.5 | `build_cna_layer.py:30` | No | Partly — `test_run_ploidy_normalised_cna.py:59` restates it; gap 8 notes it assumes diploidy |
| `DEL_THRESHOLD` | 1.5 | `build_cna_layer.py:31` | No | Same |
| `DEL_THRESHOLD` (gate) | **0.1** | `test_run_cna_gate.py:76` | **Yes** — sweep at 0.05/0.1/0.2/0.3 in `test_run_cna_gate_results.json` | Yes |
| `RHO_THRESHOLD` | **0.30** | `build_chronos_validation.py:72` | No | **Yes, extensively** — `build_chronos_validation.py:28-52` |
| `RHO_WEAK_THRESHOLD` | 0.10 | `build_chronos_validation.py:74` | No | Yes, same block |
| `RHO_THRESHOLD` (legacy 0.1) | 0.1 | `test_run_metabolomics_mirna.py:43`, `test_run_project_score_replication.py:59` | No | No — **the superseded value survives in two live scripts** |
| `P_THRESHOLD` | 0.05 | `build_chronos_validation.py:75`, +2 sites | No | No |
| `MIN_N` | 30 | `build_chronos_validation.py:76`, `test_run_project_score_replication.py:61` | No | No |
| `DRIVER_BOOST` | 0.15 | `src/Track - C/03_mutations_scoring.ipynb:62` (marked `# locked`); mirrored `test_run_signature_discount.py:66` | No | No — "locked" is asserted, not derived |
| `INFRAME_BOOST` | 0.15 | `src/Track - C/04_fusions_scoring.ipynb:68` (`# locked`) | No | No |
| `ALPHA` | 1.0 | `test_run_metabolomics_mirna.py:174` | No | No |
| Hill `(p₀,k)` — VEP | (3.0, 2.0) | `03_mutations_scoring.ipynb`; mirrored `test_run_signature_discount.py:62` | No | No |
| Hill `(p₀,k)` — pathogenicity | (0.5, 2.0) | same | No | No |
| Hill `(p₀,k)` — burden | (3.0, 1.5) | same | No | No |
| Hill `(p₀,k)` — fusion conf | (2.0, 2.0) | `04_fusions_scoring.ipynb:PARAMS_FUS` | No | Comment only |
| Hill `(p₀,k)` — fusion ffpm | (0.50, 1.5) | same | No | Comment only |
| Hill `(p₀,k)` — fusion recur | (0.50, 1.5) | same | No | Comment only |
| Hill `(p₀,k)` — scoring layer | (0.5, 2.0) | `scoring_variants.py:41` | No | **Self-declared placeholder**: "Uncalibrated — post-validation tuning would fit k, p0" |
| `rho_EP` | 0.46 | `scoring_variants.py:75` default; `test_run_abundance_gate.py:140`, `test_run_stouffer_regime_split.py:113`, `test_run_stratum_centring.py:86` | No | Partly — gap 10 flags a scale mismatch |
| `rho_EC`, `rho_PC` | 0.005 | `scoring_variants.py:75` | No | Comment: `rho_pc` "estimated equal to rho_ec" |
| `MIN_RATIO_EFFECT` | 0.10 | `test_run_detection_vs_percentile.py:160` | No | No |
| `MIN_PAUC_DELTA` | 0.01 | `test_run_detection_vs_percentile.py:161` | No | No |
| `MAX_OR_EFFECT` | 0.80 | `test_run_detection_vs_percentile.py:162` | No | No |

### A3.2 Constants **not** on the brief's list — the undocumented ones

None of the following is swept, and none carries a written justification unless noted.

| Constant | Value | File:line |
|---|---|---|
| `DEP_THRESHOLD` | **−0.5** | `test_run_abundance_gate.py:145`, `test_run_cna_gate.py:75`, `test_run_detection_vs_percentile.py:157`, `test_run_gate_mechanism.py:89`, `test_run_geo_gate_replication.py:116`, `test_run_improvement_scan.py:77`, `test_run_screen_replication.py:81`, `test_run_stratum_centring.py:88` — **8 sites, defines "dependent" for every gate result in the project** |
| `LIFT_DEFAULT` | 2.0 | `build_gene_regime.py:118` — the abundance_tracking/activation_driven class cut |
| `K_DEFAULT` | 20 | `build_gene_regime.py:117` |
| `MIN_LINES` | 100 / 150 | `test_run_abundance_gate.py:144`, `test_run_cna_gate.py:77`, `test_run_geo_gate_replication.py:117`, `test_run_screen_replication.py:82` / `test_run_improvement_scan.py:78` — **inconsistent across scripts** |
| `MIN_POS` / `MIN_NEG` | 10/10, 8/8 | 6 sites — **inconsistent** |
| `MIN_LINEAGE` | 15 | `test_run_improvement_scan.py:80` — the closest thing in-repo to `MIN_PEERS_FOR_LINEAGE` |
| `FLAT_FOLD_CEILING` | 3.0 | `evidence_state.py:74` |
| `MIN_LINES_FOR_RANKING` | 30 | `evidence_state.py:77` |
| `JOINT_FLOOR` | 0.50 | `multi_gene.py:64` — the weakest-link multi-gene cut |
| `DEFAULT_TOP_N` | 10 | `multi_gene.py:60` |
| `IDENTICAL` / `AMBIGUOUS` | 100.0 / 95.0 | `build_paralog_flags.py:44-45` |
| `VAF_LOH_THRESHOLD` | 0.75 | `test_run_tsg_two_hit.py:57` |
| `FLOOR_RATE` / `FLOOR_BASE` | 0.0139 / 0.0169 | `test_run_gate_mechanism.py:181-182` — **hardcoded results of a previous run, used as constants** |
| `HART_CUT` / `TPM_CUT` | −3.0 / 1.0 | `test_run_detection_vs_percentile.py:166-167` |
| `ABSENT_Q` | 10 | `test_run_geo_gate_replication.py:120` |
| `COVERAGE` | 0.90 | `test_run_improvement_scan.py:82` (conformal target) |
| `PARALOG_MAX_GROUP` | 30 / 10 | `test_run_detection_vs_percentile.py:202-203` |
| `MIN_CELL`, `MIN_HALF_POS/NEG`, `HALF_SAMPLE_REPS`, `KDE_BINS` | 5, 5, 5, 20, 1024 | `test_run_detection_vs_percentile.py:156,180-182,168` |
| `N_GENES` subsample | 2000 / 3000 / 6000 | `test_run_improvement_scan.py:81`, `test_run_gate_mechanism.py:90`, `test_run_screen_replication.py:85` |
| `N_BOOT` / `N_PERM` | 2000 / 10,000 | `test_run_tsg_two_hit.py:60`, `test_run_signature_discount.py:67`, `test_run_stouffer_regime_split.py:119` |
| `SEED` | 42 | 11 sites (uniform — good) |
| `BATCH` | 200 / 2000 | `build_evidence_ledger.py:43`, `build_full_predictions.py:83` (performance only) |

**The two that matter most.** `DEP_THRESHOLD = −0.5` silently defines the positive class for every gate/ranker result the project quotes, at 8 independent sites, unswept. `LIFT_DEFAULT = 2.0` decides the abundance/activation class of all 141 curated genes and therefore which routing rule each gene receives, unswept.

---

# PART B — Documents vs code

## B1.1 The Chronos sign convention — **the brief's hypothesis is not supported**

**Sites located** (all three, expression quoted verbatim):

| Site | Expression |
|---|---|
| `build_chronos_validation.py:92` (live) | `ch["chronos_pct"] = 1.0 - ch.groupby("ensg_id")["essentiality"].rank(pct=True, method="average")` |
| `discarded/code/validate_chronos.py:25` | `ch["chronos_pct"] = 1.0 - ch.groupby("ensg_id")["essentiality"].rank(pct=True, method="average")` |
| `discarded/code/build_chronos_layer.py:43-48` | `ch["chronos_pct"] = …rank(pct=True…)` then `ch["chronos_pct"] = 1.0 - ch["chronos_pct"]` |

**The sign is applied identically at all three sites.** VERIFIED. There is no layer-vs-validation discrepancy.

**Polarity of the source column, against DepMap's own reference sets** (VERIFIED):

| Gene set | median `essentiality` |
|---|---|
| DepMap `ReferenceEssentials` (n=1,043 mapped) | **−5.198** |
| DepMap `ReferenceNonEssentials` (n=558) | **+3.492** |
| `RPL*`/`RPS*`/`POLR2A`/`PSMA*`/`EIF*` (n=147) | **−10.682** |
| all genes (n=17,645) | +3.202 |

Mann–Whitney essentials vs non-essentials: **p = 3.29e-229**. Negative = essential. The `1 − rank` inversion therefore maps essential → high `chronos_pct`, which is **correct**.

**The diagnostic the brief specifies cannot work, and this is itself a finding.** `chronos_pct` is ranked *within gene across cell lines*. Mean `chronos_pct` is **0.4995 for reference essentials and 0.4995 for reference non-essentials** — identical, by construction. The transform carries **zero** between-gene information. Pan-essentiality is a between-gene property, so "pan-essentials should come out confirmed" is not a testable prediction of this metric. Any pipeline claim that reads `chronos_check` as evidence about *whether a gene is essential* is reading a quantity that cannot carry that information. VERIFIED.

**Control-set result, reported as requested** (control = 1,149 genes from DepMap essentials ∪ the symbol families, ≥30 shared lines):

| Threshold | confirmed | contradicted | % backwards |
|---|---|---|---|
| as-shipped `RHO_THRESHOLD = 0.30` | 1 | 3 | 75.0% |
| legacy `RHO_THRESHOLD = 0.10` | 34 | 391 | 92.0% |

median ρ on the control set = **−0.0655**, 79.3% negative.

**The brief's 769 / 1,258 figures are stale.** They describe the superseded `RHO_THRESHOLD = 0.10` regime, which `build_chronos_validation.py:28-52` documents as deliberately replaced. The **shipped** `chronos_validation.parquet` (16,866 genes) now reads: `none` 14,911, `weak_negative` 1,188, `weak_positive` 730, `validated` **31**, `inverted` **6**. VERIFIED.

**The directional asymmetry is real and survives the threshold change.** Of 4,332 genes significant at p<0.05: 2,468 negative vs 1,864 positive = **57.0% negative, binomial p = 4.46e-20**. In the weak band alone it is 1,188 vs 730 = 61.9%, matching the brief's 62.1%.

**The line-level confound is ruled out** (VERIFIED): Spearman(per-line screen strength over 1,043 reference essentials, per-line mean `core_score`) = **−0.019, p = 0.58, n=891**. The negative skew is not a proliferation/screen-quality artefact.

**Conclusion.** Sign correct, consistently applied. The 62% backwards figure is real but is **not** evidence of an inversion; it is evidence that within-gene abundance→dependency is weakly *negative* on average — which contradicts §C.15's implied positive relationship. That contradiction stands and is unexplained. It belongs in the risk register, not in a sign fix.

## B1.2 Does `w = 1/(1+ρ)` ever produce weights ≠ (½,½)?

**It never executes at all.** VERIFIED.

- `combine_three()` (`scoring_variants.py:75-136`) has **zero callers** across `src/`. The only imports of `scoring_variants` are `from scoring_variants import combine` in `test_run_metabolomics_lineage_confound.py:23` and `test_run_metabolomics_mirna.py:35` — the two-layer `combine()`, not `combine_three()`.
- `core_score.parquet` is written by `02_core_score.ipynb` cell 15, which **hardcodes** the two-layer rule:
  `core_score[both] = 0.5 * E[both] + 0.5 * P[both]`
  The cell's own comment states the correlation-penalised weights "degenerate to w_e = w_p = 0.5" and that it replaced a stale stored variant built before ProCan was merged.

**Distinct weight tuples actually written to `core_score.parquet`** (29,781,274 rows, VERIFIED):

| Stratum | Weights | Rows |
|---|---|---|
| `n_layers = 1` | (1.0) — pass-through of the single present layer | **24,504,210** |
| `n_layers = 2` | **(0.5, 0.5)** — hardcoded | **5,277,064** |
| `n_layers = 3` | — | **0** (no three-layer stratum exists) |

§2.4's argument is correct, but understated: the formula does not merely degenerate, it is never invoked. Had it been invoked, the E–C and P–C pairs would **not** have been (½,½) — with `rho_ep=0.46`, `rho_ec=rho_pc=0.005`, the pair weights work out to w_ec ≈ (0.448, 0.552). That asymmetry is latent in the dead code and would appear the moment a third layer is wired in.

## B1.3 Metric definitions — two quantities share one name

**Every implementation in the repo computes recall@k.** VERIFIED.

| Implementation | Expression | Actually computes |
|---|---|---|
| `discarded/code/validate_chronos.py:87-96` `hit_at_k()` | `len(top_k & sensitive) / len(sensitive)` | **recall@k** |
| `test_run_diversity_harness.py:56-60` `hit_at_k()` | `len(top_k & sensitive_model_ids) / len(sensitive_model_ids)` | **recall@k** |
| `test_run_metabolomics_mirna.py:108-114` `hit_at_k()` | same | **recall@k** |
| `build_gene_regime.py:167` `hit_at_20_train` | `hits_at_k / n_known` | **recall@k** |
| `06_held_out_eval.ipynb` cell 6 / `stage4_eval_save.py` | `hits_flat / known` | **recall@k** |

All five call it "hit-rate@20" / "hit@20" / `hit_at_20_train`. `build_gene_regime.py:84` already flags the naming: *"The column name `hit_at_20_train` suggests the original author knew…"*.

**The .docx's top-1/top-5/top-20 figures have no implementation in this repo.** `UNVERIFIED — could not locate`. Grep for `23.4`/`65.2`/`93.6` across `docs/*.md` and `src/pipeline/*.py` returns one unrelated hit (`MATH_REFERENCE.md:1651`, a stratum ratio). No code computes an any-hit@k metric anywhere.

**Reconstructing it confirms the brief's diagnosis** (VERIFIED, over all 141 curated genes):

| k | any-hit@k (recomputed) | recall@k (recomputed) | .docx claim |
|---|---|---|---|
| 1 | 17.7% | 0.1% | 23.4% |
| 5 | 53.9% | 0.4% | 65.2% |
| 20 | **89.4%** | 1.6% | 93.6% |

The .docx figures track the **any-hit** shape, not the recall shape. Two different quantities share the name "hit-rate@20": the formula the .docx *states* is recall, the numbers it *reports* are any-hit.

**The decisive point is the baseline.** With a median of 228 sensitive lines among 1,563 scored lines per gene, the hypergeometric chance baseline for any-hit@20 is **90.8%**. The observed 89.4% is *below* chance. The headline number is not evidence of anything.

## B1.4 In-sample vs held-out

**The 113 tuning genes are included** in the .docx's "141 curated genes" figure — the set is literally `train ∪ test`. VERIFIED.

Recomputed with hypergeometric baselines and a 10,000-resample paired bootstrap over genes (`*` = 95% CI on excess-over-chance excludes zero):

| Split | n | Metric | Observed | Chance | Excess | 95% CI |
|---|---|---|---|---|---|---|
| train | 113 | recall@1 | 0.0008 | 0.0006 | +0.0002 | [−0.0001, +0.0005] |
| train | 113 | any-hit@1 | 0.1858 | 0.1458 | +0.0400 | [−0.0273, +0.1125] |
| train | 113 | recall@5 | 0.0039 | 0.0032 | +0.0007 | [−0.0002, +0.0016] |
| train | 113 | any-hit@5 | 0.5310 | 0.5226 | +0.0083 | [−0.0787, +0.0946] |
| train | 113 | **recall@20** | 0.0168 | 0.0127 | **+0.0041** | **[+0.0017, +0.0067]** `*` |
| train | 113 | any-hit@20 | 0.9027 | 0.9083 | −0.0056 | [−0.0580, +0.0421] |
| **test** | **28** | recall@1 | 0.0013 | 0.0006 | +0.0007 | [−0.0005, +0.0026] |
| **test** | **28** | any-hit@1 | 0.1429 | 0.1548 | −0.0119 | [−0.1310, +0.1279] |
| **test** | **28** | recall@5 | 0.0042 | 0.0031 | +0.0010 | [−0.0011, +0.0037] |
| **test** | **28** | any-hit@5 | 0.5714 | 0.5405 | +0.0309 | [−0.1501, +0.2166] |
| **test** | **28** | **recall@20** | 0.0140 | 0.0126 | +0.0014 | [−0.0027, +0.0062] |
| **test** | **28** | any-hit@20 | 0.8571 | 0.9086 | −0.0514 | [−0.1840, +0.0640] |

**Exactly one cell in the table is above chance, and it is in-sample.** On the 28 held-out genes, every metric at every k is indistinguishable from a random ranker. This is the strongest available statement of §9 gap 15, and it is now measured rather than asserted.

The shipped `stage4_eval_summary.json` agrees: held-out `hr_flat = 0.014025`, `hr_driver = 0.017735`. The "26.45% driver lift" it reports is a lift of 0.0037 in absolute recall, with no interval attached.

## B2.1 `n_layers == 2 → "high"` — **still live**

VERIFIED at both sites:
- `build_full_predictions.py:151`: `(c.regime_source == "measured") & (c["class"] == "abundance_tracking") & (c.n_layers == 2)` → `"high"`
- `05_confidence_tiers.ipynb` cell 5: identical predicate → `"high"`

Unchanged despite §C.18 measuring the protein layer as a net negative **on global AUROC** (−0.0086 GDSC / −0.0042 Chronos). On **pAUC** — the metric that matches a top-N recommender, since that is the region of the ranking a user actually sees — the production two-layer score beats RNA alone (+0.0234 GDSC / +0.0171 Chronos; see `PROTEOMICS_DEFENCE.md` §6.2). So the tier rule is directionally justified by the metric that matches the product, not simply contradicted by measurement — but it has never been calibrated: nobody has checked whether "high"-tier rows are actually more often correct than "moderate"-tier ones. Gap 14 **open**, recharacterised.

## B2.2 Noisy-OR in Track C — **live, and user-visible**

VERIFIED.
- Produced by `03_mutations_scoring.ipynb` (`p_base = 1 − (1−p_vep)(1−p_path)(1−p_burden)`, then `+ DRIVER_BOOST·driver_flag`, clipped to 1.0) and `04_fusions_scoring.ipynb` (same form over conf/ffpm/recur, `+ INFRAME_BOOST·inframe_flag`).
- **Score path:** reaches exactly one consumer, as a boolean — `has_alteration = (p_mutation > 0.5) | (p_fusion > 0.5)` (`test_run_signature_discount.py:21` documents this; `flags_with_driver.parquet` carries it).
- **User path:** the raw float **does** reach user-visible output. `explain_pair.py:178-206` emits `"p_mutation": float(...)` and `"p_fusion": float(...)` into the pair explanation, surfaced by `cli.py:89`. They are tagged `"CONFIDENCE MODIFIER — does not affect stratum_rank"`.

So: display-only with respect to *ordering*, but a user-facing probability nonetheless — one produced by an independence assumption that has never been tested. Gap 2 **open**.

## B2.3 Multiple-testing correction — **none exists**

VERIFIED. `grep -rniE "\b(bh|fdr|benjamini|hochberg|false_discovery|multipletests|bonferroni)\b" src/` returns **zero** genuine hits. The only matches are `bh` as a bootstrap-upper-bound variable (`src/figures/compute_stats.py:495-498`) and `bh` as a bar height (`src/figures/fig_example_output.py`). Gap 3 **open**.

This matters concretely: `build_chronos_validation.py` runs **16,866** simultaneous per-gene hypothesis tests at an uncorrected p<0.05 and assigns user-facing labels from the result. At α=0.05 uncorrected, ~843 genes are expected to pass by chance alone; 4,332 pass. The 31 `validated` and 6 `inverted` labels are the only ones that survive an effect-size floor, but `weak_positive`/`weak_negative` (1,918 genes) are assigned on significance alone.

## B2.4 Composite identifiers — **still present**

VERIFIED, from `test_run_identifier_ambiguity_audit_results.json`:

| Defect | Count |
|---|---|
| ProCan accessions mapping to semicolon-joined composite symbols (`ClassA`) | **15** (e.g. `a8mya2 → CXorf49;CXorf49B`, `p0c5z0 → H2AB2;H2AB3`, `p12532 → CKMT1A;CKMT1B`, `p23610 → F8A1;F8A2;F8A3`, `p43362 → MAGEA9;MAGEA9B`) |
| CCLE/TMT composite | **1** (`q8wwn9-2 → CNK3/IPCEF1`) |
| ProCan proteins lost at the join | **1,633** (1,618 symbol-present-but-unmatched + 15 composite) |
| CCLE proteins lost at the join | **420** (402 unmatched + 17 no-symbol + 1 composite) |
| Ambiguous accessions, CCLE / ProCan | 23 / 16 (5 / 2 linked to multiple genes) |
| HGNC accessions mapping to >1 gene | 60 of 19,187 |
| Composite `model_id` dry-run rows | 400 (`test_run_composite_model_id_dryrun.parquet`) |

Semicolon-joined gene symbols are **not** fixed. The audit that detects them is read-only and no build-time assertion blocks them.

## B2.5 The protein-coding filter — **two different filters coexist**

VERIFIED. Neither is a versioned GENCODE biotype filter.

| Filter | Definition | Genes admitted | Used by |
|---|---|---|---|
| A — gene_lookup biotype | `biotype == "protein_coding" & hgnc_status == "Approved"` on `reference/gene_lookup.parquet` | **19,213** | all 8 `test_run_*` gate/scan scripts |
| B — HGNC protein-product file | membership in `data/gene_with_protein_product.txt` | **19,187** accessions / 19,273 rows in the `hgnc` warehouse table | `build_ambiguity_flags.py:36`, `test_run_identifier_ambiguity_audit.py:52,467` (via `locus_group`) |
| C — `is_protein_coding` flag | boolean on `gene_enriched.parquet` (20,163 rows) | not equal to either A or B; the table has 20,163 rows against A's 19,213 | `build_warehouse_views.py:184,212`, `metadata.py:59`, `track_c_eda.py:134` |

`gene_lookup.parquet` carries a `biotype` column but **no source-version column**; the HGNC download `data/gene_with_protein_product.txt` is dated 2026-07-10 by mtime only. Neither is pinned to a GENCODE release. Any statement of the form "the universe is N protein-coding genes" is therefore ambiguous between 19,213, 19,187 and 20,163 depending on which script is speaking. The scored universe in `core_score.parquet` is a fourth number: **19,177**.

---

## Appendix — status of the 19 known gaps

| # | Gap | Brief's priority | Status |
|---|---|---|---|
| 13 | Floor effect not separable from profiling confound (MNAR) | highest | **OPEN**, and reinforced. `test_run_gate_mechanism_results.json` verdict: `no_floor_specific_effect__profiling_confound_only`; floor Fisher p = 0.246, mid-stratum ratio 1.618, p = 0. Proteomics coverage bias measured here at Cliff's δ = +0.251 (B3.3). |
| 14 | `n_layers == 2 → high` contradicted by §C.18 | high | **OPEN** — live at `build_full_predictions.py:151` and `05_confidence_tiers.ipynb` cell 5 (B2.1) |
| 15 | Sort key unsupported by any measurement | high | **OPEN**, now quantified: ranker pAUC 0.5002/0.5015; held-out recall@20 excess CI [−0.0027, +0.0062] (B1.4) |
| 1 | No null model for hit@k | medium | **OPEN** — no baseline is computed anywhere in `src/`. Supplied here (B1.3/B1.4) |
| 2 | Noisy-OR still live in Track C | medium | **OPEN** — and reaches user output via `explain_pair.py:178-206` (B2.2) |
| 3 | No multiple-testing correction | medium | **OPEN** — zero hits repo-wide (B2.3) |
| 5 | Per-layer variances never estimated | medium | **OPEN** — no variance estimate exists; weights are hardcoded 0.5/0.5 (B1.2) |
| 6 | Confidence tiers uncalibrated | medium | **OPEN** — tiers are rule-based by construction (`05_confidence_tiers.ipynb` cell 0: "not a fabricated 0–1 score"); no calibration curve exists |
| 18 | Label reliability ceiling not reported alongside metrics | medium | **PARTLY CLOSED** — `label_fpr = 0.00928` is now computed and reported in 5 `test_run_*` result JSONs, but never alongside the headline @k figures |
| 4 | Max-aggregation extreme-value bias | low | **OPEN** — `.groupby(["ensg_id","model_id"]).max()` live in `03_mutations_scoring.ipynb` and `04_fusions_scoring.ipynb` |
| 7 / 17 | Clustering ignored in bootstrap and CV splits | low | **OPEN** — Stage 6 split is `rng.choice` over genes, unstratified (`06_held_out_eval.ipynb` cell 3) |
| 8 | CNA thresholds assume diploidy | low | **OPEN but now unblocked** — ploidy is on disk (A1.3); 55.8% of lines have ploidy > 2.5 (B3.4) |
| 9 | `product_floor` unanalysed | low | **MOOT** — `product_floor` is not used in the shipped path; `core_score` is a hardcoded mean (B1.2) |
| 10 | ρ_EP scale mismatch | low | **OPEN** — 0.46 hardcoded at 4 sites, unswept |
| 11 | Cross-stratum calibration deferred | low | **OPEN** — `stratum_rank` is computed per `(ensg_id, n_layers)`; no cross-stratum comparison exists |
| 12 | Single-anchor regression tests | low | **OPEN** — `rank_cell_lines.py:122-131` asserts on A375 only |
| 16 | Paralog buffering unaddressed | low | **DATA NO LONGER BLOCKING** — paralogy table + `build_paralog_flags.py` exist (A1.7). The *buffering analysis* itself is still absent |
| 19 | CNA residual exceeds label noise floor | low | **OPEN** — `test_run_cna_gate_results.json`: residual 0.0304 vs label FPR 0.00928, verdict `residual_above_noise__check_calls` |
