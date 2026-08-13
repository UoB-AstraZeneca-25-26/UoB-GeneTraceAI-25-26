# DECISIONS_REQUIRED

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

Questions only you can answer. Each states the options, what each costs, and what evidence bears on it. Nothing in `BUILD_SPEC.md` marked *pending decision* should be implemented until the relevant question here is closed.

---

## Q1. Should alterations demote, promote, or gate? — **highest priority**

**Why it's yours.** Both directions are defensible for different questions, and the answer depends on what a user is asking the tool, which is a product decision, not a measurement.

**Current state.** The pipeline **promotes** at **11 sites** with no demotion anywhere: `DRIVER_BOOST = 0.15`, `INFRAME_BOOST = 0.15`, and `has_driver_alteration · 1e6 + core_score` — an absolute precedence that no score can overcome (`stage4_eval_save.py:45-47`, `06_held_out_eval.ipynb` cell 5, `explain_pair.py:109-112`). Reversing this is a full reversal of a documented, load-bearing decision.

**Evidence that bears on it:**
- A global rule touches **1.597%** of scored pairs (475,695 of 29,781,274). Median gene: 1.16% of its lines. But 87 genes exceed 10% and one reaches 62%. VERIFIED.
- Existing role classification is **not sufficient to condition on**: `gene_role` has evidence for 768 COSMIC CGC genes and calls the other ~18,400 "neither" by `fillna`, not by evidence; `class` covers only the 141 curated genes and is cut at an unswept `LIFT_DEFAULT = 2.0`.
- Both aggregations use `.groupby(...).max()`, so any *continuous* score in either direction inherits the `E[max] = m/(m+1)` bias (8.7% exposure).
- The gate region is the only place this architecture has replicated above-chance evidence (pAUC 0.596–0.608 vs ranker 0.500).

**Options:**

| | What it means | Cost |
|---|---|---|
| **A. Demote globally** | Altered pairs rank lower everywhere | Reverses 11 sites. Answers "clean model of normal function". Near-no-op for the median gene, large re-ordering for ~87 genes. Loses oncogene-addiction use case entirely |
| **B. Keep promoting** | Status quo | No work. But the promotion has never been validated on held-out data — `hr_driver` 0.0177 vs `hr_flat` 0.0140, and that difference's CI includes zero |
| **C. Role-conditional** | Demote for TSG/normal-function queries, promote for oncogene-addiction | Needs a role classification that covers more than 768 genes with evidence. Currently blocked on that |
| **D. Query-time switch (recommended)** | User declares intent; the tool applies the matching rule and says which it applied | Most honest, most work. Requires a UI/API decision and an explicit default |
| **E. Gate, not tie-break (recommended form, orthogonal to A–D)** | Boolean exclusion on a pre-specified variant class rather than a continuous score shift | Avoids the max-aggregation bias entirely. Matches the only evidence class this architecture has. Precedent: CNA deletion already works this way (ratio 0.380, 69 genes with ≥20 deleted lines) |

**Recommendation:** **D + E** — a query-time use-case switch, implemented as a gate rather than a continuous demotion. If a single default is required, **E with demotion off** preserves current behaviour while removing the max-aggregation bias.

**If you pick A or C, say so and I will specify the sweep;** the parameters are listed in `BUILD_SPEC.md` C4 and must be swept, not picked.

---

## Q2. Do we quote any @k number publicly, and if so which?

**Why it's yours.** This determines what the dissertation can claim.

**Evidence:**
- The `.docx` headline **top-20 93.6%** uses an *any-hit@k* metric whose **hypergeometric chance baseline is 90.8%**. Recomputed on the current build: 89.4% observed — *below* chance. VERIFIED.
- On the 28 held-out genes, **every** metric at **every** k has a 95% CI on excess-over-chance that includes zero. VERIFIED.
- The only above-chance cell in the entire table is `recall@20` on the **113 tuning genes**: +0.0041, CI [+0.0017, +0.0067].
- The `.docx`'s stated formula is recall; its reported numbers are any-hit. Two quantities, one name.

**Options:**

| | Cost |
|---|---|
| **A. Quote nothing @k; lead with the gate** | Honest. Gate pAUC 0.596–0.608, p = 2.8e-76, replicated on two screens — this is real and defensible. Loses a headline number that looks impressive |
| **B. Quote @k with baseline and CI always attached** | Also honest, and shows the work. The numbers are unflattering: "recall@20 = 1.4% against a 1.3% baseline, CI on the difference [−0.003, +0.006]" |
| **C. Keep quoting 93.6%** | Not defensible. The figure is below its own chance baseline |

**Recommendation: A, with B in an appendix.** Standing rule 3 makes C unavailable regardless. This is not a close call on the evidence; the decision you actually own is how to frame the change in the write-up.

---

## Q3. Is the ProCan normalisation resolvable? — blocks C6

**Why it's yours.** It needs an external lookup or a contact, not a computation.

**State.** `Protein_matrix_averaged_20250211.tsv` carries **log2 intensities** (1.4% negative, median 3.44) with per-line medians spread over **0.98 log2 units ≈ 2× linear** — so between-sample differences survive and the proteomic ruler is not dead. But **what normalisation ProCan applied upstream is recorded nowhere in this repository**. A 2× residual spread is consistent with no median scaling, but also with a partial or per-batch normalisation.

**Options:** (A) find the ProCan publication/portal documentation; (B) obtain less-processed intensities; (C) proceed with a stated assumption and a prominent caveat; (D) drop C6.

**Recommendation: A first, then B.** C is acceptable only if the output is framed as an order-of-magnitude band. Note independently that the histone/total ratio has a **41% CV**, which caps C6's usefulness at a gate regardless of how Q3 resolves — it is not precise enough to be a tie-break.

---

## Q4. Can the teammate's `02_transcriptonomics.ipynb` be obtained? — blocks C5

**State.** The notebook is **not in this repository**. Searches for the file and for `silent_lineage`, `MIN_PEERS_FOR_LINEAGE`, `mad_floor`, `MAD_FLOOR` return nothing. The decisive C5 experiment — does lineage conditioning *with guards* recover an effect that plain lineage percentiles do not — cannot be run.

**What is reproducible:** the unguarded result, which is neutral, confirming §C.15: pAUC panel 0.51756 vs lineage 0.51673, Δ = **−0.0018**, p = 0.540, better in 47.8% of genes.

**Options:** (A) get the notebook; (B) get a written spec of the three guards precise enough to reimplement; (C) implement the silent-lineage guard as a **pooled-percentile validity gate**, which needs neither and is testable now; (D) drop lineage to annotation-only.

**Recommendation: C now, A/B in parallel.** C is the salvageable part — the guard's question ("is this distribution degenerate enough that a percentile is meaningless?") does not require the stratum to be a lineage, and `gene_dispersion.parquet` already holds the raw material. Lineage should be adopted as a stratification and annotation variable regardless, because gaps 7 and 17 need it.

---

## Q5. Which protein-coding universe is *the* universe?

**State.** Four numbers are in simultaneous use, and no code reconciles them:

| Filter | Count | Used by |
|---|---|---|
| `gene_lookup.biotype == 'protein_coding' & hgnc_status == 'Approved'` | **19,213** | all 8 `test_run_*` gate/scan scripts |
| HGNC `gene_with_protein_product.txt` | **19,187** accessions | `build_ambiguity_flags.py`, identifier audit |
| `gene_enriched.is_protein_coding` | table has **20,163** rows | warehouse views, `metadata.py`, `track_c_eda.py` |
| actually scored in `core_score.parquet` | **19,177** | the shipped product |

None is a versioned GENCODE biotype filter, and **neither source carries a version column** — the HGNC file is dated only by filesystem mtime (2026-07-10).

**Cost of leaving it:** any statement of the form "the universe is N protein-coding genes" is ambiguous by up to 986 genes depending on which script is speaking. Assertion A11 in `BUILD_SPEC.md` C3 will fail on day one because of it — deliberately.

**Recommendation:** pin to a **named GENCODE release**, record it in a `dataset_version` column, and make every script read the same one. Cheap, and it closes a whole class of reconciliation questions.

---

## Q6. Is `chronos_long.parquet` the label set we want?

**State.** `chronos_long.parquet` is **not DepMap Chronos.** Its values have median +3.05 and range −55.1…+29.1 — a Bayes-factor-like Project Score scale, not Chronos gene effect (which centres on 0 with essentials near −1). The genuine DepMap Chronos matrices sit **unused** at `data/DepMap_Chronos/GeneFitnessEffect_Chronos_Achilles.hdf5` (110 MB) and `…_Score.hdf5` (44 MB).

This does **not** break the sign convention — verified against DepMap's own reference sets, negative = essential, and the `1 − rank` inversion is correct. But every file, variable and doc naming this "Chronos" is mislabelled, and a second, genuinely independent label set is available and unused.

**Options:** (A) rename throughout to `project_score_*` and keep using it; (B) switch to the real Chronos HDF5s; (C) use both as independent replication — the strongest option, and `test_run_screen_replication.py` already has the shape for it.

**Recommendation: A immediately** (it is a naming fix and costs nothing), **C as the substantive improvement.** Two independent screens agreeing is much better evidence than one screen, and it directly addresses gap 18.
