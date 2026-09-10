# COMBINED_SCORE_PREREGISTRATION — Task B1

**Written before any Stage 1 measurement.** Nothing in this file was informed by
drug-response data. Commit this file before running
`src/pipeline/combined_score_b1_eligibility.py` or anything downstream.

Governing rule: every claim traces to a file path or a computed number. Where
something cannot be established, it is marked **UNVERIFIED — could not locate**.

---

## 1. The supplied A–B list, verbatim

Supplied by the user. Not derived from data, not searched for by outcome.

### Established / well-documented relationships (n = 10)

| # | A | B | stated relationship |
|---|---|---|---|
| 1 | EGFR | KRAS | KRAS mutation confers resistance to anti-EGFR therapy |
| 2 | BRAF | MAP2K1 (MEK1) | same MAPK pathway, sequential signalling |
| 3 | BRCA1 | BRCA2 | DNA repair, homologous recombination |
| 4 | TP53 | MDM2 | MDM2 negatively regulates p53 |
| 5 | PTEN | PIK3CA | both converge on PI3K/AKT; often reciprocal mutation |
| 6 | ERBB2 (HER2) | ERBB3 | heterodimerisation partners |
| 7 | ALK | EML4 | EML4–ALK fusion, oncogenic driver |
| 8 | MYC | MAX | obligate heterodimer for transcriptional activity |
| 9 | RB1 | CDK6 | RB is a substrate of the CDK4/6–cyclin D complex |
| 10 | VHL | HIF1A | VHL targets HIF1A for degradation |

### Contrast set — weak or no established direct relationship (n = 10)

| # | A | B |
|---|---|---|
| 11 | GAPDH | KRAS |
| 12 | ACTB | EGFR |
| 13 | BRCA1 | ALK |
| 14 | TP53 | EML4 |
| 15 | PTEN | MYC |
| 16 | VHL | BRAF |
| 17 | RB1 | VHL |
| 18 | ERBB2 (HER2) | TP53 |
| 19 | MAX | PTEN |
| 20 | CDK6 | VHL |

### Two resolution decisions, made now rather than later

- **#9 "CDK4/CDK6" → CDK6.** The supplied text names both. CDK6 is used, for
  consistency with #20 which names CDK6 alone. Testing both would add a second
  test per pair without a pre-specified rule for combining them.
- **#2 "MEK1/MAP2K1" → MAP2K1**; **"HER2" → ERBB2**. Symbol normalisation only.

### Citations — UNVERIFIED

The list was supplied with stated relationships but **without citations**. B1
requires a citation per pair establishing the relationship independently of this
dataset.

**Status: UNVERIFIED — could not locate.** No citation has been checked against
the literature by this analysis. The relationships in column 3 are as stated by
the user and are biologically conventional, but they are recorded here as
assertions, not as sourced claims. Before any Stage 1 result is written up, each
row needs a reference verified against the source.

This does not block the eligibility check (§3), which touches no outcome data. It
does qualify any Stage 1 result.

---

## 2. The contrast set is the strongest feature of this design

The ten control pairs are a **negative control built into the hypothesis**, and
they are better evidence than a label-randomised placebo alone:

- a label permutation tests whether the *statistic* is calibrated
- the contrast set tests whether the *biology* is doing the work

Both are reported. If the established set and the contrast set show the same
enrichment, the effect is not about the modifier relationship regardless of how
significant it looks.

**Pre-specified:** the contrast set is a control, not a second hypothesis. It is
never promoted to a positive finding if it happens to score well; that would be
selection by outcome.

---

## 3. Inclusion criteria — fixed before measurement

A pair is **eligible** only if all of the following hold. Every criterion uses
coverage or alteration data only; **none uses drug response**.

1. **A is a GDSC2 drug target** with at least one (A, drug) pair in the 294-pair
   set from [COMBINED_SCORE_PRECONDITION.md](COMBINED_SCORE_PRECONDITION.md) §1.
2. **A passes the expression validity guard** — `frac_expressed ≥ 0.20`.
3. **B is in the 4,173 varying genes** — ≥10 lines carrying a driver alteration
   (§2c of the precondition). A modifier that never varies cannot modify
   anything.
4. **Within the pair's own line set**, B has ≥10 altered *and* ≥10 unaltered
   lines, so the `NOT altered(B)` term actually partitions.
5. **≥100 lines** with expression for A, alteration status for B, and drug
   response for the (A, drug) pair.

Direction is fixed: **A is the drug target, B is the modifier.** Where both
members of a supplied pair are GDSC2 targets, only the orientation named in the
list is tested. No orientation is chosen after seeing results.

---

## 4. Kill switch

**Fewer than 10 eligible established pairs → stop.** Report the eligibility
table as the only Stage 1 deliverable and ship Task A.

Rationale for 10 rather than the B0 threshold of 20: B0's 20 counted (gene, drug)
pairs; this counts A–B *relationships*, of which only 10 established ones exist in
the supplied list. Requiring 20 would be unsatisfiable by construction.

---

## 5. Adoption rule — written before measuring

The Stage 1 statistic is:

```
delta = enrichment[ expressed(A) AND NOT altered(B) ]
      - enrichment[ expressed(A) ]
```

where enrichment is the drug-sensitivity rate in the selected set divided by the
pair's base rate, at prevalence-matched binarisation. **The comparator is
single-gene A**, not random and not tissue.

Interval: cluster bootstrap over **cell lines**, never over pairs.

| outcome | decision |
|---|---|
| 95% CI on `delta` over the **established** set **excludes 0**, **and** the contrast set's CI **includes 0**, **and** the label-randomised placebo CI **includes 0** | **ADOPT** — proceed to Stage 2 (full battery) |
| established CI includes 0 | **ABANDON** — report the negative, ship Task A |
| established CI excludes 0 **but** the contrast set also excludes 0 | **ABANDON** — the effect is not about the modifier relationship |
| established CI excludes 0 but is driven by <3 pairs (per-pair reporting) | **ABANDON** — a single dominant pair is not an interaction result |

**No combined score reaches the user interface on Stage 1 alone.** Stage 2 —
combination-rule sweep with held-out selection, cluster bootstrap, placebo in the
same table, readout reliability, second readout if available — is required before
anything ships. Task A's `combined_score` stays `None` throughout.

---

## 6. What is measured, in order

1. `combined_score_b1_eligibility.py` — eligibility only, no outcome data.
2. If the kill switch does not fire: `combined_score_b2_stage1.py` — the
   interaction test under §5.

Results: `COMBINED_SCORE_RESULT.md`.
