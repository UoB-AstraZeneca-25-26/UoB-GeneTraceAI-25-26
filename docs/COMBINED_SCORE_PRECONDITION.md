# COMBINED_SCORE_PRECONDITION — Task B0

**Verdict: kill switch does NOT fire. Task B is BLOCKED at B1 pending a
pre-registered modifier list, which is a user input I must not derive.**

Script `src/pipeline/combined_score_b0.py` → `src/pipeline/outputs/combined_score_b0_results.json`.
Task A ships regardless — see [MULTIGENE_SPEC.md](MULTIGENE_SPEC.md).

---

## 1. Gene–drug pairs

| filter | pairs | genes | drugs |
|---|---|---|---|
| raw (target gene, drug) in GDSC2 | **322** | 141 | 174 |
| + RNA coverage for A, ≥100 lines, ≥10 sensitive | 322 | 141 | — |
| + A passes the validity guard (`frac_expressed ≥ 0.20`) | **294** | **123** | — |

The middle filter removes nothing: every GDSC2 target-gene/drug pair already has
RNA coverage and clears the line and sensitivity minima. The validity guard is
the only binding constraint, costing 28 pairs and 18 genes.

Source: `validation/prepared/gdsc_scored_ready.parquet` (verified 100% GDSC2, no
GDSC1 contamination — `cell_similarity/08`). Median 640 lines per pair.

**294 usable pairs against a kill-switch threshold of 20.** Pair count is not the
constraint.

---

## 2. The marginal, reproduced on this exact pair set

Not the 141-target aggregate the gate audit used — the 294 (gene, drug) pairs
that a Stage 1 test would actually run on. Per-gene prevalence-matched
binarisation (threshold transfer across assays is invalid, `gate_audit/07`).
**Cluster bootstrap over cell lines, not pairs** — 294 pairs are built from ~640
shared lines, so a pair bootstrap would repeat the pseudo-replication fault the
gate audit measured at a design factor of 2.4×.

A working abundance gate predicts decile-1 ratio ≪ 1: low-expression lines
depleted of sensitivity.

| q (target prevalence) | pairs | decile-1 ratio | 95% CI (line cluster boot) | frac < 1 |
|---|---|---|---|---|
| **0.05** | 294 | **0.8646** | **[0.6313, 0.8882]** | 0.57 |
| 0.10 | 294 | 1.0043 | [0.8049, 1.0153] | 0.49 |
| 0.25 | 294 | 1.0015 | [0.9275, 1.0413] | 0.49 |

### This partly contradicts the premise the task was written on

The brief states *"The marginal effect is zero."* On this pair set it is **not
zero at the strictest threshold**: at q = 0.05 the decile-1 ratio is 0.865 with a
line-clustered interval excluding 1.0 (upper bound 0.888) — roughly **13.5%
depletion**. At q = 0.10 and q = 0.25 the intervals span 1.0 and the marginal is
flat, as expected.

Reasons to treat the q = 0.05 result cautiously, stated so it is not over-read:

- **One of three thresholds**, with no multiplicity correction. Three tests, one
  positive, is weak on its own.
- **q = 0.05 is the noisiest cell** — about 32 positives per pair, so the
  decile-1 estimate rests on few events.
- It differs from `gate_audit/07`'s GDSC null (0.907 at q = 0.02, n.s.) because
  that measured a *different object*: sensitivity to **any** drug annotated to
  the gene, aggregated per gene. This is per (gene, drug) pair. Both can be
  correct.

**Consequence for Stage 1, and it makes the test harder rather than easier.** The
comparator in B2 is single-gene A. If A has a real marginal at q = 0.05, then a
combined score must beat *that*, not beat zero. The brief's framing — "the only
case for a combined score is that the interaction carries signal where the
marginals do not" — needs amending at q = 0.05: there, the case would have to be
that the interaction adds to a marginal that already exists.

---

---

## 2b. Is the q = 0.05 marginal carried by one or two lineages?

`src/pipeline/combined_score_b0_lineage.py`. Two decompositions, because they
answer different questions. Cluster bootstrap over cell lines throughout.

### Leave one lineage out — clean

Recompute the pooled marginal with each lineage removed in turn. If dropping one
lineage collapses the effect, that lineage carried it.

| lineage dropped | lines | d1 ratio | shift | 95% CI |
|---|---|---|---|---|
| blood | 55 | 0.9404 | +0.0758 | [0.653, 0.940] |
| lung | 144 | 0.8064 | −0.0582 | [0.472, 0.846] |
| ovary | 34 | 0.9153 | +0.0507 | [0.630, 0.917] |
| skin | 37 | 0.9146 | +0.0500 | [0.642, 0.921] |
| colorectal | 42 | 0.9126 | +0.0480 | [0.630, 0.917] |
| breast | 47 | 0.9119 | +0.0473 | [0.632, 0.916] |
| pancreas | 28 | 0.9105 | +0.0459 | [0.632, 0.915] |

**0 of 10 lineages break it.** Every leave-one-out interval still excludes 1.0.
The largest single influence is `blood` — removing it moves the ratio from 0.865
to 0.940 — but even then the effect survives. **The marginal is not carried by
one or two lineages.**

### Within lineage — not estimable at this coverage

| lineage | lines | pairs | d1 ratio | 95% CI | |
|---|---|---|---|---|---|
| colorectal | 42 | 136 | 0.8200 | [0.335, 1.252] | estimable |
| lymphocyte | 47 | 154 | 0.9100 | [0.723, 1.451] | estimable |
| breast | 47 | 190 | 0.8800 | [0.000, 0.900] | **NOT ESTIMABLE** |
| blood | 55 | 216 | 0.9200 | [0.000, 0.951] | **NOT ESTIMABLE** |
| lung | 144 | 277 | 1.0182 | [0.000, 0.972] | **NOT ESTIMABLE** |

Three of five produce an interval that **does not contain its own point
estimate**, with a lower bound pinned at exactly 0. That is not a result, it is a
diagnostic: 711 GDSC∩expression lines spread over 28 lineages leaves at most ~144
in one lineage, so at q = 0.05 a pair carries ~7 positives and the decile-1 cell
(14 lines) is frequently empty. The per-pair ratio is then exactly 0, the
distribution is zero-inflated, and the median and its bootstrap both become
unreliable.

The script now flags this automatically rather than printing the numbers as
findings. Of the two estimable lineages, neither independently excludes 1.0.

**Answer to the question:** the marginal is **not fragile** — leave-one-out is
clean. But it is **not demonstrated as a within-lineage property either**,
because that decomposition cannot be estimated at this coverage. If Stage 1 runs,
it should report leave-one-out robustness and must not claim per-lineage effects.

---

## 2c. Inclusion criterion for the modifier list — B must come from the 4,173

The abstention measurement in [MULTIGENE_COVERAGE.md](MULTIGENE_COVERAGE.md) has
a direct consequence for B1 that was not in the original plan:

**`NOT altered(G)` is trivially satisfied for 78.2% of genes. Only 4,173 of
18,623 have ≥10 lines carrying a driver alteration.**

A modifier that never varies cannot modify anything, so B must be drawn from
those 4,173. This is a **principled restriction, not an outcome-driven one** — it
is a property of the alteration layer alone and is fixed before any drug response
is examined.

Effect on the search space:

| | combinations |
|---|---|
| naive: 294 pairs × all candidate B | ~333,000 |
| **restricted: B from the 4,173 varying genes only** | **cut by >75%** |

This must be written into the pre-registration as an inclusion criterion, with
the per-pair count reported: candidate B genes with ≥10 altered **and** ≥10
unaltered lines *within that pair's own line set* (median 1,135 per pair, already
measured in §3).

It reduces the multiplicity problem but does not remove it. At α = 0.05, even
~80,000 combinations would yield on the order of 4,000 significant by chance.
The pre-registered list remains mandatory.

---

## 3. Three-way coverage capacity

| quantity | value |
|---|---|
| pairs with ≥100 lines having expression + alteration + drug | **294** |
| median lines per pair after three-way intersection | **638** |
| median candidate modifier genes per pair (≥10 altered *and* ≥10 unaltered) | **1,135** |

Three-way coverage costs almost nothing here: 640 → 638 lines. That is because
1,424 of 1,673 lines are sequenced (85.1%) and the GDSC2 lines are
disproportionately well-characterised.

**This is capacity, not a selection.** It says the data could support an
interaction test on ~294 pairs with ~1,135 candidate modifiers each. It does not
name a modifier, and deliberately so — see below.

---

## Verdict

| gate | result |
|---|---|
| B0 kill switch (< 20 usable pairs → stop) | **does not fire** — 294 pairs |
| marginal reproduced on the exact pair set | **done** — flat at q = 0.10/0.25, **not flat at q = 0.05** |
| q = 0.05 marginal, leave-one-lineage-out | **robust** — 0 of 10 lineages break it |
| q = 0.05 marginal, within lineage | **not estimable** at this coverage — do not claim per-lineage effects |
| three-way coverage | **ample** — 638 median lines, 1,135 candidate B per pair |
| modifier inclusion criterion | **B must come from the 4,173 genes with ≥10 altered lines** |
| **B1 pre-registered modifier list** | **NOT SUPPLIED — Task B stops here** |

### Why B1 blocks, and why I did not work around it

B1: *"The user supplies the A–B pairs, with citations, before anything is run. Do
not search for modifiers that work."*

With 1,135 candidate modifiers per pair and 294 pairs, roughly **333,000**
A–B–drug combinations are available. Selecting B by outcome from that space is
the `.head(30)` error at a scale where it would be guaranteed to produce a
significant-looking result, and no subsequent battery could rescue it. So the
modifier list has to arrive as an input.

**To start Stage 1 I need, committed before any measurement:**

1. A list of A–B pairs — target gene, modifier gene, direction of the expected
   effect. **B must be drawn from the 4,173 genes with ≥10 altered lines**
   (§2c); a modifier that never varies cannot modify anything.
2. A citation for each, establishing the modifier relationship independently of
   this dataset.
3. Confirmation of the term form. The motivating example is EGFR therapy failing
   in KRAS-mutant lines, which requires `expressed(A) AND NOT altered(B)` —
   KRAS is expressed in nearly every lung line, so `NOT expressed(KRAS)` carries
   almost nothing. The `altered()` term exists in Task A (A2) for this reason.
4. The **adoption rule**, written before measuring: what result adopts the
   combined score and what result abandons it.

Item 4 is not a formality. Writing it after seeing the effect size is how a
negative becomes a positive.

### What ships now

Task A — the Kleene layer — is complete and independent of all of this. It makes
no prediction, produces no combined score, and returns `combined_score: None`
with a pointer to this file.

`COMBINED_SCORE_RESULT.md` is not written, because Stage 1 has not run.
