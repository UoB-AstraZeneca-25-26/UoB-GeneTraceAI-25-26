# COMBINED_SCORE_RESULT — Task B, Stage 1

**Verdict: the B1 kill switch fired. 2 of 10 established pairs are eligible,
against a pre-registered threshold of 10. Stage 1 does not run. Task A ships.**

Pre-registration: [COMBINED_SCORE_PREREGISTRATION.md](COMBINED_SCORE_PREREGISTRATION.md),
committed before this measurement. Precondition:
[COMBINED_SCORE_PRECONDITION.md](COMBINED_SCORE_PRECONDITION.md).
Script `src/pipeline/combined_score_b1_eligibility.py` →
`src/pipeline/outputs/combined_score_b1_eligibility.json`.

**No drug-response value was read.** Eligibility used only expression coverage,
alteration counts, and whether an (A, drug) pair exists with enough lines. No
enrichment, no sensitivity rate, no interaction statistic was computed. The
negative below cannot have been influenced by an outcome.

---

## Eligibility

### Established set (n = 10)

| A | B | eligible | drugs | lines | B altered | B unaltered | detail |
|---|---|---|---|---|---|---|---|
| **EGFR** | **KRAS** | **YES** | 8 | 699 | 110 | 589 | AZD3759 |
| BRAF | MAP2K1 | no | | | | | criterion 4/5 — no (A, drug) pair with ≥100 lines and ≥10 per B arm |
| BRCA1 | BRCA2 | no | | | | | criterion 1 — A is not a GDSC2 drug target |
| TP53 | MDM2 | no | | | | | criterion 4/5 |
| PTEN | PIK3CA | no | | | | | criterion 1 |
| **ERBB2** | **ERBB3** | **YES** | 3 | 657 | 15 | 642 | Sapitinib |
| ALK | EML4 | no | | | | | criterion 2 — A fails the validity guard (11.1% expressed) |
| MYC | MAX | no | | | | | criterion 1 |
| RB1 | CDK6 | no | | | | | criterion 1 |
| VHL | HIF1A | no | | | | | criterion 1 |

### Contrast set (n = 10)

| A | B | eligible | lines | B altered | detail |
|---|---|---|---|---|---|
| TP53 | EML4 | YES | 695 | 13 | MIRA-1 |
| ERBB2 | TP53 | YES | 657 | 448 | Sapitinib |
| CDK6 | VHL | YES | 699 | 11 | Ribociclib |
| GAPDH, ACTB, BRCA1, PTEN, VHL, RB1, MAX | | no (7) | | | criterion 1 |

**2 established eligible, 3 contrast eligible. Kill switch: fired.**

---

## Why it failed, and it is structural rather than incidental

**Six of ten established A genes are not GDSC2 drug targets:** BRCA1, PTEN, MYC,
RB1, VHL — and in the contrast set GAPDH, ACTB, MAX.

That is not a coverage accident that a larger panel would fix. The Stage 1 design
requires the readout to be **sensitivity to a drug that targets A**, so A must be
druggable. But the biology of modifier relationships is dominated by **tumour
suppressors and transcription factors as A** — BRCA1, PTEN, RB1, VHL, MYC. Those
are precisely the genes for which no direct inhibitor exists, and therefore
precisely the genes that can never appear as A in this design.

The supplied list is a fair sample of textbook modifier biology. Its failure rate
here measures a **mismatch between the question and the assay**, not a defect in
the list.

Two further failures, both informative:

- **ALK fails the expression validity guard at 11.1%.** ALK is expressed in under
  20% of the panel, so its within-gene percentile is not a meaningful relative
  position — exactly the condition the guard exists to catch. The EML4–ALK
  relationship is real; it is a *fusion*, and fusion presence is not the same
  quantity as expression percentile.
- **BRAF–MAP2K1 and TP53–MDM2 fail criterion 4/5.** A drug and ≥100 lines exist,
  but the B arms do not partition: too few lines carry a driver alteration in
  MAP2K1 or MDM2 within the pair's own line set. This is the §2c inclusion
  criterion biting exactly as intended — `NOT altered(MAP2K1)` is satisfied by
  essentially every line, so it cannot modify anything.

### The contrast set surviving better is meaningless, and worth saying so

3 contrast pairs are eligible against 2 established. That is **not** evidence
about the biology. ERBB2, TP53 and CDK6 happen to be druggable in GDSC2, so pairs
naming them survive criterion 1 regardless of whether the relationship is real.
Eligibility measures assay coverage, not truth.

It does have one consequence: with 2 established and 3 contrast pairs, the
contrast set could not have functioned as a control even if Stage 1 had run. The
design needs both arms populated.

---

## What was not done, deliberately

- **No Stage 1 statistic was computed.** The kill switch is honoured as written.
- **No modifier was substituted** to rescue the pair count. Swapping BRCA1→PARP1,
  or picking a different B for BRAF, would be choosing pairs by eligibility after
  seeing which ones failed — the same class of error as choosing by outcome, one
  step removed. If a revised list is wanted it must be written and committed as a
  new pre-registration.
- **No relaxation of the criteria.** Dropping criterion 4 to ≥5 altered lines
  would admit BRAF–MAP2K1 and TP53–MDM2 and reach 4 established pairs. Still short
  of 10, and the threshold was fixed before measurement precisely so it could not
  be moved afterwards.

---

## What would make Task B viable

Not a recommendation to proceed — an honest statement of what the obstacle is.

1. **A must be druggable.** A list drawn from kinase and receptor targets present
   in GDSC2 — the 123 valid target genes in the precondition §1 — would clear
   criterion 1 by construction. The cost is that it excludes most tumour-suppressor
   modifier biology, which is where the strongest published relationships live.
2. **B must vary.** ≥10 altered lines within the pair's own line set. Only 4,173
   of 18,623 genes qualify panel-wide, and fewer within any one pair's lines.
3. **Both constraints at once** is the binding problem: druggable A **and** a
   varying B **and** a documented relationship between them. EGFR–KRAS satisfies
   all three, which is why it is the motivating example — it may be close to the
   only clean case at this coverage.

A viable list would need ~20 pairs meeting all three, written before measurement.
Whether that many exist in published resistance biology is the question I raised
before and cannot answer from the data.

---

## Cost

Task B consumed one afternoon: the precondition (B0), the per-lineage follow-up,
the pre-registration, and this eligibility check. The staging worked as designed
— an unfavourable answer cost hours rather than the three weeks the similarity
graph took.

**Task A is unaffected and ships.** See
[MULTIGENE_SPEC.md](MULTIGENE_SPEC.md). `query()` continues to return
`combined_score: None`.
