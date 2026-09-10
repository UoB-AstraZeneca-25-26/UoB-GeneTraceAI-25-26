# MULTIGENE_GROUND_TRUTH — single-gene test against named cell lines

53 rows, 20 genes, cell lines with published molecular status. The layer either
reproduces the textbook fact or it does not. Script
`src/pipeline/test_kleene_ground_truth.py` →
`src/pipeline/outputs/kleene_ground_truth_results.json`.

**Headline after fixes: 38 PASS / 9 FAIL / 5 UNKNOWN on 52 resolved lines —
80.9% agreement.** The `mutation` category is now at **92%**.

| stage | PASS | FAIL | UNKNOWN | agreement |
|---|---|---|---|---|
| as first run | 37 | 12 | 3 | 75.5% |
| after Fix 4 (AR/LNCaP corrected) | 38 | 11 | 3 | 77.6% |
| **after Fix 1 (MODERATE → UNKNOWN)** | **38** | **9** | **5** | **80.9%** |

Fix 1 converted two confident wrong answers into honest abstentions; it did not
manufacture a single new PASS. That is the intended trade.

The remaining 9 failures split into three causes, all documented below and none
of them a logic defect.

Scope limits were declared *before* measuring: `altered()` reads mutations and
in-frame fusions, not copy number; and "null" is usually a protein-level
statement, so mRNA may persist.

---

## By category

Final, after all four fixes:

| category | n | pass | fail | unknown | agreement |
|---|---|---|---|---|---|
| amplification | 8 | 8 | 0 | 0 | **100%** |
| overexpression | 1 | 1 | 0 | 0 | **100%** |
| wild-type | 2 | 2 | 0 | 0 | **100%** |
| **mutation** | 28 | 22 | 2 | 4 | **92%** |
| fusion | 4 | 3 | 1 | 0 | 75% |
| deletion / null | 9 | 2 | 6 | 1 | **25%** |

`mutation` rose from 84% to 92% because Fix 1 moved the two EGFR rows out of FAIL
into UNKNOWN. `wild-type` reached 100% because Fix 4 recategorised AR/LNCaP as
the mutant it actually is. `deletion / null` is unchanged at 25% and is expected
to stay there — it is a property of transcript-versus-protein, addressed by the
`LOSS_IS_PROTEIN_LEVEL` warning rather than by a logic change.

Expression-based categories are near-perfect: every ERBB2-amplified line
(SKBR3 10.84, BT474 10.52, AU565 11.61, NCI-N87 11.76), both MET-amplified lines
(EBC-1 10.09, Hs746T 11.09), MYC-amplified HL-60 (8.91), AR-amplified VCaP (8.70)
and EGFR-overexpressing A431 all resolve TRUE.

Point mutations are strong: all four KRAS lines (G12S/G13D/G12V/G12D), all three
BRAF V600E lines, both NRAS Q61 lines, both PIK3CA lines, both BRCA1 lines,
BRCA2 CAPAN-1, both STK11 lines, AR 22Rv1 and APC SW480 all resolve correctly.

---

## The failures — four causes, not one

Twelve as first run, nine after the fixes. The causes are unchanged;
what changed is that two of them are now abstentions and one was an
error in the reference.

### 1. "Null" is protein-level, mRNA persists (6 rows) — expectation wrong, not the layer

| gene | line | expected | got | log2TPM |
|---|---|---|---|---|
| TP53 | SAOS2 | FALSE | TRUE | 5.09 |
| TP53 | H1299 | FALSE | TRUE | 4.65 |
| PTEN | U87MG | FALSE | TRUE | 4.93 |
| PTEN | LNCaP | FALSE | TRUE | 4.18 |
| PTEN | MDA-MB-468 | FALSE | TRUE | 4.33 |
| RB1 | SAOS2 | FALSE | TRUE | 3.40 |

This is the caveat declared up front, confirmed at scale: the `deletion_null`
category scores 25%. These lines are functionally protein-null but still
transcribe the locus — nonsense transcripts escaping decay, truncating mutations,
or partial deletions retaining the probe region. **`expressed()` is doing exactly
what it says: reporting transcript, not protein.** The lesson is for query
design — `NOT expressed(PTEN)` will not find PTEN-null lines, and `altered(PTEN)`
is the correct term for that question.

By contrast the two genuinely transcript-level losses pass: CDKN2A in U2OS
(0.30) and RB1 in WERI-Rb-1 (0.00).

### 2. In-frame indels are missed by `any_driver` (2 rows) — FIXED, now UNKNOWN

| gene | line | expected | got | mutation row exists? | `any_driver` | `max_vep_rank` |
|---|---|---|---|---|---|---|
| EGFR | PC9 | TRUE | FALSE | **yes** | **False** | 2 |
| EGFR | HCC827 | TRUE | FALSE | **yes** | **False** | 2 |

The exon-19 deletion **is recorded** in `mutations_collapsed.parquet`. It is not
flagged as a driver. The mechanism is visible in the table:

| `max_vep_rank` | rows | fraction with `any_driver` |
|---|---|---|
| 3 (HIGH impact) | 110,670 | **100%** |
| 2 (MODERATE impact) | 521,934 | **0.5%** |

VEP scores an in-frame deletion as MODERATE, so activating in-frame indels land
in rank 2 and miss the flag. EGFR exon-19 deletion is the canonical example of
this class, and it is one of the most consequential activating mutations in
oncology.

**Consequence, and it runs opposite to the fault the Kleene layer was built
for.** The layer protects against *unmeasured → negative*. This was
*measured-but-under-annotated → negative*: `NOT altered(EGFR)` **wrongly included
PC9 and HCC827** as confirmed non-altered, with full confidence and no UNKNOWN
flag. A false Match, not a false Possible.

**Fixed.** `altered()` now returns UNKNOWN when a sequenced line has no
HIGH-impact row, at least one MODERATE row, and the gene is on the pre-committed
22-gene uncertainty list (CGC Tier 1 oncogene with `"O"` in `MUTATION_TYPES`).
PC9 and HCC827 both moved FALSE → UNKNOWN, and no other ground-truth row changed.
Cost: 0% of Matches for a random modifier, median 3.05% when the modifier is on
the list. See `MULTIGENE_SPEC.md` A2.

### 3. Structural variants absent from the mutation table (3 rows)

| gene | line | lesion | mutation rows |
|---|---|---|---|
| MYC | Raji | IGH-MYC t(8;14) translocation | **none** |
| CTNNB1 | HepG2 | exon 3 deletion | **none** |
| CTNNB1 | DLD-1 | mutant | **none** |

No row at all for the gene in that line. The MYC translocation is an
enhancer-juxtaposition event, not a protein fusion, so it is invisible to both
the SNV table and the in-frame-fusion table. `altered()` correctly returns FALSE
under its own definition — the line is sequenced and no driver is recorded — but
the biological fact is a driver alteration.

### 4. The supplied annotation is wrong (1 row) — CORRECTED

| gene | line | supplied | layer says |
|---|---|---|---|
| AR | LNCaP | "wild-type, hormone-sensitive" | altered = **TRUE** |

LNCaP carries the well-known **AR T877A** (T878A) point mutation in the
ligand-binding domain — it is the classic promiscuous-AR model, hormone-sensitive
but not wild-type. **The layer was right and the ground truth was wrong.**

**Corrected 2026-08-10** in `test_kleene_ground_truth.py`: the row now reads
`("AR", "LNCAP", "T878A mutant [corrected]", "mutation")` and passes. The
original annotation and the reason for the change are kept in a comment beside
it, because a ground-truth exercise in which the system corrects the reference is
worth recording rather than quietly editing away.

---

## Would a widened `altered()` fix it? No.

Tested rather than assumed: `any_driver OR oncogene_hit OR tsg_hit`, with the
fusion contribution retained in both arms.

| rule | correct of 34 altered() rows |
|---|---|
| `any_driver` (current) | 26 |
| widened | 26 |

**0 calls change.** The failures are not a flag-threshold problem. The EGFR rows
exist but are scored MODERATE; the CTNNB1 and MYC lesions have no row at all.
`altered()` is left unchanged, because widening it would add false positives
without recovering a single true one.

*(A first version of this check omitted fusions from the widened arm and appeared
to flip the three ALK calls to FALSE. That was an artefact of the check. Fixed;
the corrected comparison shows zero changes.)*

---

## The 5 UNKNOWNs — the layer declining to call

| gene | line | term | why |
|---|---|---|---|
| VHL | RCC4 | altered | line never sequenced |
| APC | Caco-2 | altered | line never sequenced |
| CDKN2A | MIA PaCa-2 | expressed | log2TPM 0.38, inside the ±0.699 measurement band |
| **EGFR** | **PC9** | **altered** | **only MODERATE-impact variants; in-frame indels not rankable by VEP severity** |
| **EGFR** | **HCC827** | **altered** | **same** |

All five are correct abstentions. The last two are Fix 1 working: they were
confident FALSE calls that would have admitted PC9 and HCC827 to Matches. MIA PaCa-2 is the uncertainty band earning its
place on a real row: 0.38 is below the 1.0 threshold but not by more than
measurement error, so the layer declines rather than calling FALSE. RCC4 and
Caco-2 are the never-sequenced case that the whole three-valued design exists
for.

Only **P493-6** failed to resolve to DepMap — correctly, since it is an
engineered tet-inducible construct, not a catalogue line.

---

## What this changes

**Nothing in the layer's logic.** Kleene propagation, the buckets, the coverage
reporting and the abstention rules all behave as specified.

**One documented defect to carry forward:** `altered()` has a systematic
false-negative class — activating in-frame indels (VEP MODERATE) and structural
variants absent from the SNV table. This should be stated wherever `altered()` is
offered, because unlike an UNKNOWN it produces a *confident wrong answer*.

**Option 3 was adopted** — return UNKNOWN rather than FALSE where only
MODERATE-impact rows are present on a pre-committed uncertainty gene. It converts
a confident wrong answer into an honest abstention, which is the layer's own
principle applied to itself.

Still open, and declared rather than fixed:

1. **Structural variants remain invisible.** CTNNB1 exon-3 and IGH-MYC have no
   row at all, so no severity rule can reach them. Wiring the copy-number layer
   would make exon-level deletion and amplification visible; that is a
   data-integration job, not a logic change.
2. **A variant-level hotspot source** (OncoKB / COSMIC MutantCensus) would be
   stronger than the gene-level list now in use, which flags whole genes as
   untrustworthy rather than identifying activating variants.
