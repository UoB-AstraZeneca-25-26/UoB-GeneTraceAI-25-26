# DECISIONS.md

Reflective log of predictions made *before* running a diagnostic, and what
actually happened. One entry per decision point. Newest first.

---

## 2026-07-30 — Block 2: within-stratum validation gate (141-gene GDSC ground truth)

**Context:** Direct continuation of the entry above (the wide 141-gene
`gdsc_scored_ready.parquet` ground truth) and the ladder entry before it
(disease -> lineage -> full panel). Question: for how many (gene, stratum)
cells is there enough GDSC ground truth to say anything at all, and is the
within-stratum baseline easier or harder than the whole-panel one. Gate:
proceed if >=30 usable cells spanning >=15 distinct genes, and >=8 genes
with >=3 usable strata each (enough to compare a gene's ranking across
strata, not just measure it once).

**Outcome: passes by a wide margin on all three thresholds.**

| Criterion | Bar | Actual |
|---|---|---|
| Usable cells | >=30 | **2,421** |
| Distinct genes with >=1 usable cell | >=15 | **141 / 141** |
| Genes with >=3 usable strata | >=8 | **140** |

Pipeline: 1,485 scored lines -> 1,418 with a ladder assignment (825 lineage
rung, 518 disease rung, 75 fell through to full panel) -> joined against
87,618 GDSC (gene, line, sensitive) pairs after dropping 26.3% not in the
scored set -> 5,222 total (gene, stratum) cells -> 3,072 evaluable
(>=10 lines, >=2 sensitive) -> 2,421 usable after removing 651 degenerate
cells (baseline >=0.5, where any ranking looks good by construction).

**Collapse-rule choice, logged as instructed ("this choice matters"):**
GDSC has a median of 2 drugs per target gene (max 8, 76/141 genes have
more than one). Collapsing multi-drug (gene, line) pairs to a single
sensitivity label via **ANY** drug sensitive gives 33.9% sensitive overall;
via **ALL** drugs sensitive gives 18.7%. Used ANY (the script's stated
default) — not re-justified here beyond that it's the more permissive,
better-powered choice for a feasibility gate; worth revisiting if/when this
becomes a scoring input rather than a diagnostic, since ANY vs ALL changes
what "ground truth sensitive" means for genes with several drugs of
differing potency.

**The baseline-inflation check the script warns about did NOT occur, and
that's a good sign, not a null result:** within-stratum baseline (median
0.237, mean 0.251) is *lower* than the whole-panel baseline (median 0.277,
mean 0.316), not higher. Per the script's own logic ("if within-stratum
baseline is HIGHER, any hit-rate gain is partly denominator, not signal"),
this means a future within-stratum hit-rate measurement would be a
*harder*, more conservative test than the whole-panel one currently
reported — any lift found there is less likely to be a denominator
artifact, not more.

**Conclusion:** all three gates in this chain now pass —
stratum-availability (category/diseases/parent_hierarchy individually:
fail; ladder: pass), nesting (93.1%, pass), and ground-truth-density
(this entry: pass by 1-2 orders of magnitude on every threshold). Building
within-stratum percentile scoring (v2) is now measured-viable, not just
computable — the open question left is scoring-mechanism design (how a
within-stratum percentile actually gets computed and blended, if at all),
not data availability.

---

## 2026-07-30 — Stratum profiling gate: is `category`/`diseases`/`parent_hierarchy` usable for within-lineage percentiles?

**Context:** Before building within-lineage percentile scoring (a "v2"
candidate), gate it on whether `cell_line_lookup.parquet` actually has a
column with enough per-stratum coverage to compute a meaningful percentile.
Decision rule (given, not mine): viable if at least one column gives ≥15
strata at n≥20 covering ≥50% of scored lines; if the best column gives a
handful of strata covering a fifth of the panel, stop and write it up as a
negative rather than build it.

Requester's prediction: `category` usable at n≥20 for ~15-25 strata
covering roughly half the scored lines; `diseases` too granular.

**My prediction, made before running Block 1 — and it diverges from the
requester's on `category` specifically:**

I already looked at `category` in this session (checking it as a lineage
candidate for the metabolomics lineage-confound diagnostic, several turns
back) and it is **not a tissue/lineage field at all** — it's cell-line
*type*, with exactly 7 distinct values, heavily skewed:
`cancer cell line` 1699, `transformed cell line` 76, `undefined cell line
type` 38, `None` 22, `spontaneously immortalized cell line` 2,
`telomerase immortalized cell line` 2, `finite cell line` 1 (n=1840 total,
pre-scoring-filter). I predict **`category` fails the decision rule
outright** — at most 2-3 strata clear n≥20 (`cancer cell line`,
`transformed cell line`, maybe `undefined cell line type`), and one bucket
alone (`cancer cell line`) covers ~90%+ of all lines, which defeats the
purpose of stratification (nearly everything lands in one stratum, not 15-25
meaningfully-sized ones).

`diseases`: I agree with the requester here — from the earlier lineage
work I sampled it directly (free-text NCIT/ORDO terms like "gastric
adenocarcinoma", "ewing sarcoma"), and it looks like hundreds of distinct
values, most represented by 1-2 lines. I predict very few `diseases` strata
clear n≥20, and it fails the ≥50%-coverage bar even if some individual
strata pass n≥20.

`parent_hierarchy`: neither of us predicted this one. From the sample
values I saw earlier (`'cvcl_0531 ! sk-n-sh'`, `'cvcl_d317 ! 1181n1'`), this
looks like a per-cell-line CVCL/name pair, near-unique per row — I predict
it is not a grouping column at all (distinct-value count close to row
count), and fails trivially.

**Net prediction: none of the three columns tested in Block 1 pass the
gate.** Separately (not part of Block 1, flagging for the record): I
already know of a column that *would* likely satisfy the requester's
original pattern — `lineage` in `data/parquet/data_clean/sample_info_clean.parquet`
(a different file, not `cell_line_lookup.parquet`), which I profiled in the
earlier lineage-confound entries: 30 distinct lineages, e.g. lung=274,
blood=141, skin=120, lymphocyte=110, CNS=109, breast=86, colorectal=85,
bone=79, upper_aerodigestive=79, ovary=74 (n=1840 total). That shape (30
strata, several well above n=20, no single stratum dominating) looks like
what the requester expected from `category`. I'm not substituting it into
Block 1 — running the given code unmodified first — but flagging it now so
it's on the record before the result comes back, not after.

**Outcome:** My prediction was **directionally right on all three columns**;
the requester's prediction was **wrong on which column the numbers belong
to** — the shape they described for `category` is almost exactly what
`diseases` produced instead.

Scored lines: 1,418 of 1,840 lookup rows (1,485 distinct scored `model_id`,
67 not present in the lookup at all).

| Column | Distinct strata | n≥20: strata / coverage | Verdict at n≥20 |
|---|---|---|---|
| `category` | 4 | **3 strata / 99.0%** | FAILS (strata count) — 1,346/1,418 lines (94.9%) sit in one bucket, `cancer cell line`. Only `undefined cell line type` (35) and `transformed cell line` (23) clear n≥20 alongside it. |
| `diseases` | 247 | **16 strata / 36.5%** | FAILS (coverage) — clears the ≥15-strata bar by one, but falls well short of ≥50% coverage. At n≥15: 25 strata / 47.2% (still short). At n≥10: 37 strata / **57.4%** (clears both bars). |
| `parent_hierarchy` | 65 | **0 strata / 0.0%** | FAILS outright — median stratum size 1, max 3. Not a grouping column (95.1% null on top of that). |

Applying the decision rule exactly as given (≥15 strata at n≥20, ≥50%
coverage): **none of the three columns pass.** `category` fails because it
isn't stratifying anything (one dominant bucket, exactly as I predicted from
already having seen its value_counts). `diseases` fails on coverage at
n≥20, but is close, and passes cleanly if the threshold relaxes to n≥10.
`parent_hierarchy` fails trivially, exactly as predicted.

The requester's predicted shape — "~15-25 strata at n≥20, roughly half the
panel" — doesn't land on `category` (way off: 3 strata, 99% coverage, all
in one bucket) but comes close to landing on `diseases` (16 strata at n≥20,
36.5% coverage; 25 strata / 47.2% at n≥15). It looks like the right ballpark
number was predicted, attributed to the wrong column.

Not run in Block 1 (different file, flagged in the prediction above but not
substituted in): `lineage` from `sample_info_clean.parquet` — 30 lineages,
several well above n=20 with no single dominant bucket, which was always
the more likely candidate for the pattern being tested. Worth running the
same profiling against it before concluding v2 is dead, rather than stopping
on `diseases` alone at a threshold it narrowly misses.

**Decision needed before Block 2:** the gate as literally specified did not
pass, so per "run only if Block 1 passes," Block 2 was not run. Options:
(1) treat this as the documented negative and stop, (2) relax to n≥10 where
`diseases` clears both bars, or (3) profile `lineage` from
`sample_info_clean.parquet` first, since it wasn't tested and looks like the
strongest candidate. Reported back to the requester rather than picking one.

**Addendum (requester's own note, logged verbatim because it's the correct
attribution and matters more than the numbers):** the requester chose
option 3 and was explicit that this was *their* error, not a data finding —
"My decision rule was a gate on those three columns, and I picked the
columns wrong... The gate was testing 'is there a usable lineage field,'
and I pointed you at files that don't contain one." My prediction that
`category` would fail (right) and the requester's prediction that it would
pass at "~15-25 strata, roughly half the panel" (wrong) reflects that split
of error precisely: I had already seen `category`'s value_counts earlier
this session, the requester was reasoning from the column name alone. The
magnitude match between the requester's predicted numbers and `diseases`'s
actual numbers (16 strata / 36.5% at n≥20) was **coincidence, not
insight** — `diseases` was never the column being reasoned about, `category`
was, and `category`'s real shape (3 strata, one bucket holding 94.9%) looks
nothing like either prediction.

---

## 2026-07-30 — Finding the wide gene-level GDSC ground truth (resolves the diversity-harness n=9 loose end)

**Context:** Block 2 (GDSC coverage per stratum) needs a ground-truth file
with a gene column, a per-line sensitivity label, and wide gene coverage —
`curated_validation_pairs.parquet` (9 genes) can't carry this gate alone:
stratifying 9 genes across ~22 lineage strata leaves most gene x stratum
cells with a handful of lines and often zero sensitive ones, close to
unfalsifiable. Printed schemas for the 3 other `gdsc_*` files plus
`curated_validation_pairs.parquet`, and searched the repo for anything wider.

**Found it:** `validation/prepared/gdsc_scored_ready.parquet` — 255,142
rows, **141 distinct genes** (`target_ensg`), 962 distinct `model_id`, with
a boolean `sensitive` column (63,976 True / 191,166 False). This is the
141-curated-gene ground truth the architecture doc references. It wasn't
missed for lacking "gdsc" in the name (it has it) — it was missed because
the glob's *first alphabetical match* (`gdsc_ach.parquet`) got used instead
without checking the other three.

Confirmed shape of the pipeline that produces these files:
`gdsc_long.parquet` (229,934 rows: drug response only, no gene column) +
`gdsc_pairs.parquet` (687 rows: drug->target map, 192 distinct target
genes, no response data) join on `drug_id` -> `gdsc_scored_ready.parquet`
(141 genes survive the join — fewer than 192 because not every mapped
target has actual screened response data). `curated_validation_pairs.parquet`
is a further hand-picked 9-gene subset of that (with `is_curated` /
`ground_truth` columns added on top), not a separate or wider source.

**Resolves the open question from the "Diversity-aware top-20" entry
above:** that harness ran on exactly 9 genes (EGFR, MAP2K1, MAP2K2, CDK4,
CDK6, BRAF, MET, ERBB2, BCL2) because `curated_validation_pairs.parquet`'s
`target_symbol` has exactly 9 distinct values — that file *is* its entire
gene universe, not an arbitrarily-small sample of a larger available set.
Nine wasn't a scoping choice made for that harness; it was inherited
unnoticed from which ground-truth file got loaded. The wider option
(`gdsc_scored_ready.parquet`, 141 genes) was sitting in the same directory
the whole time.

**What this means going forward:** per the stated rule — "a file with 100+
genes and a per-line sensitivity label -> run Block 2 against it, and the
gate is meaningful" — `gdsc_scored_ready.parquet` clears that bar (141
genes). The GDSC-coverage-per-stratum gate should be run against this file,
not the 9-gene one. Also worth flagging: the diversity-harness result
itself (mean delta-hit@20 ~= +0.0015, reported as "a wash") was measured on
the narrowest ground truth available in the repo, not a deliberately chosen
small-N check — that context belongs alongside the original number if it's
cited again.

---

## 2026-07-30 — Fallback ladder (disease → lineage → full panel) for within-stratum percentiles

**Context:** Direct continuation of the entry above. `lineage` (from
`sample_info_clean.parquet`, not `cell_line_lookup.parquet`) profiled
separately, then tested as a fallback rung underneath `diseases` rather than
as a standalone replacement for it — per the requester's point that
`diseases` (247 values) and `lineage` (30 values) are the fine and coarse
levels of one hierarchy, not two competing candidates. Ladder: disease
stratum if n≥20, else parent lineage if n≥20, else full panel. Revised gate:
passes if ≥60% of scored lines land on a real rung (not the full panel) and
≥90% of disease terms nest inside exactly one lineage.

My prediction from the entry above, restated for the record: `lineage`
would be the strongest of the four columns now tested — 30 categories, no
single dominant bucket. Not re-predicting the ladder's pass/fail here since
the requester specified the exact mechanism and gate; this entry records
the measurement, not a forecast.

**Outcome:** **Ladder passes both parts of the revised gate, clearly.**

- `lineage` alone: 30 distinct strata, n=1,418 scored lines, 100% have a
  lineage value (0 missing — better than either `category` or `diseases`).
  At n≥20: **22 strata / 94.7% coverage**. Largest stratum (`lung`, 212) is
  14.9% of scored lines — nothing like `category`'s 94.9%-in-one-bucket
  problem. This confirms the prediction: `lineage` is the strongest single
  column of the four tested (`category`, `diseases`, `parent_hierarchy`,
  `lineage`).
- Nesting check: **230/247 disease terms (93.1%) map to exactly one
  lineage** — clears the ≥90% bar. 17 violations exist (e.g. `melanoma`
  spans `skin`/`central_nervous_system`; `colon adenocarcinoma` spans
  `colorectal`/`lung` — almost certainly single mis-annotated lines rather
  than genuine biological ambiguity, worth a spot-check before v2 build but
  not disqualifying at this rate).
- Ladder simulation (MIN_N=20): **54.4% of scored lines resolve at the
  lineage rung, 40.8% at the disease rung, only 4.7% fall through to the
  full panel** → **95.3% land on a real stratum**, clearing the ≥60% gate
  by a wide margin. 39 distinct strata actually used across both rungs.
  Median stratum size: disease rung 42, lineage rung 65 — both comfortably
  above MIN_N=20, so no stratum is being kept alive right at the edge.

**Conclusion: v2 (within-stratum percentile scoring) is viable as a
disease→lineage→full-panel ladder.** This reverses the tentative negative
implied by the first entry — the negative was about three specific columns,
not about the underlying question, exactly as the requester flagged before
running anything. Block 2 (GDSC coverage per stratum) is the next gate per
the requester's stated sequencing — a stratum that's computable isn't
necessarily validatable.

---

## 2026-07-30 — Diversity-aware top-20 (best-line-per-cluster-then-fill) vs raw top-20

**Context:** Per-review guidance: similarity-within-the-returned-list should
be grouping, not re-ranking (every prior attempt to use similarity to
*reorder* failed because it was lineage-dominated -- here that's the
feature, not the bug). Built `cell_line_grouping.py`: metadata grouping by
`diseases` (cheap, explainable) first, expression-profile hierarchical
clustering second (only kept if it adds something metadata grouping
doesn't), then a diversity re-rank (`diversify_top_n` -- best line per
cluster, then fill) that never touches `core_score`. Testing it through the
same hit@20 harness used throughout, on the 9 genes with real GDSC ground
truth (EGFR, MAP2K1, MAP2K2, CDK4, CDK6, BRAF, MET, ERBB2, BCL2).

**Prediction, made before running:**

The requester's own prior was "a modest improvement." I agree directionally
less than the requester does, for a structural reason: hit@20 here is
defined as *recall* of sensitive lines within the top-20
(`hits / n_sensitive`, not precision) -- see `hit_at_k` in
`validate_chronos.py`. Diversifying necessarily swaps some raw-top-20
members for lower-ranked-but-cluster-distinct members drawn from the top-50
pool. That swap can only helps recall if a sensitive line that was *outside*
the raw top-20 but *inside* the top-50 pool gets pulled in by
diversification -- otherwise recall stays flat or drops. My prediction:
**flat-to-slightly-negative on average**, with the exception being
genes where sensitive lines are lineage-clustered *differently* from the
raw top-20's dominant cluster (i.e. the raw top-20 is genuinely redundant
in a way that's actively hiding a sensitive line just past position 20).
BCL2 is the one gene in this set where I'd bet on that exception showing up,
given its lineage concentration (blood/lymphoid) already surfaced in the
earlier lineage-confound entries.

**Outcome:** Overall magnitude prediction was **right**; the specific case
I called out was **wrong**.

| Gene | ARI (meta vs profile) | Raw hit@20 | Diversified hit@20 | Δ |
|---|---|---|---|---|
| EGFR | 0.228 | 0.0035 | 0.0140 | +0.0105 |
| MAP2K1 | 0.238 | 0.0111 | 0.0148 | +0.0037 |
| MAP2K2 | 0.556 | 0.0074 | 0.0074 | 0.0000 |
| CDK4 | 0.310 | 0.0137 | 0.0103 | -0.0034 |
| CDK6 | 0.484 | 0.0034 | 0.0103 | +0.0068 |
| BRAF | 0.400 | 0.0107 | 0.0107 | 0.0000 |
| MET | 0.210 | 0.0070 | 0.0070 | 0.0000 |
| ERBB2 | 0.173 | 0.0148 | 0.0111 | -0.0037 |
| BCL2 | 0.396 | 0.0000 | 0.0000 | 0.0000 |

Mean Δ ≈ **+0.0015** — essentially flat, 3 genes improved, 2 worsened, 4
unchanged, all small in magnitude. This matches my "flat-to-slightly"
prediction in direction and size, and lands closer to flat than to the
requester's "modest improvement." The specific exception I predicted --
BCL2 benefiting because its sensitive lines are lineage-clustered away from
the raw top-20 -- **did not materialize**: BCL2 shows zero change in either
direction, because its true GDSC-sensitive lines apparently don't fall
within the top-50 pool at all under the current 2-layer core_score (this
tracks with BCL2's weak/backwards baseline correlation already documented
in the earlier metabolomics-confound entry).

**ARI check (does profile clustering just recover metadata?):** No — ARI
ranged 0.17-0.56 across the 9 genes, moderate but well short of 1.0. Profile
clustering is not redundant with disease labels, but it's also not a
strong, universal improvement on them. Per the stated decision rule
("if clusters just recover disease labels, use the labels"), the ARI
values don't force that conclusion either way — recommend keeping metadata
grouping as the default (cheaper, fully explainable) and treating
profile-based clustering as a validated-but-optional enrichment, not a
required upgrade.

**What this means for the diversity re-rank:** it's a legitimate,
explainable presentation option (never touches `core_score`), but the
measured effect on hit@20 is a wash, not a win. Recommend presenting it as
an optional view ("show diverse contexts" toggle) rather than replacing the
raw ranked list by default -- the number doesn't support making it the
default behaviour, which is itself the reportable result the requester
asked for ("either result is reportable").

**Context:** Same question as the metabolomics entry below, applied to
miRNA: a per-cell-line "global miRNA burden" scalar (mean z-score across
all 734 miRNAs), added as a constant term rather than the curated
target-specific inhibitor already tested (`test_run_metabolomics_mirna.py`
Part A — mixed/weak result there: MET improved slightly, BCL2's already-odd
baseline direction shrank toward zero).

**Prediction, made before running:**

1. **Lineage ANOVA:** I expect miRNA to show a **larger** lineage effect
   than metabolomics did (R²=0.042). Reasoning: miRNA profiling is used
   clinically as a tissue-of-origin classifier precisely because several
   miRNAs are near-binary lineage markers (miR-142/144/451 for
   blood/haematopoietic, miR-1/133 for muscle, miR-124 for neural) — much
   more discretely tissue-restricted than most metabolite pools, which
   overlap continuously across tissues. Prior: **R² in the 0.08–0.20 range**,
   i.e. above metabolomics's 4.2% but still not dominant, since averaging
   across 734 miRNAs from many different tissue programs should still dilute
   any single lineage's signal.

2. **Constant-term blend (BCL2, MET):** Given blood-specific miRNAs
   (miR-142/144/451) are typically expressed at much larger fold-changes in
   haematopoietic lines than most metabolites are, I predict **higher risk
   of a real confound-driven inflation for BCL2 specifically** than the
   metabolomics version showed (which degraded BCL2 instead of inflating
   it). If blood/lymphocyte lines pull the mean miRNA z-score up sharply
   and those same lines are disproportionately venetoclax-sensitive, BCL2's
   rho could rise here even though metabolomics's rho fell. For MET, no
   comparably strong lineage-restricted miRNA prior exists, so I predict
   **little change, plausibly still a mild degradation** as with metabolomics.

**Outcome:** Prediction was **half right**.

1. **Lineage ANOVA:** R² = **0.072** (7.2%), F=2.975, p=2.5e-06 across 25
   lineage groups (n=948 lines). Direction correct — higher than
   metabolomics's 4.2%, as predicted — but landed just below my stated
   8–20% range. Also wrong on *which* lineages carry it: blood (+0.0415) and
   lymphocyte (-0.1018) sit near the middle of the distribution, not at the
   extremes. The real outliers are `prostate` (+0.263, n=7 — small group,
   noisy), `fibroblast` (+0.209), and `upper_aerodigestive` (-0.196) — none
   of which have an obvious BCL2 story. The "blood-specific miRNA markers
   dominate the mean" mechanism I proposed doesn't show up in the per-lineage
   breakdown.

2. **Constant-term blend (BCL2, MET):** **Wrong on the risky part of the
   prediction.** No confound-driven inflation for BCL2 — instead the same
   degradation pattern as metabolomics: rho 0.3034 → 0.2427 (Δ=-0.0607, 95%
   CI [-0.1850, +0.0782], not significant). MET also degraded in its
   meaningful direction: rho -0.2253 → -0.1591 (Δ=+0.0662, 95% CI
   [-0.0245, +0.1550], not significant but the lower bound sits right at
   zero) — this part matched the "mild degradation" prediction.

**What this means, combined with the metabolomics entry above:** two
independent modalities (metabolomics, miRNA), tested the same way, produced
the same qualitative result — a small-but-real lineage effect (4–7% of
variance, not dominant) and *no* confound-driven inflation on either
testable gene. The specific trap I was watching for (numbers improve for
the wrong reason) didn't materialize in either case. What both actually show
is that a single averaged, gene-agnostic scalar built by collapsing hundreds
of measurements into one number is closer to noise than to a usable global
signal — it drags correlation down slightly rather than up. That's a more
useful and more general negative result than either test alone: **the
global-modifier approach doesn't work here, and it's specifically because
of dilution-by-averaging, not because it's secretly a lineage proxy.**
Recommendation is unchanged — don't fold either into `core_score` — but the
justification is now measured rather than assumed for both modalities.

---

## 2026-07-30 — Metabolomics global per-line modifier: lineage confound check

**Context:** Tested whether a curated gene-level metabolomics layer helps or
hurts known-gene predictions (see `src/pipeline/test_run_metabolomics_mirna.py`,
Part B) — result was a null (no gene cleared rho>0.1, p<0.05 against Chronos
essentiality). Before abandoning metabolomics entirely, checking a cheaper
alternative: a single per-cell-line "global metabolic activity" scalar
(mean z-score across all 225 metabolites), added as a constant term rather
than a gene-specific one.

**The concern:** this global scalar has no gene axis by construction, so any
predictive power it shows is suspect — it could just be a lineage proxy
riding on GDSC's well-known lineage-structured drug response, not a real
metabolic signal.

**Prediction, made before running step 3 (the harness Δρ check):**

1. **Step 2 (lineage ANOVA):** I expect lineage to explain a *moderate but
   not dominant* share of variance in mean metabolite z-score — my prior is
   **R² in the 0.15–0.30 range**, p << 0.001 given n=928. Reasoning: CCLE
   metabolomics (Li et al. 2019, the source of this exact dataset) reports
   lineage/tissue-of-origin as one of the dominant axes of variation, but
   averaging across 225 metabolites spanning very different pathways (amino
   acids, TCA intermediates, lipids, nucleotides) should partially cancel
   lineage-specific signatures that point in different directions across
   metabolite classes, capping how much a single averaged scalar can carry.

2. **Step 3 (Δρ with global term added, BCL2 & MET):** I predict **BCL2 is
   at high risk of exactly the lineage-confound trap** described — BCL2
   inhibitors (venetoclax-class) are one of the most clinically
   well-established lineage-concentrated drug responses in oncology
   (dramatically more effective in blood/lymphoid lineage lines than solid
   tumours). If blood/lymphoid lines have a systematically different mean
   metabolic z-score from solid-tumour lines (plausible — blood-cell
   metabolism differs substantially from solid epithelial tumour
   metabolism), adding the global term should **inflate BCL2's rho against
   GDSC ground truth for the wrong reason** (lineage proxy, not BCL2
   biology). Predicted direction: rho improves for BCL2.
   For MET, I have no comparably strong lineage-concentration prior (MET-driven
   cancers span lung/gastric/hepatocellular, not one dominant lineage), so I
   predict **little to no change, possibly a small decrease**.

**If the prediction is right:** this is not evidence the layer works — it's
a measured demonstration of the trap. Per the reasoning already given: (a) it
would make the existing across-line coverage confounding strictly worse,
now confounded with lineage and `n_layers` too; (b) Stage 7's explanation
layer cannot produce an honest reason string for it — "this line scores
highly because it has generally high metabolite levels" explains nothing
about the queried gene. The correct action is to report the lineage
component separately and labelled, never fold it into `core_score`.

**Outcome:** Prediction was **wrong on both parts**.

1. **Step 2 (lineage ANOVA):** R² = **0.042** (4.2% of variance), F=1.713,
   p=0.0198 across 24 lineage groups (n=928 lines). Statistically
   significant only because n is large — the *effect size* is small, well
   below my 0.15–0.30 prior. Lineage does **not** dominate this particular
   metric. Per-lineage means also don't show the pattern I expected:
   `lymphocyte` (-0.0387) and `blood` (-0.0205) are the *lowest*-scoring
   lineages, not an extreme outlier group that would obviously carry a
   BCL2-relevant signal — `bone` (+0.0396) and `thyroid` (+0.0347) are the
   highest.

2. **Step 3 (constant-term blend, paired bootstrap):**
   - BCL2: rho **0.3040 → 0.1490** (Δ = **-0.1550**, 95% CI [-0.3172, +0.0100]
     — not quite significant at 95%, but the interval's upper edge sits right
     at zero). The global term made BCL2's correlation with real GDSC
     response *worse*, not better as predicted.
   - MET: rho **-0.2326 → -0.1395** (Δ = **+0.0931**, 95% CI [-0.0043, +0.1840]
     — also not significant, again right on the boundary). Since MET's
     established "good" direction is more-negative (oncogene-addiction:
     high score should track high sensitivity), this is also a move in the
     *wrong* direction, i.e. also a degradation — same conclusion as BCL2,
     just via a different sign.

**What this means:** I predicted a lineage-confound rescue (numbers look
better, for the wrong reason). What actually happened is closer to the
global term acting as noise-to-mild-antagonist against the existing 2-layer
score for both testable genes — not a confound that flatters the number,
just a drag on it. Lineage explaining only 4% of the metric's variance is
the likely reason: averaging across 225 metabolites spanning unrelated
pathways (amino acids, TCA cycle, lipids, nucleotides) washes out whatever
lineage-coherent signal exists in any one pathway, leaving a summary scalar
that's closer to noise than to a usable lineage proxy — let alone a gene
proxy. Both interpretations converge on the same recommendation as before:
**do not fold this into `core_score`.** Not because it's a masked lineage
confound (it measurably isn't, at 4.2%), but because it isn't carrying
enough signal of any kind to justify the coverage/confound cost or the lost
explainability (Stage 7 still can't write an honest reason string for it).
Being wrong about *why* it doesn't work is itself worth keeping in the
methods section — the a priori concern (lineage confound) and the actual
failure mode (near-noise after averaging) are different findings, and only
measuring it surfaced the real one.
