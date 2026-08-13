# TRANSCRIPTOMICS_VALIDATION_REPAIR — what the tests measure, and what they should

Parts B and C of the transcriptomics brief. All numbers `VERIFIED` (computed here) unless
tagged otherwise. Scripts: `docs/audit_scripts/b1_null.py`, `b2_b3_recovery.py`, `c1_chain.py`.

Reference gene for the runs that need a single gene: **TSPAN6**, 1,580 lines, 914 of them
with ≥2 sources. Standing rule 3 governs everything here: *a test that cannot fail is not
evidence*, and every replacement below can fail.

---

## Summary

| Test | Currently | After repair | Verdict on the current test |
|---|---|---|---|
| **B1** null calibration | FPR 4.64%, "calibrated ✓" | true rate **11.59%** at nominal α=0.05 (**2.32× nominal**) | cannot fail — the shuffle enforces the assumption under test |
| **B2** recovery | 100% sensitivity at 2 and 3 MAD | **85%** at 3 MAD, **65%** at 2 MAD | optimistic at high effect, *pessimistic* at low effect |
| **B3** source-count | median z drop **+1.273** | **−0.043** on a random sample | the entire effect is a selection artefact |
| **C1** significance chain | 115 lines significant | **0** | the p-value is not a p-value |

---

# PART B — The two validation tests

## B1. The null calibration test

### What it currently measures

`validate_null_calibration` (cell 12) permutes values **independently within each
`(source, lineage)`**:

```python
d["value"] = (d.groupby(["source","lineage"], dropna=False)["value"]
                .transform(lambda s: rng.permutation(s.values)))
```

A line that was high in every source receives three independent draws. **The shuffle
destroys the cross-source correlation**, so the false-positive rate measured is the rate
under independence — which is the assumption being tested. Finer strata make it worse,
because the shuffle destroys the correlation more thoroughly within each.

### Measured, as-is — VERIFIED

30 permutations, TSPAN6:

| | Value |
|---|---|
| mean FPR | **0.0464** |
| p95 | 0.0551 |
| max | 0.0557 |
| `calibrated` (≤ 0.07) | **True** |

It passes.

### Why a joint permutation is not sufficient either — VERIFIED

The brief proposes shuffling the cell-line label once per lineage and applying the same
shuffle across every source. Implemented and run:

| | Significant rate |
|---|---|
| observed data (no shuffle) | 0.0791 |
| joint permutation | **0.0785** |

A joint permutation is very nearly a **pure relabelling**: for lines appearing in the same
set of sources it simply renames them, leaving the multiset of per-line z-vectors — and
therefore the score distribution — unchanged. The 0.0006 difference comes only from lines
swapping between different source-coverage patterns. It returns the observed rate, so it
cannot serve as a null.

This is worth stating in the methods rather than quietly substituting something else: the
brief's proposed fix is directionally right (preserve the cross-source pattern) but is
degenerate as a permutation scheme.

### The corrected null used here

**Specification.** Draw each line's per-source z-vector from a multivariate normal with
the **measured correlation matrix R** and zero mean — no signal by construction — while
preserving each line's real source-presence pattern and real `√n_samples` weights. Then run
the notebook's own combination and BH correction.

This has (i) no signal, (ii) the real dependence structure, and (iii) can fail.

```
R =  DepMap  1.000  0.863  0.495
     HPA     0.863  1.000  0.530
     GEO     0.495  0.530  1.000
```

Source-presence patterns preserved across 1,453 lines: DepMap+HPA 478, all three 383,
DepMap only 366, GEO only 145, DepMap+GEO 32, HPA only 28, HPA+GEO 21.

### Before and after — VERIFIED

2,000 simulations each:

| Quantity | Under independence (R = I) | Under measured R | Ratio |
|---|---|---|---|
| **rejection rate at nominal α = 0.05, before BH** | **0.0498** | **0.1159** | **2.32×** |
| FPR after BH at q < 0.05 | 0.0000 | 0.0066 | 158× |

**The headline number: the test operates at α = 0.116 while claiming α = 0.05.**

Two things to note, and the second is the point:

- The as-is null's measured 0.0464 sits within sampling error of the independence
  rejection rate 0.0498. **The as-is test reproduces the independence rate because its
  shuffle makes independence true.** That is the mechanical proof that it cannot fail for
  the reason it exists.
- After BH the absolute FPR is small either way (0.66% vs 0.00%), because BH across ~1,450
  lines is severe. The *relative* inflation is 158×. Report both — the absolute number is
  reassuring and the relative number is the finding.

**Treat the rise as a finding, not a failure.** The design is not broken; its stated
operating characteristic is wrong by a factor of 2.3, and it is wrong in the direction that
matters.

---

## B2. The recovery test

### Two faults, both confirmed

**Fault 1 — planted signal is near-perfectly correlated across sources by construction.**
`validate_recovery` (cell 12) plants `med + eff*sd + rng.normal(0, .15)` into *every*
source. Measured, the noise term `0.15` on the raw value corresponds in z units to:

| Source | 0.15 / MAD |
|---|---|
| depmap_expr | 0.165 |
| geo_expr | 0.147 |
| hpa_rna | 0.143 |

Against the **real** cross-source disagreement, measured as the per-line SD of z across
sources on lines with ≥2 sources (n = 914):

| | Value |
|---|---|
| median | **0.299** |
| IQR | 0.156 – 0.499 |

So planted lines agree roughly **twice as tightly** as real lines do — the best possible
case for a combination rule that rewards agreement.

**Fault 2 — planted lines are appended, not substituted.** `d = pd.concat([d, ...])`.
Lung has **236** real lines; the as-is test inflates it to **256**, and every planted line
sits at `+eff` above the median, which shifts the lineage median up and inflates its MAD.

### Corrected specification

1. **Plant with realistic disagreement.** Per-source effect `eff + N(0, 0.299)` in z units,
   with 0.299 measured from the data as above — the direct manifestation of the A2
   correlation.
2. **Substitute, do not append.** Overwrite `n_planted` randomly chosen existing lines in
   the lineage, so its size and centre are unchanged.

### Detection curve, before and after — VERIFIED

20 planted lines per effect size, lineage `lung`:

| Effect (MAD) | As-is sensitivity | As-is median z | **Corrected sensitivity** | Corrected median z |
|---|---|---|---|---|
| 0.5 | 0.00 | 0.76 | **0.30** | 1.01 |
| 1.0 | 0.00 | 1.29 | **0.40** | 1.33 |
| 2.0 | **1.00** | 2.33 | **0.65** | 2.27 |
| 3.0 | **1.00** | 3.39 | **0.85** | 3.54 |

Lineage size: as-is 256 (inflated from 236), corrected 236.

### What this shows

- **The as-is test is optimistic at high effect.** It reports 100% sensitivity at 2 and 3
  MAD. With realistic disagreement it is **65% and 85%**. The notebook's own reporting
  template (cell 37) states *"planted deviations of 3 MAD are recovered with 100%
  sensitivity"* — that claim is **not supported**; the honest figure is 85%.
- **The as-is test is also *pessimistic* at low effect**, which is the less obvious half.
  Appending 20 lines at `+eff` to a 236-line lineage shifts the reference median up and
  inflates the MAD, suppressing the planted lines' own z. That is exactly the
  "changing its size and centre" fault, and it drags the as-is curve to 0% at 0.5 and 1.0
  MAD where the corrected curve reaches 30% and 40%.
- So the as-is curve is not uniformly biased in one direction; it is **a step function
  where the truth is a gradual curve**. A step function is what makes "3 MAD works, 2 MAD
  does not" sound like a clean operating characteristic. It is not one.

### The honest minimum detectable effect

**There is no effect size at which this method reliably detects a planted deviation.** At
3 MAD — a very large deviation — sensitivity is 85%. Stated for the methods section:

> With realistic cross-source disagreement, sensitivity is 30% at 0.5 MAD, 40% at 1.0 MAD,
> 65% at 2.0 MAD and 85% at 3.0 MAD. Deviations below ~2 MAD should be described as
> suggestive; even at 3 MAD roughly one in seven is missed.

**Free parameter.** The disagreement SD (0.299) is measured, not chosen. If the method
ships, sweep it over the measured IQR {0.156, 0.299, 0.499} and report the curve at each,
so the sensitivity claim carries its own sensitivity analysis.

---

## B3. The source-count test

### Is the sample random or outcome-selected? — **outcome-selected** — VERIFIED

`validate_source_count` (cell 12):

```python
full,_ = score_gene_lineage(df, ...)                       # returns sorted by z_shrunk DESC
three = full[full.n_sources==3].model_id.head(30).tolist()
```

`score_gene_lineage` returns `res.sort_values("z_shrunk", ascending=False)`. So `.head(30)`
takes the **30 highest-scoring** three-source lines.

| | z_shrunk median | range |
|---|---|---|
| the 30 selected | **1.827** | 1.597 – 4.457 |
| all 383 three-source lines | **0.037** | −8.346 – 4.457 |

The sample is the extreme upper tail of the distribution it is meant to represent.

### Effect of hiding two of three sources — VERIFIED

| Sample | n | Median z_shrunk drop | Mean | IQR |
|---|---|---|---|---|
| top-30 (notebook) | 30 | **+1.273** | +1.319 | +1.142 / +1.494 |
| **random-30 (corrected)** | 30 | **−0.043** | +0.005 | −0.383 / +0.313 |

**Absolute overstatement: 1.316 z units.** (A percentage is not quotable here — the
corrected denominator is indistinguishable from zero.)

### What this shows

The notebook's conclusion — *"a well-behaved method should shrink the single-source version
harder"*, and it appears to — is **entirely an artefact of selecting on the outcome**. On a
random sample of three-source lines the effect is **−0.043 with an IQR straddling zero**:
removing two of three sources changes the score by nothing on average.

That is mechanically unsurprising once A4 is in hand. The shrinkage factor `k/(k+1)` moves
from 0.75 to 0.50 when k drops from 3 to 1, but the combined Z itself *rises* when the
redundant sources are removed, and the two effects cancel. Selecting the top 30 breaks the
cancellation because those lines were selected for having high Z with three sources.

### Corrected sampling

Draw the k-source lines **uniformly at random** with a fixed seed, and report the drop with
a paired bootstrap interval over lines. Report the drop for k=3→1 *and* k=2→1 separately;
they are different questions and the k=2 case is the common one (478 of 1,106 lines).

---

# PART C — The significance chain

## C1. The wrong quantity is being tested

### The problem, measured

Stouffer's method combines **standardised test statistics**. The notebook feeds it
`robust_z` — an effect size in MAD units. Under the null this is not N(0,1):

| Quantity fed to `stouffer()` | n | mean | **sd** | skew | kurtosis | **P(\|z\| > 1.96)** |
|---|---|---|---|---|---|---|
| `robust_z` (current) | 2,750 | −0.117 | **1.606** | −4.54 | 78.1 | **0.1055** |
| N(0,1) predicts | | 0 | 1 | 0 | 0 | 0.0500 |

**The tails are 2.1× heavier than normal before any combination happens.** `stouffer()`
then computes `p = 2*(1-norm.cdf(|Z|))`, reading a p-value off a distribution the statistic
does not follow. That p-value is not a p-value, and `false_discovery_control` inherits the
fault.

Note this is a **second, independent** source of anti-conservatism, on top of the 1.365–1.503×
correlation inflation in A4. They compound.

### The fix already exists in the notebook

`empirical_p(x, v)` is computed for every `(line, source)` — a rank-based, distribution-free
p-value — and is used only to build `fisher_p`, which is itself never read. Converting it to
a standardised statistic makes the whole chain valid:

```python
z_std = sign(value - lineage_median) * norm.isf(p_emp / 2)
```

Measured on the converted statistic: sd **0.818**, P(|z| > 1.96) = **0.0079**. Conservative
rather than anti-conservative — because of the resolution floor below.

### Effect on the verdicts — VERIFIED

TSPAN6, 1,580 lines, current chain vs empirical-p chain:

| | Value |
|---|---|
| lines changing verdict | **232 (14.7%)** |
| significant at q < 0.05 | **115 → 0** |

Direction of change — **uniformly downward**, no line moves up in strength:

| From | To | n |
|---|---|---|
| Low (weak) | Typical | 62 |
| LOW | Low (weak) | 52 |
| High (Weak) | Typical | 44 |
| HIGH | High (Weak) | 18 |
| LOW (Single Source) | Typical | 17 |
| HIGH (Single Source) | Typical | 13 |
| Conflicting Sources | Typical | 9 |
| Conflicting Sources | Low (weak) | 6 |
| LOW (Small Lineage) | Low (weak) | 4 |
| LOW (Small Lineage) | Typical | 3 |
| Typical | High (Weak) | 3 |
| Conflicting Sources | High (Weak) | 1 |

**Every significant call for this gene disappears under a valid chain.** That is not a
reason to prefer the invalid one; it is the measurement of how much of the current output
is an artefact of the distributional assumption.

## C1c. The resolution floor — a real limit that belongs in the methods

`empirical_p` is two-tailed: `p = 2 · min(hi, lo)` with `hi = ((v >= x).sum()+1)/(n+1)`.
Its **minimum attainable value is `2/(n+1)`** for a lineage of n peers.

| Lineage size n | p_min | Can reach p < 0.05? |
|---|---|---|
| 3 | 0.5000 | no |
| **15** (`MIN_PEERS_FOR_LINEAGE`) | **0.1250** | **no — 2.5× the threshold** |
| 20 | 0.0952 | no |
| 39 | 0.0500 | no |
| **40** | 0.0488 | **yes, marginally** |
| 100 | 0.0198 | yes |
| 200 | 0.0100 | yes |

**The smallest lineage in which any line can reach p < 0.05 at all is n = 39** — before any
BH correction, and before combining sources.

### How much of the panel is below that floor — VERIFIED

DepMap coverage, 31 lineages, 1,479 lines:

| Lineage | n | | Lineage | n |
|---|---|---|---|---|
| lung | 210 | | urinary_tract | 37 |
| blood | 104 | | peripheral_nervous_system | 35 |
| central_nervous_system | 86 | | esophagus | 32 |
| skin | 86 | | plasma_cell | 30 |
| lymphocyte | 84 | | liver | 24 |
| colorectal | 73 | | cervix | 19 |
| (no lineage recorded) | 67 | | thyroid | 18 |
| ovary | 65 | | eye | 17 |
| breast | 63 | | prostate | 12 |
| soft_tissue | 60 | | embryo | 3 |
| upper_aerodigestive | 57 | | unknown | 3 |
| pancreas | 53 | | epidermoid_carcinoma | 1 |
| gastric / bile_duct | 41 / 41 | | adrenal_cortex | 1 |
| bone / uterus | 40 / 40 | | | |
| kidney | 39 | | | |
| fibroblast | 38 | | | |

**14 of 31 lineages are below n = 39. They contain 270 of 1,479 cell lines (18.3%).**
No line in any of them can reach significance through the empirical route, at any effect
size, ever.

**This is a real limit on the design and belongs in the methods, not in a workaround.**
Stated for the write-up:

> The rank-based empirical p-value has a resolution floor of 2/(n+1) set by lineage size.
> For 14 of 31 lineages (18.3% of cell lines) that floor exceeds 0.05, so no line in those
> lineages can be called significant regardless of effect size. Results for those lineages
> are reported as effect sizes with abstention on significance.

Do not paper over it by switching those lineages back to the parametric route — that is the
route C1 has just shown to be invalid.

## C2. The two unused routes — VERIFIED

| Quantity | Computed for | Used by |
|---|---|---|
| `fisher_p` | 1,580 lines (all) | **nothing** — absent from `verdict()`, `q_value`, and the sort key |
| `p_emp` | every (line, source) | only to build `fisher_p`, which is itself unused |
| `lineage_silent_any` | every line | **nothing** — never read |

`corr(fisher_p, p_combined) = 0.826` — they agree substantially but are not redundant.

**Why it was discarded:** `UNVERIFIED — could not locate`. No comment, markdown cell or
commit message explains it. The most likely reading is that it was built as a cross-check,
found to agree, and left in place.

**Recommendation: use it, do not remove it.** Fisher's method on the empirical p-values is
the *more assumption-light* of the two routes — it needs no normality of the effect size,
only valid per-source p-values, which `empirical_p` provides by construction. A reviewer
will ask why the assumption-light quantity was the one thrown away, and there is no good
answer. Concretely:

- Promote the empirical-p route to primary (C1), with Fisher or Stouffer-of-standardised-
  empirical-p as the combiner.
- Keep both and report agreement; where they disagree, that disagreement is itself
  informative and cheaper to explain than a single unexplained number.
- Note that Fisher's method **also assumes independence** and needs the same R correction
  (Brown's method, 1975 — the Fisher analogue of Strube/Hartung, `PROVEN`). Switching
  routes does not escape Part A.
- Delete `lineage_silent_any`, or wire it in. A computed-and-unread column is a maintenance
  trap.

## C3. Multiple testing across lines

### Is the dependence structure what BH assumes? — no — VERIFIED

Lines within a lineage are scored against a **shared** median and MAD, so their z's are
mechanically constrained — the deviations about the median sum to approximately zero.
Measured mean z per lineage (should be ≈0 by construction):

| | Value |
|---|---|
| median across lineages | **−0.093** |
| range | −0.718 to +0.306 |

This is **negative dependence**, not the positive-regression-dependence (PRDS) structure
Benjamini–Hochberg assumes.

### Does it matter in practice? — **no, and this is the one place the news is good**

BH under negative dependence is generally **conservative**, not anti-conservative. It
controls FDR at or below the nominal level rather than above it. So this route does not
inflate the false-discovery rate.

Given the B1 measurement — the combination operates at α = 0.116 while claiming 0.05 —
**the correlation inflation is the binding problem and this is not.** Fixing C3 first would
be optimising the wrong term.

**Recommendation: leave BH in place, and state the dependence structure in the methods.**
If a belt-and-braces guarantee is later wanted, Benjamini–Yekutieli controls FDR under
*arbitrary* dependence at the cost of a `log(m) ≈ 7.4` factor at m = 1,580 — which, given
that a valid chain already yields zero significant calls for TSPAN6 (C1), would be
purely decorative. Do not adopt it.
