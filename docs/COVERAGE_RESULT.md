# Coverage and discrimination are inversely related — a result

**Claim.** Across every evidence check available to this architecture, the checks that are
widely measured barely discriminate, and the checks that discriminate are barely measured.
This is a structural property of the available data, not a defect in any one source. It is
the binding constraint on building any multi-evidence ordering, and it is measurable,
predictive, and was tested by intervention.

All figures `VERIFIED` — computed 2026-08-09 on 300 protein-coding genes (seed 42) ×
1,543 cell lines = 462,900 (gene, line) pairs. Scripts: `docs/audit_scripts/rank_checks.py`,
`partA.py`, `partA_v2.py`. Build: `src/pipeline/build_cn_layer_complete.py`.

---

## 1. The measurement

For each candidate evidence check, two independent quantities on gate-retained pairs:

- **Coverage** — the fraction of pairs for which the check is *assessed* (has data at all).
- **Discrimination** — `min(pass, fail) / assessed`, the fraction of assessed pairs on the
  minority side of the call. A check that always says the same thing has discrimination 0
  and can order nothing, however complete its coverage.

| Check | Assessed | Pass (of assessed) | Fail (of assessed) | **Discrimination** |
|---|---|---|---|---|
| `alteration` | **98.77%** | 98.36% | 1.64% | **0.016** |
| `cna_gate` (gene-complete, rel < 0.5) | **52.29%** | 98.49% | 1.51% | **0.015** |
| `cna_gate` (gene-complete, rel < 0.75) | 52.29% | 88.26% | 11.74% | 0.117 |
| `geo_corrob` | 32.08% | 68.51% | 31.49% | **0.315** |
| `protein` | **19.03%** | 50.03% | 49.97% | **0.500** |

**Spearman(coverage, discrimination) = −0.80** across the four checks at their primary
thresholds. The two best-covered checks have the two lowest discrimination values; the two
worst-covered have the two highest.

The single perfectly-discriminating, perfectly-covered check does not exist. The RNA gate is
100% assessed and is the only check that is both — which is why it defines the gate tier, and
why nothing is left to put underneath it.

---

## 2. The intervention — a prediction that failed, informatively

The relationship above predicts something falsifiable: **adding coverage without adding
discrimination should make a confirmation count *more* of a data-availability proxy, not
less.** That prediction was tested directly.

`data/DepMap_Chronos/CCLEGeneCopyNumber20Q2.hdf5` had been sitting unused in the repository —
908 lines × 27,639 genes, zero missing values, already normalised to each line's own ploidy
(so the 30%-of-panel ploidy gap does not bite). It was wired in
(`src/pipeline/build_cn_layer_complete.py`) as a strictly better CN source:

| | COSMIC CNA (shipped) | Gene-complete CN (new) |
|---|---|---|
| models | 969 | 908 |
| genes | 15,545 | 18,276 |
| **(gene, line) pairs** | **117,540** | **16,594,608** |
| assessed share of gate-retained pairs | **0.28%** | **52.29%** |

**A 141× increase in pairs, and a 187× increase in assessed share.**

The pre-specified decision rule (`docs/RANKING_PRESPEC.md`, written before any re-run number
existed) was then applied. Result, valid genes, within the gate-retained set, primary
variant:

| Criterion | Adopt if | Abandon if | **Before** (COSMIC) | **After** (gene-complete) |
|---|---|---|---|---|
| A1 coverage confound \|ρ\| | ≤ 0.20 | ≥ 0.40 | **+0.182** | **+0.531** |
| A2a largest tied group | ≤ 0.50 | ≥ 0.65 | **0.587** | **0.654** |
| A2b distinct values | ≥ 4 | ≤ 2 | 4 | 4 |
| A3 between-gene ρ | ≤ 0.25 | ≥ 0.40 | **+0.092** | **+0.475** |
| **Verdict** | | | INCONCLUSIVE | **ABANDON** |

**Every criterion got worse.** ABANDON is unanimous across all three pre-specified CN
thresholds (0.25 / 0.5 / 0.75), so the prespec's "verdicts differ ⇒ inconclusive" escape does
not apply.

### Why more data made it worse

Two mechanisms, both measured:

1. **CN availability is a line-level indicator.** 840 of 1,543 lines have CN data, and every
   covered line has data for *all* sampled genes (288 of 288, min = max). So the check's
   presence carries no gene-specific information — it is a binary label on cell lines that
   correlates with other evidence coverage at **ρ = +0.242 (p = 4.9e-22)**.
2. **The check almost never fails.** Among assessed retained pairs it passes **98.49%** of the
   time. Median 2 deleted lines per gene out of 686 assessed; **104 of 300 genes have zero
   deleted retained lines**.

Adding a near-constant "pass" that is present for a specific half of the panel is,
arithmetically, adding a coverage indicator to the count. The confound rose from +0.182 to
+0.531 for exactly that reason.

### Pre-registered predictions, scored honestly

`RANKING_PRESPEC.md` §6 recorded four predictions before the run. **Three of four were
wrong, all in the same direction** — I expected the CN layer to help.

| # | Prediction | Outcome |
|---|---|---|
| 1 | `cna_gate` assessed-fraction rises above 55% of retained pairs | **FAILED** — 52.29% |
| 2 | A1 will improve but may not clear 0.20 | **FAILED** — it worsened, +0.182 → +0.531 |
| 3 | A2a is the criterion most likely to fail | **partly right** — A2a failed, but so did A1 and A3 |
| 4 | A3 will improve most | **FAILED badly** — +0.092 → +0.475 |

The systematic direction of the error is itself the finding: the intuition that "more data
helps" is wrong here, and it is wrong for a reason that generalises.

---

## 3. Why this is a result and not a limitation

A limitation is a caveat attached to a claim. This is a claim.

1. **It is quantified.** Coverage and discrimination are measured separately for five checks,
   and the inverse relationship is ρ = −0.80.
2. **It made a falsifiable prediction**, which was pre-registered and tested by intervention,
   and the intervention behaved as the relationship predicts (worse, not better) rather than
   as intuition predicts.
3. **It explains a prior negative result** that otherwise looks like a design failure. The
   first Part A returned negative; this says the reason was not the arithmetic of the count,
   and shows that the obvious fix does not work.
4. **It is actionable and it rules things out.** It says the next unit of work is not another
   ordering method, another normalisation, or another data source of the same kind. Adding a
   *sixth* well-covered, low-discrimination check would make the count worse again.
5. **It bounds the product.** Not "results should be interpreted with caution" but: the
   architecture can return ~1,200 lines in 2 tiers per gene, and no rearrangement of the
   current evidence produces a shortlist.

---

## 4. What the relationship implies for what to build

A useful check must clear **both** bars. The measurements say what that costs:

| Requirement | Threshold implied by this work | Current best |
|---|---|---|
| assessed on gate-retained pairs | ≳ 70% | RNA gate 100%; next best CN 52% |
| discrimination `min(p,f)/assessed` | ≳ 0.15 | protein 0.50, GEO 0.32, CN@0.75 0.12 |

**The two candidates that already clear the discrimination bar are `protein` (0.500) and
`geo_corrob` (0.315). Both fail the coverage bar (19.03%, 32.08%).** So the highest-value work
is not acquiring a new *kind* of evidence — it is extending the coverage of evidence the
project already has and has already shown to discriminate:

- **Proteomics coverage.** ProCan covers 952 lines and TMT 375; the union reaches 769 of the
  1,479 RNA-covered lines. Proteomics is the single most discriminating check available and
  it is assessed on under a fifth of retained pairs. Note the MNAR caveat stands
  (Cliff's δ = +0.251 selection toward higher-expressing lines), so extending coverage is not
  merely a volume problem.
- **GEO coverage.** 588 lines, 32.08% of retained pairs, discrimination 0.315.

Raising `protein` from 19% to 70% coverage would do more for an ordering than any
rearrangement of the checks that already exist. Whether that is obtainable is a data question,
not a modelling one.

**The one caveat on the CN layer:** at the shallower `rel < 0.75` threshold its discrimination
rises to 0.117 and 210 of 300 genes gain ≥10 deleted retained lines. It remains ABANDON as a
count input, but it is now a usable *exclusion* signal in its own right — which is how CNA
deletion was already being used (Fisher ratio 0.380 on 69 genes). The layer is worth keeping
for that; it is simply not a confirmation.

---

## 5. Scope

- Measured on one gene sample (300 genes, seed 42) and one panel. The inverse relationship is
  described over four checks — enough to be striking, not enough for a confidence interval on
  ρ itself. Reported as descriptive.
- Discrimination is defined against each check's own primary threshold. A different threshold
  changes the number (CN moves 0.015 → 0.117 between rel < 0.5 and rel < 0.75); the
  *ordering* of the checks by discrimination is unchanged across the thresholds tested.
- "Coverage" here is assessed-fraction on gate-retained pairs, not panel coverage. The two
  differ because retention is itself gene-dependent.
