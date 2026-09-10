"""
combined_score_b0_lineage.py
----------------------------
B0 follow-up: is the q = 0.05 marginal carried by one or two lineages?

If it is, both the marginal and any interaction built on top of it are fragile,
and that is worth knowing BEFORE Stage 1 rather than during it.

Two decompositions, because they answer different questions:

  A. LEAVE-ONE-LINEAGE-OUT. Recompute the pooled marginal with each lineage
     removed in turn. If dropping one lineage collapses the effect, that lineage
     was carrying it. This is the fragility test.

  B. WITHIN-LINEAGE. Compute the marginal separately inside each lineage with
     enough lines. This asks whether the effect is a general property or a
     between-lineage artefact -- low-expression lines could be depleted of
     sensitivity simply because whole lineages differ in both.

Cluster bootstrap over CELL LINES throughout, per the standing rule.

    python src/pipeline/combined_score_b0_lineage.py
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(Path(__file__).resolve().parent))
import multi_gene_kleene as K  # noqa: E402

GDSC = ROOT / "validation" / "prepared" / "gdsc_scored_ready.parquet"
OUT = ROOT / "src" / "pipeline" / "outputs"
Q = 0.05
MIN_LINES_PER_PAIR = 100
RNG = np.random.default_rng(42)
res = {"q": Q}

print("=" * 78)
print(f"B0 FOLLOW-UP -- IS THE q={Q} MARGINAL CARRIED BY ONE OR TWO LINEAGES?")
print("=" * 78)

E = K.expression()
fe = K.gene_validity()
lin = K.lineage()

gd = pd.read_parquet(GDSC, columns=["model_id", "drug_id", "target_ensg",
                                    "target_symbol", "AUC"])
gd["model_id"] = gd.model_id.astype(str).str.lower()
gd["target_ensg"] = (gd.target_ensg.astype("string").str.split(".").str[0]
                     .str.lower())
expr_lines = set(E.index)
gd = gd[gd.model_id.isin(expr_lines) & gd.target_ensg.isin(E.columns)]

# rebuild the same 294-pair set
pairs = []
for (g, did), sub in gd.groupby(["target_ensg", "drug_id"]):
    if fe.get(g, 0) < K.SILENT_FRAC:
        continue
    lines = sorted(set(sub.model_id))
    if len(lines) < MIN_LINES_PER_PAIR:
        continue
    x = E[g].reindex(lines).values
    auc = sub.set_index("model_id")["AUC"].reindex(lines).values
    ok = ~np.isnan(x) & ~np.isnan(auc)
    if ok.sum() < MIN_LINES_PER_PAIR:
        continue
    ls = np.array(lines)[ok]
    pairs.append((ls, x[ok], auc[ok]))
print(f"\n  pairs rebuilt: {len(pairs)}")

all_lines = sorted({m for ls, _, _ in pairs for m in ls})
line_ix = {m: i for i, m in enumerate(all_lines)}
lv = np.array([str(lin.lineage.get(m, "unknown")) for m in all_lines])
print(f"  lines: {len(all_lines):,}   lineages: {len(set(lv))}")


def d1_median(keep_mask=None, weights=None, min_lines=None):
    """Median per-pair decile-1 ratio at q=Q.
    keep_mask: boolean over all_lines, restricting which lines are used."""
    ml = MIN_LINES_PER_PAIR if min_lines is None else min_lines
    vals = []
    for ls, x, auc in pairs:
        idx = np.array([line_ix[m] for m in ls])
        if keep_mask is not None:
            sel = keep_mask[idx]
            if sel.sum() < ml:
                continue
            x2, a2, idx = x[sel], auc[sel], idx[sel]
        else:
            x2, a2 = x, auc
        if weights is not None:
            w = weights[idx]
            if w.sum() < ml:
                continue
            x2 = np.repeat(x2, w)
            a2 = np.repeat(a2, w)
        n = len(x2)
        k = max(10, int(round(n * Q)))
        if n - k < 10:
            continue
        cut = np.partition(a2, k - 1)[k - 1]
        pos = a2 <= cut
        base = pos.mean()
        if base <= 0:
            continue
        order = np.argsort(x2, kind="mergesort")
        d1 = order[:max(5, n // 10)]
        vals.append(pos[d1].mean() / base)
    return (float(np.median(vals)) if vals else np.nan), len(vals)


def boot_ci(keep_mask=None, n_boot=400, min_lines=None):
    """Cluster bootstrap over cell lines.

    When a keep_mask restricts the analysis to a subset (one lineage), the
    RESAMPLING must be restricted to that subset too. Resampling all 711 lines
    and then masking draws only ~n_subset lines by chance, most pairs then fail
    the per-pair minimum, and the statistic is silently computed on a handful of
    survivors. On the first attempt that produced an interval which did not
    contain its own point estimate -- lung 1.0143 with CI [0.0000, 0.9641].
    """
    b = []
    n = len(all_lines)
    pool = np.where(keep_mask)[0] if keep_mask is not None else np.arange(n)
    for _ in range(n_boot):
        draw = RNG.choice(pool, len(pool), replace=True)
        w = np.bincount(draw, minlength=n)
        v, _ = d1_median(keep_mask, w, min_lines=min_lines)
        if np.isfinite(v):
            b.append(v)
    return (np.percentile(b, [2.5, 97.5]) if len(b) > 20 else (np.nan, np.nan))


base_val, base_n = d1_median()
lo, hi = boot_ci()
print(f"\n  POOLED (all lineages): {base_val:.4f}  95% CI [{lo:.4f}, {hi:.4f}]  "
      f"({base_n} pairs)")
res["pooled"] = {"d1_ratio": base_val, "ci": [float(lo), float(hi)],
                 "n_pairs": base_n}

# ---------------------------------------------------------------- A
print("\nA -- LEAVE ONE LINEAGE OUT")
print("   If removing one lineage collapses the effect, that lineage carried it.")
print()
counts = pd.Series(lv).value_counts()
big = [l for l in counts.index if counts[l] >= 25]
print(f"  {'lineage dropped':<26} {'lines':>6} {'d1 ratio':>10} {'shift':>9} "
      f"{'95% CI':>22}")
rows = []
for l in big:
    mask = lv != l
    v, npair = d1_median(mask)
    l2, h2 = boot_ci(mask, n_boot=250)
    rows.append({"lineage": l, "n_lines": int(counts[l]), "d1_ratio": v,
                 "shift": v - base_val, "ci": [float(l2), float(h2)],
                 "still_excludes_1": bool(h2 < 1.0)})
rows.sort(key=lambda r: -abs(r["shift"]))
for r in rows[:10]:
    flag = "" if r["still_excludes_1"] else "   <- effect no longer excludes 1.0"
    ci = "[%.4f, %.4f]" % (r["ci"][0], r["ci"][1])
    print(f"  {r['lineage']:<26} {r['n_lines']:>6} {r['d1_ratio']:>10.4f} "
          f"{r['shift']:>+9.4f} {ci:>22}{flag}")
res["leave_one_out"] = rows
n_breaks = sum(1 for r in rows if not r["still_excludes_1"])
print(f"\n  lineages whose removal makes the interval span 1.0: {n_breaks} of "
      f"{len(rows)}")

# ---------------------------------------------------------------- B
WITHIN_MIN_LINES = 40
print(f"\nB -- WITHIN LINEAGE (lineages with >= {WITHIN_MIN_LINES} lines)")
print("   Is the effect general, or a between-lineage artefact?")
print(f"   Per-pair minimum relaxed to {WITHIN_MIN_LINES} lines -- the pooled")
print(f"   minimum of {MIN_LINES_PER_PAIR} is unreachable inside one lineage.")
print()
print(f"  {'lineage':<26} {'lines':>6} {'pairs':>6} {'d1 ratio':>10} {'95% CI':>22}")
wrows = []
for l in [x for x in counts.index if counts[x] >= WITHIN_MIN_LINES]:
    mask = lv == l
    v, npair = d1_median(mask, min_lines=WITHIN_MIN_LINES)
    if not np.isfinite(v) or npair < 20:
        continue
    l2, h2 = boot_ci(mask, n_boot=250, min_lines=WITHIN_MIN_LINES)
    wrows.append({"lineage": l, "n_lines": int(counts[l]), "n_pairs": npair,
                  "d1_ratio": v, "ci": [float(l2), float(h2)],
                  "excludes_1": bool(h2 < 1.0)})
wrows.sort(key=lambda r: r["d1_ratio"])

# ESTIMABILITY GUARD. Inside a lineage there are at most ~144 lines, so at
# q = 0.05 a pair carries ~7 positives and the decile-1 cell (14 lines) very
# often contains ZERO of them. The per-pair ratio is then exactly 0, the
# distribution is zero-inflated, and both the median and its bootstrap interval
# become unreliable. The diagnostic is unmissable: an interval that does not
# contain its own point estimate, or a lower bound pinned at exactly 0.
for r in wrows:
    r["ci_contains_point"] = bool(r["ci"][0] <= r["d1_ratio"] <= r["ci"][1])
    r["lower_bound_pinned_at_zero"] = bool(r["ci"][0] <= 1e-9)
    r["estimable"] = bool(r["ci_contains_point"] and not r["lower_bound_pinned_at_zero"])

for r in wrows:
    ci = "[%.4f, %.4f]" % (r["ci"][0], r["ci"][1])
    note = "" if r["estimable"] else "   NOT ESTIMABLE"
    print(f"  {r['lineage']:<26} {r['n_lines']:>6} {r['n_pairs']:>6} "
          f"{r['d1_ratio']:>10.4f} {ci:>22}{note}")
res["within_lineage"] = wrows

n_est = sum(1 for r in wrows if r["estimable"])
n_sig = sum(1 for r in wrows if r["estimable"] and r["excludes_1"])
print()
print(f"  estimable lineages: {n_est} of {len(wrows)}")
if n_est < len(wrows):
    print(f"  {len(wrows) - n_est} lineage(s) produced an interval that does not")
    print(f"  contain its own point estimate, or a lower bound pinned at 0. Those")
    print(f"  are zero-inflation artefacts of ~7 positives per pair inside a")
    print(f"  lineage, NOT results, and are not counted.")
print(f"  of the estimable, independently exclude 1.0: {n_sig}")
res["within_lineage_estimable"] = int(n_est)

# ---------------------------------------------------------------- verdict
print("\n" + "=" * 78)
print("VERDICT")
print("=" * 78)
frag = n_breaks > 0
gen = n_est >= 3 and n_sig >= max(2, n_est // 3)
print(f"  pooled {base_val:.4f} [{lo:.4f}, {hi:.4f}]")
print(f"  removal of any single lineage breaks it : "
      f"{'YES -- ' + ', '.join(r['lineage'] for r in rows if not r['still_excludes_1']) if frag else 'no'}")
print(f"  within-lineage decomposition            : {n_est} of {len(wrows)} "
      f"lineages estimable; {n_sig} of those exclude 1.0")
print()
if frag:
    verdict = "FRAGILE__carried_by_specific_lineages"
    print("  FRAGILE. The marginal depends on particular lineages, so an")
    print("  interaction built on it would inherit that fragility. Stage 1 would")
    print("  need per-lineage reporting from the start, and the pooled number")
    print("  should not be quoted alone.")
elif gen:
    verdict = "GENERAL__survives_leave_one_out_and_appears_within_lineages"
    print("  GENERAL. The marginal survives removal of every single lineage and")
    print("  appears independently inside multiple lineages. It is a property of")
    print("  the relationship, not of the panel composition.")
else:
    verdict = "NOT_FRAGILE__within_lineage_not_estimable_at_this_coverage"
    print("  NOT FRAGILE, and the within-lineage question is UNANSWERABLE here.")
    print("  Leave-one-out is clean: no single lineage carries the effect.")
    print("  But the within-lineage decomposition is not estimable at this")
    print("  coverage -- 711 GDSC-and-expression lines over 28 lineages leaves at")
    print("  most ~144 per lineage, which at q = 0.05 is ~7 positives per pair and")
    print("  a decile-1 cell that is frequently empty. Report the leave-one-out")
    print("  result; do not report per-lineage point estimates.")
print(f"\n  VERDICT: {verdict}")
res["verdict"] = verdict

p = OUT / "combined_score_b0_lineage_results.json"
p.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {p}")
