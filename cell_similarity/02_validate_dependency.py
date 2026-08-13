"""
cell_similarity/02_validate_dependency.py
-----------------------------------------
THE TEST
--------
Take pairs the graph calls close and ask whether their dependency profiles
correlate more than random pairs do. This is a pair-to-pair comparison, not a
prediction of a value, so it tolerates the messiness in the dependency data —
noise attenuates both sides equally and the contrast survives.

Dependency = DepMap 24Q4 CRISPR gene effect, restricted to SELECTIVE genes
(sd > 0.20, not pan-essential). Pan-essential genes sit near -1 in every line;
including them would make every pair of lines correlate highly for reasons that
have nothing to do with the graph, and the test would pass no matter what.

THE CONFOUND THAT DECIDES WHETHER ANY OF THIS IS WORTH ANYTHING
----------------------------------------------------------------
Cell-line similarity is dominated by lineage, and same-lineage lines share
dependencies for reasons a tissue label already captures. So "close pairs
correlate more than random pairs" is nearly guaranteed to pass, and passing it
would tell you almost nothing: a lookup table of tissue labels would pass it too.

Every comparison here is therefore run twice:

  POOLED          close pairs vs all pairs. The stated test. Expect it to pass.
  WITHIN LINEAGE  close pairs vs random pairs FROM THE SAME LINEAGE. This is the
                  one that matters. It asks whether the graph resolves anything
                  once tissue is held fixed — i.e. whether it beats the free
                  baseline of "same tissue".

An axis that passes pooled and fails within-lineage is re-deriving the tissue
label. That is not useless, but it is not worth a graph, and it must not be
sold as more than it is.
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, spearmanr

import common as C

TOP_FRACS = (0.001, 0.005, 0.01, 0.05)
res = {}

C.banner("WEEK 2, TEST 1 -- DO SIMILAR LINES BEHAVE SIMILARLY?")

# ------------------------------------------------------------ dependency
print("\n  loading dependency profiles")
W = C.load_dependency(selective_only=True, verbose=True)
lin = C.load_lineage()

# Dependency-profile correlation between lines, computed once.
Wz = W.sub(W.mean(axis=1), axis=0)
M = Wz.values
M = np.nan_to_num(M, nan=0.0)
norm = np.linalg.norm(M, axis=1, keepdims=True)
norm[norm == 0] = 1.0
Mn = M / norm
DEP = pd.DataFrame(Mn @ Mn.T, index=W.index, columns=W.index)
np.fill_diagonal(DEP.values, np.nan)
print(f"    dependency profile correlation: {DEP.shape[0]:,} lines")
res["dependency"] = {"n_lines": int(DEP.shape[0]), "n_selective_genes": int(W.shape[1])}

rng = np.random.default_rng(C.SEED)


def pair_vectors(S, dep, lines):
    """Aligned similarity and dependency-correlation vectors over unique pairs,
    plus a same-lineage mask and the lineage of each pair."""
    S = S.loc[lines, lines]
    D = dep.loc[lines, lines]
    n = len(lines)
    iu, ju = np.triu_indices(n, k=1)
    s = S.values[iu, ju]
    d = D.values[iu, ju]
    lv = lin.reindex(lines).lineage.values
    same = lv[iu] == lv[ju]
    ok = ~np.isnan(s) & ~np.isnan(d)
    return s[ok], d[ok], same[ok], lv[iu][ok], lv[ju][ok]


# ------------------------------------------------------------ 2A
C.banner("2A -- POOLED:  close pairs vs all pairs   (the stated test)")
print("  dependency-profile correlation, median")
print()
print(f"  {'axis':<8} {'lines':>6} {'pairs':>10} " +
      "".join(f"{'top ' + f'{f:.1%}':>11}" for f in TOP_FRACS) +
      f"{'ALL pairs':>12} {'rho(sim,dep)':>14}")
res["pooled"] = {}
store = {}
for axis in C.available_axes():
    S = C.load_axis_similarity(axis)
    lines = sorted(set(S.index) & set(DEP.index) & set(lin.index))
    if len(lines) < 100:
        print(f"  {axis:<8} too few shared lines ({len(lines)})")
        continue
    s, d, same, la, lb = pair_vectors(S, DEP, lines)
    store[axis] = (s, d, same, la, lb, lines)
    order = np.argsort(-s)
    cells = []
    rec = {"n_lines": len(lines), "n_pairs": int(len(s)),
           "all_pairs_median": float(np.median(d))}
    for f in TOP_FRACS:
        k = max(10, int(len(s) * f))
        m = float(np.median(d[order[:k]]))
        cells.append(m)
        rec[f"top_{f}"] = m
    rho = spearmanr(s, d).statistic
    rec["spearman_sim_vs_dep"] = float(rho)
    print(f"  {axis:<8} {len(lines):>6,} {len(s):>10,} " +
          "".join(f"{c:>11.4f}" for c in cells) +
          f"{np.median(d):>12.4f} {rho:>14.4f}")
    res["pooled"][axis] = rec
print()
print("  Every axis will beat the all-pairs baseline here. That is expected and")
print("  it is not evidence of much -- see 2B.")

# ------------------------------------------------------------ 2B
C.banner("2B -- WITHIN LINEAGE:  the test that actually matters")
print("  Close pairs vs random pairs DRAWN FROM THE SAME LINEAGE. If an axis")
print("  loses its advantage here, it was re-deriving the tissue label.")
print()
print(f"  {'axis':<8} {'same-lin pairs':>15} {'top 1% (same-lin)':>19} "
      f"{'all same-lin':>14} {'delta':>9} {'MWU p':>11} {'rho':>9}")
res["within_lineage"] = {}
for axis, (s, d, same, la, lb, lines) in store.items():
    ss, dd = s[same], d[same]
    if len(ss) < 200:
        print(f"  {axis:<8} too few same-lineage pairs ({len(ss)})")
        continue
    order = np.argsort(-ss)
    k = max(10, int(len(ss) * 0.01))
    top = dd[order[:k]]
    rest = dd[order[k:]]
    p = mannwhitneyu(top, rest, alternative="greater").pvalue
    rho = spearmanr(ss, dd).statistic
    print(f"  {axis:<8} {len(ss):>15,} {np.median(top):>19.4f} "
          f"{np.median(dd):>14.4f} {np.median(top)-np.median(dd):>+9.4f} "
          f"{p:>11.3g} {rho:>9.4f}")
    res["within_lineage"][axis] = {
        "n_same_lineage_pairs": int(len(ss)),
        "top1pct_median": float(np.median(top)),
        "all_same_lineage_median": float(np.median(dd)),
        "delta": float(np.median(top) - np.median(dd)),
        "mwu_p": float(p), "spearman_sim_vs_dep": float(rho)}

# ------------------------------------------------------------ 2C
C.banner("2C -- HOW MUCH IS THE GRAPH, HOW MUCH IS THE TISSUE LABEL?")
print("  Three baselines on the same pairs, same axis, same lines:")
print("    random pair                  -- no information")
print("    same lineage, random within  -- what a free tissue lookup gives you")
print("    top 1% by similarity         -- what the graph gives you")
print()
print(f"  {'axis':<8} {'random':>10} {'same lineage':>14} {'graph top 1%':>14} "
      f"{'graph - tissue':>16}")
res["baselines"] = {}
for axis, (s, d, same, la, lb, lines) in store.items():
    order = np.argsort(-s)
    k = max(10, int(len(s) * 0.01))
    rnd = float(np.median(d))
    tis = float(np.median(d[same])) if same.sum() > 50 else np.nan
    gph = float(np.median(d[order[:k]]))
    print(f"  {axis:<8} {rnd:>10.4f} {tis:>14.4f} {gph:>14.4f} "
          f"{gph - tis:>+16.4f}")
    res["baselines"][axis] = {"random": rnd, "same_lineage": tis,
                              "graph_top1pct": gph, "graph_minus_tissue": gph - tis}
print()
print("  'graph - tissue' is the honest headline for this test. It is what the")
print("  similarity graph adds over knowing the tissue, which you already know.")

# ------------------------------------------------------------ 2D
C.banner("2D -- IS THE TOP-1% ADVANTAGE JUST MORE SAME-LINEAGE PAIRS?")
print(f"  {'axis':<8} {'same-lin share, top 1%':>24} {'same-lin share, all':>21}")
res["top_pair_composition"] = {}
for axis, (s, d, same, la, lb, lines) in store.items():
    order = np.argsort(-s)
    k = max(10, int(len(s) * 0.01))
    a = float(same[order[:k]].mean())
    b = float(same.mean())
    print(f"  {axis:<8} {a:>23.1%} {b:>20.1%}")
    res["top_pair_composition"][axis] = {"top1pct_same_lineage": a,
                                         "all_same_lineage": b}
print()
print("  A top-1% that is overwhelmingly same-lineage while the panel is not")
print("  explains a pooled result entirely, and is why 2B is the real test.")

# ------------------------------------------------------------ verdict
C.banner("VERDICT")
res["verdict"] = {}
for axis in store:
    wl = res["within_lineage"].get(axis)
    bl = res["baselines"].get(axis, {})
    if not wl:
        continue
    beats_tissue = bl.get("graph_minus_tissue", 0) > 0.01
    within_sig = wl["mwu_p"] < 0.05 and wl["delta"] > 0.01
    if within_sig and beats_tissue:
        v = "adds_information_beyond_tissue"
    elif within_sig:
        v = "marginal__within_lineage_signal_but_barely_beats_tissue"
    else:
        v = "re_derives_the_tissue_label__no_added_resolution"
    print(f"  {axis:<8} within-lineage delta {wl['delta']:+.4f} "
          f"(p {wl['mwu_p']:.3g})   over tissue "
          f"{bl.get('graph_minus_tissue', float('nan')):+.4f}   -> {v}")
    res["verdict"][axis] = v

out = C.OUT / "02_validate_dependency_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
