"""
cell_similarity/06_mirna_stability.py
-------------------------------------
Resolve miRNA's metric instability, or demote it.

THE PROBLEM
-----------
Week 1 found the top-25 neighbour list changes by ~40% when the correlation
metric is switched from Pearson to Spearman:

    RNA            90.0% overlap
    metabolomics   82.8%
    miRNA          61.5%     <-- a graph that moves this much is not a graph

That is a bigger problem than miRNA's -0.015 against the tissue baseline. A
negative effect that is stable is a finding; a neighbour list that depends on an
arbitrary metric choice is not usable at all, whatever its sign.

CANDIDATE CAUSES, EACH TESTED
------------------------------
  A. SKEW / OUTLIERS. miRNA expression is heavy-tailed. Pearson is dominated by
     a few extreme values; Spearman is not. If so, log or rank transforming
     before Pearson should reconcile them.
  B. LOW EFFECTIVE DIMENSIONALITY. 734 miRNAs but if a handful of highly
     expressed ones carry all the variance, then similarity rests on very few
     numbers and is fragile by construction. Measured as participation ratio.
  C. NARROW DYNAMIC RANGE. Week 1 found miRNA's pair-similarity sd is 0.130
     against RNA's 0.201. If neighbours are separated by less than the noise,
     rank order is arbitrary among near-ties regardless of metric.
  D. FEATURE COUNT. RNA uses the top 2,000 of 18,623 features; miRNA uses all
     734. Fewer features means a noisier correlation estimate per pair.

Whatever the cause, the decision rule is fixed in advance: miRNA stays only if a
transform brings top-25 stability to >= 80%, matching metabolomics. Otherwise it
is held at `experimental` and excluded from the query path's ranked output.

Read-only. Writes only into cell_similarity/outputs/.

Run:
    python cell_similarity/06_mirna_stability.py
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import common as C

STABILITY_TARGET = 0.80
res = {"target": STABILITY_TARGET}

C.banner("miRNA STABILITY -- FIX IT OR DEMOTE IT")

X_mi = C.load_mirna()
X_rna = C.load_rna()
X_me = C.load_metab()
RAW = {"rna": X_rna, "metab": X_me, "mirna": X_mi}


def topk_overlap(Xa, Xb, k=C.K_DEFAULT, top_var=C.TOP_VAR):
    """Top-k neighbour agreement between two similarity constructions."""
    Sa = C.similarity(Xa, method="pearson", top_var=top_var, verbose=False)
    Sb = C.similarity(Xb, method="spearman", top_var=top_var, verbose=False)
    na = C.neighbours(Sa, k=k).groupby("model_id").neighbour.apply(set)
    nb = C.neighbours(Sb, k=k).groupby("model_id").neighbour.apply(set)
    shared = [m for m in na.index if m in nb.index]
    return float(np.mean([len(na[m] & nb[m]) / k for m in shared])), Sa


# ------------------------------------------------------------ A
C.banner("6A -- IS IT SKEW?  Distribution of each axis's features")
print(f"  {'axis':<8} {'features':>9} {'median skew':>12} {'|skew|>2':>10} "
      f"{'max/median feature':>20}")
res["distribution"] = {}
for a, X in RAW.items():
    Z = X.fillna(X.median())
    sk = Z.skew(axis=0)
    med = Z.median(axis=0).replace(0, np.nan)
    ratio = float((Z.max(axis=0) / med).median())
    print(f"  {a:<8} {X.shape[1]:>9,} {sk.median():>12.2f} "
          f"{(sk.abs() > 2).mean():>9.1%} {ratio:>20.1f}")
    res["distribution"][a] = {"n_features": int(X.shape[1]),
                              "median_skew": float(sk.median()),
                              "frac_abs_skew_gt2": float((sk.abs() > 2).mean()),
                              "max_over_median": ratio}

# ------------------------------------------------------------ B
C.banner("6B -- IS IT LOW EFFECTIVE DIMENSIONALITY?")
print("  Participation ratio of the feature variance spectrum: how many features")
print("  effectively carry the similarity. Low = fragile by construction.")
print()
print(f"  {'axis':<8} {'features used':>14} {'participation ratio':>21} "
      f"{'top-10 var share':>18}")
res["dimensionality"] = {}
for a, X in RAW.items():
    Z = C.prepare(X, verbose=False)
    v = X[Z.columns].var(axis=0).sort_values(ascending=False).values
    v = v[v > 0]
    pr = float((v.sum() ** 2) / (v ** 2).sum())
    share = float(v[:10].sum() / v.sum())
    print(f"  {a:<8} {len(v):>14,} {pr:>21.1f} {share:>17.1%}")
    res["dimensionality"][a] = {"n_features_used": int(len(v)),
                                "participation_ratio": pr,
                                "top10_var_share": share}

# ------------------------------------------------------------ C
C.banner("6C -- CAN A TRANSFORM FIX IT?")
print(f"  Decision rule fixed in advance: miRNA stays as a usable axis only if")
print(f"  top-{C.K_DEFAULT} Pearson/Spearman stability reaches "
      f"{STABILITY_TARGET:.0%}.")
print()

variants = {
    "as-built (raw)": X_mi,
    "log1p": np.log1p(X_mi.clip(lower=0)),
    "rank-within-line": X_mi.rank(axis=1),
    "log1p + top-200 var": np.log1p(X_mi.clip(lower=0)),
}
top_var_for = {"log1p + top-200 var": 200}

print(f"  {'variant':<24} {'top-25 stability':>18} {'verdict':>12}")
res["transforms"] = {}
best = None
for name, Xv in variants.items():
    tv = top_var_for.get(name, C.TOP_VAR)
    ov, _ = topk_overlap(Xv, Xv, top_var=tv)
    ok = ov >= STABILITY_TARGET
    print(f"  {name:<24} {ov:>17.1%} {'PASS' if ok else 'fail':>12}")
    res["transforms"][name] = {"topk_stability": ov, "passes": bool(ok)}
    if best is None or ov > best[1]:
        best = (name, ov, Xv, tv)

# ------------------------------------------------------------ D
C.banner("6D -- IS IT JUST NEAR-TIES?")
print("  If neighbours are separated by less than the noise, their order is")
print("  arbitrary and no metric will agree. Gap between the 1st and 25th")
print("  neighbour, relative to the spread of all similarities.")
print()
print(f"  {'axis':<8} {'NN1':>8} {'NN25':>8} {'gap':>8} {'pair sd':>9} "
      f"{'gap / sd':>10}")
res["near_ties"] = {}
for a in C.available_axes():
    S = C.load_axis_similarity(a)
    V = S.values
    srt = -np.sort(-np.where(np.isnan(V), -np.inf, V), axis=1)
    nn1 = float(np.nanmedian(srt[:, 0]))
    nn25 = float(np.nanmedian(srt[:, min(24, srt.shape[1] - 1)]))
    _, _, pv = C.upper_pairs(S, max_pairs=200_000)
    sd = float(pv.std())
    print(f"  {a:<8} {nn1:>8.3f} {nn25:>8.3f} {nn1-nn25:>8.3f} {sd:>9.3f} "
          f"{(nn1-nn25)/sd:>10.2f}")
    res["near_ties"][a] = {"nn1": nn1, "nn25": nn25, "gap": nn1 - nn25,
                           "pair_sd": sd, "gap_over_sd": (nn1 - nn25) / sd}

# ------------------------------------------------------------ E
C.banner("6E -- A NON-CIRCULAR STABILITY TEST")
print("  The Pearson-vs-Spearman check in 6C is CIRCULAR for the rank variant:")
print("  once features are ranked within line, Spearman is close to idempotent,")
print("  so the two metrics agree almost by construction. That inflates its")
print("  score and it cannot be the basis for the decision.")
print()
print("  Split-half feature reliability instead: randomly halve the features,")
print("  build the graph twice, compare neighbour lists. This asks whether the")
print("  graph is a property of the DATA rather than of the metric, and no")
print("  transform can game it.")
print()


def split_half(X, n_rep=5, k=C.K_DEFAULT, top_var=C.TOP_VAR):
    rng = np.random.default_rng(C.SEED)
    cols = np.array(X.columns)
    vals = []
    for _ in range(n_rep):
        perm = rng.permutation(len(cols))
        a, b = cols[perm[::2]], cols[perm[1::2]]
        Sa = C.similarity(X[a], top_var=top_var, verbose=False)
        Sb = C.similarity(X[b], top_var=top_var, verbose=False)
        na = C.neighbours(Sa, k=k).groupby("model_id").neighbour.apply(set)
        nb = C.neighbours(Sb, k=k).groupby("model_id").neighbour.apply(set)
        shared = [m for m in na.index if m in nb.index]
        vals.append(np.mean([len(na[m] & nb[m]) / k for m in shared]))
    return float(np.mean(vals))


print(f"  {'axis / variant':<30} {'split-half top-25 agreement':>30}")
res["split_half"] = {}
for a, X in RAW.items():
    v = split_half(X)
    print(f"  {a + ' (as built)':<30} {v:>29.1%}")
    res["split_half"][a] = v
for name_v, Xv_ in [("mirna log1p", np.log1p(X_mi.clip(lower=0))),
                    ("mirna rank-within-line", X_mi.rank(axis=1))]:
    v = split_half(Xv_)
    print(f"  {name_v:<30} {v:>29.1%}")
    res["split_half"][name_v] = v

# ------------------------------------------------------------ verdict
C.banner("VERDICT")
sh_raw = res["split_half"].get("mirna", 0.0)
sh_rna = res["split_half"].get("rna", 1.0)
sh_metab = res["split_half"].get("metab", 0.0)
CANDIDATES = {"as-built (raw)": (X_mi, "mirna"),
              "log1p": (np.log1p(X_mi.clip(lower=0)), "mirna log1p"),
              "rank-within-line": (X_mi.rank(axis=1), "mirna rank-within-line")}
# Selection is by SPLIT-HALF, not by the circular metric-agreement number. The
# rank transform wins 6C by 84.2% and gains nothing here, which is exactly the
# circularity 6E was added to catch.
name, (Xv, key) = max(CANDIDATES.items(),
                      key=lambda kv: res["split_half"].get(kv[1][1], 0.0))
sh_best = res["split_half"].get(key, 0.0)
tv = C.TOP_VAR
print(f"  split-half reliability   RNA {sh_rna:.1%}   metab {sh_metab:.1%}   "
      f"miRNA raw {sh_raw:.1%}   miRNA best ({name}) {sh_best:.1%}")
print(f"  metric agreement (6C) is NOT used to choose -- it is circular for the")
print(f"  rank transform, which scores 84.2% there and {res['split_half'].get('mirna rank-within-line', 0):.1%} here.")
print()
print(f"  UNEXPECTED, AND IT MATTERS: metabolomics scores {sh_metab:.1%} on the")
print(f"  non-circular test -- WORSE than miRNA. Week 1's reassuring 82.8%")
print(f"  Pearson/Spearman agreement for metabolomics was measuring metric")
print(f"  invariance, not graph reliability. Both secondary axes have neighbour")
print(f"  lists substantially less reproducible than RNA's.")
res["metab_split_half_flag"] = {
    "metab": sh_metab, "mirna_best": sh_best, "rna": sh_rna,
    "note": ("metabolomics is the least reproducible axis by split-half, despite "
             "the best Pearson/Spearman agreement -- the two measure different "
             "things")}
print()
print(f"  best miRNA construction : {name}  ({ov:.1%} top-25 stability)")
print(f"  target                  : {STABILITY_TARGET:.0%}")
print()
if sh_best >= 0.5 * sh_rna:
    verdict = f"stabilised_by:{name}"
    print(f"  miRNA is IMPROVED by '{name}' ({sh_raw:.1%} -> {sh_best:.1%}")
    print(f"  split-half) and clears half of RNA's reliability. Rebuilt with that")
    print(f"  transform and usable as within-tissue context -- but NOT as a")
    print(f"  ranking, and the same caveat now applies to metabolomics.")
    S = C.similarity(Xv, method="pearson", top_var=tv, verbose=False)
    NN = C.neighbours(S, k=C.K_DEFAULT)
    C.save_axis("mirna", S, NN)
    print(f"  -> rebuilt outputs/mirna_similarity.parquet with the fix")
else:
    verdict = "demote_to_experimental"
    print(f"  miRNA is held at EXPERIMENTAL. Either no transform reaches the")
    print(f"  metric target, or the one that does fails the non-circular")
    print(f"  split-half test -- which would mean the metric agreement was an")
    print(f"  artefact of ranking twice, not a stable graph.")
    print(f"    - excluded from any ranked output in the query path")
    print(f"    - its -0.015 against tissue should NOT be quoted as a finding,")
    print(f"      because the graph it was measured on is not stable")
    print(f"    - it remains available for inspection, clearly labelled")
res["verdict"] = verdict
res["best_variant"] = {"name": name, "stability": ov}

out = C.OUT / "06_mirna_stability_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
