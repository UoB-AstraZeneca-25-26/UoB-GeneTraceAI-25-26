"""
cell_similarity/01_build_axes.py
--------------------------------
Build the similarity graph — one axis at a time, kept separate.

RNA first, then metabolomics, then miRNA. Each axis gets its own line x line
similarity matrix and its own k-nearest-neighbour list. Nothing is fused.

---------------------------------------------
  A. COVERAGE. How many lines each axis has, and how the three overlap. This is
     the first thing that constrains everything downstream: an axis covering 928
     lines cannot answer questions about the other 745, and a fused score would
     hide exactly that.

  B. THE SIMILARITY DISTRIBUTION per axis. If an axis puts every pair at 0.98 it
     has no resolving power regardless of what it correlates with. Reported as
     the spread of off-diagonal similarity and the gap between a line's nearest
     neighbour and its median neighbour — that gap is what a "most similar
     lines" query actually sells.

  C. WHAT EACH AXIS CONTRIBUTES THAT THE OTHERS DO NOT. On the lines all three
     share, how much does axis A's similarity agree with axis B's? Two axes
     correlating at 0.9 are one axis. Two correlating at 0.2 are genuinely
     different views and worth keeping apart. This is the measurement that
     tells you whether fusing is even worth considering later.

  D. A ROBUSTNESS CHECK on the metric itself: Pearson versus Spearman. The
     pipeline median-imputes missing features, and Spearman is much less
     sensitive to that. If the two metrics disagree about who is close to whom,
     the imputation is driving the graph and that has to be fixed before any of
     it is used.


"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import common as C

ap = argparse.ArgumentParser()
ap.add_argument("--k", type=int, default=C.K_DEFAULT)
ap.add_argument("--top-var", type=int, default=C.TOP_VAR)
ap.add_argument("--max-pairs", type=int, default=200_000,
                help="pair subsample for the cross-axis and distribution stats")
args = ap.parse_args()

res = {"params": {"k": args.k, "top_var": args.top_var,
                  "max_feature_missing": C.MAX_FEATURE_MISSING}}

C.banner("WEEK 1 -- BUILD THE SIMILARITY GRAPH, ONE AXIS AT A TIME")

lin = C.load_lineage()
raw, sims = {}, {}

for axis in C.AXES:
    print(f"\n  [{axis}]  {C.AXIS_LABEL[axis]}")
    X = C.LOADERS[axis]()
    X = X[~X.index.duplicated()]
    print(f"    lines {X.shape[0]:,}   features {X.shape[1]:,}")
    S = C.similarity(X, method="pearson", top_var=args.top_var)
    NN = C.neighbours(S, k=args.k)
    C.save_axis(axis, S, NN)
    raw[axis] = X
    sims[axis] = S
    print(f"    wrote {C.axis_path(axis, 'similarity').name} and "
          f"{C.axis_path(axis, 'neighbours').name}")

# ------------------------------------------------------------ A
C.banner("1A -- COVERAGE, AND HOW THE AXES OVERLAP")
sets = {a: set(sims[a].index) for a in C.AXES}
print(f"  {'axis':<8} {'lines':>7} {'with lineage':>14}")
res["coverage"] = {}
for a in C.AXES:
    withlin = len(sets[a] & set(lin.index))
    print(f"  {a:<8} {len(sets[a]):>7,} {withlin:>14,}")
    res["coverage"][a] = {"n_lines": len(sets[a]), "n_with_lineage": withlin,
                          "n_features_raw": int(raw[a].shape[1])}
print()
print(f"  {'':<8}" + "".join(f"{b:>10}" for b in C.AXES))
for a in C.AXES:
    print(f"  {a:<8}" + "".join(f"{len(sets[a] & sets[b]):>10,}" for b in C.AXES))
all3 = sets["rna"] & sets["metab"] & sets["mirna"]
print(f"\n  lines on ALL THREE axes: {len(all3):,}")
print(f"  lines on RNA only       : {len(sets['rna'] - sets['metab'] - sets['mirna']):,}")
res["overlap"] = {f"{a}_{b}": len(sets[a] & sets[b]) for a in C.AXES for b in C.AXES}
res["n_all_three"] = len(all3)
print()
print("  This is why the axes stay separate. Any fused score would be RNA alone")
print("  for most of the panel and something else for a minority, with no way to")
print("  tell afterwards which lines got which.")

# ------------------------------------------------------------ B
C.banner("1B -- DOES EACH AXIS ACTUALLY RESOLVE ANYTHING?")
print("  A similarity that is 0.98 for every pair ranks nothing. The useful")
print("  quantity is the SPREAD, and the gap between a line's best neighbour and")
print("  its typical one.")
print()
print(f"  {'axis':<8} {'median sim':>11} {'p01':>8} {'p99':>8} {'sd':>8} "
      f"{'NN1 - median':>14}")
res["resolution"] = {}
for a in C.AXES:
    _, _, v = C.upper_pairs(sims[a], max_pairs=args.max_pairs)
    S = sims[a]
    best = np.nanmax(S.values, axis=1)
    med = np.nanmedian(S.values, axis=1)
    gap = float(np.nanmedian(best - med))
    q = np.percentile(v, [1, 50, 99])
    print(f"  {a:<8} {q[1]:>11.4f} {q[0]:>8.4f} {q[2]:>8.4f} {v.std():>8.4f} "
          f"{gap:>14.4f}")
    res["resolution"][a] = {"median": float(q[1]), "p01": float(q[0]),
                            "p99": float(q[2]), "sd": float(v.std()),
                            "nn1_minus_median_gap": gap}

# ------------------------------------------------------------ C
C.banner("1C -- WHAT DOES EACH AXIS CONTRIBUTE THAT THE OTHERS DO NOT?")
print("  Agreement between axes, on the lines they share. Computed on the SAME")
print("  pairs, so this is a like-for-like comparison.")
print()
common_all = sorted(all3)
res["axis_agreement"] = {}
if len(common_all) >= 50:
    print(f"  on the {len(common_all):,} lines shared by all three axes:")
    print(f"  {'pair':<18} {'Spearman':>10} {'Pearson':>9}")
    sub = {a: sims[a].loc[common_all, common_all] for a in C.AXES}
    n = len(common_all)
    iu, ju = np.triu_indices(n, k=1)
    if len(iu) > args.max_pairs:
        rng = np.random.default_rng(C.SEED)
        sel = rng.choice(len(iu), args.max_pairs, replace=False)
        iu, ju = iu[sel], ju[sel]
    vecs = {a: sub[a].values[iu, ju] for a in C.AXES}
    for i, a in enumerate(C.AXES):
        for b in C.AXES[i + 1:]:
            m = ~np.isnan(vecs[a]) & ~np.isnan(vecs[b])
            rho = spearmanr(vecs[a][m], vecs[b][m]).statistic
            r = float(np.corrcoef(vecs[a][m], vecs[b][m])[0, 1])
            print(f"  {a + ' vs ' + b:<18} {rho:>10.4f} {r:>9.4f}")
            res["axis_agreement"][f"{a}_vs_{b}"] = {"spearman": float(rho),
                                                    "pearson": r}
    print()
    print("  Reading this: rho near 1 means the two axes are one axis wearing two")
    print("  names. rho near 0 means they are genuinely different views and")
    print("  keeping them separate is buying you something.")
else:
    print("  too few shared lines for a stable comparison")

# also: do the axes agree on NEIGHBOURS, which is what a user actually sees
print()
print("  Neighbour-level agreement (what a user actually sees):")
print(f"  {'pair':<18} {'shared in top-' + str(args.k):>18}")
res["neighbour_overlap"] = {}
nbrs = {a: C.load_axis_neighbours(a) for a in C.AXES}
nsets = {a: nbrs[a].groupby("model_id").neighbour.apply(set) for a in C.AXES}
for i, a in enumerate(C.AXES):
    for b in C.AXES[i + 1:]:
        shared = [m for m in nsets[a].index if m in nsets[b].index]
        if len(shared) < 20:
            continue
        ov = np.mean([len(nsets[a][m] & nsets[b][m]) / args.k for m in shared])
        print(f"  {a + ' vs ' + b:<18} {ov:>17.1%}")
        res["neighbour_overlap"][f"{a}_vs_{b}"] = float(ov)

# ------------------------------------------------------------ D
C.banner("1D -- IS THE GRAPH AN ARTEFACT OF THE IMPUTATION?")
print("  Pearson (used) vs Spearman (rank-based, far less sensitive to the")
print("  median-imputed cells). If these disagree, the imputation is driving it.")
print()
print(f"  {'axis':<8} {'pair-level rho':>15} {'top-' + str(args.k) + ' overlap':>18}")
res["metric_robustness"] = {}
for a in C.AXES:
    Sp = C.similarity(raw[a], method="spearman", top_var=args.top_var)
    n = Sp.shape[0]
    iu, ju = np.triu_indices(n, k=1)
    if len(iu) > args.max_pairs:
        rng = np.random.default_rng(C.SEED)
        sel = rng.choice(len(iu), args.max_pairs, replace=False)
        iu, ju = iu[sel], ju[sel]
    v1 = sims[a].values[iu, ju]
    v2 = Sp.values[iu, ju]
    m = ~np.isnan(v1) & ~np.isnan(v2)
    rho = spearmanr(v1[m], v2[m]).statistic
    NNs = C.neighbours(Sp, k=args.k)
    s2 = NNs.groupby("model_id").neighbour.apply(set)
    s1 = nsets[a]
    shared = [x for x in s1.index if x in s2.index]
    ov = np.mean([len(s1[x] & s2[x]) / args.k for x in shared]) if shared else np.nan
    print(f"  {a:<8} {rho:>15.4f} {ov:>17.1%}")
    res["metric_robustness"][a] = {"pair_spearman": float(rho),
                                   "topk_overlap": float(ov)}


print("  Built, separately:")
for a in C.AXES:
    print(f"    {a:<8} {res['coverage'][a]['n_lines']:>6,} lines   "
          f"-> outputs/{a}_similarity.parquet, outputs/{a}_neighbours.parquet")
print()
print("  Nothing is fused. Week 2 tests what each one is worth.")

out = C.OUT / "01_build_axes_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
