"""
cell_similarity/03_validate_tissue.py
-------------------------------------
WEEK 2, test 2.  The boring test with a known answer.

Do lines from the same tissue come out close together? If the graph puts two
lung lines next to each other, good. If it puts a lung line next to a blood
line, something is wrong. Cheap, and the answer is known in advance, which is
exactly what makes it worth running — a sanity check you could fail is the only
kind worth having.

Lineage comes from `cleaned_track_data/sample_info.parquet` (30 clean values,
a controlled vocabulary), not from parsing free-text disease strings.

WHAT IS REPORTED
----------------
  A. NEIGHBOUR PURITY. Of a line's k nearest neighbours, what fraction share its
     lineage, against the base rate of that lineage in the panel. Enrichment is
     the ratio. This is the direct form of the question.
  B. SEPARATION. AUROC of similarity for same-lineage versus different-lineage
     pairs. One number per axis, and it answers "could you tell tissue from this
     graph alone".
  C. PER LINEAGE. Which tissues the graph recovers and which it does not. A
     pooled number hides that blood is trivially separable and, say, upper
     aerodigestive is not.
  D. THE FAILURE CASES, NAMED. The specific lines whose nearest neighbour is a
     different lineage, with what it was confused for. A confusion table of the
     worst offenders — because "94% pure" is not actionable and "these 40 lines
     are wrong, and they are mostly X mistaken for Y" is.

Note on interpretation, which matters for how 02 is read: a HIGH score here is
not straightforwardly good. An axis that reproduces the tissue label perfectly
and adds nothing beyond it has scored 1.0 here and is worth nothing over a
lookup table. Read this together with `02`'s within-lineage result.

Read-only. Writes only into cell_similarity/outputs/.

Run:
    python cell_similarity/03_validate_tissue.py
"""
import json

import numpy as np
import pandas as pd

import common as C

res = {}
C.banner("WEEK 2, TEST 2 -- DO SAME-TISSUE LINES COME OUT CLOSE TOGETHER?")

lin = C.load_lineage()
print(f"\n  lineage vocabulary: {lin.lineage.nunique()} values over "
      f"{len(lin):,} lines")

# ------------------------------------------------------------ A
C.banner("3A -- NEIGHBOUR PURITY")
print("  Of a line's k nearest neighbours, what share are the same lineage?")
print("  'base' is what you would get by picking neighbours at random.")
print()
print(f"  {'axis':<8} {'lines':>7} {'k':>4} {'purity':>9} {'base':>8} "
      f"{'enrichment':>12} {'NN1 correct':>13}")
res["purity"] = {}
store = {}
for axis in C.available_axes():
    NN = C.load_axis_neighbours(axis)
    NN = NN.join(lin.lineage.rename("lin_self"), on="model_id")
    NN = NN.join(lin.lineage.rename("lin_nbr"), on="neighbour")
    NN = NN.dropna(subset=["lin_self", "lin_nbr"])
    if not len(NN):
        continue
    NN["same"] = NN.lin_self == NN.lin_nbr
    k = int(NN["rank"].max())
    lines = NN.model_id.unique()
    # base rate: probability a random other line shares your lineage
    counts = lin.loc[lin.index.isin(set(NN.model_id) | set(NN.neighbour))
                     ].lineage.value_counts()
    tot = counts.sum()
    base_by_line = lin.reindex(lines).lineage.map(
        lambda l: (counts.get(l, 0) - 1) / (tot - 1) if pd.notna(l) else np.nan)
    purity = NN.groupby("model_id")["same"].mean()
    nn1 = NN[NN["rank"] == 1].set_index("model_id")["same"]
    base = float(base_by_line.mean())
    p = float(purity.mean())
    print(f"  {axis:<8} {len(lines):>7,} {k:>4} {p:>9.1%} {base:>8.1%} "
          f"{p/base:>11.1f}x {nn1.mean():>12.1%}")
    res["purity"][axis] = {"n_lines": int(len(lines)), "k": k,
                           "purity": p, "base_rate": base,
                           "enrichment": float(p / base),
                           "nn1_same_lineage": float(nn1.mean())}
    store[axis] = NN

# ------------------------------------------------------------ B
C.banner("3B -- SEPARATION:  could you recover tissue from the graph alone?")
print("  AUROC of similarity for same-lineage vs different-lineage pairs.")
print("  0.5 = no tissue information. 1.0 = tissue is perfectly recoverable.")
print()
print(f"  {'axis':<8} {'pairs':>11} {'same-lin':>10} {'AUROC':>9}")
res["separation"] = {}
for axis in C.available_axes():
    S = C.load_axis_similarity(axis)
    lines = sorted(set(S.index) & set(lin.index))
    Sx = S.loc[lines, lines]
    lv = lin.reindex(lines).lineage.values
    n = len(lines)
    iu, ju = np.triu_indices(n, k=1)
    v = Sx.values[iu, ju]
    same = lv[iu] == lv[ju]
    ok = ~np.isnan(v) & pd.notna(lv[iu]) & pd.notna(lv[ju])
    v, same = v[ok], same[ok]
    from scipy.stats import rankdata
    r = rankdata(v)
    n1, n0 = int(same.sum()), int((~same).sum())
    auc = (r[same].sum() - n1 * (n1 + 1) / 2) / (n1 * n0)
    print(f"  {axis:<8} {len(v):>11,} {n1:>10,} {auc:>9.4f}")
    res["separation"][axis] = {"n_pairs": int(len(v)), "n_same": n1,
                               "auroc": float(auc)}

# ------------------------------------------------------------ C
C.banner("3C -- PER LINEAGE  (RNA axis).  Which tissues does it actually recover?")
axis = "rna" if "rna" in store else C.available_axes()[0]
NN = store[axis]
per = NN.groupby(["model_id", "lin_self"])["same"].mean().reset_index()
agg = (per.groupby("lin_self")
       .agg(lines=("model_id", "nunique"), purity=("same", "mean"))
       .sort_values("purity", ascending=False))
agg = agg[agg.lines >= 5]
print(f"  {'lineage':<28} {'lines':>6} {'purity':>9}")
for l, r in agg.iterrows():
    bar = "#" * int(round(r.purity * 30))
    print(f"  {str(l):<28} {int(r.lines):>6} {r.purity:>8.1%}  {bar}")
res["per_lineage_rna"] = {str(l): {"lines": int(r.lines), "purity": float(r.purity)}
                          for l, r in agg.iterrows()}

# ------------------------------------------------------------ D
C.banner("3D -- THE FAILURE CASES, NAMED")
print(f"  Lines on the {axis} axis whose NEAREST neighbour is a different lineage.")
print()
nn1 = NN[NN["rank"] == 1]
bad = nn1[~nn1.same]
print(f"  {len(bad):,} of {len(nn1):,} lines ({len(bad)/len(nn1):.1%}) have a "
      f"nearest neighbour from another lineage")
print()
conf = (bad.groupby(["lin_self", "lin_nbr"]).size()
        .sort_values(ascending=False).head(15))
print(f"  most common confusions:")
print(f"  {'actual':<26} {'confused with':<26} {'n':>4}")
for (a, b), n in conf.items():
    print(f"  {str(a):<26} {str(b):<26} {n:>4}")
res["confusions"] = [{"actual": str(a), "confused_with": str(b), "n": int(n)}
                     for (a, b), n in conf.items()]
res["nn1_wrong_frac"] = float(len(bad) / len(nn1))

print()
print("  Worst individual cases (lowest-similarity nearest neighbour that is")
print("  also a lineage mismatch -- these are lines the graph has no real")
print("  neighbour for, and a query should say so rather than return one):")
worst = bad.nsmallest(10, "similarity")
print(f"  {'line':<16} {'lineage':<24} {'returned':<16} {'its lineage':<22} {'sim':>7}")
for _, r in worst.iterrows():
    print(f"  {r.model_id:<16} {str(r.lin_self):<24} {r.neighbour:<16} "
          f"{str(r.lin_nbr):<22} {r.similarity:>7.3f}")
res["worst_cases"] = worst[["model_id", "lin_self", "neighbour", "lin_nbr",
                            "similarity"]].to_dict("records")

# ------------------------------------------------------------ verdict
C.banner("VERDICT")
for axis in res["purity"]:
    p = res["purity"][axis]
    a = res["separation"][axis]["auroc"]
    if p["enrichment"] > 3 and a > 0.75:
        v = "recovers_tissue_strongly"
    elif p["enrichment"] > 1.5 and a > 0.6:
        v = "recovers_tissue_weakly"
    else:
        v = "does_NOT_recover_tissue__investigate_before_using"
    print(f"  {axis:<8} purity {p['purity']:.1%} ({p['enrichment']:.1f}x base), "
          f"AUROC {a:.3f}  -> {v}")
    res.setdefault("verdict", {})[axis] = v
print()
print("  Remember this cuts both ways. A perfect score here means the axis has")
print("  reproduced the tissue label and nothing more. Read with 02's")
print("  within-lineage column, which is where the added value lives.")

out = C.OUT / "03_validate_tissue_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
