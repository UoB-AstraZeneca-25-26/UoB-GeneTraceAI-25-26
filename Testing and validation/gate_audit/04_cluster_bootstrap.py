"""
gate_audit/04_cluster_bootstrap.py
----------------------------------
AUDIT 4 of 6.  "The p-value is inflated, certainly."

THE PROBLEM
-----------
The published gate-region p is 2.8e-76 (Wilcoxon signed-rank of per-gene pAUC
against 0.5). It is computed as though the units were independent. Two
dependence structures are ignored, and they act in opposite places:

  ACROSS LINES.  Cell lines cluster by lineage. Lines from the same lineage share
      expression programmes AND dependency profiles, so the effective number of
      independent lines is far below the nominal count. This is documented as
      gaps 7 and 17 in the project's own risk register.

  ACROSS GENES.  The Wilcoxon is taken over genes, and genes are not independent
      either — co-expressed and paralogous genes carry near-identical abundance
      vectors. n_genes is the Wilcoxon's sample size, so this inflates it
      directly.

WHAT THIS SCRIPT DOES
---------------------
Three nulls, increasingly honest, all on the same statistic:

  1. NOMINAL          the published test. Reported for comparison only.
  2. LINEAGE CLUSTER BOOTSTRAP
        Resample LINEAGES with replacement, rebuild the panel from the lines in
        the drawn lineages, recompute the statistic. This is the standard
        cluster bootstrap (Field & Welsh 2007); it propagates within-lineage
        correlation into the sampling distribution instead of assuming it away.
  3. TWO-WAY (lineage x gene) CLUSTER BOOTSTRAP
        Resample lineages AND genes. This is the interval to quote.

A design-effect number is reported alongside — the ratio of the bootstrap
standard error to the nominal one — because that single figure is what makes the
inflation legible to a statistician on the panel.

WHAT TO QUOTE
-------------
Not 2.8e-76. Quote the two-way bootstrap CI on the effect and say the effect is
significant at that interval. An overwhelming effect stays overwhelming under an
honest interval; the difference is that the honest one survives scrutiny.

Read-only. Writes only into gate_audit/outputs/.

Run:
    python gate_audit/04_cluster_bootstrap.py
    python gate_audit/04_cluster_bootstrap.py --n-boot 5000 --valid-only
"""
import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import wilcoxon

import common as C

ap = argparse.ArgumentParser()
ap.add_argument("--n-boot", type=int, default=600)
ap.add_argument("--n-genes", type=int, default=800,
                help="genes carried into the bootstrap. The statistic is a MEDIAN "
                     "over genes, so a random subsample estimates it with an "
                     "error far below the cluster SE this audit exists to "
                     "measure; the full 3.9M pAUC evaluations buy nothing and "
                     "cost hours. Raise it to check.")
ap.add_argument("--valid-only", action="store_true",
                help="restrict to genes passing the silence guard (audit 2)")
ap.add_argument("--min-lineage-n", type=int, default=5,
                help="lineages with fewer lines than this are pooled into 'other'")
args = ap.parse_args()

panel = C.read_panel()
G = pd.read_parquet(C.OUT / "genes.parquet")
if args.valid_only:
    keep = set(G.loc[G.is_valid, "ensg_id"])
    panel = panel[panel.ensg_id.isin(keep)]
_all_genes = np.array(sorted(panel.ensg_id.unique()))
if args.n_genes < len(_all_genes):
    _sub = np.random.default_rng(C.SEED).choice(_all_genes, args.n_genes,
                                                replace=False)
    panel = panel[panel.ensg_id.isin(set(_sub))]

res = {"n_boot": args.n_boot, "valid_only": args.valid_only}

C.banner("AUDIT 4 -- CLUSTER BOOTSTRAP FOR AN HONEST p")

# ------------------------------------------------------------ clusters
counts = panel.drop_duplicates("model_id").lineage.value_counts()
small = set(counts[counts < args.min_lineage_n].index)
panel = panel.copy()
panel["cluster"] = panel.lineage.where(~panel.lineage.isin(small), "other")
line_cluster = panel.drop_duplicates("model_id").set_index("model_id")["cluster"]
clusters = sorted(line_cluster.unique())
cl_lines = {c: list(line_cluster.index[line_cluster == c]) for c in clusters}
n_lines = line_cluster.shape[0]

print(f"  genes                      : {panel.ensg_id.nunique():,}")
print(f"  lines                      : {n_lines:,}")
print(f"  lineage clusters           : {len(clusters)} "
      f"(pooled {len(small)} lineages with <{args.min_lineage_n} lines)")
sizes = pd.Series({c: len(v) for c, v in cl_lines.items()})
print(f"  cluster size               : median {sizes.median():.0f}  "
      f"max {sizes.max()}  ('{sizes.idxmax()}')")
# Kish effective sample size under exchangeable within-cluster correlation
kish = float(sizes.sum() ** 2 / (sizes ** 2).sum())
print(f"  Kish effective n (clusters): {kish:.1f} of {len(clusters)} clusters")
print(f"     -> for line-level claims the effective n is nearer {kish:.0f} than "
      f"{n_lines:,}")
res["clusters"] = {"n_clusters": len(clusters), "n_lines": int(n_lines),
                   "kish_effective_n": kish,
                   "largest": str(sizes.idxmax()), "largest_n": int(sizes.max())}

# ------------------------------------------------------------ statistic
# Pre-index once. model_id is mapped to a dense integer code so a bootstrap
# replicate's line multiplicities are an ARRAY lookup, not 2.9M dict lookups per
# replicate -- that difference is the whole runtime of this audit.
line_codes = {m: i for i, m in enumerate(line_cluster.index)}
panel = panel[panel.model_id.isin(line_codes)]
by_gene = {}
for g, d in panel.groupby("ensg_id", sort=False):
    by_gene[g] = (d.Z.values.astype(float), d.dep.values.astype(bool),
                  np.fromiter((line_codes[m] for m in d.model_id.values),
                              np.int32, len(d)))
genes = sorted(by_gene)
cl_codes = {c: np.fromiter((line_codes[m] for m in v), np.int32, len(v))
            for c, v in cl_lines.items()}
n_codes = len(line_codes)
print(f"  bootstrap replicates       : {args.n_boot:,}")
print(f"  genes carried              : {len(genes):,}")


def stat_from(gene_subset, mult_arr):
    """Median per-gene gate-region pAUC. `mult_arr[code]` is how many times that
    line was drawn, which is how a cluster bootstrap resamples lines without
    materialising duplicated frames."""
    vals = []
    for g in gene_subset:
        s, p, codes = by_gene[g]
        w = mult_arr[codes]
        sel = w > 0
        if sel.sum() < C.MIN_LINES:
            continue
        ss = np.repeat(s[sel], w[sel])
        pp = np.repeat(p[sel], w[sel])
        if pp.sum() < C.MIN_POS or (~pp).sum() < C.MIN_NEG:
            continue
        v = C.partial_auc(ss, pp, *C.GATE_REGION)
        if np.isfinite(v):
            vals.append(v)
    return float(np.median(vals)) if vals else np.nan


# ------------------------------------------------------------ 1. nominal
C.banner("4A -- NOMINAL (the published test)")
pa = []
for g in genes:
    s, p, _ = by_gene[g]
    v = C.partial_auc(s, p, *C.GATE_REGION)
    if np.isfinite(v):
        pa.append(v)
pa = np.array(pa)
obs = float(np.median(pa))
w_nom = wilcoxon(pa - 0.5)
se_nom = float(np.std(pa, ddof=1) / np.sqrt(len(pa)))
print(f"  gate-region pAUC (median over {len(pa):,} genes) : {obs:.4f}")
print(f"  Wilcoxon signed-rank vs 0.5                     : p = {w_nom.pvalue:.4g}")
print(f"  nominal SE of the mean                          : {se_nom:.5f}")
print(f"  -> treats {len(pa):,} genes x {n_lines:,} lines as independent units")
res["nominal"] = {"pauc": obs, "wilcoxon_p": float(w_nom.pvalue),
                  "se": se_nom, "n_genes": int(len(pa))}

rng = np.random.default_rng(C.SEED)

# ------------------------------------------------------------ 2. lineage
C.banner("4B -- LINEAGE CLUSTER BOOTSTRAP")
print(f"  resampling {len(clusters)} lineage clusters with replacement, "
      f"{args.n_boot:,} times")
boot_l = []
for b in range(args.n_boot):
    drawn = rng.choice(clusters, size=len(clusters), replace=True)
    mult = np.zeros(n_codes, np.int32)
    np.add.at(mult, np.concatenate([cl_codes[c] for c in drawn]), 1)
    v = stat_from(genes, mult)
    if np.isfinite(v):
        boot_l.append(v)
    if (b + 1) % max(1, args.n_boot // 5) == 0:
        print(f"    {b+1:,}/{args.n_boot:,}")
boot_l = np.array(boot_l)
lo_l, hi_l = np.percentile(boot_l, [2.5, 97.5])
se_l = float(boot_l.std(ddof=1))
p_l = 2 * min((boot_l <= 0.5).mean(), (boot_l >= 0.5).mean())
print(f"\n  pAUC {obs:.4f}   95% CI [{lo_l:.4f}, {hi_l:.4f}]   SE {se_l:.5f}")
print(f"  bootstrap p (share of replicates at or across 0.5): "
      f"{'< ' + format(1/len(boot_l), '.2g') if p_l == 0 else format(p_l, '.4g')}")
print(f"  DESIGN EFFECT vs nominal SE: {se_l/se_nom:.1f}x")
res["lineage_bootstrap"] = {"ci": [float(lo_l), float(hi_l)], "se": se_l,
                            "p": float(p_l), "design_effect": float(se_l / se_nom),
                            "n_replicates": int(len(boot_l))}

# ------------------------------------------------------------ 3. two-way
C.banner("4C -- TWO-WAY (lineage x gene) CLUSTER BOOTSTRAP  <- quote this one")
boot_2 = []
gi = np.arange(len(genes))
for b in range(args.n_boot):
    drawn = rng.choice(clusters, size=len(clusters), replace=True)
    mult = np.zeros(n_codes, np.int32)
    np.add.at(mult, np.concatenate([cl_codes[c] for c in drawn]), 1)
    gsub = [genes[j] for j in rng.choice(gi, size=len(gi), replace=True)]
    v = stat_from(gsub, mult)
    if np.isfinite(v):
        boot_2.append(v)
    if (b + 1) % max(1, args.n_boot // 5) == 0:
        print(f"    {b+1:,}/{args.n_boot:,}")
boot_2 = np.array(boot_2)
lo2, hi2 = np.percentile(boot_2, [2.5, 97.5])
se2 = float(boot_2.std(ddof=1))
p2 = 2 * min((boot_2 <= 0.5).mean(), (boot_2 >= 0.5).mean())
print(f"\n  pAUC {obs:.4f}   95% CI [{lo2:.4f}, {hi2:.4f}]   SE {se2:.5f}")
print(f"  bootstrap p: "
      f"{'< ' + format(1/len(boot_2), '.2g') if p2 == 0 else format(p2, '.4g')}")
print(f"  DESIGN EFFECT vs nominal SE: {se2/se_nom:.1f}x")
res["twoway_bootstrap"] = {"ci": [float(lo2), float(hi2)], "se": se2,
                           "p": float(p2), "design_effect": float(se2 / se_nom),
                           "n_replicates": int(len(boot_2))}

np.save(C.OUT / "04_bootstrap_lineage.npy", boot_l)
np.save(C.OUT / "04_bootstrap_twoway.npy", boot_2)

# ------------------------------------------------------------ verdict
C.banner("VERDICT -- WHAT TO WRITE IN THE THESIS")
survives = lo2 > 0.5
print(f"  DO NOT quote p = {w_nom.pvalue:.3g}. It assumes independent cell lines")
print(f"  and independent genes, and both assumptions fail here.")
print()
print(f"  QUOTE INSTEAD:")
print(f"    gate-region pAUC = {obs:.3f}, 95% CI [{lo2:.3f}, {hi2:.3f}]")
print(f"    (two-way cluster bootstrap over {len(clusters)} lineages and "
      f"{len(genes):,} genes, {len(boot_2):,} replicates)")
print(f"    design effect {se2/se_nom:.1f}x the nominal standard error")
print()
if survives:
    verdict = "effect_survives_clustering__quote_the_bootstrap_interval"
    print("  The effect survives clustering: the interval excludes chance. The")
    print("  honest version is still decisive, and it is defensible.")
else:
    verdict = "effect_does_not_survive_clustering"
    print("  The interval includes chance once clustering is respected. The")
    print("  nominal p was carrying the claim.")
print(f"\n  VERDICT: {verdict}")
res["verdict"] = verdict
res["quote"] = (f"gate-region pAUC = {obs:.3f}, 95% CI [{lo2:.3f}, {hi2:.3f}] "
                f"(two-way cluster bootstrap, {len(clusters)} lineage clusters, "
                f"{len(genes)} genes, {len(boot_2)} replicates)")

out = C.OUT / "04_cluster_bootstrap_results.json"
out.write_text(json.dumps(res, indent=2, default=float), encoding="utf-8")
print(f"\nwrote {out}")
