"""B1: the null calibration test, as-is and corrected.

AS-IS (notebook cell 12): values permuted INDEPENDENTLY within each
(source, lineage). That destroys the cross-source correlation, so the FPR it
measures is the FPR under independence -- the assumption being tested.

WHY A JOINT PERMUTATION IS NOT ENOUGH (proved here, not asserted): permuting the
cell-line label once per lineage and applying it identically to every source is a
pure RELABELLING of lines that appear in the same source set. The multiset of
per-line z-vectors is unchanged, so the score distribution -- and hence the
number called significant -- is exactly the observed one. Demonstrated below.

THE CORRECTED NULL used here: draw per-line source z-vectors from a multivariate
normal with the MEASURED correlation matrix R and zero mean (no signal), keeping
each line's real source-presence pattern and real sqrt(n_samples) weights. Then
run the notebook's own combination + BH. This has no signal by construction, has
the real dependence structure, and can fail.
"""
import os, json, numpy as np, pandas as pd
from scipy import stats
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))
import tx_lib as T

A2 = json.load(open(os.path.join(W, "a2_results.json")))
NAME = {"depmap_expr": "DepMap", "hpa_rna": "HPA", "geo_expr": "GEO"}
SRC = ["depmap_expr", "hpa_rna", "geo_expr"]
R = np.array([[A2["R_zscore"][NAME[a]][NAME[b]] for b in SRC] for a in SRC])
print("measured R:\n", pd.DataFrame(R, index=list(NAME.values()),
                                    columns=list(NAME.values())).round(4).to_string())

GENE = "TSPAN6"
con = T.connect()
df, gid = T.fetch_gene_expression(con, GENE, verbose=False)

# ---------------- 1. as-is null (notebook), reduced n_perm for runtime -----------
def null_asis(df, n_perm=30, seed=0):
    rng = np.random.default_rng(seed)
    rates = []
    for _ in range(n_perm):
        d = df.copy()
        d["value"] = (d.groupby(["source", "lineage"], dropna=False)["value"]
                        .transform(lambda s: rng.permutation(s.values)))
        r, _ = T.score_gene_lineage(d, verbose=False)
        q = r.q_value.dropna()
        rates.append((q < .05).mean() if len(q) else 0.0)
    return np.array(rates)


print(f"\n[1] AS-IS null ({GENE}, 30 permutations) ...")
r_asis = null_asis(df, n_perm=30)
print(f"    mean FPR = {r_asis.mean():.4f}   p95 = {np.percentile(r_asis,95):.4f}   "
      f"max = {r_asis.max():.4f}   'calibrated' (<=0.07) = {r_asis.mean() <= 0.07}")

# ---------------- 2. observed significant rate (= joint-permutation rate) --------
res_obs, per_obs = T.score_gene_lineage(df, verbose=False)
q_obs = res_obs.q_value.dropna()
print(f"\n[2] OBSERVED significant rate on real data = {(q_obs < .05).mean():.4f}")
print("    A joint (same-shuffle-across-sources) permutation is a pure relabelling")
print("    and returns exactly this number -- it cannot be a null. Proof below.")

# empirical proof of the relabelling invariance
rng = np.random.default_rng(0)
d2 = df.copy()
mapping = {}
for lg, g in df.groupby("lineage", dropna=False):
    lines = g.model_id.unique()
    mapping.update(dict(zip(lines, rng.permutation(lines))))
d2["model_id"] = d2["model_id"].map(lambda m: mapping.get(m, m))
res_joint, _ = T.score_gene_lineage(d2, verbose=False)
qj = res_joint.q_value.dropna()
print(f"    joint-permutation significant rate       = {(qj < .05).mean():.4f}  "
      f"(identical: {abs((qj<.05).mean() - (q_obs<.05).mean()) < 1e-12})")

# ---------------- 3. corrected parametric null under measured R ------------------
pres = (per_obs.pivot_table(index="model_id", columns="source", values="z",
                            aggfunc="first").notna())
for s in SRC:
    if s not in pres.columns:
        pres[s] = False
pres = pres[SRC]
wt = (per_obs.pivot_table(index="model_id", columns="source", values="n_samples",
                          aggfunc="first")
      .reindex(index=pres.index, columns=SRC).fillna(1.0))
print(f"\n[3] CORRECTED null: {len(pres):,} lines, source-presence pattern preserved")
print("    pattern counts:",
      pres.apply(lambda r: "+".join(NAME[s] for s in SRC if r[s]), axis=1)
          .value_counts().to_dict())

L = np.linalg.cholesky(R)
Wm = np.sqrt(wt.to_numpy(float))
P = pres.to_numpy()


def null_corrected(n_perm=2000, seed=0, Rm=None):
    Lm = np.linalg.cholesky(Rm) if Rm is not None else L
    rng = np.random.default_rng(seed)
    n = len(P)
    rates = []
    for _ in range(n_perm):
        Zs = rng.standard_normal((n, 3)) @ Lm.T           # correlated, zero-mean
        Zs = np.where(P, Zs, np.nan)
        w = np.where(P, Wm, np.nan)
        num = np.nansum(w * Zs, axis=1)
        den = np.sqrt(np.nansum(w ** 2, axis=1))
        k = P.sum(axis=1)
        ok = k > 0
        Z = np.full(n, np.nan); Z[ok] = num[ok] / den[ok]
        p = 2 * (1 - stats.norm.cdf(np.abs(Z[ok])))
        q = stats.false_discovery_control(p)
        rates.append((q < .05).mean())
    return np.array(rates)


r_corr = null_corrected(n_perm=2000, seed=0, Rm=R)
r_indep = null_corrected(n_perm=2000, seed=0, Rm=np.eye(3))
print(f"\n    FPR under INDEPENDENCE (R = I)   : mean {r_indep.mean():.4f}  "
      f"p95 {np.percentile(r_indep,95):.4f}  max {r_indep.max():.4f}")
print(f"    FPR under MEASURED R             : mean {r_corr.mean():.4f}  "
      f"p95 {np.percentile(r_corr,95):.4f}  max {r_corr.max():.4f}")
print(f"    inflation in FPR                 : {r_corr.mean()/max(r_indep.mean(),1e-9):.1f}x")

# uncorrected per-line significance rate (before BH) is the cleaner diagnostic
def raw_alpha(n_perm=2000, seed=1, Rm=R):
    Lm = np.linalg.cholesky(Rm)
    rng = np.random.default_rng(seed); n = len(P); out = []
    for _ in range(n_perm):
        Zs = rng.standard_normal((n, 3)) @ Lm.T
        Zs = np.where(P, Zs, np.nan); w = np.where(P, Wm, np.nan)
        Z = np.nansum(w * Zs, axis=1) / np.sqrt(np.nansum(w ** 2, axis=1))
        out.append((2 * (1 - stats.norm.cdf(np.abs(Z))) < 0.05).mean())
    return np.array(out)


a_corr, a_indep = raw_alpha(Rm=R), raw_alpha(Rm=np.eye(3))
print(f"\n    nominal alpha=0.05, BEFORE BH:")
print(f"      under independence : {a_indep.mean():.4f}")
print(f"      under measured R   : {a_corr.mean():.4f}   "
      f"({a_corr.mean()/0.05:.2f}x nominal)")

json.dump({"gene": GENE,
           "asis_mean_fpr": float(r_asis.mean()), "asis_p95": float(np.percentile(r_asis, 95)),
           "observed_sig_rate": float((q_obs < .05).mean()),
           "joint_perm_sig_rate": float((qj < .05).mean()),
           "corrected_fpr_indep": float(r_indep.mean()),
           "corrected_fpr_measuredR": float(r_corr.mean()),
           "alpha_indep": float(a_indep.mean()), "alpha_measuredR": float(a_corr.mean())},
          open(os.path.join(W, "b1_results.json"), "w"), indent=1)
print("\nwritten b1_results.json")
