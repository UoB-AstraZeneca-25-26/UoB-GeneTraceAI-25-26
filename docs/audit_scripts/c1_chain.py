"""C1/C2/C3: the significance chain.

C1  Stouffer combines standardised test statistics. The notebook feeds it robust
    z -- an EFFECT SIZE in MAD units, which is not N(0,1) under the null. The
    rank-based empirical p is already computed (per['p_emp']) and discarded.
    Convert it to a standardised statistic and combine those instead.
C2  fisher_p is computed and never used.
C3  BH across lines that share a lineage median/MAD.
"""
import os, json, numpy as np, pandas as pd
from scipy import stats
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))
import tx_lib as T

GENE = "TSPAN6"
con = T.connect()
df, gid = T.fetch_gene_expression(con, GENE, verbose=False)

# ---------------- is robust z actually N(0,1) under the null? -------------------
res0, per0 = T.score_gene_lineage(df, verbose=False)
z = per0["z"].dropna()
print("=== C1a: distribution of the quantity being fed to Stouffer ===")
print(f"  robust z: n={len(z):,} mean={z.mean():.4f} sd={z.std():.4f} "
      f"skew={stats.skew(z):.3f} kurtosis={stats.kurtosis(z):.3f}")
print(f"  P(|z|>1.96) observed = {(z.abs()>1.96).mean():.4f}   (N(0,1) predicts 0.0500)")
print(f"  -> sd {z.std():.3f} != 1 and tails "
      f"{'heavier' if (z.abs()>1.96).mean()>0.05 else 'lighter'} than normal, so a")
print("     p-value read off the normal CDF is not a p-value.")

# the empirical p, converted to a standardised statistic
sgn = np.sign(per0["value"] - per0["lineage_median"]).replace(0, 1.0)
pe = per0["p_emp"].clip(lower=1e-12, upper=1 - 1e-12)
zemp = (sgn * stats.norm.isf(pe / 2.0)).replace([np.inf, -np.inf], np.nan).dropna()
print(f"\n  empirical-p-derived z: n={len(zemp):,} mean={zemp.mean():.4f} "
      f"sd={zemp.std():.4f}  P(|z|>1.96)={(zemp.abs()>1.96).mean():.4f}")

# ---------------- C1b: effect on verdicts ---------------------------------------
res_emp, _ = T.score_gene_lineage(df, verbose=False, use_empirical_z=True)
m = (res0[["model_id", "z_shrunk", "q_value", "verdict"]]
     .merge(res_emp[["model_id", "z_shrunk", "q_value", "verdict"]],
            on="model_id", suffixes=("_cur", "_emp")))
print("\n=== C1b: verdict changes, current chain vs empirical-p chain ===")
print(f"  lines compared            : {len(m):,}")
changed = m[m.verdict_cur != m.verdict_emp]
print(f"  lines changing verdict    : {len(changed):,} ({100*len(changed)/len(m):.1f}%)")
sig_cur = (m.q_value_cur < .05).sum(); sig_emp = (m.q_value_emp < .05).sum()
print(f"  significant at q<0.05     : {sig_cur:,} -> {sig_emp:,} "
      f"({sig_emp-sig_cur:+,})")
print("\n  direction of change (top transitions):")
print(changed.groupby(["verdict_cur", "verdict_emp"]).size()
      .sort_values(ascending=False).head(12).to_string())

# ---------------- C1c: the resolution floor -------------------------------------
print("\n=== C1c: resolution floor of the empirical p ===")
print("  empirical_p is two-tailed: p_min = 2 * 1/(n+1) for a lineage of n peers.")
for n in (3, 15, 20, 39, 40, 50, 100, 200):
    print(f"    n={n:4d} -> p_min = {2/(n+1):.4f}"
          + ("   <- cannot reach p<0.05" if 2/(n+1) >= 0.05 else ""))
n_need = int(np.ceil(2/0.05 - 1))
print(f"  smallest lineage that can reach p<0.05 at all: n = {n_need}")

lin = T.lineage_map(con)
sizes = df[df.source == "depmap_expr"].assign(
    lg=lambda d: d.model_id.map(lin).fillna(T.NO_LINEAGE)).groupby("lg").model_id.nunique()
sizes = sizes.sort_values(ascending=False)
print(f"\n  lineage sizes (DepMap coverage, {len(sizes)} lineages):")
print(sizes.to_string())
too_small = sizes[sizes < n_need]
print(f"\n  lineages where NO line can reach p<0.05 via the empirical route: "
      f"{len(too_small)} of {len(sizes)}")
print(f"  cell lines in them: {int(too_small.sum()):,} of {int(sizes.sum()):,} "
      f"({100*too_small.sum()/sizes.sum():.1f}%)")
print(f"  at MIN_PEERS_FOR_LINEAGE={T.MIN_PEERS_FOR_LINEAGE}: p_min = "
      f"{2/(T.MIN_PEERS_FOR_LINEAGE+1):.4f} -- {2/(T.MIN_PEERS_FOR_LINEAGE+1)/0.05:.1f}x "
      f"the 0.05 threshold, before any BH correction")

# ---------------- C2: the unused routes -----------------------------------------
print("\n=== C2: quantities computed and discarded ===")
fp = res0["fisher_p"].dropna()
print(f"  fisher_p    : computed for {len(fp):,} lines, referenced nowhere in the")
print(f"                verdict / q_value / ranking chain. Range "
      f"{fp.min():.2e}-{fp.max():.3f}")
agree = np.corrcoef(res0.loc[res0.fisher_p.notna() & res0.p_combined.notna(), "fisher_p"],
                    res0.loc[res0.fisher_p.notna() & res0.p_combined.notna(), "p_combined"])[0, 1]
print(f"  correlation(fisher_p, p_combined) = {agree:.4f}")
print(f"  p_emp       : computed per (line, source); used only to build fisher_p")
print(f"  lineage_silent_any : computed in the notebook, never read")

# ---------------- C3: dependence among lines within a lineage -------------------
print("\n=== C3: are within-lineage line scores dependent? ===")
sub = per0[(per0.source == "depmap_expr") & per0.z.notna()]
cors = []
for lg, g in sub.groupby("lineage"):
    if len(g) >= 20:
        cors.append({"lineage": lg, "n": len(g), "mean_z": g.z.mean(),
                     "sum_z": g.z.sum()})
c = pd.DataFrame(cors)
print("  lines in a lineage are scored against a SHARED median and MAD, so their")
print("  z's are constrained: sum of deviations about the median is near zero.")
print(f"  observed mean z per lineage (should be ~0 by construction):")
print(f"    median across lineages = {c.mean_z.median():+.4f}  "
      f"range {c.mean_z.min():+.4f} to {c.mean_z.max():+.4f}")
print("  -> negative dependence, which is NOT the PRDS structure BH assumes.")
print("     BH under negative dependence is generally conservative, so the FDR")
print("     is unlikely to be violated by this route -- the B1 correlation")
print("     inflation is the binding problem, not this.")

json.dump({"robust_z_sd": float(z.std()), "robust_z_tail": float((z.abs() > 1.96).mean()),
           "emp_z_sd": float(zemp.std()), "emp_z_tail": float((zemp.abs() > 1.96).mean()),
           "n_lines": int(len(m)), "n_verdict_changed": int(len(changed)),
           "sig_current": int(sig_cur), "sig_empirical": int(sig_emp),
           "p_min_at_min_peers": 2 / (T.MIN_PEERS_FOR_LINEAGE + 1),
           "smallest_lineage_reaching_p05": n_need,
           "lineage_sizes": {k: int(v) for k, v in sizes.items()},
           "n_lineages_too_small": int(len(too_small)),
           "n_lines_in_too_small": int(too_small.sum())},
          open(os.path.join(W, "c1_results.json"), "w"), indent=1)
print("\nwritten c1_results.json")
