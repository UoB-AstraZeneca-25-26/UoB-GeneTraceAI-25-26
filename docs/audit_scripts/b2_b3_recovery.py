"""B2 recovery test (as-is vs corrected) and B3 source-count sampling."""
import os, json, numpy as np, pandas as pd
from scipy import stats
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))
import tx_lib as T

GENE, LINEAGE = "TSPAN6", "lung"
con = T.connect()
df, gid = T.fetch_gene_expression(con, GENE, verbose=False)
res0, per0 = T.score_gene_lineage(df, verbose=False)

# ---------- how much do real sources actually disagree, in z units? -------------
zsd = per0.dropna(subset=["z"]).groupby("model_id")["z"].agg(["std", "size"])
zsd = zsd[zsd["size"] >= 2]["std"].dropna()
REAL_DISAGREE = float(zsd.median())
print(f"real cross-source disagreement (per-line SD of z, n>=2 sources):")
print(f"  n={len(zsd):,}  median={REAL_DISAGREE:.4f}  IQR="
      f"{zsd.quantile(.25):.4f}-{zsd.quantile(.75):.4f}")
print(f"the as-is planting uses N(0, 0.15) on the RAW value -- in z units that is")
base = df[df.lineage == LINEAGE]
sd_by_src = base.groupby("source")["value"].apply(
    lambda s: stats.median_abs_deviation(s, scale="normal"))
print(f"  ~0.15 / MAD = {(0.15/sd_by_src).round(3).to_dict()}  -> far below {REAL_DISAGREE:.3f}")


def recovery(df, corrected, n_planted=20, effect_sizes=(0.5, 1.0, 2.0, 3.0), seed=0):
    rng = np.random.default_rng(seed)
    base = df[df.lineage == LINEAGE]
    sd = base.groupby("source")["value"].apply(
        lambda s: stats.median_abs_deviation(s, scale="normal")).mean()
    med = base.groupby("source")["value"].median().mean()
    srcs = sorted(df.source.unique())
    rows = []
    for eff in effect_sizes:
        d = df.copy()
        if corrected:
            # SUBSTITUTE: overwrite existing lines in the lineage, so the
            # lineage size and centre are unchanged.
            pool = sorted(base.model_id.unique())
            planted = list(rng.choice(pool, size=min(n_planted, len(pool)), replace=False))
            for mid in planted:
                # per-source effect with REALISTIC disagreement (z units)
                e = eff + rng.normal(0, REAL_DISAGREE, size=len(srcs))
                for s, es in zip(srcs, e):
                    m = (d.model_id == mid) & (d.source == s)
                    if m.any():
                        d.loc[m, "value"] = med + es * sd
        else:
            # AS-IS: append new lines, near-identical value in every source
            planted = []
            add = []
            for i in range(n_planted):
                mid = f"PLANT_{eff}_{i}"; planted.append(mid)
                for s in srcs:
                    add.append({"model_id": mid, "source": s,
                                "value": med + eff * sd + rng.normal(0, .15),
                                "lineage": LINEAGE, "n_samples": 1})
            d = pd.concat([d, pd.DataFrame(add)], ignore_index=True)
        r, _ = T.score_gene_lineage(d, verbose=False)
        got = r[r.model_id.isin(planted)]
        rows.append({"effect_size_mad": eff, "n_planted": len(planted),
                     "detected_q05": int((got.q_value < .05).sum()),
                     "sensitivity": round(float((got.q_value < .05).mean()), 3),
                     "median_z_shrunk": round(float(got.z_shrunk.median()), 2),
                     "lineage_n": int((d.lineage == LINEAGE).groupby(
                         d[d.lineage == LINEAGE].model_id).ngroups)
                     if False else int(d[d.lineage == LINEAGE].model_id.nunique())})
    return pd.DataFrame(rows)


print(f"\n=== B2 AS-IS recovery (append, N(0,0.15) across sources) ===")
as_is = recovery(df, corrected=False)
print(as_is.to_string(index=False))
print(f"\n=== B2 CORRECTED recovery (substitute, disagreement SD={REAL_DISAGREE:.3f} z) ===")
corr = recovery(df, corrected=True)
print(corr.to_string(index=False))

n_lung_real = df[df.lineage == LINEAGE].model_id.nunique()
print(f"\nlineage '{LINEAGE}' real size = {n_lung_real} lines; as-is planting inflates it "
      f"to {as_is.lineage_n.iloc[0]}, corrected keeps it at {corr.lineage_n.iloc[0]}")

# ---------------- B3: source-count sampling ------------------------------------
print("\n=== B3 source-count test: is the sample random or outcome-selected? ===")
full = res0                                  # already sorted by z_shrunk DESC
three_all = full[full.n_sources == 3]
sel_head = three_all.model_id.head(30).tolist()          # notebook's choice
print(f"  3-source lines available            : {len(three_all):,}")
print(f"  notebook takes .head(30) of a frame sorted by z_shrunk DESCENDING")
print(f"  z_shrunk of the 30 selected         : median "
      f"{three_all.set_index('model_id').loc[sel_head,'z_shrunk'].median():.3f}  "
      f"range {three_all.set_index('model_id').loc[sel_head,'z_shrunk'].min():.3f} to "
      f"{three_all.set_index('model_id').loc[sel_head,'z_shrunk'].max():.3f}")
print(f"  z_shrunk of ALL 3-source lines      : median {three_all.z_shrunk.median():.3f}  "
      f"range {three_all.z_shrunk.min():.3f} to {three_all.z_shrunk.max():.3f}")
print(f"  -> the sample is the TOP of the distribution, selected on the outcome.")


def source_count(df, lines, label):
    d = df[~((df.model_id.isin(lines)) & (df.source != "depmap_expr"))]
    one, _ = T.score_gene_lineage(d, verbose=False)
    m = (res0[res0.model_id.isin(lines)][["model_id", "z_shrunk", "q_value"]]
         .merge(one[one.model_id.isin(lines)][["model_id", "z_shrunk", "q_value"]],
                on="model_id", suffixes=("_3src", "_1src")))
    m["drop"] = m.z_shrunk_3src - m.z_shrunk_1src
    print(f"  {label:22s} n={len(m):3d}  median drop={m['drop'].median():+.3f}  "
          f"mean={m['drop'].mean():+.3f}  IQR={m['drop'].quantile(.25):+.3f}"
          f"/{m['drop'].quantile(.75):+.3f}")
    return m


print("\n  effect of hiding 2 of 3 sources on z_shrunk:")
m_top = source_count(df, sel_head, "top-30 (notebook)")
rng = np.random.default_rng(42)
sel_rand = list(rng.choice(three_all.model_id.to_numpy(), size=30, replace=False))
m_rand = source_count(df, sel_rand, "random-30 (corrected)")
print(f"\n  overstatement from selecting on the outcome: "
      f"{m_top['drop'].median() - m_rand['drop'].median():+.3f} z units "
      f"({100*(m_top['drop'].median()/max(m_rand['drop'].median(),1e-9)-1):+.0f}%)")

json.dump({"real_disagreement_z_sd": REAL_DISAGREE,
           "as_is": as_is.to_dict("records"), "corrected": corr.to_dict("records"),
           "lineage_real_n": int(n_lung_real),
           "b3_top30_median_drop": float(m_top["drop"].median()),
           "b3_random30_median_drop": float(m_rand["drop"].median())},
          open(os.path.join(W, "b2_b3_results.json"), "w"), indent=1)
print("\nwritten b2_b3_results.json")
