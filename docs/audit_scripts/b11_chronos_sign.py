"""B1.1 diagnostic: is the Chronos/essentiality sign convention correct?

Test 1: polarity of the source column, using pan-essential vs non-essential
        reference gene sets shipped with DepMap (ReferenceEssentials.csv /
        ReferenceNonEssentials.csv) and by symbol families (RPL*/RPS*/POLR2A/
        PSMA*/EIF*).
Test 2: recompute the confirmed/contradicted classification restricted to the
        pan-essential control set.
"""
import os, re, numpy as np, pandas as pd
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")

gl = pd.read_parquet("reference/gene_lookup.parquet", columns=["ensg_id", "hgnc_symbol"])
gl["ensg_id"] = gl["ensg_id"].str.upper()
sym2ensg = dict(zip(gl.hgnc_symbol, gl.ensg_id))

ess = pd.read_csv("data/DepMap_Chronos/ReferenceEssentials.csv")
non = pd.read_csv("data/DepMap_Chronos/ReferenceNonEssentials.csv")
print("ReferenceEssentials cols:", list(ess.columns)[:5], len(ess))
print("ReferenceNonEssentials cols:", list(non.columns)[:5], len(non))


def to_symbols(df):
    col = df.columns[0] if df.columns[0] != "Unnamed: 0" else df.columns[1]
    s = df[col].astype(str)
    return set(s.str.split(" ").str[0].str.strip())


ess_sym, non_sym = to_symbols(ess), to_symbols(non)
ess_ensg = {sym2ensg[s] for s in ess_sym if s in sym2ensg}
non_ensg = {sym2ensg[s] for s in non_sym if s in sym2ensg}
print(f"mapped essentials {len(ess_ensg)}/{len(ess_sym)}, non-essentials {len(non_ensg)}/{len(non_sym)}")

fam = gl[gl.hgnc_symbol.str.match(r"^(RPL\d|RPS\d|POLR2A$|PSMA\d|EIF\d)", na=False)]
fam_ensg = set(fam.ensg_id)
print(f"symbol-family pan-essential set: {len(fam_ensg)} genes")

ch = pd.read_parquet("validation/prepared/chronos_long.parquet")
ch["ensg_id"] = ch["ensg_id"].str.upper()

med = ch.groupby("ensg_id")["essentiality"].median()
print("\n--- TEST 1: source-column polarity ---")
print(f"all genes            median-of-medians = {med.median():.3f}   n={len(med):,}")
for name, s in [("DepMap ReferenceEssentials", ess_ensg),
                ("DepMap ReferenceNonEssentials", non_ensg),
                ("RPL*/RPS*/POLR2A/PSMA*/EIF*", fam_ensg)]:
    sub = med[med.index.isin(s)]
    print(f"{name:32s} median = {sub.median():.3f}   n={len(sub):,}")

# ---------------- TEST 2: confirmed/contradicted on the control set -------------
ib = pd.read_parquet("validation/prepared/id_bridge.parquet")
ch2 = ch.merge(ib, on="sanger_model_id", how="inner")
ch2["model_id"] = ch2["model_id"].str.lower()
ch2["ensg_l"] = ch2["ensg_id"].str.lower()

# EXACTLY the pipeline's expression (build_chronos_validation.py:92)
ch2["chronos_pct"] = 1.0 - ch2.groupby("ensg_l")["essentiality"].rank(pct=True, method="average")

pred = pd.read_parquet("src/pipeline/outputs/predictions_with_confidence.parquet",
                       columns=["model_id", "ensg_id", "core_score"])
pred["model_id"] = pred["model_id"].str.lower()

ctrl = {e.lower() for e in (fam_ensg | ess_ensg)}
m = pred[pred.ensg_id.isin(ctrl)].merge(
    ch2[["model_id", "ensg_l", "chronos_pct"]].rename(columns={"ensg_l": "ensg_id"}),
    on=["model_id", "ensg_id"], how="inner")
print(f"\n--- TEST 2: control-set classification --- rows={len(m):,} genes={m.ensg_id.nunique():,}")

from scipy.stats import t as tdist
m["r1"] = m.groupby("ensg_id")["core_score"].rank()
m["r2"] = m.groupby("ensg_id")["chronos_pct"].rank()
g = m.groupby("ensg_id")
m["dx"] = m.r1 - g.r1.transform("mean"); m["dy"] = m.r2 - g.r2.transform("mean")
agg = m.assign(dxdy=m.dx*m.dy, dx2=m.dx**2, dy2=m.dy**2).groupby("ensg_id").agg(
    n=("r1", "size"), sxy=("dxdy", "sum"), sxx=("dx2", "sum"), syy=("dy2", "sum"))
agg["rho"] = agg.sxy/np.sqrt(agg.sxx*agg.syy)
agg = agg[agg.n >= 30]
tst = agg.rho*np.sqrt((agg.n-2)/(1-agg.rho**2))
agg["p"] = 2*tdist.sf(np.abs(tst), df=agg.n-2)

sig = agg.p < 0.05
for lbl, thr in [("as-shipped RHO_THRESHOLD=0.30", 0.30), ("legacy RHO_THRESHOLD=0.10", 0.10)]:
    conf = int((sig & (agg.rho >= thr)).sum()); contra = int((sig & (agg.rho <= -thr)).sum())
    tot = conf+contra
    print(f"{lbl:32s} confirmed={conf:4d}  contradicted={contra:4d}"
          + (f"  -> {100*contra/tot:.1f}% backwards" if tot else ""))
print(f"median rho on control set = {agg.rho.median():.4f}  (n_genes={len(agg)}, "
      f"{100*(agg.rho<0).mean():.1f}% negative)")

# Same again with the sign NOT inverted, to show which convention is right
m["chronos_pct_noinv"] = 1.0 - m["chronos_pct"]
m["r2b"] = m.groupby("ensg_id")["chronos_pct_noinv"].rank()
g2 = m.groupby("ensg_id"); m["dyb"] = m.r2b - g2.r2b.transform("mean")
agg2 = m.assign(dxdy=m.dx*m.dyb, dx2=m.dx**2, dy2=m.dyb**2).groupby("ensg_id").agg(
    n=("r1","size"), sxy=("dxdy","sum"), sxx=("dx2","sum"), syy=("dy2","sum"))
agg2["rho"] = agg2.sxy/np.sqrt(agg2.sxx*agg2.syy); agg2 = agg2[agg2.n >= 30]
print(f"median rho on control set WITHOUT the 1-rank inversion = {agg2.rho.median():.4f} "
      f"({100*(agg2.rho>0).mean():.1f}% positive)")
