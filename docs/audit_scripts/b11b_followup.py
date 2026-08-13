"""B1.1 follow-up: separate the two questions the control test conflates.

Q1 (polarity, BETWEEN genes): after the 1-rank transform, do essential genes
    carry higher chronos_pct than non-essential genes?  This is the actual
    sign test.
Q2 (WITHIN gene): the pipeline's rho is a within-gene rank correlation across
    cell lines. Characterise it on pan-essentials vs a matched random set, and
    test the proliferation/screen-quality confound (per-LINE mean essentiality).
Q3 what the SHIPPED chronos_validation.parquet actually contains now.
"""
import os, numpy as np, pandas as pd
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
from scipy.stats import mannwhitneyu, spearmanr

gl = pd.read_parquet("reference/gene_lookup.parquet", columns=["ensg_id", "hgnc_symbol"])
gl["ensg_id"] = gl["ensg_id"].str.upper()
sym2ensg = dict(zip(gl.hgnc_symbol, gl.ensg_id))
ess = pd.read_csv("data/DepMap_Chronos/ReferenceEssentials.csv")
non = pd.read_csv("data/DepMap_Chronos/ReferenceNonEssentials.csv")
f = lambda d: {s.split(" ")[0].strip() for s in d[d.columns[0]].astype(str)}
ess_e = {sym2ensg[s] for s in f(ess) if s in sym2ensg}
non_e = {sym2ensg[s] for s in f(non) if s in sym2ensg}

ch = pd.read_parquet("validation/prepared/chronos_long.parquet")
ch["ensg_id"] = ch["ensg_id"].str.upper()
ch["chronos_pct"] = 1.0 - ch.groupby("ensg_id")["essentiality"].rank(pct=True, method="average")

print("=== Q1: BETWEEN-gene polarity after the 1-rank transform ===")
print("NOTE: chronos_pct is ranked WITHIN gene, so its mean is ~0.5 for every")
print("gene by construction -- it carries no between-gene information at all.")
gm = ch.groupby("ensg_id")["chronos_pct"].mean()
print(f"  mean chronos_pct, essentials     = {gm[gm.index.isin(ess_e)].mean():.4f}")
print(f"  mean chronos_pct, non-essentials = {gm[gm.index.isin(non_e)].mean():.4f}")
raw = ch.groupby("ensg_id")["essentiality"].median()
u, p = mannwhitneyu(raw[raw.index.isin(ess_e)].dropna(), raw[raw.index.isin(non_e)].dropna())
print(f"  RAW essentiality  essentials median={raw[raw.index.isin(ess_e)].median():.3f} "
      f"vs non-essentials={raw[raw.index.isin(non_e)].median():.3f}  MWU p={p:.3e}")

print("\n=== Q2: the within-gene confound ===")
# per-line mean essentiality across pan-essential genes = screen strength / proliferation proxy
sub = ch[ch.ensg_id.isin(ess_e)]
line_q = sub.groupby("sanger_model_id")["essentiality"].mean()
print(f"  per-line mean essentiality over {len(ess_e)} reference essentials:")
print(f"    min={line_q.min():.2f} median={line_q.median():.2f} max={line_q.max():.2f} "
      f"IQR={line_q.quantile(.25):.2f}-{line_q.quantile(.75):.2f}  n_lines={len(line_q):,}")

ib = pd.read_parquet("validation/prepared/id_bridge.parquet")
lq = line_q.rename("line_screen_strength").reset_index().merge(ib, on="sanger_model_id")
lq["model_id"] = lq["model_id"].str.lower()

pred = pd.read_parquet("src/pipeline/outputs/predictions_with_confidence.parquet",
                       columns=["model_id", "ensg_id", "core_score"])
pred["model_id"] = pred["model_id"].str.lower()
line_expr = pred.groupby("model_id")["core_score"].mean().rename("line_mean_core")
j = lq.merge(line_expr, on="model_id", how="inner")
r, pv = spearmanr(j.line_screen_strength, j.line_mean_core)
print(f"  Spearman(per-line screen strength, per-line mean core_score) = {r:.3f} "
      f"p={pv:.2e}  n_lines={len(j)}")
print("  -> if materially non-zero, the within-gene rho is partly a LINE-level")
print("     artefact, not a gene-level abundance->dependency relationship.")

print("\n=== Q3: what chronos_validation.parquet actually ships ===")
cv = pd.read_parquet("src/pipeline/outputs/chronos_validation.parquet")
print(cv.chronos_check.value_counts().to_string())
print(cv.groupby("chronos_check")["rho"].agg(["count", "min", "median", "max"]).round(3).to_string())
print(f"  total genes classified: {len(cv):,}")
