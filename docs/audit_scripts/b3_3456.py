"""B3.3 proteomics coverage bias, B3.4 ploidy, B3.5 lineage sizes, B3.6 alteration prevalence."""
import os, re, json, numpy as np, pandas as pd
from scipy.stats import mannwhitneyu
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))
OUT = "src/pipeline/outputs"
R = {}

cm = pd.read_parquet(f"{OUT}/coverage_matrix_enriched.parquet")
print(f"coverage matrix: {len(cm):,} models")

# ---------------------------------------------------------------- B3.3
import pyarrow.parquet as pq
DEP = "data/parquet/2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.parquet"
gl = pd.read_parquet("reference/gene_lookup.parquet", columns=["ensg_id", "biotype", "hgnc_status"])
uni = set(gl[(gl.biotype == "protein_coding") & (gl.hgnc_status == "Approved")].ensg_id.str.lower())
cols = pq.ParquetFile(DEP).schema_arrow.names
c2e = {c: m.group(1).lower() for c in cols[1:] if (m := re.search(r"(ENSG\d+)", str(c)))}
keep = [cols[0]] + [c for c, e in c2e.items() if e in uni]
dep = pd.read_parquet(DEP, columns=keep).reset_index()
dep = dep.rename(columns={dep.columns[0]: "profileid"})
import ast
h = pd.read_parquet(f"{OUT}/harmonised_enriched.parquet", columns=["model_id", "profile_ids"])
p2m = {}
for mid, v in zip(h.model_id, h.profile_ids):
    try:
        for x in (ast.literal_eval(v) if isinstance(v, str) else list(v or [])):
            p2m[str(x).lower()] = mid
    except Exception:
        pass
dep["model_id"] = dep["profileid"].astype(str).str.lower().map(p2m)
dep = dep.dropna(subset=["model_id"])
mean_rna = dep.drop(columns=["profileid"]).groupby("model_id").mean().mean(axis=1).rename("mean_rna")
print(f"\nB3.3: mean log2(TPM+1) per line computed for {len(mean_rna):,} lines")

cov = cm.set_index("model_id")
res33 = {}
for plat, col in [("ProCan DIA", "procan_proteomics"), ("CCLE/Gygi TMT", "proteomics"),
                  ("either platform", None)]:
    has = (cov["procan_proteomics"] | cov["proteomics"]) if col is None else cov[col]
    j = pd.concat([mean_rna, has.rename("has")], axis=1).dropna()
    a, b = j[j.has].mean_rna, j[~j.has].mean_rna
    if len(a) < 5 or len(b) < 5:
        continue
    u, p = mannwhitneyu(a, b)
    cliff = 2 * u / (len(a) * len(b)) - 1
    d = (a.mean() - b.mean()) / np.sqrt(((len(a)-1)*a.var() + (len(b)-1)*b.var()) / (len(a)+len(b)-2))
    res33[plat] = {"n_with": len(a), "n_without": len(b),
                   "mean_with": round(a.mean(), 4), "mean_without": round(b.mean(), 4),
                   "median_with": round(a.median(), 4), "median_without": round(b.median(), 4),
                   "cliffs_delta": round(cliff, 4), "cohens_d": round(d, 4), "mwu_p": float(p)}
    print(f"  {plat:18s} with={len(a):4d} without={len(b):4d}  "
          f"mean {a.mean():.4f} vs {b.mean():.4f}  Cliff's d={cliff:+.4f}  "
          f"Cohen's d={d:+.4f}  MWU p={p:.3e}")
R["B3.3_proteomics_coverage_bias"] = res33

# ---------------------------------------------------------------- B3.4 ploidy
g = pd.read_parquet(f"{OUT}/gdsc_models.parquet",
                    columns=["model_id", "ploidy_wes", "ploidy_wgs", "tissue", "cancer_type"])
print("\nB3.4 PLOIDY (source: GDSC/Sanger model_list -> gdsc_models.parquet)")
pl = {}
for c in ["ploidy_wes", "ploidy_wgs"]:
    s = g[c].dropna()
    pl[c] = {"n": int(len(s)), "min": round(s.min(), 3), "q25": round(s.quantile(.25), 3),
             "median": round(s.median(), 3), "q75": round(s.quantile(.75), 3),
             "max": round(s.max(), 3), "frac_gt_2_5": round(float((s > 2.5).mean()), 4),
             "frac_within_1_9_2_1": round(float(s.between(1.9, 2.1).mean()), 4)}
    print(f"  {c}: n={len(s):,}  median={s.median():.3f}  IQR={s.quantile(.25):.3f}-{s.quantile(.75):.3f}"
          f"  range={s.min():.2f}-{s.max():.2f}")
    print(f"     frac >2.5 = {(s>2.5).mean():.1%}   frac near-diploid (1.9-2.1) = {s.between(1.9,2.1).mean():.1%}")
comb = g["ploidy_wgs"].fillna(g["ploidy_wes"]).dropna()
pl["combined_any"] = {"n": int(len(comb)), "median": round(comb.median(), 3),
                      "frac_near_diploid": round(float(comb.between(1.9, 2.1).mean()), 4)}
print(f"  any ploidy estimate: {len(comb):,} of {len(cm):,} panel models "
      f"({100*len(comb)/len(cm):.1f}%)  median={comb.median():.3f}")
R["B3.4_ploidy"] = pl

# ---------------------------------------------------------------- B3.5 lineage
print("\nB3.5 LINEAGE SIZE DISTRIBUTION (source: gdsc_models.tissue)")
t = g.dropna(subset=["tissue"]).groupby("tissue").size().sort_values(ascending=False)
print(t.to_string())
MINP = 15
big = t[t >= MINP]
print(f"\n  distinct lineages: {len(t)}")
print(f"  lineages with >= {MINP} lines: {len(big)} of {len(t)}")
print(f"  models in those lineages: {big.sum():,} of {t.sum():,} labelled ({100*big.sum()/t.sum():.1f}%)")
print(f"  models with NO tissue label: {len(cm) - t.sum():,} of {len(cm):,} panel models")
R["B3.5_lineage"] = {"n_lineages": int(len(t)), "min_peers": MINP,
                     "n_lineages_ge_min": int(len(big)),
                     "models_in_large_lineages": int(big.sum()),
                     "models_labelled": int(t.sum()), "panel_models": int(len(cm)),
                     "models_unlabelled": int(len(cm) - t.sum()),
                     "frac_panel_covered": round(float(big.sum()/len(cm)), 4),
                     "distribution": {k: int(v) for k, v in t.items()}}

# ---------------------------------------------------------------- B3.6 alterations
print("\nB3.6 ALTERATION PREVALENCE")
fl = pd.read_parquet(f"{OUT}/flags_with_driver.parquet",
                     columns=["model_id", "ensg_id", "p_mutation", "p_fusion",
                              "has_alteration", "has_driver_alteration", "has_cna_alteration",
                              "has_cna_context", "cna_type"])
core = pd.read_parquet(f"{OUT}/core_score.parquet", columns=["model_id", "ensg_id"])
core["model_id"] = core["model_id"].str.lower()
n_pairs = len(core)
n_genes = core.ensg_id.nunique(); n_lines = core.model_id.nunique()
print(f"  scored (gene,line) pairs in core_score: {n_pairs:,} "
      f"({n_genes:,} genes x {n_lines:,} lines)")
fl["model_id"] = fl["model_id"].str.lower()
key = core.merge(fl, on=["model_id", "ensg_id"], how="left")
a6 = {}
for label, mask in [("mutation (p>0.5)", key.p_mutation > 0.5),
                    ("fusion (p>0.5)", key.p_fusion > 0.5),
                    ("any alteration", key.has_alteration.fillna(False)),
                    ("driver alteration", key.has_driver_alteration.fillna(False)),
                    ("CNA alteration", key.has_cna_alteration.fillna(False))]:
    n = int(mask.sum())
    a6[label] = {"n_pairs": n, "frac_of_scored_pairs": round(n / n_pairs, 6)}
    print(f"  {label:22s} {n:>10,} pairs  = {100*n/n_pairs:6.3f}% of scored pairs")

pg = key.assign(alt=key.has_alteration.fillna(False)).groupby("ensg_id")["alt"].mean()
print(f"\n  per-gene fraction of lines altered: median={pg.median():.4f} "
      f"IQR={pg.quantile(.25):.4f}-{pg.quantile(.75):.4f} max={pg.max():.4f}")
print(f"  genes with >10% of lines altered: {int((pg>0.10).sum()):,} of {len(pg):,}")
print(f"  genes with  0% of lines altered: {int((pg==0).sum()):,} of {len(pg):,}")
a6["per_gene"] = {"median": round(float(pg.median()), 5), "q25": round(float(pg.quantile(.25)), 5),
                  "q75": round(float(pg.quantile(.75)), 5), "max": round(float(pg.max()), 5),
                  "n_genes_gt_10pct": int((pg > 0.10).sum()), "n_genes_zero": int((pg == 0).sum()),
                  "n_genes": int(len(pg))}
# CNA coverage caveat
print(f"\n  CNA context available on {int(key.has_cna_context.fillna(False).sum()):,} pairs "
      f"({100*key.has_cna_context.fillna(False).mean():.2f}%) -- absence elsewhere is NOT 'no CNA'")
a6["cna_context_pairs"] = int(key.has_cna_context.fillna(False).sum())
R["B3.6_alteration_prevalence"] = a6

json.dump(R, open(os.path.join(W, "b3_3456_results.json"), "w"), indent=1)
print("\nwritten b3_3456_results.json")
