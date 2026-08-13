"""PART D: the disagreement (I^2) flag.  PART F2: the silence guard applied to
the MAIN pipeline's own pooled scores."""
import os, json, numpy as np, pandas as pd
import pyarrow.parquet as pq
from scipy import stats
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))
import tx_lib as T

# =============================== PART D ========================================
GENES = ["TSPAN6", "EGFR", "TP53", "MYC", "BRAF", "PTEN", "KRAS", "CDKN2A"]
con = T.connect()
rows, alld = [], []
for g in GENES:
    try:
        df, gid = T.fetch_gene_expression(con, g, verbose=False)
    except ValueError:
        print(f"  {g}: not found"); continue
    res, per = T.score_gene_lineage(df, verbose=False)
    res["gene"] = g
    alld.append(res)
    conf = res[res.verdict == "Conflicting Sources"]
    scoreable = res[res.n_sources >= 2]
    rows.append({
        "gene": g, "n_lines": len(res),
        "n_ge2_sources": int(len(scoreable)),
        "n_conflicting": int(len(conf)),
        "pct_of_scoreable": round(100 * len(conf) / max(len(scoreable), 1), 2),
        "n_conf_with_2_sources": int((conf.n_sources == 2).sum()),
        "pct_conf_on_2_sources": round(100 * (conf.n_sources == 2).mean(), 1) if len(conf) else np.nan,
        "n_conf_with_3_sources": int((conf.n_sources == 3).sum()),
    })
D = pd.DataFrame(rows)
print("=== PART D: 'Conflicting Sources' verdict incidence ===")
print(D.to_string(index=False))
A = pd.concat(alld, ignore_index=True)
conf_all = A[A.verdict == "Conflicting Sources"]
sc_all = A[A.n_sources >= 2]
print(f"\n  pooled over {len(GENES)} genes: {len(conf_all):,} conflicting of "
      f"{len(sc_all):,} scoreable ({100*len(conf_all)/len(sc_all):.2f}%)")
print(f"  of those, on only TWO sources: {int((conf_all.n_sources==2).sum()):,} "
      f"({100*(conf_all.n_sources==2).mean():.1f}%)")
print(f"                    THREE sources: {int((conf_all.n_sources==3).sum()):,} "
      f"({100*(conf_all.n_sources==3).mean():.1f}%)")

# --- how biased is I^2 at k=2 and k=3? simulate under NO true heterogeneity ----
print("\n=== D: I^2 under NO true heterogeneity (k=2, k=3), 200,000 draws ===")
rng = np.random.default_rng(42)
for k in (2, 3, 5, 7):
    z = rng.standard_normal((200000, k))
    Q = ((z - z.mean(axis=1, keepdims=True)) ** 2).sum(axis=1)
    I2 = np.where(Q > 0, np.maximum(0, (Q - (k - 1)) / Q) * 100, 0.0)
    print(f"  k={k}: mean I2 = {I2.mean():5.1f}%   median = {np.median(I2):5.1f}%   "
          f"P(I2 > {T.I2_CUT}) = {(I2 > T.I2_CUT).mean():.4f}")
print(f"  -> with NO real disagreement, the I2>{T.I2_CUT} flag fires on a large")
print("     fraction of lines. At k=2 the statistic has 1 degree of freedom.")

# =============================== PART F2 =======================================
print("\n\n=== PART F2: silence guard applied to the MAIN pipeline's pooled scores ===")
DATA_CLEAN = "data/parquet/data_clean"
names = pq.ParquetFile(f"{DATA_CLEAN}/depmap_expr_clean.parquet").schema_arrow.names
col = {c.split(".")[0].lower(): c for c in names if c != "index"}
gl = pd.read_parquet("reference/gene_lookup.parquet",
                     columns=["ensg_id", "biotype", "hgnc_status"])
uni = set(gl[(gl.biotype == "protein_coding") & (gl.hgnc_status == "Approved")]
          .ensg_id.str.lower())
genes = sorted(g for g in col if g in uni)
print(f"  protein-coding genes with DepMap expression: {len(genes):,}")

prof = pd.read_parquet("reference/depmap_profiles.parquet")
rp = prof[prof.datatype == "rna"][["profileid", "modelid"]].copy()
rp["model_id"] = rp.modelid.str.lower(); rp["profileid"] = rp.profileid.astype(str)

CH = 2000
frac_expr, n_lines_tot = {}, None
for i in range(0, len(genes), CH):
    blk = genes[i:i + CH]
    E = pq.read_table(f"{DATA_CLEAN}/depmap_expr_clean.parquet",
                      columns=["index"] + [col[g] for g in blk]).to_pandas()
    if "index" not in E.columns:
        E = E.reset_index()
    E["profileid"] = E["index"].astype(str)
    E = (E.drop(columns=["index"]).merge(rp[["profileid", "model_id"]],
                                         on="profileid", how="inner")
           .drop(columns=["profileid"]).groupby("model_id").mean())
    E.columns = [c.split(".")[0].lower() for c in E.columns]
    n_lines_tot = len(E)
    fe = (E > T.EXPRESSED_MIN).sum() / E.notna().sum()
    frac_expr.update(fe.to_dict())
    print(f"    ...{min(i+CH, len(genes)):,}/{len(genes):,}", end="\r")
fe = pd.Series(frac_expr).dropna()
print(f"\n  panel lines: {n_lines_tot:,}")
print(f"  fraction of lines expressing each gene (log2(TPM+1) > {T.EXPRESSED_MIN}):")
print(f"    median {fe.median():.3f}   IQR {fe.quantile(.25):.3f}-{fe.quantile(.75):.3f}")

silent = fe[fe < T.SILENT_FRAC]
print(f"\n  GENES failing the silence guard (expressed in < {100*T.SILENT_FRAC:.0f}% of the panel):")
print(f"    {len(silent):,} of {len(fe):,} genes ({100*len(silent)/len(fe):.1f}%)")
print(f"    gene x line combinations affected: {len(silent)*n_lines_tot:,} "
      f"of {len(fe)*n_lines_tot:,} ({100*len(silent)/len(fe):.1f}%)")
for cut in (0.05, 0.10, 0.20, 0.30):
    s = fe[fe < cut]
    print(f"      at SILENT_FRAC={cut:.2f}: {len(s):,} genes  "
          f"({100*len(s)/len(fe):.1f}%),  {len(s)*n_lines_tot:,} gene x line pairs")

# does the pipeline's own dispersion table already catch these?
gd = pd.read_parquet("src/pipeline/outputs/gene_dispersion.parquet")
gd["ensg_id"] = gd.ensg_id.str.lower()
j = gd.set_index("ensg_id").join(fe.rename("frac_expressed"), how="inner")
print(f"\n  cross-check against gene_dispersion.parquet ({len(j):,} genes matched):")
print(f"    of the {int((j.frac_expressed < T.SILENT_FRAC).sum()):,} silence-guard "
      f"failures, signal_spread says:")
print(j[j.frac_expressed < T.SILENT_FRAC].signal_spread.value_counts().to_string())
print(f"\n    MAD-floor-style check alone -- dynamic_range on those genes:")
sub = j[j.frac_expressed < T.SILENT_FRAC]
print(f"      median dynamic_range {sub.dynamic_range.median():.3f}  "
      f"IQR {sub.dynamic_range.quantile(.25):.3f}-{sub.dynamic_range.quantile(.75):.3f}")
print(f"      of these, {int((sub.dynamic_range > 3.0).sum()):,} exceed "
      f"FLAT_FOLD_CEILING=3.0 and so would NOT be caught by a spread floor")

json.dump({"part_d": D.to_dict("records"),
           "pooled_conflicting": int(len(conf_all)),
           "pooled_scoreable": int(len(sc_all)),
           "pct_conf_two_sources": float(100 * (conf_all.n_sources == 2).mean()),
           "f2_n_genes": int(len(fe)), "f2_panel_lines": int(n_lines_tot),
           "f2_silent_genes": int(len(silent)),
           "f2_silent_pairs": int(len(silent) * n_lines_tot),
           "f2_by_cut": {str(c): int((fe < c).sum()) for c in (0.05, 0.10, 0.20, 0.30)}},
          open(os.path.join(W, "d_f2_results.json"), "w"), indent=1)
fe.rename("frac_expressed").to_frame().to_parquet(os.path.join(W, "f2_frac_expressed.parquet"))
print("\nwritten d_f2_results.json")
