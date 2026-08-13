"""B3.1 source correlation matrix + B3.2 line overlap.

Builds (model_id x ensg_id) matrices for DepMap, HPA and GEO from the RAW
sources, mapped through src/pipeline/outputs/harmonised_enriched.parquet.
Per-gene Spearman WITHIN gene ACROSS shared cell lines, summarised across genes.
"""
import os, json, ast, numpy as np, pandas as pd
import pyarrow.parquet as pq
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")
W = os.path.dirname(os.path.abspath(__file__))

h = pd.read_parquet("src/pipeline/outputs/harmonised_enriched.parquet",
                    columns=["model_id", "profile_ids", "geo_accessions", "hpa_cell_lines"])


def explode_map(col):
    """model_id -> list of source ids, as a source_id -> model_id dict."""
    out = {}
    for mid, v in zip(h.model_id, h[col]):
        if v is None:
            continue
        try:
            lst = ast.literal_eval(v) if isinstance(v, str) else list(v)
        except Exception:
            continue
        for x in lst:
            if x:
                out[str(x).lower()] = mid
    return out


prof2m = explode_map("profile_ids")
gsm2m = explode_map("geo_accessions")
hpa2m = explode_map("hpa_cell_lines")
print(f"maps: depmap profiles={len(prof2m):,}  geo gsm={len(gsm2m):,}  hpa names={len(hpa2m):,}")

# protein-coding universe, same filter the test_run_* scripts use
gl = pd.read_parquet("reference/gene_lookup.parquet", columns=["ensg_id", "biotype", "hgnc_status"])
uni = set(gl[(gl.biotype == "protein_coding") & (gl.hgnc_status == "Approved")].ensg_id.str.lower())
print(f"protein-coding Approved universe: {len(uni):,} genes")

# ---------------------------------------------------------------- DepMap
import re
DEP_PATH = "data/parquet/2_DepMap_OmicsExpressionAllGenesTPMLogp1Profile.parquet"
f = pq.ParquetFile(DEP_PATH)
dcols = f.schema_arrow.names
idcol = dcols[0]
# headers are 'SYMBOL (ENSG00000000003)' -- pull the ENSG out
col2ensg = {}
for c in dcols[1:]:
    m = re.search(r"(ENSG\d+)", str(c))
    if m:
        col2ensg[c] = m.group(1).lower()
keep = [idcol] + [c for c, e in col2ensg.items() if e in uni]
print(f"DepMap: {len(dcols):,} cols, {len(col2ensg):,} carry an ENSG "
      f"-> keeping {len(keep)-1:,} protein-coding")
dep = pd.read_parquet(DEP_PATH, columns=keep).reset_index()
if idcol not in dep.columns:
    dep = dep.rename(columns={dep.columns[0]: idcol})
dep[idcol] = dep[idcol].astype(str).str.lower()
dep["model_id"] = dep[idcol].map(prof2m)
print(f"  profiles {len(dep):,} -> mapped to model {dep.model_id.notna().sum():,}")
dep = dep.dropna(subset=["model_id"]).drop(columns=[idcol])
dep.columns = [col2ensg.get(c, c) for c in dep.columns]
DEP = dep.groupby("model_id").mean()
print(f"  DepMap matrix: {DEP.shape[0]:,} lines x {DEP.shape[1]:,} genes  (scale: log2(TPM+1))")
print(f"  value range {np.nanmin(DEP.values):.2f} .. {np.nanmax(DEP.values):.2f}")

# ---------------------------------------------------------------- HPA (long)
hpa = pd.read_parquet("data/parquet/1_4_hpa_rna_celline.parquet")
print(f"HPA raw cols: {list(hpa.columns)}")
gcol = [c for c in hpa.columns if c.lower() in ("gene",)][0]
ccol = [c for c in hpa.columns if "cell" in c.lower()][0]
lower2real = {c.lower(): c for c in hpa.columns}
vcol = lower2real.get("ntpm") or lower2real["tpm"]   # nTPM = HPA's cross-sample normalised value
hpa["model_id"] = hpa[ccol].astype(str).str.lower().map(hpa2m)
hpa[gcol] = hpa[gcol].astype(str).str.lower()
hpa = hpa[hpa.model_id.notna() & hpa[gcol].isin(uni)]
print(f"  HPA rows kept {len(hpa):,}; value column = {vcol}")
HPA = hpa.groupby(["model_id", gcol])[vcol].mean().unstack()
print(f"  HPA matrix: {HPA.shape[0]:,} lines x {HPA.shape[1]:,} genes  (scale: {vcol})")
print(f"  value range {np.nanmin(HPA.values):.2f} .. {np.nanmax(HPA.values):.2f}")
del hpa

# ---------------------------------------------------------------- GEO (genes x samples)
gf = pq.ParquetFile("data/parquet/3_GEOexpression.parquet")
gcols = gf.schema_arrow.names
gid = gcols[0]
gsm_keep = [c for c in gcols[1:] if c.lower() in gsm2m]
print(f"GEO: {len(gcols)-1:,} sample cols -> {len(gsm_keep):,} map to a model_id")
geo = pd.read_parquet("data/parquet/3_GEOexpression.parquet", columns=[gid] + gsm_keep)
geo[gid] = geo[gid].astype(str).str.lower()
geo = geo[geo[gid].isin(uni)].set_index(gid)
print(f"  GEO value range {np.nanmin(geo.values):.2f} .. {np.nanmax(geo.values):.2f}")
GEO = geo.T
GEO.index = [gsm2m[str(i).lower()] for i in GEO.index]
GEO = GEO.groupby(level=0).mean()
print(f"  GEO matrix: {GEO.shape[0]:,} lines x {GEO.shape[1]:,} genes")
del geo

# ---------------------------------------------------------------- B3.2 Venn
D, H, G = set(DEP.index), set(HPA.index), set(GEO.index)
venn = {
    "depmap_only": len(D - H - G), "hpa_only": len(H - D - G), "geo_only": len(G - D - H),
    "depmap_hpa_not_geo": len((D & H) - G), "depmap_geo_not_hpa": len((D & G) - H),
    "hpa_geo_not_depmap": len((H & G) - D), "all_three": len(D & H & G),
    "n_depmap": len(D), "n_hpa": len(H), "n_geo": len(G),
    "depmap_hpa": len(D & H), "depmap_geo": len(D & G), "hpa_geo": len(H & G),
    "union": len(D | H | G),
}
print("\n=== B3.2 LINE OVERLAP (Venn, model_id) ===")
for k, v in venn.items():
    print(f"  {k:22s} {v:5d}")

# ---------------------------------------------------------------- B3.1 correlations
def per_gene_spearman(A, B, nameA, nameB, min_lines=20):
    lines = sorted(set(A.index) & set(B.index))
    genes = sorted(set(A.columns) & set(B.columns))
    if len(lines) < min_lines or not genes:
        return None
    a = A.loc[lines, genes]; b = B.loc[lines, genes]
    ra = a.rank(axis=0, na_option="keep"); rb = b.rank(axis=0, na_option="keep")
    mask = ra.notna() & rb.notna()
    ra = ra.where(mask); rb = rb.where(mask)
    ra = ra - ra.mean(); rb = rb - rb.mean()
    num = (ra * rb).sum()
    den = np.sqrt((ra ** 2).sum() * (rb ** 2).sum())
    rho = (num / den).replace([np.inf, -np.inf], np.nan)
    n = mask.sum()
    rho = rho[n >= min_lines]
    return {"pair": f"{nameA}-{nameB}", "n_lines_shared": len(lines),
            "n_genes": int(len(rho)), "median": float(rho.median()),
            "q25": float(rho.quantile(.25)), "q75": float(rho.quantile(.75)),
            "mean": float(rho.mean()), "frac_gt_0": float((rho > 0).mean()),
            "_rho": rho}


res = []
for A, B, na, nb in [(DEP, HPA, "DepMap", "HPA"), (DEP, GEO, "DepMap", "GEO"), (HPA, GEO, "HPA", "GEO")]:
    r = per_gene_spearman(A, B, na, nb)
    if r:
        res.append(r)
        print(f"\n{r['pair']}: shared lines={r['n_lines_shared']}, genes={r['n_genes']:,}")
        print(f"  Spearman rho within-gene across-lines: median={r['median']:.3f} "
              f"IQR={r['q25']:.3f}-{r['q75']:.3f}  mean={r['mean']:.3f}  "
              f"frac>0={r['frac_gt_0']:.3f}")

print("\n=== B3.1 3x3 SOURCE CORRELATION MATRIX (median per-gene Spearman) ===")
M = pd.DataFrame(1.0, index=["DepMap", "HPA", "GEO"], columns=["DepMap", "HPA", "GEO"])
for r in res:
    a, b = r["pair"].split("-")
    M.loc[a, b] = M.loc[b, a] = round(r["median"], 3)
print(M.to_string())

json.dump({"venn": venn,
           "pairs": [{k: v for k, v in r.items() if k != "_rho"} for r in res],
           "matrix": M.to_dict()},
          open(os.path.join(W, "b3_1_2_results.json"), "w"), indent=1)
for r in res:
    r["_rho"].rename("rho").to_frame().to_parquet(
        os.path.join(W, f"b31_rho_{r['pair'].replace('-', '_')}.parquet"))
print("\nwritten b3_1_2_results.json")
