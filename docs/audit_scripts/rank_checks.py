"""Shared builder for the five candidate confirmation checks.

Each check is THREE-STATE per (gene, line):
    1.0  = passed        0.0 = failed        NaN = not assessed (no data)

Nothing here is production code; it exists so Part A can be measured.
"""
import os, numpy as np, pandas as pd, duckdb
import pyarrow.parquet as pq
os.chdir(r"C:\Disertation\UoB-GeneTraceAI-25-26")

DB = "src/pipeline/outputs/celllineselector.db"
OUT = "src/pipeline/outputs"
EXPRESSED_MIN = 1.0          # log2(TPM+1) > 1  <=>  TPM > 1   (HPA detection convention)
GEO_TOL = 0.25               # percentile-units tolerance for GEO corroboration
SILENT_FRAC = 0.20           # validity: gene expressed in >=20% of the panel

CHECKS = ["rna_gate", "cna_gate", "protein", "alteration", "geo_corrob"]
# gene-independent evidence-type availability per line, used for the A1 confound
COVERAGE_COLS = {
    "rna_gate":   ["depmap_expr", "hpa_rna"],
    "cna_gate":   ["cosmic_cna"],
    "protein":    ["procan_proteomics", "proteomics"],
    "alteration": ["mutations", "fusions"],
    "geo_corrob": ["geo_expr"],
}


def connect():
    return duckdb.connect(DB, read_only=True)


def sample_genes(con, n=300, seed=42):
    dep = set(x for x in [r[0] for r in con.execute('DESCRIBE depmap_expr').fetchall()]
              if x.startswith("ensg"))
    gl = pd.read_parquet("reference/gene_lookup.parquet",
                         columns=["ensg_id", "hgnc_symbol", "biotype", "hgnc_status"])
    gl["ensg_id"] = gl.ensg_id.str.lower()
    uni = gl[(gl.biotype == "protein_coding") & (gl.hgnc_status == "Approved")]
    cand = sorted(set(uni.ensg_id) & dep)
    rng = np.random.default_rng(seed)
    genes = sorted(rng.choice(cand, size=min(n, len(cand)), replace=False))
    sym = dict(zip(uni.ensg_id, uni.hgnc_symbol.astype(str).str.lower()))
    return genes, {g: sym.get(g) for g in genes}


def build(con, genes, sym):
    """Returns dict of (line x gene) DataFrames, one per check, three-state."""
    qg = ",".join(f'"{g}"' for g in genes)

    # ---- RNA: DepMap + HPA consensus on log2(TPM+1)
    D = con.execute(f"SELECT model_id,{qg} FROM depmap_expr").df()
    D = D.set_index("model_id").groupby(level=0).median()
    hp = con.execute(
        f"SELECT model_id, lower(gene) AS gene, avg(ntpm) AS v FROM hpa_rna "
        f"WHERE lower(gene) IN ({','.join(chr(39)+g+chr(39) for g in genes)}) "
        f"GROUP BY 1,2").df()
    H = np.log2(hp.pivot(index="model_id", columns="gene", values="v").clip(lower=0) + 1)
    lines = sorted(l for l in (set(D.index) | set(H.index)) if isinstance(l, str))
    D = D.reindex(index=lines, columns=genes)
    H = H.reindex(index=lines, columns=genes)
    RNA = pd.concat([D, H]).groupby(level=0).mean().reindex(index=lines, columns=genes)

    # ---- GEO
    geo_have = [g for g in genes
                if g in {r[0] for r in con.execute("DESCRIBE geo_expr").fetchall()}]
    if geo_have:
        G = con.execute(
            f"SELECT model_id,{','.join(f'median(\"{g}\") AS \"{g}\"' for g in geo_have)} "
            f"FROM geo_expr WHERE model_id IS NOT NULL GROUP BY model_id").df()
        G = np.log2(G.set_index("model_id").clip(lower=0) + 1).reindex(columns=genes)
    else:
        G = pd.DataFrame(np.nan, index=lines, columns=genes)
    print(f"  GEO covers {len(geo_have)} of {len(genes)} sampled genes")

    # ---- Protein: ProCan, uniprot -> symbol, mapped to model_id via gdsc_models
    pm = con.execute("SELECT uniprot_id, gene_symbol FROM procan_protein_map").df()
    s2u = {}
    for u, s in zip(pm.uniprot_id, pm.gene_symbol):
        s2u.setdefault(str(s).lower(), []).append(str(u).lower())
    pro = pd.read_parquet(f"{OUT}/procan_proteomics.parquet")
    gm = pd.read_parquet(f"{OUT}/gdsc_models.parquet",
                         columns=["sanger_model_id", "model_id"]).dropna()
    s2m = dict(zip(gm.sanger_model_id, gm.model_id))
    pro["model_id"] = pro["sanger_model_id"].map(s2m)
    pro = pro.dropna(subset=["model_id"]).set_index("model_id")
    P = pd.DataFrame(np.nan, index=pro.index, columns=genes)
    for g in genes:
        for u in s2u.get(sym.get(g) or "", []):
            if u in pro.columns:
                P[g] = pro[u].astype(float)
                break
    P = P.groupby(level=0).mean()

    # ---- CNA + alteration flags
    cna = pd.read_parquet(f"{OUT}/cna_flags.parquet")
    cna = cna[cna.ensg_id.isin(genes)]
    fl = pd.read_parquet(f"{OUT}/flags_with_driver.parquet",
                         columns=["model_id", "ensg_id", "has_alteration"])
    fl = fl[fl.ensg_id.isin(genes)]

    cov = pd.read_parquet(f"{OUT}/coverage_matrix_enriched.parquet").set_index("model_id")

    idx = pd.Index(lines, name="model_id")

    def wide(df, col, fill=np.nan):
        w = df.pivot_table(index="model_id", columns="ensg_id", values=col, aggfunc="first")
        return w.reindex(index=idx, columns=genes).astype(float)

    cna_ctx = wide(cna.assign(v=cna.has_cna_context.astype(float)), "v")
    cna_loss = wide(cna.assign(v=cna.cna_type.isin(["loss", "both"]).astype(float)), "v")
    alt = wide(fl.assign(v=fl.has_alteration.astype(float)), "v")

    # =============== the five checks, three-state ===============
    C = {}

    # 1. RNA gate: assessed where RNA covers; pass where expressed
    C["rna_gate"] = pd.DataFrame(np.where(RNA.notna(), (RNA > EXPRESSED_MIN).astype(float), np.nan),
                                 index=idx, columns=genes)

    # 2. CNA gate: assessed where COSMIC CNA context exists; FAIL where deleted
    C["cna_gate"] = pd.DataFrame(np.where(cna_ctx.fillna(0) > 0, 1.0 - cna_loss.fillna(0), np.nan),
                                 index=idx, columns=genes)

    # 3. Protein: assessed where quantified; pass above the protein's panel median
    pmed = P.median(axis=0)
    C["protein"] = pd.DataFrame(np.where(P.notna(), (P > pmed).astype(float), np.nan),
                                index=P.index, columns=genes).reindex(index=idx)

    # 4. Alteration: assessed where the LINE has mutation/fusion coverage.
    #    Default direction = alterations DEMOTE, so unaltered passes. See D1.
    has_mut_cov = cov.reindex(idx)[["mutations", "fusions"]].any(axis=1)
    a = alt.fillna(0.0)
    C["alteration"] = pd.DataFrame(
        np.where(has_mut_cov.to_numpy()[:, None], 1.0 - a.to_numpy(), np.nan),
        index=idx, columns=genes)

    # 5. GEO corroboration: assessed where GEO covers; pass where its within-gene
    #    percentile agrees with the RNA consensus percentile within GEO_TOL
    rp = RNA.rank(axis=0, pct=True, na_option="keep")
    gp = G.reindex(index=idx, columns=genes).rank(axis=0, pct=True, na_option="keep")
    both = rp.notna() & gp.notna()
    C["geo_corrob"] = pd.DataFrame(
        np.where(both, ((rp - gp).abs() <= GEO_TOL).astype(float), np.nan),
        index=idx, columns=genes)

    # ---- gene-independent coverage count per line (A1)
    covcount = pd.Series(0, index=idx, dtype=float)
    for chk, cols in COVERAGE_COLS.items():
        have = cov.reindex(idx)[cols].any(axis=1).fillna(False)
        covcount += have.astype(float)

    # ---- validity: gene expressed in >= SILENT_FRAC of the panel (F2 guard)
    frac_expr = (RNA > EXPRESSED_MIN).sum(axis=0) / RNA.notna().sum(axis=0)
    valid_gene = (frac_expr >= SILENT_FRAC)

    return {"checks": C, "RNA": RNA, "cov_count": covcount, "lines": idx,
            "genes": genes, "valid_gene": valid_gene, "frac_expr": frac_expr,
            "cov": cov.reindex(idx)}


def confirmation_count(C, genes, idx):
    """n_pass (count of checks passed) and n_assessed (count with data)."""
    npass = pd.DataFrame(0.0, index=idx, columns=genes)
    nass = pd.DataFrame(0.0, index=idx, columns=genes)
    for chk in CHECKS:
        v = C[chk].reindex(index=idx, columns=genes)
        npass += v.fillna(0.0)
        nass += v.notna().astype(float)
    return npass, nass
