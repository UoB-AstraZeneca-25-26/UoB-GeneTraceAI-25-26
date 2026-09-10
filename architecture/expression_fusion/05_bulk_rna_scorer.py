"""
expression_fusion/05_bulk_rna_scorer.py
-----------------------------------------
Bulk-computes lineage-conditioned RNA z-scores for every gene in the universe,
using all three sources (DepMap + HPA + warehouse GEO), then Stouffer-combines
them into a single rna_z_t per (gene, model_id).

OUTPUT
------
expression_fusion/outputs/bulk_rna_z.parquet
    columns: gene_id, model_id, rna_z_t, n_sources
    one row per (gene, model_id) with a valid combined z-score.

HOW IT FEEDS INTO 02_core_score.ipynb
--------------------------------------
    rna_z = pd.read_parquet("expression_fusion/outputs/bulk_rna_z.parquet")
    # convert lineage-conditioned z to global percentile per gene (preserves
    # the lineage-relative ordering the z-score encodes)
    expr_pct = (rna_z.pivot(index="model_id", columns="gene_id", values="rna_z_t")
                      .rank(pct=True, method="min"))
    # use expr_pct in place of the DepMap-only percentile in Noisy-OR

SCALE DECISIONS
---------------
- DepMap:       already log2(TPM+1). Left as-is.
- HPA nTPM:     linear. log2(x+1) applied.
- GEO warehouse: linear (median ~489, max ~8239, confirmed by diagnostic).
                 log2(x+1) applied. Same source the shared notebook validated
                 against; warehouse geo_expr covers 590 model_ids vs 484 for
                 geo_expr_v2, so we use the warehouse for broader coverage here.

SPEED STRATEGY
--------------
Wide-table column fetching in per-gene batches costs ~60s per 200 genes (100+
sequential SQL round-trips). Instead each source is loaded ONCE in full, held
in memory, and the per-gene z-score is computed via vectorized numpy. Peak RAM:
DepMap ~250MB wide + HPA ~180MB long + GEO ~200MB wide = ~630MB total.
"""
from __future__ import annotations

import math
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats as scipy_stats

import common as C

# ============================================================== constants
MAD_FLOOR  = 0.384   # matches shared notebook's calibrated default
MIN_PEERS  = 5       # minimum lines per (source, lineage) stratum to score
N_SOURCES  = 3       # DepMap, HPA, GEO -- used in shrinkage denominator


def _robust_z_matrix(X: np.ndarray, min_peers=MIN_PEERS, mad_floor=MAD_FLOOR) -> np.ndarray:
    """
    Vectorised robust z-score across rows for each column of X.
    X shape: (n_lines, n_genes). Returns same shape; columns with < min_peers
    non-NaN rows get all-NaN.
    """
    n_valid = np.sum(np.isfinite(X), axis=0)            # (n_genes,)
    med     = np.nanmedian(X, axis=0)                   # (n_genes,)
    abs_dev = np.abs(X - med)
    mad     = np.nanmedian(abs_dev, axis=0) * 1.4826    # scale='normal'
    mad     = np.maximum(mad, mad_floor)
    Z       = (X - med) / mad                           # (n_lines, n_genes)
    Z[:, n_valid < min_peers] = np.nan
    return Z


def _score_source(wide: pd.DataFrame,
                  gene_ids: list[str],
                  lineage_map: pd.Series,
                  source_name: str) -> pd.DataFrame:
    """
    wide   : model_id (index) × gene_id (columns), values already log2-scaled.
    Returns long DataFrame: gene_id, model_id, source, z
    """
    wide = wide.loc[wide.index.isin(lineage_map.index), gene_ids].copy()
    wide["lineage"] = lineage_map.reindex(wide.index)
    wide = wide.dropna(subset=["lineage"])

    rows = []
    for lin, grp in wide.groupby("lineage"):
        X = grp[gene_ids].values.astype(float)       # (n_lines, n_genes)
        Z = _robust_z_matrix(X)                      # same shape
        model_ids = grp.index.tolist()
        z_wide = pd.DataFrame(Z, index=model_ids, columns=gene_ids)
        # future_stack=True suppresses FutureWarning in pandas ≥2.1
        try:
            z_long = z_wide.stack(future_stack=True).dropna().reset_index()
        except TypeError:
            z_long = z_wide.stack().dropna().reset_index()
        z_long.columns = ["model_id", "gene_id", "z"]
        rows.append(z_long)

    if not rows:
        return pd.DataFrame(columns=["gene_id", "model_id", "source", "z"])
    out = pd.concat(rows, ignore_index=True)
    out["source"] = source_name
    return out


def _stouffer_combine(z_long: pd.DataFrame) -> pd.DataFrame:
    """
    Vectorised Stouffer combination per (gene_id, model_id).
    Weight  = sqrt(n_lines this source covers for this gene).
    Shrinkage: divide z_raw by sqrt(N_SOURCES / n_available) when n_available < N_SOURCES.
    """
    valid = z_long.dropna(subset=["z"]).copy()

    # per (gene, source): count distinct model_ids that have a valid z
    src_n = (valid.groupby(["gene_id", "source"])["model_id"]
             .nunique().reset_index(name="src_n"))
    valid = valid.merge(src_n, on=["gene_id", "source"], how="left")
    valid["w"]  = np.sqrt(valid["src_n"])
    valid["wz"] = valid["w"] * valid["z"]
    valid["w2"] = valid["w"] ** 2

    agg = valid.groupby(["gene_id", "model_id"]).agg(
        wz_sum    = ("wz", "sum"),
        w2_sum    = ("w2", "sum"),
        n_sources = ("source", "nunique"),
    ).reset_index()

    denom = np.sqrt(agg["w2_sum"].values)
    safe  = denom > 0
    z_raw = np.where(safe, agg["wz_sum"].values / np.where(safe, denom, 1.0), np.nan)
    shrink = np.where(agg["n_sources"] < N_SOURCES,
                      np.sqrt(N_SOURCES / agg["n_sources"]), 1.0)
    agg["rna_z_t"]  = z_raw / shrink
    agg["n_sources"] = agg["n_sources"].astype("int8")
    return agg[["gene_id", "model_id", "rna_z_t", "n_sources"]]


# ============================================================== source loaders
def _load_depmap(con) -> tuple[pd.DataFrame, list[str]]:
    print("Loading DepMap wide table...", flush=True)
    t0 = time.time()
    prof = con.execute(
        "SELECT profileid, model_id FROM main.depmap_profiles WHERE datatype = 'rna'"
    ).df().drop_duplicates("profileid").set_index("profileid")["model_id"]

    wide = con.execute("SELECT * FROM main.depmap_expr").df()
    wide["model_id"] = wide["profileid"].map(prof)
    wide = wide.dropna(subset=["model_id"]).drop(columns=["profileid"])
    wide = wide.groupby("model_id").median()  # collapse dup profiles

    gene_ids = [c for c in wide.columns]
    print(f"  DepMap: {wide.shape} ({time.time()-t0:.0f}s)")
    return wide, gene_ids


def _load_hpa(con, gene_ids: set) -> pd.DataFrame:
    print("Loading HPA...", flush=True)
    t0 = time.time()
    long = con.execute(
        "SELECT gene AS gene_id, model_id, ntpm AS value "
        "FROM main.hpa_rna WHERE is_ambiguous = FALSE"
    ).df().dropna()
    long = long[long.gene_id.isin(gene_ids)]
    long["value"] = np.log2(np.clip(long["value"].values, 0, None) + 1)
    long = long.groupby(["gene_id", "model_id"], as_index=False)["value"].mean()
    wide = long.pivot(index="model_id", columns="gene_id", values="value")
    print(f"  HPA: {wide.shape} ({time.time()-t0:.0f}s)")
    return wide


def _load_geo(con, gene_ids: set) -> pd.DataFrame:
    print("Loading warehouse GEO (linear -> log2)...", flush=True)
    t0 = time.time()
    geo_gene_cols = (set(con.execute("DESCRIBE main.geo_expr").df()["column_name"])
                     - {"sample", "model_id", "is_ambiguous"})
    available = list(geo_gene_cols & gene_ids)
    if not available:
        return pd.DataFrame()
    cols = ", ".join(f'"{g}"' for g in available)
    wide = con.execute(
        f'SELECT model_id, {cols} FROM main.geo_expr '
        f'WHERE model_id IS NOT NULL AND is_ambiguous = FALSE'
    ).df()
    wide = wide.groupby("model_id").median()   # index = model_id, columns = gene_ids
    # warehouse values are linear; apply log2(x+1) globally
    wide = wide.apply(lambda col: np.log2(np.clip(col.values, 0, None) + 1), axis=0)
    print(f"  GEO:   {wide.shape} ({time.time()-t0:.0f}s)")
    return wide


# ============================================================== main
def main():
    C.banner("Bulk RNA z-scorer: DepMap + HPA + warehouse GEO -> rna_z_t")

    con = C.connect()
    lineage_map = C.load_lineage(con)
    print(f"Lineage map: {len(lineage_map):,} model_ids, {lineage_map.nunique()} lineages\n")

    # --- load all sources ---
    dm_wide, dm_genes = _load_depmap(con)

    gene_set = set(dm_genes)
    hpa_wide = _load_hpa(con, gene_set)
    geo_wide  = _load_geo(con, gene_set)

    # union of gene ids across all sources
    all_genes = sorted(
        set(dm_wide.columns)
        | set(hpa_wide.columns if not hpa_wide.empty else [])
        | set(geo_wide.columns if not geo_wide.empty else [])
    )
    print(f"\nTotal genes in union: {len(all_genes):,}")

    # --- per-source lineage z-scores ---
    print("\nComputing DepMap z-scores per lineage...", flush=True)
    t0 = time.time()
    dm_genes_list = [g for g in all_genes if g in dm_wide.columns]
    dm_z = _score_source(dm_wide, dm_genes_list, lineage_map, "depmap_expr")
    print(f"  {len(dm_z):,} rows ({time.time()-t0:.0f}s)")

    print("Computing HPA z-scores per lineage...", flush=True)
    t0 = time.time()
    hpa_genes_list = [g for g in all_genes if g in hpa_wide.columns]
    hpa_z = _score_source(hpa_wide, hpa_genes_list, lineage_map, "hpa_rna")
    print(f"  {len(hpa_z):,} rows ({time.time()-t0:.0f}s)")

    print("Computing GEO z-scores per lineage...", flush=True)
    t0 = time.time()
    if not geo_wide.empty:
        geo_genes_list = [g for g in all_genes if g in geo_wide.columns]
        geo_z = _score_source(geo_wide, geo_genes_list, lineage_map, "geo_expr")
        print(f"  {len(geo_z):,} rows ({time.time()-t0:.0f}s)")
    else:
        geo_z = pd.DataFrame(columns=["gene_id", "model_id", "source", "z"])

    # --- combine ---
    print("\nStoreing Stouffer combination...", flush=True)
    t0 = time.time()
    all_z = pd.concat([dm_z, hpa_z, geo_z], ignore_index=True)
    combined = _stouffer_combine(all_z)
    combined["rna_z_t"] = combined["rna_z_t"].astype("float32")
    combined["n_sources"] = combined["n_sources"].astype("int8")
    print(f"  {len(combined):,} (gene, model_id) pairs ({time.time()-t0:.0f}s)")

    # --- write ---
    out = C.OUT / "bulk_rna_z.parquet"
    combined.to_parquet(out, index=False)

    print(f"\nWrote {len(combined):,} rows -> {out}")
    print(f"  genes:     {combined.gene_id.nunique():,}")
    print(f"  model_ids: {combined.model_id.nunique():,}")
    print(f"  n_sources distribution:\n{combined.n_sources.value_counts().sort_index().to_string()}")
    print(f"  rna_z_t range: [{combined.rna_z_t.min():.2f}, {combined.rna_z_t.max():.2f}]")
    con.close()


if __name__ == "__main__":
    main()
