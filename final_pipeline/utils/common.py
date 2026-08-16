"""
utils/common.py
---------------
Shared loaders and z-score primitives for the GeneTraceAI pipeline.
Reads from the Stage-0 DuckDB warehouse. Writes nothing — callers own outputs.
"""
from __future__ import annotations
import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from config import DB

SEED          = 42
EXPRESSED_MIN = 1.0    # log2(TPM+1) > 1.0 <=> TPM > 1 (HPA/Uhlen 2015)
SILENT_FRAC   = 0.20   # if < 20% of lineage lines express a gene, z-scores are uninformative
MAD_FLOOR     = 0.384  # calibrated against the shared transcriptomics notebook
MIN_PEERS     = 5      # minimum within-(source, lineage) peers for a z-score


def connect(read_only: bool = True) -> duckdb.DuckDBPyConnection:
    if not DB.exists():
        raise SystemExit(f"Warehouse not found: {DB}\nRun 00_harmonisation first.")
    return duckdb.connect(str(DB), read_only=read_only)


def gene_universe(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    dm_cols  = {c for c in con.execute("DESCRIBE main.depmap_expr").df()["column_name"] if c != "profileid"}
    geo_cols = {c for c in con.execute("DESCRIBE main.geo_expr").df()["column_name"] if c != "sample"}
    hpa_genes = set(con.execute("SELECT DISTINCT gene FROM main.hpa_rna").df()["gene"])
    gene = con.execute("SELECT gene_id, hugo_symbol FROM main.gene").df()
    common = dm_cols & geo_cols & hpa_genes & set(gene.gene_id)
    return gene[gene.gene_id.isin(common)].reset_index(drop=True)


def load_lineage(con) -> pd.Series:
    d = con.execute(
        "SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL"
    ).df().drop_duplicates("model_id")
    return d.set_index("model_id")["lineage"]


def detect_scale(v, log_max=25.0, log_min=-10.0) -> str:
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if not len(v):
        return "unknown"
    if v.min() >= log_min and np.nanpercentile(v, 99) < log_max and v.max() <= log_max * 2:
        return "log2"
    return "linear"


def _to_log2(v: np.ndarray) -> np.ndarray:
    return np.log2(np.clip(v, 0, None) + 1)


def fetch_depmap(con, gene_ids: list[str]) -> pd.DataFrame:
    cols = ", ".join(f'"{g}"' for g in gene_ids)
    wide = con.execute(f'SELECT profileid, {cols} FROM main.depmap_expr').df()
    prof = con.execute(
        "SELECT profileid, model_id FROM main.depmap_profiles WHERE datatype = 'rna'"
    ).df().drop_duplicates("profileid")
    wide = wide.merge(prof, on="profileid", how="inner").drop(columns=["profileid"])
    long = wide.melt(id_vars="model_id", var_name="gene_id", value_name="value").dropna()
    if detect_scale(long["value"]) == "linear":
        long["value"] = _to_log2(long["value"].values)
    long = long.groupby(["gene_id", "model_id"], as_index=False)["value"].median()
    long["source"] = "depmap_expr"
    return long


def fetch_hpa(con, gene_ids: list[str]) -> pd.DataFrame:
    ph = ", ".join(f"'{g}'" for g in gene_ids)
    long = con.execute(
        f"SELECT gene AS gene_id, model_id, ntpm AS value FROM main.hpa_rna "
        f"WHERE gene IN ({ph}) AND is_ambiguous = FALSE"
    ).df().dropna()
    if detect_scale(long["value"]) == "linear":
        long["value"] = _to_log2(long["value"].values)
    long = long.groupby(["gene_id", "model_id"], as_index=False)["value"].mean()
    long["source"] = "hpa_rna"
    return long


def fetch_geo(con, gene_ids: list[str]) -> pd.DataFrame:
    cols = ", ".join(f'"{g}"' for g in gene_ids)
    wide = con.execute(f'SELECT sample, {cols} FROM main.geo_expr').df()
    gi   = con.execute(
        "SELECT geo_accession, model_id FROM main.geo_info WHERE model_id IS NOT NULL"
    ).df().drop_duplicates("geo_accession")
    wide = wide.merge(gi, left_on="sample", right_on="geo_accession", how="inner").drop(
        columns=["sample", "geo_accession"])
    long = wide.melt(id_vars="model_id", var_name="gene_id", value_name="value").dropna()
    if detect_scale(long["value"]) == "linear":
        long["value"] = _to_log2(long["value"].values)
    agg = long.groupby(["gene_id", "model_id"]).agg(
        value=("value", "median"), n_samples=("value", "size")).reset_index()
    agg["source"] = "geo_expr"
    return agg


def robust_z_matrix(X: np.ndarray, min_peers: int = MIN_PEERS,
                    mad_floor: float = MAD_FLOOR) -> np.ndarray:
    n_valid = np.sum(np.isfinite(X), axis=0)
    with np.errstate(all="ignore"):
        med = np.nanmedian(X, axis=0)
        mad = np.nanmedian(np.abs(X - med), axis=0) * 1.4826
    mad = np.maximum(np.where(np.isfinite(mad), mad, mad_floor), mad_floor)
    Z = (X - med) / mad
    Z[:, n_valid < min_peers] = np.nan
    return Z


def score_source_lineage(wide: pd.DataFrame, gene_ids: list,
                         lineage_map: pd.Series, source_name: str) -> pd.DataFrame:
    wide = wide.loc[wide.index.isin(lineage_map.index), gene_ids].copy()
    wide["lineage"] = lineage_map.reindex(wide.index)
    wide = wide.dropna(subset=["lineage"])
    rows = []
    for _, grp in wide.groupby("lineage"):
        X = grp[gene_ids].values.astype(float)
        Z = robust_z_matrix(X)
        # SILENCE GUARD — added 2026-08-15.
        # This restores the third guard from the source notebook
        # (02_transcriptomics.ipynb, robust_z), which was dropped during the
        # adaptation into robust_z_matrix. SILENT_FRAC = 0.20 and
        # EXPRESSED_MIN = 1.0 are the notebook's own values.
        # Placed here because raw log2 expression is not reachable from
        # core_score.py without a full warehouse batch-fetch loop (the same
        # loop rna_scorer.py already runs). This is a change inside the
        # adapted RNA utility layer; it does not modify 02_transcriptomics.ipynb.
        frac_expressed = np.mean(X > EXPRESSED_MIN, axis=0)
        Z[:, frac_expressed < SILENT_FRAC] = np.nan
        z_wide = pd.DataFrame(Z, index=grp.index.tolist(), columns=gene_ids)
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


def stouffer_combine(z_long: pd.DataFrame, n_sources_total: int) -> pd.DataFrame:
    valid = z_long.dropna(subset=["z"]).copy()
    src_n = valid.groupby(["gene_id", "source"])["model_id"].nunique().reset_index(name="src_n")
    valid = valid.merge(src_n, on=["gene_id", "source"], how="left")
    valid["w"]  = np.sqrt(valid["src_n"])
    valid["wz"] = valid["w"] * valid["z"]
    valid["w2"] = valid["w"] ** 2
    agg = valid.groupby(["gene_id", "model_id"]).agg(
        wz_sum=("wz", "sum"), w2_sum=("w2", "sum"), n_sources=("source", "nunique")
    ).reset_index()
    denom  = np.sqrt(agg["w2_sum"].values)
    safe   = denom > 0
    z_raw  = np.where(safe, agg["wz_sum"].values / np.where(safe, denom, 1.0), np.nan)
    # Conservative penalty for sparse-source pairs: divide by √(n_total/n_observed),
    # equivalent to multiplying by √(n_observed/n_total) < 1. This SHRINKS z toward
    # zero beyond what the Stouffer weighting already captures. At n_observed=1 of
    # n_total=3, z_t = z_raw × √(1/3) ≈ 0.577 × z_raw. Documented as a deliberate
    # epistemic discount for low-corroboration measurements; not yet theoretically derived.
    shrink = np.where(agg["n_sources"] < n_sources_total,
                      np.sqrt(n_sources_total / agg["n_sources"]), 1.0)
    agg["z_t"]       = (z_raw / shrink).astype("float32")
    agg["n_sources"] = agg["n_sources"].astype("int8")
    return agg[["gene_id", "model_id", "z_t", "n_sources"]]
