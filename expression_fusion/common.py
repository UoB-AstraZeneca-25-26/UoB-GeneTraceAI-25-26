"""
expression_fusion/common.py
----------------------------
Loaders and primitives for the multi-source (DepMap + HPA + GEO) RNA
combination work. Read-only against `src/pipeline/outputs/celllineselector.db`
(the Stage-0 warehouse, see docs/DUCKDB_LAYER.md section 6). Writes only into
`expression_fusion/outputs/`.

WHY A NEW MODULE RATHER THAN REUSING THE SHARED NOTEBOOK
----------------------------------------------------------
`02_transcriptonomics_proteins.ipynb` (not part of this repo's pipeline,
supplied separately) assumes a long-format `celllineselector.db` with a
`source`/`value`/`model_id`/`lineage` shape. This repo's actual warehouse
stores `depmap_expr` and `geo_expr` WIDE (one row per profile/sample, one
column per gene) and `hpa_rna` LONG. The row counts match the notebook's
printed output almost exactly (1,495 DepMap profiles, ~591 GEO lines after
collapse, ~1,105 HPA lines) -- which means the notebook's own
`fetch_gene_expression` must pivot wide -> long per gene at query time. This
module does that pivot explicitly against the real schema, batched over many
genes at once (wide-table single-column selects are cheap; per-gene Python
loops over hundreds of genes are not).

The scale-detection and MAD-floor primitives below are deliberately copied
from the shared notebook's cells 6 and 8, not reinvented -- they are the
defensible part of that notebook and reuse beats a second, drifting copy.
"""
from __future__ import annotations

import sys
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy import stats

ROOT = Path(__file__).resolve().parent.parent
HERE = Path(__file__).resolve().parent
OUT = HERE / "outputs"
OUT.mkdir(exist_ok=True)

DB_PATH = ROOT / "src" / "pipeline" / "outputs" / "celllineselector.db"

SEED = 42                 # matches cell_similarity/ and gate_audit/ convention
EXPRESSED_MIN = 1.0        # log2(TPM+1) > 1.0  <=>  TPM > 1  (HPA/Uhlen 2015 convention)
FLOOR_FRAC_LOW = 0.02

SOURCES = ("depmap_expr", "hpa_rna", "geo_expr")


def banner(title: str):
    print("=" * 78)
    print(title)
    print("=" * 78)


def connect(read_only: bool = True) -> duckdb.DuckDBPyConnection:
    if not DB_PATH.exists():
        raise SystemExit(f"{DB_PATH} not found.")
    return duckdb.connect(str(DB_PATH), read_only=read_only)


# =============================================================== scale handling
# Copied verbatim (semantics) from 02_transcriptonomics_proteins.ipynb cell 6.
def detect_scale(v, log_max=25.0, log_min=-10.0) -> str:
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if not len(v):
        return "unknown"
    if v.min() >= log_min and np.nanpercentile(v, 99) < log_max and v.max() <= log_max * 2:
        return "log2"
    return "linear"


def maybe_log(v, name="", verbose=True) -> np.ndarray:
    v = np.asarray(v, float)
    scale = detect_scale(v)
    do = scale == "linear"
    if verbose:
        print(f"  {name:12s} {scale:6s} (min {np.nanmin(v):8.2f}, max {np.nanmax(v):10.2f})"
              f" -> {'log2(x+1) applied' if do else 'left as-is'}")
    return np.log2(np.clip(v, 0, None) + 1) if do else v


# =============================================================== gene universe
def gene_universe(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """Genes present as a queryable column in BOTH wide tables, as a value in
    hpa_rna, AND in the gene dimension. 16,060 of ~19-20k on this warehouse,
    measured 2026-08-12."""
    dm_cols = {c for c in con.execute("DESCRIBE main.depmap_expr").df()["column_name"]
              if c != "profileid"}
    geo_cols = {c for c in con.execute("DESCRIBE main.geo_expr").df()["column_name"]
               if c != "sample"}
    hpa_genes = set(con.execute("SELECT DISTINCT gene FROM main.hpa_rna").df()["gene"])
    gene = con.execute("SELECT gene_id, hugo_symbol FROM main.gene").df()
    common = dm_cols & geo_cols & hpa_genes & set(gene.gene_id)
    return gene[gene.gene_id.isin(common)].reset_index(drop=True)


def sample_genes(con: duckdb.DuckDBPyConnection, n: int, seed: int = SEED) -> pd.DataFrame:
    g = gene_universe(con)
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(g), size=min(n, len(g)), replace=False)
    return g.iloc[sorted(idx)].reset_index(drop=True)


# =============================================================== batched fetch
def fetch_depmap(con, gene_ids: list[str], verbose=True) -> pd.DataFrame:
    """long: gene_id, model_id, value (log2(TPM+1), DepMap's native scale)."""
    cols = ", ".join(f'"{g}"' for g in gene_ids)
    wide = con.execute(f'SELECT profileid, {cols} FROM main.depmap_expr').df()
    prof = con.execute(
        "SELECT profileid, model_id FROM main.depmap_profiles WHERE datatype = 'rna'"
    ).df().drop_duplicates("profileid")
    wide = wide.merge(prof, on="profileid", how="inner").drop(columns=["profileid"])
    long = wide.melt(id_vars="model_id", var_name="gene_id", value_name="value").dropna()
    scale = detect_scale(long["value"])
    if verbose:
        print(f"  depmap_expr  {scale:6s} (min {long.value.min():.2f}, max {long.value.max():.2f})"
              f" -> {'log2(x+1) applied' if scale == 'linear' else 'left as-is'}")
    if scale == "linear":
        long["value"] = np.log2(np.clip(long["value"], 0, None) + 1)
    long = long.groupby(["gene_id", "model_id"], as_index=False)["value"].median()
    long["source"] = "depmap_expr"
    return long


def fetch_geo(con, gene_ids: list[str], verbose=True) -> pd.DataFrame:
    """long: gene_id, model_id, value, n_samples. Collapsed GSM -> model_id by
    median, count as the replicate-precision weight -- the notebook's own
    stated principle (TRANSCRIPTOMICS_DECISIONS.md Q1), applied here too."""
    cols = ", ".join(f'"{g}"' for g in gene_ids)
    wide = con.execute(f'SELECT sample, {cols} FROM main.geo_expr').df()
    gi = con.execute(
        "SELECT geo_accession, model_id FROM main.geo_info WHERE model_id IS NOT NULL"
    ).df().drop_duplicates("geo_accession")
    wide = wide.merge(gi, left_on="sample", right_on="geo_accession", how="inner")
    wide = wide.drop(columns=["sample", "geo_accession"])
    long = wide.melt(id_vars="model_id", var_name="gene_id", value_name="value").dropna()
    scale = detect_scale(long["value"])
    if verbose:
        print(f"  geo_expr     {scale:6s} (min {long.value.min():.2f}, max {long.value.max():.2f})"
              f" -> {'log2(x+1) applied' if scale == 'linear' else 'left as-is'}")
    if scale == "linear":
        long["value"] = np.log2(np.clip(long["value"], 0, None) + 1)
    agg = long.groupby(["gene_id", "model_id"]).agg(
        value=("value", "median"), n_samples=("value", "size")).reset_index()
    agg["source"] = "geo_expr"
    return agg


def fetch_hpa(con, gene_ids: list[str], verbose=True) -> pd.DataFrame:
    """long: gene_id, model_id, value (nTPM -> log2(nTPM+1))."""
    placeholders = ", ".join(f"'{g}'" for g in gene_ids)
    long = con.execute(
        f"SELECT gene AS gene_id, model_id, ntpm AS value FROM main.hpa_rna "
        f"WHERE gene IN ({placeholders}) AND is_ambiguous = FALSE"
    ).df().dropna()
    scale = detect_scale(long["value"])
    if verbose:
        print(f"  hpa_rna      {scale:6s} (min {long.value.min():.2f}, max {long.value.max():.2f})"
              f" -> {'log2(x+1) applied' if scale == 'linear' else 'left as-is'}")
    if scale == "linear":
        long["value"] = np.log2(np.clip(long["value"], 0, None) + 1)
    long = long.groupby(["gene_id", "model_id"], as_index=False)["value"].mean()
    long["source"] = "hpa_rna"
    return long


def fetch_all(con, gene_ids: list[str], verbose=True) -> pd.DataFrame:
    """One long frame: gene_id, model_id, source, value[, n_samples]."""
    if verbose:
        print(f"fetching {len(gene_ids):,} genes x 3 sources")
    dm = fetch_depmap(con, gene_ids, verbose=verbose)
    hp = fetch_hpa(con, gene_ids, verbose=verbose)
    ge = fetch_geo(con, gene_ids, verbose=verbose)
    df = pd.concat([dm, hp, ge], ignore_index=True)
    if verbose:
        n = df.groupby("source").size()
        print(f"  rows: {dict(n)}")
    return df


def fetch_geo_v2(gene_ids: list[str], verbose=True) -> pd.DataFrame:
    """long: gene_id, model_id, value, n_samples -- from the per-series
    scale-corrected rebuild (02_rebuild_geo_scale_corrected.py), not the
    warehouse's pooled `geo_expr` table. Values are already on a single,
    correct log2 scale; no further transform here."""
    path = OUT / "geo_expr_v2.parquet"
    if not path.exists():
        raise SystemExit(f"{path} not found. Run 02_rebuild_geo_scale_corrected.py first.")
    df = pd.read_parquet(path)
    df = df[df.gene_id.isin(gene_ids)].copy()
    df["source"] = "geo_expr_v2"
    if verbose:
        print(f"  geo_expr_v2  (pre-scaled per series)              "
              f"{df.model_id.nunique()} model_ids, {len(df)} rows")
    return df


def load_lineage(con) -> pd.Series:
    d = con.execute(
        "SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL"
    ).df().drop_duplicates("model_id")
    return d.set_index("model_id")["lineage"]


# ============================================================= z-score utilities
# Shared by 05_bulk_rna_scorer.py and 06_bulk_protein_scorer.py.
MAD_FLOOR = 0.384   # calibrated floor matching the shared transcriptomics notebook
MIN_PEERS = 5       # minimum within-(source, lineage) peers to produce a z-score


def robust_z_matrix(X: np.ndarray,
                    min_peers: int = MIN_PEERS,
                    mad_floor: float = MAD_FLOOR) -> np.ndarray:
    """
    Vectorised robust z-score across rows (cell lines) for each column (gene).
    X shape: (n_lines, n_genes). Returns same shape; gene columns with < min_peers
    non-NaN rows are set all-NaN.
    """
    n_valid = np.sum(np.isfinite(X), axis=0)
    with np.errstate(all="ignore"):   # suppress all-NaN slice warnings for empty lineages
        med     = np.nanmedian(X, axis=0)
        abs_dev = np.abs(X - med)
        mad     = np.nanmedian(abs_dev, axis=0) * 1.4826   # scale='normal'
    mad = np.maximum(mad, mad_floor)
    mad = np.where(np.isfinite(mad), mad, mad_floor)
    Z   = (X - med) / mad
    Z[:, n_valid < min_peers] = np.nan
    return Z


def score_source_lineage(wide: pd.DataFrame,
                         gene_ids: list,
                         lineage_map: pd.Series,
                         source_name: str) -> pd.DataFrame:
    """
    Apply robust_z_matrix within each lineage stratum.
    wide   : model_id (index) × gene_id (columns), values already on log scale.
    Returns long DataFrame: gene_id, model_id, source, z
    """
    wide = wide.loc[wide.index.isin(lineage_map.index), gene_ids].copy()
    wide["lineage"] = lineage_map.reindex(wide.index)
    wide = wide.dropna(subset=["lineage"])

    rows = []
    for _, grp in wide.groupby("lineage"):
        X         = grp[gene_ids].values.astype(float)
        Z         = robust_z_matrix(X)
        model_ids = grp.index.tolist()
        z_wide    = pd.DataFrame(Z, index=model_ids, columns=gene_ids)
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
    """
    Vectorised Stouffer combination per (gene_id, model_id).
    Weight  = sqrt(n_lines this source covers for this gene).
    Shrinkage: divide z_raw by sqrt(N_SOURCES_TOTAL / n_available) when fewer
    than n_sources_total sources contributed — mirrors the RNA scorer's
    single-source penalty.
    Returns: gene_id, model_id, z_t (float32), n_sources (int8).
    """
    valid = z_long.dropna(subset=["z"]).copy()
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

    denom  = np.sqrt(agg["w2_sum"].values)
    safe   = denom > 0
    z_raw  = np.where(safe, agg["wz_sum"].values / np.where(safe, denom, 1.0), np.nan)
    shrink = np.where(agg["n_sources"] < n_sources_total,
                      np.sqrt(n_sources_total / agg["n_sources"]), 1.0)
    agg["z_t"]      = (z_raw / shrink).astype("float32")
    agg["n_sources"] = agg["n_sources"].astype("int8")
    return agg[["gene_id", "model_id", "z_t", "n_sources"]]
