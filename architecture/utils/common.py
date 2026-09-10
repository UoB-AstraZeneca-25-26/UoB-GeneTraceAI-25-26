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
    """
    Open a connection to the Stage-0 DuckDB warehouse.

    Parameters
    ----------
    read_only : bool, default True
        Passed straight to ``duckdb.connect``. Every pipeline stage after
        harmonisation should leave this at the default — only
        ``00_harmonisation`` writes to the warehouse.

    Returns
    -------
    duckdb.DuckDBPyConnection
        Open connection to ``config.DB``. Caller is responsible for
        closing it.

    Notes
    -----
    Raises ``SystemExit`` if the warehouse file does not exist yet,
    instructing the caller to run ``00_harmonisation`` first.
    """
    if not DB.exists():
        raise SystemExit(f"Warehouse not found: {DB}\nRun 00_harmonisation first.")
    return duckdb.connect(str(DB), read_only=read_only)


def gene_universe(con: duckdb.DuckDBPyConnection) -> pd.DataFrame:
    """
    Return the genes present in every one of the three RNA sources.

    A gene only enters the scored universe if it has a column in
    ``depmap_expr``, a column in ``geo_expr``, AND at least one row in
    ``hpa_rna`` — the intersection, not the union, of the three sources'
    gene sets.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the Stage-0 warehouse.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_id``, ``hugo_symbol`` — one row per gene in the
        three-way intersection, in ``main.gene``'s row order (index reset).
    """
    dm_cols  = {c for c in con.execute("DESCRIBE main.depmap_expr").df()["column_name"] if c != "profile_id"}
    geo_cols = {c for c in con.execute("DESCRIBE main.geo_expr").df()["column_name"] if c != "gsm_id"}
    hpa_genes = set(con.execute("SELECT DISTINCT gene_id FROM main.hpa_rna").df()["gene_id"])
    gene = con.execute("SELECT gene_id, hugo_symbol FROM main.gene").df()
    common = dm_cols & geo_cols & hpa_genes & set(gene.gene_id)
    return gene[gene.gene_id.isin(common)].reset_index(drop=True)


def load_lineage(con) -> pd.Series:
    """
    Load the model_id -> lineage mapping used to stratify scores.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the Stage-0 warehouse.

    Returns
    -------
    pandas.Series
        Indexed by ``model_id``, values are the lineage string. Rows
        where ``sample_info.lineage`` is NULL are excluded; if a
        ``model_id`` has more than one row in ``sample_info`` only the
        first is kept (``drop_duplicates``).
    """
    d = con.execute(
        "SELECT model_id, lineage FROM main.sample_info WHERE lineage IS NOT NULL"
    ).df().drop_duplicates("model_id")
    return d.set_index("model_id")["lineage"]


def detect_scale(v, log_max=25.0, log_min=-10.0) -> str:
    """
    Guess whether a vector of expression values is already log2-scaled.

    Parameters
    ----------
    v : array-like
        Expression values to inspect. Non-finite entries are dropped
        before the check.
    log_max : float, default 25.0
        A vector is classified "log2" only if its max does not exceed
        ``log_max * 2`` and its 99th percentile is below ``log_max``.
    log_min : float, default -10.0
        A vector is classified "log2" only if its min is at or above
        this floor (log2 values can be legitimately negative for
        near-zero linear expression; linear (TPM/counts) values cannot).

    Returns
    -------
    str
        ``"log2"``, ``"linear"``, or ``"unknown"`` (empty/all-non-finite
        input).

    Notes
    -----
    Heuristic, not a guarantee — calibrated against this project's three
    RNA sources, not a general-purpose scale detector.
    """
    v = np.asarray(v, float)
    v = v[np.isfinite(v)]
    if not len(v):
        return "unknown"
    if v.min() >= log_min and np.nanpercentile(v, 99) < log_max and v.max() <= log_max * 2:
        return "log2"
    return "linear"


def _to_log2(v: np.ndarray) -> np.ndarray:
    """
    Convert linear-scale expression values to log2(x + 1).

    Parameters
    ----------
    v : numpy.ndarray
        Linear-scale expression values (e.g. TPM). Negative inputs are
        clipped to 0 before the transform, since negative expression is
        not meaningful.

    Returns
    -------
    numpy.ndarray
        ``log2(clip(v, 0, None) + 1)``, same shape as ``v``.
    """
    return np.log2(np.clip(v, 0, None) + 1)


def fetch_depmap(con, gene_ids: list[str]) -> pd.DataFrame:
    """
    Load DepMap RNA expression for a set of genes, reshaped to long form.

    Joins ``depmap_expr`` (one column per gene, one row per profile_id)
    against ``depmap_profiles`` restricted to ``datatype = 'rna'`` to
    resolve each profile to a ``model_id``, converts to log2 if the raw
    values are detected as linear-scale (see :func:`detect_scale`), and
    collapses duplicate (gene, model) measurements with the median.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the Stage-0 warehouse.
    gene_ids : list of str
        Gene column names to pull from ``main.depmap_expr``.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_id``, ``model_id``, ``value`` (log2 scale),
        ``source`` (always ``"depmap_expr"``). One row per (gene, model)
        pair with at least one non-null measurement.
    """
    cols = ", ".join(f'"{g}"' for g in gene_ids)
    wide = con.execute(f'SELECT profile_id, {cols} FROM main.depmap_expr').df()
    prof = con.execute(
        "SELECT profile_id, model_id FROM main.depmap_profiles WHERE datatype = 'rna'"
    ).df().drop_duplicates("profile_id")
    wide = wide.merge(prof, on="profile_id", how="inner").drop(columns=["profile_id"])
    long = wide.melt(id_vars="model_id", var_name="gene_id", value_name="value").dropna()
    if detect_scale(long["value"]) == "linear":
        long["value"] = _to_log2(long["value"].values)
    long = long.groupby(["gene_id", "model_id"], as_index=False)["value"].median()
    long["source"] = "depmap_expr"
    return long


def fetch_hpa(con, gene_ids: list[str]) -> pd.DataFrame:
    """
    Load Human Protein Atlas RNA expression for a set of genes.

    Reads ``ntpm`` directly from ``main.hpa_rna`` (already one row per
    (gene, model)), excludes rows flagged ``is_ambiguous``, converts to
    log2 if detected as linear-scale, and collapses duplicates with the
    mean.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the Stage-0 warehouse.
    gene_ids : list of str
        Gene IDs to filter ``main.hpa_rna`` on.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_id``, ``model_id``, ``value`` (log2 scale),
        ``source`` (always ``"hpa_rna"``).
    """
    ph = ", ".join(f"'{g}'" for g in gene_ids)
    long = con.execute(
        f"SELECT gene_id, model_id, ntpm AS value FROM main.hpa_rna "
        f"WHERE gene_id IN ({ph}) AND is_ambiguous = FALSE"
    ).df().dropna()
    if detect_scale(long["value"]) == "linear":
        long["value"] = _to_log2(long["value"].values)
    long = long.groupby(["gene_id", "model_id"], as_index=False)["value"].mean()
    long["source"] = "hpa_rna"
    return long


def fetch_geo(con, gene_ids: list[str]) -> pd.DataFrame:
    """
    Load GEO RNA expression for a set of genes, reshaped to long form.

    Joins ``geo_expr`` (one column per gene, one row per GSM) against
    ``geo_info`` (restricted to rows with a resolved ``model_id``) to map
    each sample to a cell line, converts to log2 if detected as
    linear-scale, and collapses duplicate (gene, model) measurements
    with the median — tracking how many raw samples contributed via
    ``n_samples``.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the Stage-0 warehouse.
    gene_ids : list of str
        Gene column names to pull from ``main.geo_expr``.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_id``, ``model_id``, ``value`` (log2 scale, median
        of contributing samples), ``n_samples`` (count of raw GSM rows
        collapsed into this value), ``source`` (always ``"geo_expr"``).
    """
    cols = ", ".join(f'"{g}"' for g in gene_ids)
    wide = con.execute(f'SELECT gsm_id, {cols} FROM main.geo_expr').df()
    gi   = con.execute(
        "SELECT geo_accession, model_id FROM main.geo_info WHERE model_id IS NOT NULL"
    ).df().drop_duplicates("geo_accession")
    wide = wide.merge(gi, left_on="gsm_id", right_on="geo_accession", how="inner").drop(
        columns=["gsm_id", "geo_accession"])
    long = wide.melt(id_vars="model_id", var_name="gene_id", value_name="value").dropna()
    if detect_scale(long["value"]) == "linear":
        long["value"] = _to_log2(long["value"].values)
    agg = long.groupby(["gene_id", "model_id"]).agg(
        value=("value", "median"), n_samples=("value", "size")).reset_index()
    agg["source"] = "geo_expr"
    return agg


def robust_z_matrix(X: np.ndarray, min_peers: int = MIN_PEERS,
                    mad_floor: float = MAD_FLOOR) -> np.ndarray:
    """
    Column-wise robust z-score a matrix using median and MAD.

    Each column of ``X`` (typically one gene, rows = cell lines within a
    single lineage) is standardized as ``(x - median) / (MAD * 1.4826)``,
    the constant that makes MAD a consistent estimator of the standard
    deviation under normality. Columns with fewer than ``min_peers``
    finite values are entirely blanked to NaN — a z-score computed from
    too few peers is not treated as meaningful.

    Parameters
    ----------
    X : numpy.ndarray
        2-D array, rows = samples, columns = genes. Non-finite entries
        are ignored in the median/MAD computation via ``nanmedian``.
    min_peers : int, default MIN_PEERS (5)
        Minimum number of finite values a column must have; columns
        below this are set to all-NaN in the output.
    mad_floor : float, default MAD_FLOOR (0.384)
        Lower bound on the scale estimate, applied after the 1.4826
        rescale. Prevents division by a near-zero or non-finite MAD
        (e.g. when most values in a column are identical) from
        producing exploding or undefined z-scores.

    Returns
    -------
    numpy.ndarray
        Same shape as ``X``. ``NaN`` wherever the input was NaN, or
        wherever the column had fewer than ``min_peers`` finite values.
    """
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
    """
    Robust-z score one RNA source, independently within each lineage.

    For each lineage, builds a (lines x genes) matrix restricted to that
    lineage's lines, robust-z scores it column-wise (see
    :func:`robust_z_matrix`), then applies a silence guard: any gene
    expressed (``value > EXPRESSED_MIN``) in fewer than ``SILENT_FRAC``
    of that lineage's lines has its z-scores blanked, since a z-score
    computed mostly from "not really expressed" values is not a
    meaningful signal. Results across all lineages are stacked to long
    form.

    Parameters
    ----------
    wide : pandas.DataFrame
        Index = model_id, columns include every id in ``gene_ids`` plus
        whatever else the caller passed through; only the gene columns
        and index are used here.
    gene_ids : list
        Gene columns to score.
    lineage_map : pandas.Series
        Indexed by model_id, values are lineage strings (see
        :func:`load_lineage`). Rows of ``wide`` whose model_id is not in
        this index, or whose lineage is NaN, are dropped before scoring.
    source_name : str
        Written verbatim into the output's ``source`` column (e.g.
        ``"depmap_expr"``).

    Returns
    -------
    pandas.DataFrame
        Columns ``model_id``, ``gene_id``, ``z``, ``source``. One row
        per (gene, model) pair that survived both the peer-count floor
        in :func:`robust_z_matrix` and the silence guard here. Empty
        DataFrame with the same columns if no lineage yielded any rows.
    """
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
    """
    Combine per-source z-scores into one Stouffer-weighted z per (gene, model).

    Each source's contribution is weighted by ``sqrt(src_n)``, where
    ``src_n`` is the number of distinct lines that source scored for
    that gene (more corroborating lines -> more weight for that
    source's z on that gene). A pair scored by fewer than
    ``n_sources_total`` sources is additionally shrunk toward zero by
    ``sqrt(n_observed / n_sources_total)`` — a deliberate, not yet
    theoretically derived, discount for low-corroboration measurements
    (e.g. a pair seen by 1 of 3 sources is discounted to about 0.577x
    the Stouffer-combined value it would otherwise get).

    Parameters
    ----------
    z_long : pandas.DataFrame
        Long-form per-source z-scores with columns ``gene_id``,
        ``model_id``, ``source``, ``z`` — typically the concatenation of
        several :func:`score_source_lineage` outputs. Rows with NaN
        ``z`` are dropped before combining.
    n_sources_total : int
        The total number of sources in play for this run (e.g. 3 for
        RNA's depmap_expr/hpa_rna/geo_expr). Used only to decide how
        much to shrink pairs seen by fewer than this many sources — it
        does not filter which sources are combined.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_id``, ``model_id``, ``z_t`` (float32, the final
        shrunk Stouffer z), ``n_sources`` (int8, how many distinct
        sources actually contributed to that pair). One row per
        (gene, model) pair with at least one valid source.
    """
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
