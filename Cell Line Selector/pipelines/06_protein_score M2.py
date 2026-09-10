"""
Cell-line z-scores from two proteomics sources (ProCAN + CCLE proteomics).

Three changes over the base version:
  1. Isoforms are handled explicitly. When several protein IDs map to one gene,
     we first check whether they agree across cell lines (Spearman). If they do,
     the gene value for a line is the MAX abundance across isoforms; if they do
     not, we fall back to the mean so an unrelated isoform can't dominate.
  2. ProCAN is NOT preferred over CCLE. Both platforms are combined with
     coverage-only (symmetric) weights, and platform disagreement lowers our
     confidence in the *merged gene* instead of deleting one platform's data.
  3. No domain threshold is hardcoded. Isoform concordance, platform-agreement
     bounds and the minimum shared-line count are all derived from the data
     (null / matched distributions + statistical power). The only remaining
     numbers are standard statistical conventions (significance level, the
     percentiles that define "null" and "clear agreement"), surfaced as
     documented parameters rather than buried constants.

Output table
------------
    gene_id | model_id | lineage | z_score | z_raw | platform_confidence
            | n_sources | n_isoforms

``z_raw`` is the coverage-weighted combination; ``z_score`` is that value
scaled by ``platform_confidence``, so a gene the two platforms disagree
about is shrunk toward zero rather than having one platform's data
discarded. A QC table of isoform rollup decisions is written alongside.

Usage
-----
    python <this_file>.py
    python <this_file>.py --save-params results/protein/thresholds.json
    python <this_file>.py --params results/protein/thresholds.json
"""

from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import spearmanr, t as t_dist


ROOT = Path(__file__).resolve().parent.parent
DEFAULT_DB = ROOT / "db" / "celllineselector.duckdb"

# Method 1 source names.
PROCAN_SOURCE = "protein_matrix_averaged_20250211"
CCLE_SOURCE = "proteomics"

# Statistical CONVENTIONS (not domain-tuned values). Exposed via CLI so they are
# visible and overridable rather than hidden. 1.4826 elsewhere is the fixed
# MAD -> sigma constant for a normal, a mathematical constant, not a tunable.
DEFAULT_ALPHA = 0.05        # significance level for the power-based min overlap
DEFAULT_NULL_PCTILE = 95.0  # "more correlated than X% of unrelated pairs" = agree
DEFAULT_MATCHED_PCTILE = 60.0  # where matched agreement is treated as full trust
DEFAULT_MIN_GAP = 0.10      # minimum width of the confidence ramp (low -> high)


# ---------------------------------------------------------------------
# Identifier helpers (unchanged from base)
# ---------------------------------------------------------------------

def canonical_ensg(value) -> str | None:
    """
    Normalise a value to a bare, version-stripped Ensembl gene ID.

    Accepts ``ENSG00000141510`` or ``ENSG00000141510.17`` in any case.
    Anything not matching the ENSG pattern returns None, so callers can
    drop unmapped rows with ``dropna``.

    Parameters
    ----------
    value : any
        Candidate identifier; NaN, None and empty strings are tolerated.

    Returns
    -------
    str or None
        Uppercase ENSG ID without its version suffix, or None.
    """
    if value is None or pd.isna(value):
        return None
    text = str(value).strip().upper()
    if not text:
        return None
    match = re.fullmatch(r"ENSG\d+(?:\.\d+)?", text)
    return text.split(".", 1)[0] if match else None


def split_id_list(value) -> list[str]:
    """
    Split a delimited identifier string into clean uppercase tokens.

    Roster columns pack several identifiers into one cell with no single
    consistent delimiter, so semicolons, commas, pipes and both slashes
    are all treated as separators. String forms of missing values
    (``nan``, ``none``, ``null``) are discarded.

    Parameters
    ----------
    value : any
        Delimited identifier string; NaN and None are tolerated.

    Returns
    -------
    list of str
        Uppercased tokens, empty if nothing usable was found.
    """
    if value is None or pd.isna(value):
        return []
    text = str(value).strip()
    if not text:
        return []
    out = []
    for p in re.split(r"[;,|/\\]", text):
        token = str(p).strip().upper()
        if token and token.lower() not in {"nan", "none", "null"}:
            out.append(token)
    return out


def build_uniprot_to_ensg(con) -> dict[str, str]:
    """
    Build a lowercase UniProt accession -> ENSG lookup from ``gene_roster``.

    Reads both accession columns, ``uniprot_id(s)`` and
    ``procan_uniprot_id``, splitting each multi-value cell.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the project database.

    Returns
    -------
    dict
        Lowercased accession to canonical ENSG ID. Where an accession
        appears in both columns, ``procan_uniprot_id`` wins, since it is
        processed second.

    Notes
    -----
    A column absent from the roster is skipped rather than raising, so
    the map reflects whatever columns exist.
    """
    mapping: dict[str, str] = {}
    roster = con.execute(
        'SELECT gene_id, "uniprot_id(s)", procan_uniprot_id FROM gene_roster'
    ).fetchdf()
    for col in ["uniprot_id(s)", "procan_uniprot_id"]:
        if col not in roster.columns:
            continue
        for _, row in roster[["gene_id", col]].dropna().iterrows():
            gene = canonical_ensg(row["gene_id"])
            if gene is None:
                continue
            for token in split_id_list(row[col]):
                mapping[token.lower()] = gene
    return mapping


def build_symbol_to_ensg(con) -> dict[str, str]:
    """
    Build a lowercase gene-symbol -> ENSG lookup from ``gene_roster``.

    Reads ``gene_name(s)``, splitting each multi-value cell, so aliases
    map to the same gene as the primary symbol. This is the fallback
    identifier space for platforms whose columns are symbols rather than
    UniProt accessions.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the project database.

    Returns
    -------
    dict
        Lowercased symbol to canonical ENSG ID. Empty if the roster has
        no ``gene_name(s)`` column.
    """
    mapping: dict[str, str] = {}
    roster = con.execute('SELECT gene_id, "gene_name(s)" FROM gene_roster').fetchdf()
    if "gene_name(s)" not in roster.columns:
        return mapping
    for _, row in roster[["gene_id", "gene_name(s)"]].dropna().iterrows():
        gene = canonical_ensg(row["gene_id"])
        if gene is None:
            continue
        for token in split_id_list(row["gene_name(s)"]):
            mapping[token.lower()] = gene
    return mapping


def normalize_name(value) -> str:
    """
    Reduce a cell-line name to a punctuation-free comparison key.

    ``MCF-7``, ``MCF 7`` and ``mcf7`` all collapse to ``mcf7``, which is
    what lets a protein matrix's ``ccle_name`` match ``sample_info``
    entries that punctuate the same line differently.

    Parameters
    ----------
    value : any
        Cell-line name; NaN and None are tolerated.

    Returns
    -------
    str
        Alphanumeric-only lowercase key, empty if the value is missing.
    """
    if value is None or pd.isna(value):
        return ""
    return re.sub(r"[^a-z0-9]+", "", str(value).strip().lower())


# ---------------------------------------------------------------------
# Loading / reshaping (unchanged from base)
# ---------------------------------------------------------------------

def attach_model_ids(con, df: pd.DataFrame) -> pd.DataFrame:
    """
    Attach ``model_id`` to a protein table keyed only by ``ccle_name``.

    Builds a name lookup from ``sample_info`` across ``ccle_id``,
    ``cell_line_name`` and ``stripped_cell_line_name``, normalising each
    through :func:`normalize_name`. Columns are tried in that order and
    ``setdefault`` keeps the first hit, so ``ccle_id`` takes precedence
    on collision.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the project database.
    df : pandas.DataFrame
        Protein table. Returned unchanged if it already has ``model_id``.

    Returns
    -------
    pandas.DataFrame
        Copy of ``df`` with ``model_id`` attached; rows whose name
        matched nothing are dropped.

    Raises
    ------
    ValueError
        If ``df`` has neither ``model_id`` nor ``ccle_name``, or if
        ``sample_info`` yields no usable mappings — in both cases nothing
        downstream could be keyed, so failing loudly beats returning an
        empty frame.
    """
    if "model_id" in df.columns:
        return df
    if "ccle_name" not in df.columns:
        raise ValueError("Protein matrix missing both model_id and ccle_name.")

    meta = con.execute(
        """
        SELECT model_id, cell_line_name, stripped_cell_line_name, ccle_id
        FROM sample_info WHERE model_id IS NOT NULL
        """
    ).fetchdf()
    if meta.empty:
        raise ValueError("sample_info has no usable cell-line mappings.")

    lookup: dict[str, str] = {}
    for col in ["ccle_id", "cell_line_name", "stripped_cell_line_name"]:
        if col not in meta.columns:
            continue
        for _, row in meta[["model_id", col]].dropna().iterrows():
            key = normalize_name(row[col])
            if key:
                lookup.setdefault(key, str(row["model_id"]))

    out = df.copy()
    out["model_id"] = out["ccle_name"].map(lambda v: lookup.get(normalize_name(v), None))
    return out.dropna(subset=["model_id"]).copy()


def wide_to_long(df: pd.DataFrame, id_col: str) -> pd.DataFrame:
    """
    Melt a wide protein matrix into long form.

    Every column other than ``id_col`` is treated as a protein
    measurement. Missing measurements are dropped rather than carried as
    NaN, so the result holds observations only.

    Parameters
    ----------
    df : pandas.DataFrame
        Wide table: one row per cell line, one column per protein.
    id_col : str
        Column identifying the cell line, typically ``model_id``.

    Returns
    -------
    pandas.DataFrame
        Columns ``[id_col, "protein", "value"]`` with a reset index.

    Notes
    -----
    ``stack`` is called with ``future_stack=True`` where supported,
    falling back on ``TypeError``, which keeps this working across the
    pandas 2.x stack deprecation.
    """
    protein_cols = [c for c in df.columns if c != id_col]
    sub = df[[id_col] + protein_cols].set_index(id_col)
    try:
        long = sub.stack(future_stack=True).reset_index()
    except TypeError:
        long = sub.stack().reset_index()
    long.columns = [id_col, "protein", "value"]
    return long.dropna(subset=["value"]).reset_index(drop=True)


def table_to_long(con, df, source_name, protein_map) -> pd.DataFrame:
    """
    Convert one wide protein source table into a mapped long frame.

    Ensures the table is keyed by ``model_id`` (attaching it from
    ``ccle_name`` if needed), melts it, lowercases the protein
    identifiers to match the lookup convention, and maps each to a
    canonical ENSG ID. Both the protein column and the mapped gene ID
    are retained — the protein label is what the isoform collapse later
    groups on.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection, used only if ``model_id`` must be attached.
    df : pandas.DataFrame
        Wide protein table as read from DuckDB.
    source_name : str
        Platform name, recorded in ``source`` and used in error messages.
    protein_map : dict
        Lowercased UniProt accession or gene symbol to ENSG ID.

    Returns
    -------
    pandas.DataFrame
        Columns ``model_id``, ``protein``, ``value``, ``gene_id``,
        ``source``.

    Raises
    ------
    ValueError
        If the table has neither ``model_id`` nor ``ccle_name``.

    Notes
    -----
    Unmapped proteins are dropped silently, so a platform whose column
    identifiers are absent from ``protein_map`` yields an empty frame
    rather than an error. Worth checking per-source row counts when a
    platform contributes nothing.
    """
    work = df.copy()
    if "model_id" not in work.columns and "ccle_name" in work.columns:
        work = attach_model_ids(con, work)
    if "model_id" not in work.columns:
        raise ValueError(f"Table '{source_name}' has neither model_id nor ccle_name.")

    long = wide_to_long(work, "model_id")
    long["protein"] = long["protein"].astype(str).str.strip().str.lower()
    long["gene_id"] = long["protein"].map(protein_map)
    long = long.dropna(subset=["gene_id"]).copy()
    long["gene_id"] = long["gene_id"].map(canonical_ensg)
    long = long.dropna(subset=["gene_id"]).copy()
    long["source"] = source_name
    return long


def load_lineage_map(con) -> pd.Series:
    """
    Load the ``model_id`` -> lineage mapping from ``sample_info``.

    Lineage is what the z-scores are conditioned on, so a protein is
    scored against cell lines of the same tissue of origin rather than
    against the whole panel.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the project database.

    Returns
    -------
    pandas.Series
        Lineage indexed by ``model_id``, missing values filled as
        ``"unknown"``. Empty if ``sample_info`` yields no rows.

    Notes
    -----
    Duplicates are collapsed by keeping the first row per ``model_id``,
    so conflicting lineage labels resolve arbitrarily but
    deterministically.
    """
    df = con.execute(
        "SELECT DISTINCT model_id, lineage FROM sample_info WHERE model_id IS NOT NULL"
    ).fetchdf()
    if df.empty:
        return pd.Series(dtype=object)
    return (
        df.drop_duplicates(subset=["model_id"]).set_index("model_id")["lineage"]
        .fillna("unknown")
    )


def load_sources(con, table_names, protein_map) -> dict[str, pd.DataFrame]:
    """
    Load and convert every available protein source table.

    Each named table is read, converted to a mapped long frame, and kept
    only if it yielded rows. Tables that are missing, unreadable, empty
    or wholly unmapped are skipped, so the pipeline runs on whichever
    platforms are actually present.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the project database.
    table_names : iterable of str
        Source table names to attempt.
    protein_map : dict
        Lowercased identifier to ENSG ID.

    Returns
    -------
    dict
        Source name to its long frame. Empty if nothing loaded.

    Notes
    -----
    The ``except Exception: continue`` swallows every read error,
    including genuine SQL or permission problems, not just a missing
    table. A platform silently absent from the result will show up as
    ``n_sources = 1`` downstream rather than as a failure.
    """
    raw: dict[str, pd.DataFrame] = {}
    for name in table_names:
        try:
            table = con.execute(f'SELECT * FROM "{name}"').fetchdf()
        except Exception:
            continue
        if table.empty:
            continue
        long = table_to_long(con, table, name, protein_map)
        if not long.empty:
            raw[name] = long
    return raw


# ---------------------------------------------------------------------
# CHANGE 1: isoform check -> max if concordant, else mean
# ---------------------------------------------------------------------

def _pairwise_rho(a: pd.Series, b: pd.Series, min_overlap: int) -> float:
    """
    Spearman correlation between two series on their shared observations.

    Returns NaN rather than a number in the three cases where a
    correlation would be meaningless: too few shared observations, or
    either side being constant across them.

    Parameters
    ----------
    a, b : pandas.Series
        Aligned series, typically two isoform columns or the same gene on
        two platforms.
    min_overlap : int
        Minimum shared non-missing observations required.

    Returns
    -------
    float
        Spearman rho, or NaN where it could not be computed.
    """
    both = a.notna() & b.notna()
    if int(both.sum()) < min_overlap:
        return float("nan")
    x, y = a[both], b[both]
    if x.nunique() < 2 or y.nunique() < 2:
        return float("nan")
    with np.errstate(all="ignore"):
        rho, _ = spearmanr(x, y)
    return float(rho) if np.isfinite(rho) else float("nan")


def _min_pairwise_spearman(mat: pd.DataFrame, min_overlap: int) -> tuple[float, int]:
    """Min pairwise rho over isoform columns (min => all isoforms must agree).

    Taking the minimum rather than the mean makes the concordance test
    strict: one isoform behaving unlike the others drags the statistic
    down and sends the gene to the mean rollup, which is the conservative
    outcome.

    Parameters
    ----------
    mat : pandas.DataFrame
        Cell line x isoform matrix for a single gene.
    min_overlap : int
        Minimum shared cell lines per pair, passed to
        :func:`_pairwise_rho`.

    Returns
    -------
    min_rho : float
        Smallest finite pairwise rho, or NaN if no pair was assessable.
    n_pairs : int
        Number of pairs that yielded a finite correlation.
    """
    cols = list(mat.columns)
    rhos = []
    for i in range(len(cols)):
        for j in range(i + 1, len(cols)):
            r = _pairwise_rho(mat[cols[i]], mat[cols[j]], min_overlap)
            if np.isfinite(r):
                rhos.append(r)
    if not rhos:
        return float("nan"), 0
    return min(rhos), len(rhos)


def collapse_isoforms(long, source_name, *, min_overlap, corr_threshold):
    """
    One value per (model_id, gene_id).
      1 isoform                       -> pass through.
      >=2 concordant (min rho >= thr) -> MAX abundance per line.
      >=2 discordant / unassessable   -> MEAN per line (symmetric fallback).

    The max rollup assumes the isoforms are measuring the same thing, so
    the highest reading is the best-detected one. The mean fallback is
    used whenever that assumption fails — including when concordance
    could not be assessed at all — so an unrelated isoform cannot
    dominate the gene's value.

    Parameters
    ----------
    long : pandas.DataFrame
        Long frame for one platform, with ``model_id``, ``protein``,
        ``value``, ``gene_id``.
    source_name : str
        Platform name, recorded on both returned frames.
    min_overlap : int
        Minimum shared cell lines for an isoform pair to be assessed.
    corr_threshold : float
        Minimum pairwise rho for isoforms to count as concordant.

    Returns
    -------
    collapsed : pandas.DataFrame
        Columns ``model_id``, ``gene_id``, ``value``, ``n_isoforms``,
        ``source`` — one row per cell line and gene.
    qc : pandas.DataFrame
        One row per multi-isoform gene: ``source``, ``gene_id``,
        ``n_isoforms``, ``n_pairs_assessed``, ``min_pair_rho``,
        ``rollup``. Lets you see afterwards which genes took which path.

    Notes
    -----
    Single-isoform genes go through ``mean`` too, but over one value per
    cell line, so the aggregation is a no-op that also collapses any
    duplicate rows.
    """
    iso_counts = long.groupby("gene_id")["protein"].nunique()
    single = iso_counts[iso_counts <= 1].index
    multi = iso_counts[iso_counts > 1].index

    frames, qc = [], []

    if len(single):
        s = long[long["gene_id"].isin(single)]
        s1 = s.groupby(["model_id", "gene_id"])["value"].mean().reset_index()
        s1["n_isoforms"] = 1
        frames.append(s1[["model_id", "gene_id", "value", "n_isoforms"]])

    for gene in multi:
        g = long[long["gene_id"] == gene]
        mat = g.pivot_table(index="model_id", columns="protein", values="value", aggfunc="mean")
        min_rho, n_pairs = _min_pairwise_spearman(mat, min_overlap)
        concordant = bool(np.isfinite(min_rho) and (min_rho >= corr_threshold))
        val = mat.max(axis=1, skipna=True) if concordant else mat.mean(axis=1, skipna=True)

        gg = val.dropna().reset_index()
        gg.columns = ["model_id", "value"]
        gg["gene_id"] = gene
        gg["n_isoforms"] = int(mat.shape[1])
        frames.append(gg[["model_id", "gene_id", "value", "n_isoforms"]])
        qc.append({"source": source_name, "gene_id": gene, "n_isoforms": int(mat.shape[1]),
                   "n_pairs_assessed": int(n_pairs),
                   "min_pair_rho": float(min_rho) if np.isfinite(min_rho) else np.nan,
                   "rollup": "max" if concordant else "mean"})

    collapsed = (pd.concat(frames, ignore_index=True) if frames
                 else pd.DataFrame(columns=["model_id", "gene_id", "value", "n_isoforms"]))
    collapsed["source"] = source_name
    return collapsed, pd.DataFrame(qc, columns=["source", "gene_id", "n_isoforms",
                                                "n_pairs_assessed", "min_pair_rho", "rollup"])


# ---------------------------------------------------------------------
# Lineage-conditioned robust z (unchanged from base)
# ---------------------------------------------------------------------

def robust_z_matrix(X: np.ndarray) -> np.ndarray:
    """
    Column-wise robust z-scores using median and MAD.

    Centres each column on its median and scales by ``1.4826 * MAD``, the
    MAD-to-sigma constant for a normal distribution. Median/MAD replaces
    mean/SD because protein abundance across cell lines carries genuine
    extreme values that would otherwise inflate the denominator and
    flatten the signal being looked for.

    Parameters
    ----------
    X : numpy.ndarray
        2-D array, rows are cell lines and columns are genes.

    Returns
    -------
    numpy.ndarray
        Same shape, z-scored per column. An empty input is returned
        unchanged.

    Notes
    -----
    Two degenerate cases are handled so no column returns infinities: a
    zero MAD falls back to the standard deviation, and a zero scale after
    that becomes 1.0, leaving plain deviations from the median.
    """
    X = np.asarray(X, dtype=float)
    if X.size == 0:
        return X
    med = np.median(X, axis=0, keepdims=True)
    mad = np.median(np.abs(X - med), axis=0, keepdims=True)
    scale = np.where(np.isclose(mad, 0.0), np.std(X, axis=0, keepdims=True), 1.4826 * mad)
    scale = np.where(np.isclose(scale, 0.0), 1.0, scale)
    return (X - med) / scale


def zscore_platform(df, lineage_map, source_name) -> pd.DataFrame:
    """
    Compute lineage-conditioned robust z-scores for one platform.

    Pivots to a cell-line x gene matrix, attaches lineage, and z-scores
    each lineage group independently. Scoring within lineage means a
    value is judged against comparable cell lines, so lineage-typical
    abundance does not read as gene-specific signal.

    Parameters
    ----------
    df : pandas.DataFrame
        Collapsed long frame for one platform, with ``model_id``,
        ``gene_id``, ``value``.
    lineage_map : pandas.Series
        Lineage indexed by ``model_id``.
    source_name : str
        Platform name, recorded in ``source``.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_id``, ``model_id``, ``lineage``, ``z_score``,
        ``source``. Empty (with those columns) when nothing survives the
        lineage join.

    Notes
    -----
    Cell lines absent from ``lineage_map`` are dropped rather than pooled
    — the map's ``unknown`` fill only covers models it already knows
    about. Duplicate (model, gene) pairs resolve by ``aggfunc="first"``.
    """
    wide = df.pivot_table(index="model_id", columns="gene_id", values="value", aggfunc="first")
    wide["lineage"] = lineage_map.reindex(wide.index)
    wide = wide.dropna(subset=["lineage"])
    empty = pd.DataFrame(columns=["gene_id", "model_id", "lineage", "z_score", "source"])
    if wide.empty:
        return empty

    rows = []
    for lineage, grp in wide.groupby("lineage", sort=True):
        gene_cols = [c for c in grp.columns if c != "lineage"]
        if not gene_cols:
            continue
        Z = robust_z_matrix(grp[gene_cols].to_numpy(dtype=float))
        zdf = pd.DataFrame(Z, index=grp.index, columns=gene_cols).stack().reset_index()
        zdf.columns = ["model_id", "gene_id", "z_score"]
        zdf["lineage"] = lineage
        zdf["source"] = source_name
        rows.append(zdf)

    if not rows:
        return empty
    out = pd.concat(rows, ignore_index=True)
    return out[["gene_id", "model_id", "lineage", "z_score", "source"]]


# ---------------------------------------------------------------------
# CHANGE 2: symmetric platform agreement (no ProCAN preference)
# ---------------------------------------------------------------------

def gene_agreement(procan_long, ccle_long, min_shared) -> pd.DataFrame:
    """Per-gene ProCAN/CCLE Spearman on shared cell lines (symmetric).

    Both platforms are first reduced to one value per (cell line, gene),
    then merged on gene ID rather than on the raw protein label — so
    agreement can be measured even where the two platforms name the same
    protein differently.

    Parameters
    ----------
    procan_long, ccle_long : pandas.DataFrame
        Collapsed long frames, with ``model_id``, ``gene_id``, ``value``.
    min_shared : int
        Minimum shared cell lines before a gene's correlation is kept.
        Derived from statistical power, not hardcoded — see
        :func:`derive_thresholds`.

    Returns
    -------
    pandas.DataFrame
        Columns ``gene_id``, ``platform_rho``, ``n_shared``. Empty (with
        those columns) if the platforms share no observations.

    Notes
    -----
    Nothing here drops data. The correlation only feeds
    :func:`confidence_from_rho`, which shrinks the merged score —
    the symmetric alternative to the base version's CCLE removal.
    """
    p = procan_long.groupby(["model_id", "gene_id"])["value"].mean().rename("v_p").reset_index()
    c = ccle_long.groupby(["model_id", "gene_id"])["value"].mean().rename("v_c").reset_index()
    merged = p.merge(c, on=["model_id", "gene_id"], how="inner")
    if merged.empty:
        return pd.DataFrame(columns=["gene_id", "platform_rho", "n_shared"])

    counts = merged.groupby("gene_id").size()
    merged = merged[merged["gene_id"].isin(counts[counts >= min_shared].index)]

    rows = []
    for gene, grp in merged.groupby("gene_id"):
        with np.errstate(all="ignore"):
            rho, _ = spearmanr(grp["v_p"], grp["v_c"])
        if np.isfinite(rho):
            rows.append({"gene_id": gene, "platform_rho": float(rho), "n_shared": int(len(grp))})
    return pd.DataFrame(rows, columns=["gene_id", "platform_rho", "n_shared"])


def confidence_from_rho(rho, low, high) -> float:
    """rho -> [0,1] confidence. Symmetric: disagreement lowers trust in the
    merged gene, it does not pick a platform.

    A linear ramp: at or below ``low`` (the null level, where agreement is
    indistinguishable from unrelated genes) confidence is 0; at or above
    ``high`` it is 1; in between it interpolates.

    Parameters
    ----------
    rho : float
        Platform-agreement correlation for one gene. NaN means agreement
        could not be assessed.
    low, high : float
        Ramp bounds, derived by :func:`derive_thresholds`.

    Returns
    -------
    float
        Confidence in [0, 1].

    Notes
    -----
    An unassessed gene (NaN rho) returns 1.0, not 0.0 — a gene measured
    on one platform only is not penalised for the absence of a second
    opinion. A degenerate ramp (``high <= low``) falls back to a hard
    step at ``high``.
    """
    if not np.isfinite(rho):
        return 1.0
    if high <= low:
        return 1.0 if rho >= high else 0.0
    return float(np.clip((rho - low) / (high - low), 0.0, 1.0))


# ---------------------------------------------------------------------
# CHANGE 3: derive every threshold from the data
# ---------------------------------------------------------------------

def _protein_matrix(raw: pd.DataFrame):
    """
    Pivot a raw long frame to cell line x protein, with a protein -> gene map.

    Parameters
    ----------
    raw : pandas.DataFrame
        Uncollapsed long frame, with ``model_id``, ``protein``,
        ``value``, ``gene_id``.

    Returns
    -------
    mat : pandas.DataFrame
        Cell line x protein matrix, duplicates averaged.
    gene_of : dict
        Protein label to its ENSG ID — used to exclude same-gene pairs
        when sampling the null.
    """
    mat = raw.pivot_table(index="model_id", columns="protein", values="value", aggfunc="mean")
    gene_of = raw.drop_duplicates("protein").set_index("protein")["gene_id"].to_dict()
    return mat, gene_of


def _gene_matrix(raw: pd.DataFrame) -> pd.DataFrame:
    """
    Pivot a raw long frame to a cell line x gene matrix.

    Isoforms are averaged into one value per gene, which is enough for
    threshold derivation — the full concordance-aware rollup happens
    later in :func:`collapse_isoforms`.

    Parameters
    ----------
    raw : pandas.DataFrame
        Uncollapsed long frame.

    Returns
    -------
    pandas.DataFrame
        Cell line x gene matrix.
    """
    lvl = raw.groupby(["model_id", "gene_id"])["value"].mean().reset_index()
    return lvl.pivot_table(index="model_id", columns="gene_id", values="value")


def _power_min_n(target_r: float, alpha: float, nmax: int = 300) -> int:
    """Smallest n at which |rho| = target_r is significant at alpha (two-sided).

    Walks n upward until the critical rho at that sample size falls to or
    below ``target_r``. Used to turn a correlation the pipeline cares
    about into a minimum number of shared cell lines needed to detect it.

    Parameters
    ----------
    target_r : float
        Correlation to be detectable. Clipped into (0, 1).
    alpha : float
        Two-sided significance level.
    nmax : int, optional
        Search ceiling. Default 300.

    Returns
    -------
    int
        Smallest sufficient n, or ``nmax`` if the target is not reachable
        within the ceiling.
    """
    target_r = float(np.clip(abs(target_r), 1e-3, 0.999))
    for n in range(4, nmax + 1):
        tc = t_dist.ppf(1 - alpha / 2, n - 2)
        rc = tc / np.sqrt(n - 2 + tc ** 2)         # critical Spearman rho at this n
        if rc <= target_r:
            return n
    return nmax


def _sample_null_protein_rho(raw_longs, seed_overlap, rng, n_pairs):
    """
    Sample correlations between protein pairs belonging to different genes.

    This is the null distribution: how correlated two proteins look when
    there is no reason for them to track each other. Same-gene pairs are
    excluded, since those are exactly the signal the isoform test is
    trying to detect.

    Parameters
    ----------
    raw_longs : dict
        Source name to uncollapsed long frame. Sampled across all
        sources.
    seed_overlap : int
        Minimum shared cell lines per pair.
    rng : numpy.random.Generator
        Seeded generator, so derivation is reproducible.
    n_pairs : int
        Draws attempted per source. Fewer are kept, since same-gene and
        self pairs are skipped and some correlations are not finite.

    Returns
    -------
    numpy.ndarray
        Finite null correlations. May be short of ``n_pairs``.
    """
    null = []
    for raw in raw_longs.values():
        mat, gene_of = _protein_matrix(raw)
        cols = list(mat.columns)
        if len(cols) < 4:
            continue
        genes = np.array([gene_of.get(c) for c in cols], dtype=object)
        for _ in range(n_pairs):
            i, j = rng.integers(0, len(cols), size=2)
            if i == j or genes[i] == genes[j]:
                continue
            r = _pairwise_rho(mat[cols[i]], mat[cols[j]], seed_overlap)
            if np.isfinite(r):
                null.append(r)
    return np.asarray(null, dtype=float)


def derive_thresholds(raw_longs, lineage_map=None, *, alpha, null_pctile,
                      matched_pctile, min_gap, seed=0, n_pairs=1500, verbose=True):
    """
    Derive every threshold the pipeline needs from the data itself.

    Nothing here is a domain-tuned constant. Three quantities are
    derived, each against an empirical reference distribution:

    * **Isoform concordance** (``isoform_corr``) — the ``null_pctile``
      percentile of correlations between proteins of *different* genes.
      Isoforms must agree better than that many unrelated pairs before
      the max rollup is used.
    * **Platform-agreement ramp** (``agree_low``, ``agree_high``) —
      ``low`` from the null of mismatched gene pairs across the two
      platforms, ``high`` from the ``matched_pctile`` percentile of
      same-gene pairs. ``min_gap`` guarantees the ramp has width, so
      confidence never collapses to a hard step.
    * **Minimum shared lines** (``min_shared``) — the sample size at
      which a correlation of ``low`` would be significant at ``alpha``.

    Parameters
    ----------
    raw_longs : dict
        Source name to uncollapsed long frame.
    lineage_map : pandas.Series, optional
        Accepted but unused; thresholds are derived across all lineages.
    alpha : float
        Two-sided significance level.
    null_pctile : float
        Percentile of the null distribution treated as "agrees".
    matched_pctile : float
        Percentile of matched agreement treated as full trust.
    min_gap : float
        Minimum width of the confidence ramp.
    seed : int, optional
        RNG seed for the null sampling. Default 0.
    n_pairs : int, optional
        Random pairs attempted per null. Default 1500.
    verbose : bool, optional
        Print the derived values. Default True.

    Returns
    -------
    dict
        ``isoform_corr``, ``isoform_min_overlap``, ``agree_low``,
        ``agree_high``, ``min_shared``, ``seed_overlap``.

    Notes
    -----
    Each derivation has a significance-floor fallback via :func:`t_crit`,
    used when too few pairs were assessable — fewer than 30 null samples
    for the isoform threshold, or fewer than 10 matched / 30 null pairs
    for the ramp. A run that falls back is not wrong, but its thresholds
    are conventional rather than empirical; the printed values are worth
    checking against a run that did not fall back.
    """
    rng = np.random.default_rng(seed)
    # smallest statistically valid overlap: the fewest lines at which even a very
    # strong correlation clears significance. Used only to build the null samples.
    seed_overlap = _power_min_n(0.9, alpha)

    # --- isoform concordance threshold: null of cross-gene protein correlations
    null_iso = _sample_null_protein_rho(raw_longs, seed_overlap, rng, n_pairs)
    if null_iso.size >= 30:
        iso_thr = float(np.clip(np.percentile(null_iso, null_pctile), 0.0, 0.95))
    else:
        iso_thr = float(t_crit(seed_overlap, alpha))   # significance floor fallback

    # --- platform-agreement bounds: null (mismatched genes) vs matched (same gene)
    low = high = float("nan")
    if PROCAN_SOURCE in raw_longs and CCLE_SOURCE in raw_longs:
        P, C = _gene_matrix(raw_longs[PROCAN_SOURCE]), _gene_matrix(raw_longs[CCLE_SOURCE])
        shared = [g for g in P.columns if g in C.columns]
        matched = [r for g in shared if np.isfinite(r := _pairwise_rho(P[g], C[g], seed_overlap))]
        null_x = []
        pc, cc = list(P.columns), list(C.columns)
        for _ in range(n_pairs):
            g, h = pc[rng.integers(0, len(pc))], cc[rng.integers(0, len(cc))]
            if g == h:
                continue
            r = _pairwise_rho(P[g], C[h], seed_overlap)
            if np.isfinite(r):
                null_x.append(r)
        if len(matched) >= 10 and len(null_x) >= 30:
            low = float(np.clip(np.percentile(null_x, null_pctile), 0.0, 0.9))
            high = float(np.clip(max(low + min_gap, np.percentile(matched, matched_pctile)),
                                 low + min_gap, 0.95))
    if not (np.isfinite(low) and np.isfinite(high)):
        low = float(t_crit(seed_overlap, alpha))     # significance-floor fallback
        high = float(min(low + max(min_gap, low), 0.9))

    # --- minimum shared lines to trust an agreement estimate: power at `low`
    min_shared = int(_power_min_n(low, alpha))

    thr = {
        "isoform_corr": iso_thr,
        "isoform_min_overlap": int(seed_overlap),
        "agree_low": low,
        "agree_high": high,
        "min_shared": min_shared,
        "seed_overlap": int(seed_overlap),
    }
    if verbose:
        print("Derived thresholds (data-driven; conventions: "
              f"alpha={alpha}, null_pctile={null_pctile}, matched_pctile={matched_pctile}):")
        for k, v in thr.items():
            print(f"  {k:<20} = {v}")
    return thr


def t_crit(n: int, alpha: float) -> float:
    """Critical Spearman rho at sample size n and level alpha (two-sided).

    The smallest correlation that would be significant at ``n``
    observations. Used as the significance-floor fallback when a
    threshold cannot be derived empirically.

    Parameters
    ----------
    n : int
        Sample size. Floored at 4.
    alpha : float
        Two-sided significance level.

    Returns
    -------
    float
        Critical rho.
    """
    n = max(int(n), 4)
    tc = t_dist.ppf(1 - alpha / 2, n - 2)
    return tc / np.sqrt(n - 2 + tc ** 2)


# ---------------------------------------------------------------------
# Core computation
# ---------------------------------------------------------------------

def compute_from_sources(raw_longs, lineage_map, thr) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Score every gene and cell line from the loaded sources, given thresholds.

    Five steps:

    1. Collapse isoforms per source, keeping the QC record of which genes
       took the max path and which the mean.
    2. Measure symmetric per-gene platform agreement and convert it to a
       confidence in [0, 1]. No rows are dropped.
    3. Z-score each platform independently, conditioned on lineage.
    4. Combine platforms with coverage weights ``w = sqrt(n_models)``, as
       ``sum(w*z) / sqrt(sum(w^2))``, which preserves the z-scale. Both
       platforms are treated identically — neither is preferred.
    5. Scale the combined score by platform confidence, so disagreement
       shrinks the merged gene toward zero instead of deleting a
       platform's data. Both the scaled ``z_score`` and the unscaled
       ``z_raw`` are returned, so the effect is auditable.

    Parameters
    ----------
    raw_longs : dict
        Source name to uncollapsed long frame.
    lineage_map : pandas.Series
        Lineage indexed by ``model_id``.
    thr : dict
        Thresholds from :func:`derive_thresholds` or loaded from JSON.

    Returns
    -------
    result : pandas.DataFrame
        Columns ``gene_id``, ``model_id``, ``lineage``, ``z_score``,
        ``z_raw``, ``platform_confidence``, ``n_sources``,
        ``n_isoforms``, sorted by gene then model with IDs uppercased.
    isoform_qc : pandas.DataFrame
        Concatenated isoform rollup decisions across sources.

    Raises
    ------
    ValueError
        If no rows survive the isoform collapse, or if nothing could be
        scored.

    Notes
    -----
    Genes with no agreement estimate — measured on one platform, or too
    few shared lines — get confidence 1.0, so they are neither penalised
    nor rewarded for the missing second opinion. ``n_isoforms`` is the
    max across sources, so it reflects the most fragmented platform.
    """
    # 1. isoform collapse per source
    source_longs, qc_frames, iso_counts = {}, [], []
    for name, raw in raw_longs.items():
        collapsed, qc = collapse_isoforms(
            raw, name, min_overlap=thr["isoform_min_overlap"], corr_threshold=thr["isoform_corr"]
        )
        if collapsed.empty:
            continue
        source_longs[name] = collapsed
        iso_counts.append(collapsed[["gene_id", "n_isoforms"]])
        if not qc.empty:
            qc_frames.append(qc)
    if not source_longs:
        raise ValueError("No protein rows after isoform collapse.")

    isoform_qc = pd.concat(qc_frames, ignore_index=True) if qc_frames else pd.DataFrame()
    iso_by_gene = (pd.concat(iso_counts, ignore_index=True).groupby("gene_id")["n_isoforms"].max()
                   if iso_counts else pd.Series(dtype=int))

    # 2. symmetric per-gene agreement -> confidence (no platform preferred)
    if PROCAN_SOURCE in source_longs and CCLE_SOURCE in source_longs:
        agree = gene_agreement(source_longs[PROCAN_SOURCE], source_longs[CCLE_SOURCE], thr["min_shared"])
    else:
        agree = pd.DataFrame(columns=["gene_id", "platform_rho", "n_shared"])
    conf_by_gene = (
        agree.assign(confidence=agree["platform_rho"].map(
            lambda r: confidence_from_rho(r, thr["agree_low"], thr["agree_high"])))
        .set_index("gene_id")["confidence"]
        if not agree.empty else pd.Series(dtype=float)
    )

    # 3. lineage-conditioned robust z per platform
    blocks = [z for name, long in source_longs.items()
              if not (z := zscore_platform(long, lineage_map, name)).empty]
    if not blocks:
        raise ValueError("No protein rows could be scored.")
    all_z = pd.concat(blocks, ignore_index=True)

    # 4. SYMMETRIC coverage-weighted combination (identical treatment of both platforms)
    src_n = all_z.groupby(["gene_id", "source"])["model_id"].nunique().rename("src_n").reset_index()
    all_z = all_z.merge(src_n, on=["gene_id", "source"], how="left")
    all_z["w"] = np.sqrt(all_z["src_n"].astype(float))
    all_z["wz"] = all_z["w"] * all_z["z_score"]
    all_z["w2"] = all_z["w"] ** 2

    agg = all_z.groupby(["gene_id", "model_id", "lineage"], as_index=False).agg(
        wz_sum=("wz", "sum"), w2_sum=("w2", "sum"), n_sources=("source", "nunique"))
    denom = np.sqrt(agg["w2_sum"].to_numpy(dtype=float))
    z_comb = np.where(denom > 0, agg["wz_sum"].to_numpy(dtype=float) / denom, 0.0)

    result = agg[["gene_id", "model_id", "lineage"]].copy()
    result["z_raw"] = z_comb.astype(float)
    result["n_sources"] = agg["n_sources"].astype(int)

    # 5. disagreement lowers confidence in the merged gene (symmetric, not a drop)
    result["platform_confidence"] = result["gene_id"].map(conf_by_gene).fillna(1.0).astype(float)
    result["z_score"] = result["z_raw"] * result["platform_confidence"]
    result["n_isoforms"] = result["gene_id"].map(iso_by_gene).fillna(1).astype(int)

    result = result.dropna(subset=["gene_id", "model_id", "lineage", "z_score"]).copy()
    result["gene_id"] = result["gene_id"].astype(str).str.upper()
    result["model_id"] = result["model_id"].astype(str).str.upper()
    result["lineage"] = result["lineage"].astype(str)
    result = result.sort_values(["gene_id", "model_id"]).reset_index(drop=True)
    result = result[["gene_id", "model_id", "lineage", "z_score", "z_raw",
                     "platform_confidence", "n_sources", "n_isoforms"]]
    return result, isoform_qc


def compute_protein_z(con, table_names=None, *, alpha=DEFAULT_ALPHA,
                      null_pctile=DEFAULT_NULL_PCTILE, matched_pctile=DEFAULT_MATCHED_PCTILE,
                      min_gap=DEFAULT_MIN_GAP, seed=0, thresholds=None, verbose=True):
    """
    Run the full protein scoring, end to end.

    Builds the identifier maps, loads the lineage map and every available
    source table, derives the thresholds (unless supplied), then scores.
    Passing ``thresholds`` skips derivation entirely, which is what makes
    a run reproducible from a saved JSON.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection to the project database.
    table_names : list of str, optional
        Source tables. Defaults to ``[CCLE_SOURCE, PROCAN_SOURCE]``.
    alpha, null_pctile, matched_pctile, min_gap : float, optional
        Statistical conventions passed through to
        :func:`derive_thresholds`. Ignored when ``thresholds`` is given.
    seed : int, optional
        RNG seed for threshold derivation. Default 0.
    thresholds : dict, optional
        Pre-derived thresholds; skips derivation when supplied.
    verbose : bool, optional
        Print the derived thresholds. Default True.

    Returns
    -------
    result : pandas.DataFrame
        Scored rows, as from :func:`compute_from_sources`.
    isoform_qc : pandas.DataFrame
        Isoform rollup decisions.
    thr : dict
        The thresholds used, derived or supplied — return them so a
        caller can persist them for a reproducible rerun.

    Raises
    ------
    ValueError
        If no identifier mapping could be built from ``gene_roster``, or
        if no source table could be loaded.

    Notes
    -----
    Symbols are merged first and UniProt accessions over them, so an
    accession wins where a token appears in both identifier spaces.
    """
    if table_names is None:
        table_names = [CCLE_SOURCE, PROCAN_SOURCE]
    protein_map = {**build_symbol_to_ensg(con), **build_uniprot_to_ensg(con)}
    if not protein_map:
        raise ValueError("No UniProt/symbol -> ENSG mapping from gene_roster.")
    lineage_map = load_lineage_map(con)
    raw_longs = load_sources(con, table_names, protein_map)
    if not raw_longs:
        raise ValueError("No protein source tables could be loaded.")

    thr = thresholds or derive_thresholds(
        raw_longs, lineage_map, alpha=alpha, null_pctile=null_pctile,
        matched_pctile=matched_pctile, min_gap=min_gap, seed=seed, verbose=verbose)
    result, isoform_qc = compute_from_sources(raw_longs, lineage_map, thr)
    result.attrs_thr = thr  # for callers that want the derived values
    return result, isoform_qc, thr


# ---------------------------------------------------------------------
# Persistence + CLI
# ---------------------------------------------------------------------

def _write(con, df, table):
    """
    Write a frame to DuckDB, replacing any existing table of that name.

    Registers the frame as a temporary view and materialises it. A None
    or empty frame is a silent no-op, so an absent QC table does not
    leave an empty artefact behind.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open, writable connection.
    df : pandas.DataFrame or None
        Frame to write.
    table : str
        Destination table name.

    Returns
    -------
    None
    """
    if df is None or df.empty:
        return
    con.register(f"{table}_tmp", df.copy())
    con.execute(f'CREATE OR REPLACE TABLE "{table}" AS SELECT * FROM "{table}_tmp"')


def write_protein_z(con, df, output_table="protein_z"):
    """
    Write the scored results to DuckDB with a fixed column order.

    Selecting an explicit column list means the output schema stays
    stable regardless of extra columns the caller's frame carries.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open, writable connection.
    df : pandas.DataFrame
        Result frame from :func:`compute_protein_z`.
    output_table : str, optional
        Destination table. Defaults to ``"protein_z"``.

    Returns
    -------
    None

    Notes
    -----
    Replaces the table outright — there is no append path and no
    versioning, so save the thresholds if you need to reproduce a run.
    """
    cols = ["gene_id", "model_id", "lineage", "z_score", "z_raw",
            "platform_confidence", "n_sources", "n_isoforms"]
    _write(con, df[cols], output_table)


def parse_args():
    """
    Define and parse the command-line interface.

    Returns
    -------
    argparse.Namespace
        Parsed arguments:

        * ``db`` (str) — DuckDB file path;
        * ``tables`` (list of str) — protein source tables;
        * ``output_table``, ``isoform_qc_table`` (str) — destinations;
        * ``params`` (str or None) — JSON of pre-derived thresholds, to
          skip derivation;
        * ``save_params`` (str or None) — write the thresholds used to
          this JSON;
        * ``seed`` (int) — RNG seed for derivation;
        * ``alpha``, ``null_pctile``, ``matched_pctile``, ``min_gap``
          (float) — the statistical conventions, exposed rather than
          buried so they can be seen and overridden.

    Notes
    -----
    ``--params`` and ``--save-params`` are the reproducibility pair:
    save from one run, pass into the next to score against identical
    thresholds.
    """
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default=str(DEFAULT_DB))
    ap.add_argument("--tables", nargs="*", default=[CCLE_SOURCE, PROCAN_SOURCE])
    ap.add_argument("--output-table", default="protein_z")
    ap.add_argument("--isoform-qc-table", default="protein_isoform_qc")
    ap.add_argument("--params", default=None, help="JSON of pre-derived thresholds (skip derivation)")
    ap.add_argument("--save-params", default=None, help="write derived thresholds to this JSON")
    ap.add_argument("--seed", type=int, default=0)
    # statistical conventions (not domain values) -- exposed, not buried
    ap.add_argument("--alpha", type=float, default=DEFAULT_ALPHA)
    ap.add_argument("--null-pctile", type=float, default=DEFAULT_NULL_PCTILE)
    ap.add_argument("--matched-pctile", type=float, default=DEFAULT_MATCHED_PCTILE)
    ap.add_argument("--min-gap", type=float, default=DEFAULT_MIN_GAP)
    return ap.parse_args()


def main():
    """
    Run the protein scoring from the command line.

    Opens the database, loads thresholds from ``--params`` if given,
    scores the requested sources, writes the results and the isoform QC
    table, optionally saves the thresholds used, then prints a summary
    read back from the written table.

    The summary covers row, gene, model and lineage counts, how many rows
    rest on one platform versus two, how many involve multi-isoform
    genes, and the mean platform confidence — which together show how
    much of the output was affected by the two headline changes.

    Returns
    -------
    None
        The output tables, optional threshold JSON and console summary
        are the side effects.

    Raises
    ------
    ValueError
        Propagated from :func:`compute_protein_z` when no mapping or no
        source rows are available.

    Notes
    -----
    The connection is closed in a ``finally`` block, so it is released
    even when scoring raises. The summary is read back from the database
    rather than from the in-memory frame, which verifies the write landed.
    """
    args = parse_args()
    con = duckdb.connect(args.db, read_only=False)
    try:
        thresholds = json.loads(Path(args.params).read_text()) if args.params else None
        result, isoform_qc, thr = compute_protein_z(
            con, table_names=list(args.tables), alpha=args.alpha,
            null_pctile=args.null_pctile, matched_pctile=args.matched_pctile,
            min_gap=args.min_gap, seed=args.seed, thresholds=thresholds)

        write_protein_z(con, result, args.output_table)
        _write(con, isoform_qc, args.isoform_qc_table)
        if args.save_params:
            Path(args.save_params).write_text(json.dumps(thr, indent=2))

        summary = con.execute(
            f"""
            SELECT COUNT(*) AS n_rows, COUNT(DISTINCT gene_id) AS n_genes,
                   COUNT(DISTINCT model_id) AS n_models, COUNT(DISTINCT lineage) AS n_lineages,
                   COUNT(*) FILTER (WHERE n_sources = 1) AS n_sources_1,
                   COUNT(*) FILTER (WHERE n_sources = 2) AS n_sources_2,
                   COUNT(*) FILTER (WHERE n_isoforms > 1) AS n_multi_isoform_rows,
                   ROUND(AVG(platform_confidence), 3) AS avg_confidence
            FROM "{args.output_table}"
            """
        ).fetchdf()
        print(summary.to_string(index=False))
    finally:
        con.close()


if __name__ == "__main__":
    main()