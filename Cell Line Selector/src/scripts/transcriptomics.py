#!/usr/bin/env python3
"""
transcriptomics.py

Function library for the lineage-stratified transcriptomics z-score layer.

Sources / methods
    hpa_rna     (long)  RNA-seq, TPM                                -> 1 method
    depmap_expr (wide)  RNA-seq, log2(TPM+1)                        -> 1 method
    geo_expr    (wide)  arrays; method = platform_id x proc_family  -> N methods

Design notes
    - Lineage is the reference population for the z, not an output key.
      One cell line has one lineage.
    - GEO methods are split on platform AND processing family: two samples on
      the same GPL processed by MAS5 vs RMA sit on different scales and are
      not commensurable.
    - Censoring cutoffs, the lineage-shrinkage constant n0 and the array floor
      quantile are selected empirically (see calibrate_params), not hardcoded.
    - Flag don't drop: unmapped samples get their own method and a status
      flag; below-detection cells get NaN + a status, never a fabricated 0.
    - Output is gene-major (all cell lines for gene 1, then gene 2, ...) so
      DuckDB zone maps can skip row groups on `ensg`.

Output goes to DuckDB only. At ~30M rows CSV is not a viable target.

No side effects at import. All orchestration lives in
05_transcriptomics_stat_layer.py.

Module layout
-------------
1. GEO SOFT metadata -> method assignment
2. Schema discovery
3. Axes: method rows -> model_id -> lineage
4. Chunked fetch
5. Scale detection and detection floors
6. Z-score core
7. Sensitivity analysis
8. Assembling and writing output
9. Plots
"""

from __future__ import annotations

import itertools
import logging
import re
from collections import defaultdict
from dataclasses import dataclass, field, asdict
from pathlib import Path
import warnings

import numpy as np
import pandas as pd
from scipy.stats import norm, rankdata, spearmanr  # type: ignore

log = logging.getLogger("transcriptomics")

# ----------------------------------------------------------------------
# Constants
# ----------------------------------------------------------------------

TPM_FLOOR = 1.0          # HPA's own >=1 TPM detection criterion (Uhlen 2015)
EPS = 1e-9
GENE_CHUNK = 500         # default chunk size for full-power runs
LOW_MEMORY_GENE_CHUNK = 128  # safer default for laptops / warm CPUs
N_CALIB_GENES = 1200
MIN_CALIB_CELLS = 200    # min comparable (line, gene) pairs to score a combo

#: Integer status code -> label, written to the output's ``status`` column.
#: 0 parametric z; 1 rank inverse-normal under censoring; 2 below the
#: detection floor (NaN, never 0); 3 nothing measured.
STATUS = {0: "ok", 1: "censored_ranked", 2: "below_detection", 3: "no_data"}

#: Integer scale-source code -> label, written to ``scale_source``. Records
#: which population supplied the z denominator: the lineage alone, the
#: lineage shrunk toward the method-global scale, or the global fallback
#: used when the lineage was too thin to trust.
SCALE_SRC = {0: "none", 1: "lineage", 2: "shrunk", 3: "fallback_global"}

#: Parameter values swept by the sensitivity analysis. Combinations where
#: ``lower_cut >= upper_cut`` are discarded in :func:`grid_combinations`.
PARAM_GRID = dict(
    floor_q=[0.02, 0.05, 0.10, 0.20],
    n0=[5, 10, 25, 50, 100],
    upper_cut=[0.30, 0.40, 0.50, 0.60, 0.70, 0.80],
    lower_cut=[0.02, 0.05, 0.10, 0.20],
)

# Reference params for the held-out source during calibration. Held fixed so
# the reference does not move with the parameters being scored (circularity).
REF_PARAMS_RAW = dict(floor_q=0.10, n0=25, upper_cut=0.50, lower_cut=0.10)

# First match wins. gcRMA before RMA; MAS5 cues include the GCOS/global-scaling
# language that Affymetrix pipelines emit.
PROC_RULES = [
    ("gcrma",    r"\bgc[\s\-]?rma\b",                                   "log"),
    ("rma",      r"\brma\b|robust\s+multi[\s\-]?array",                 "log"),
    ("mas5",     r"\bmas\s*5|\bgcos\b|microarray\s+suite\s*5|"
                 r"trimmed\s+mean.*target|global\s+scaling",            "linear"),
    ("plier",    r"\bplier\b",                                          "log"),
    ("dchip",    r"\bdchip\b|model[\s\-]based\s+expression",            "linear"),
    ("mas4",     r"\bmas\s*4|average\s+difference",                     "linear"),
    ("vsn",      r"\bvsn\b|variance\s+stabilis?z?ing",                  "log"),
    ("quantile", r"quantile\s+normal",                                  "unknown"),
    ("lowess",   r"\blowess\b|\bloess\b",                               "log"),
    ("rnaseq",   r"\brsem\b|\btpm\b|\bfpkm\b|\brpkm\b|\bdeseq\b",       "linear"),
]

_SOFT_KEY = re.compile(r"^!Sample_([A-Za-z0-9_]+)\s*=\s*(.*)$")
_MULTI_KEYS = {"characteristics_ch1", "data_processing", "series_id",
               "supplementary_file", "relation", "contact_name"}
_ENSG = re.compile(r"(?i)^(ENSG\d+)(\.\d+)?$")


# ----------------------------------------------------------------------
# Data containers (no behaviour)
# ----------------------------------------------------------------------

@dataclass
class Params:
    """
    The four empirically selected scoring parameters.

    Attributes
    ----------
    floor_q : float
        Quantile of a method's own value distribution used as its detection
        floor. Applies to arrays only; RNA-seq sources use the fixed TPM
        criterion instead.
    n0 : int
        Lineage shrinkage constant. A lineage of ``n`` cell lines gets weight
        ``n / (n + n0)`` on its own scale, the remainder on the method-global
        scale, so small lineages borrow strength rather than producing a
        noisy MAD.
    upper_cut : float
        Detection fraction at or above which the parametric robust z is used.
    lower_cut : float
        Detection fraction below which a cell is marked below_detection and
        returns NaN. Between the two cuts, the rank inverse-normal applies.

    Notes
    -----
    Chosen by :func:`calibrate_params`, not set by hand. Serialised to
    ``chosen_params.json`` so a scoring run can be reproduced exactly.
    """
    floor_q: float
    n0: int
    upper_cut: float
    lower_cut: float


@dataclass
class MethodAxis:
    """Row axis for one method: raw sample ids, their model_id, replicate map.

    One method is one commensurable measurement scale — a whole source for
    RNA-seq, or a single platform x processing combination for GEO.

    Attributes
    ----------
    name : str
        Unique axis name, e.g. ``geo_GPL570|rma``.
    source : str
        Parent source: ``hpa_rna``, ``depmap`` or ``geo``. Determines the
        vote weight at the final combination stage.
    table, id_col : str
        Physical DuckDB table and the column holding its raw sample IDs.
    units : str
        Declared units — from the table for RNA-seq, from SOFT metadata for
        GEO. A hint only; :func:`profile_method` trusts the data over it.
    raw_ids, model_ids : numpy.ndarray
        Parallel arrays, one entry per raw sample row.
    models : numpy.ndarray
        Unique ``model_id`` values, the axis's collapsed row order.
    inv : numpy.ndarray
        Index from each raw row into ``models``; drives replicate collapse.
    n_samples : numpy.ndarray
        Raw rows per model, used as the sqrt-n weight within GEO.
    floor : float
        Detection floor in log2 space. Set by :func:`profile_method`.
    scale_detected : str
        ``"log"``, ``"linear"`` or ``"unknown"``, set from the data by
        :func:`profile_method`.

    Notes
    -----
    Fields from ``models`` down are populated by :func:`make_axis` and
    :func:`profile_method`, not at construction.
    """
    name: str
    source: str
    table: str
    id_col: str
    units: str
    raw_ids: np.ndarray
    model_ids: np.ndarray
    models: np.ndarray = field(default=None)
    inv: np.ndarray = field(default=None)
    n_samples: np.ndarray = field(default=None)
    floor: float = None
    scale_detected: str = None


# ======================================================================
# 1. GEO SOFT metadata -> method assignment
# ======================================================================

def parse_soft_file(path: Path) -> dict:
    """Parse one GSMxxxx.txt SOFT sample record into a flat dict.

    Keys in :data:`_MULTI_KEYS` legitimately repeat within a record and are
    joined with a pipe; every other key keeps its first value.

    Parameters
    ----------
    path : Path
        SOFT sample file. Decoding errors are replaced rather than raised,
        since a stray byte should not lose a whole record.

    Returns
    -------
    dict
        ``!Sample_<key>`` names mapped to their values. If no
        ``geo_accession`` line was present, one is recovered from the
        filename where possible.
    """
    rec = defaultdict(list)
    with open(path, errors="replace") as fh:
        for line in fh:
            m = _SOFT_KEY.match(line.rstrip("\n"))
            if m:
                rec[m.group(1)].append(m.group(2).strip())
    out = {k: (" | ".join(v) if k in _MULTI_KEYS else v[0]) for k, v in rec.items()}
    if "geo_accession" not in out:
        m = re.search(r"(GSM\d+)", path.stem, re.I)
        if m:
            out["geo_accession"] = m.group(1).upper()
    return out


def classify_processing(text: str) -> tuple[str, str]:
    """!Sample_data_processing -> (proc_family, scale_hint).

    Matches :data:`PROC_RULES` in order, first match wins — gcRMA is tested
    before RMA so the more specific label is not swallowed by the general
    one.

    Parameters
    ----------
    text : str
        Free-text processing description from the SOFT record.

    Returns
    -------
    proc_family : str
        Pipeline family, or ``"unknown"`` if nothing matched.
    scale_hint : str
        ``"log"``, ``"linear"`` or ``"unknown"``. A hint only — the scale
        actually used is detected from the values in
        :func:`profile_method`.
    """
    t = (text or "").lower()
    for name, pat, scale in PROC_RULES:
        if re.search(pat, t):
            return name, scale
    return "unknown", "unknown"


def build_geo_platform_frame(meta_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse geo_meta/*.txt -> (records, parse_failures).

    Each sample's method is ``platform_id|proc_family``, the unit of
    commensurability for GEO.

    Parameters
    ----------
    meta_dir : Path
        Directory of SOFT files. ``GSM*.txt`` is tried first, falling back
        to all ``*.txt``.

    Returns
    -------
    records : pandas.DataFrame
        One row per sample: ``gsm_id``, ``platform_id``, ``proc_family``,
        ``scale_hint``, ``method``, plus series and description fields.
        Duplicate accessions are collapsed to the first occurrence.
    failures : pandas.DataFrame
        Columns ``file``, ``reason``, for records that could not be used.

    Raises
    ------
    FileNotFoundError
        If no ``.txt`` files are found under ``meta_dir``.

    Notes
    -----
    Flag don't drop: a record missing ``platform_id`` is still kept, under
    ``UNKNOWN_GPL``, and also listed in ``failures``. Only a missing
    accession is fatal to the record, since without it nothing can be
    joined.
    """
    files = sorted(meta_dir.glob("GSM*.txt")) or sorted(meta_dir.glob("*.txt"))
    if not files:
        raise FileNotFoundError(f"no SOFT .txt files under {meta_dir}")
    log.info("parsing %d SOFT files from %s", len(files), meta_dir)

    rows, bad = [], []
    for i, f in enumerate(files, 1):
        r = parse_soft_file(f)
        gsm = r.get("geo_accession", "").upper()
        plat = r.get("platform_id", "").strip().upper()
        if not gsm:
            bad.append((f.name, "no accession"))
            continue
        if not plat:
            bad.append((f.name, "no platform_id"))
            plat = "UNKNOWN_GPL"
        proc, scale = classify_processing(r.get("data_processing", ""))
        rows.append(dict(
            gsm_id=gsm, platform_id=plat, proc_family=proc, scale_hint=scale,
            method=f"{plat}|{proc}",
            gse_id=r.get("series_id", ""), title=r.get("title", ""),
            source_name=r.get("source_name_ch1", ""),
            molecule=r.get("molecule_ch1", ""),
            sample_type=r.get("type", ""),
            channel_count=r.get("channel_count", ""),
            characteristics=r.get("characteristics_ch1", ""),
            data_processing=(r.get("data_processing", "") or "")[:2000],
            src_file=f.name))
        if i % 2000 == 0:
            log.info("  %d/%d", i, len(files))

    df = pd.DataFrame(rows)
    dup = int(df.gsm_id.duplicated().sum()) if len(df) else 0
    if dup:
        log.warning("%d duplicate GSM records -> keeping first", dup)
        df = df.drop_duplicates("gsm_id", keep="first")
    return df, pd.DataFrame(bad, columns=["file", "reason"])


def write_geo_platform(con, df: pd.DataFrame, table: str = "geo_platform") -> None:
    """
    Write the GEO method assignments to DuckDB, replacing any existing table.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open, writable connection.
    df : pandas.DataFrame
        Records frame from :func:`build_geo_platform_frame`.
    table : str, optional
        Destination table. Defaults to ``"geo_platform"``.

    Returns
    -------
    None

    Notes
    -----
    The temporary registration is dropped in a ``finally`` block, so a
    failed write does not leave a stale view bound to the connection.
    """
    con.register("_geo_platform_tmp", df)
    try:
        con.execute(
            f'CREATE OR REPLACE TABLE "{table}" AS SELECT * FROM _geo_platform_tmp')
    finally:
        con.unregister("_geo_platform_tmp")
    log.info("wrote %s (%d rows)", table, len(df))


def audit_geo_platform(con, df: pd.DataFrame) -> pd.DataFrame:
    """Before/after coverage audit. Returns the method summary table.

    Reports three things worth seeing before committing to an assignment:
    how many samples each method holds, which processing strings failed to
    classify, and how many ``geo_expr`` samples have no metadata record at
    all.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection, used for the coverage check only.
    df : pandas.DataFrame
        Records frame from :func:`build_geo_platform_frame`.

    Returns
    -------
    pandas.DataFrame
        Sample counts by ``platform_id``, ``proc_family`` and
        ``scale_hint``, descending.

    Notes
    -----
    The ``geo_expr`` coverage check is best-effort: if the table is absent
    the check is logged as skipped rather than raising, so an audit can run
    before the expression tables are loaded.
    """
    summary = (df.groupby(["platform_id", "proc_family", "scale_hint"])
                 .size().rename("n_gsm").reset_index()
                 .sort_values("n_gsm", ascending=False))
    log.info("methods (platform x processing):\n%s", summary.to_string(index=False))

    unk = df[df.proc_family == "unknown"]
    if len(unk):
        log.warning("%d GSMs with unclassified processing; top strings:", len(unk))
        for s, n in unk.data_processing.value_counts().head(8).items():
            log.warning("   [%5d] %s", n, s[:140])

    try:
        have = con.execute(
            "SELECT DISTINCT upper(trim(gsm_id)) g FROM geo_expr").fetchdf().g
        s = set(df.gsm_id)
        log.info("geo_expr samples: %d | with metadata: %d | WITHOUT: %d",
                 len(have), int(have.isin(s).sum()), int((~have.isin(s)).sum()))
    except Exception as e:
        log.warning("geo_expr coverage check skipped: %s", e)
    return summary


# ======================================================================
# 2. Schema discovery
# ======================================================================

def describe_tables(con) -> pd.DataFrame:
    """
    Summarise every table in the database.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.

    Returns
    -------
    pandas.DataFrame
        Columns ``table``, ``rows``, ``n_cols`` and ``head`` — the first ten
        column names, enough to recognise a wide expression matrix without
        printing thousands of gene columns.
    """
    rows = []
    for t in con.execute("SHOW TABLES").fetchdf().iloc[:, 0]:
        n = con.execute(f'SELECT count(*) FROM "{t}"').fetchone()[0]
        cols = con.execute(f"PRAGMA table_info('{t}')").fetchdf()["name"].tolist()
        rows.append(dict(table=t, rows=n, n_cols=len(cols),
                         head=", ".join(cols[:10])))
    return pd.DataFrame(rows)


def ensg_columns(con, table: str) -> dict[str, str]:
    """{normalised ENSG -> actual column name} for a wide expression table.

    depmap_expr uses lowercase headers, geo_expr uppercase; (?i) catches both
    and .upper() puts them on one axis. Version suffixes are stripped.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Wide expression table to inspect.

    Returns
    -------
    dict
        Canonical uppercase ENSG ID to the column's actual name, so the
        original header can be quoted back in SQL.

    Notes
    -----
    Non-ENSG columns are skipped silently up to a count of three (the
    expected ID and metadata columns); beyond that a warning fires, since a
    large skip count usually means the regex needs widening rather than
    that the table genuinely has that much metadata.
    """
    cols = con.execute(f"PRAGMA table_info('{table}')").fetchdf()["name"].tolist()
    out, skipped = {}, []
    for c in cols:
        m = _ENSG.fullmatch(c.strip())
        if m:
            out[m.group(1).upper()] = c
        else:
            skipped.append(c)
    if len(skipped) > 3:
        log.warning("%s: %d non-ENSG columns skipped, e.g. %s -- widen _ENSG "
                    "if these are genes", table, len(skipped), skipped[:6])
    return out


# ======================================================================
# 3. Axes: method rows -> model_id -> lineage
# ======================================================================

def make_axis(name, source, table, id_col, raw_ids, model_ids, units) -> MethodAxis:
    """
    Construct a :class:`MethodAxis` and derive its replicate structure.

    Computes the unique model list, the inverse index from raw rows to
    models, and the raw-sample count per model — the three things
    :func:`collapse_replicates` and the GEO weighting need.

    Parameters
    ----------
    name, source, table, id_col, units : str
        Axis metadata; see :class:`MethodAxis`.
    raw_ids, model_ids : iterable
        Parallel sequences, one entry per raw sample row.

    Returns
    -------
    MethodAxis
        With ``models``, ``inv`` and ``n_samples`` populated. ``floor`` and
        ``scale_detected`` remain unset until :func:`profile_method` runs.
    """
    raw = np.asarray(list(raw_ids), dtype=object)
    mod = np.asarray(list(model_ids), dtype=object)
    models, inv = np.unique(mod, return_inverse=True)
    ax = MethodAxis(name=name, source=source, table=table, id_col=id_col,
                    units=units, raw_ids=raw, model_ids=mod)
    ax.models = models
    ax.inv = inv
    ax.n_samples = np.bincount(inv, minlength=len(models)).astype(np.float32)
    return ax


def build_method_axes(con, geo_platform_table: str = "geo_platform"
                      ) -> list[MethodAxis]:
    """
    Build one :class:`MethodAxis` per commensurable measurement scale.

    Three sources contribute: ``hpa_rna`` and ``depmap_expr`` give one axis
    each, while ``geo_expr`` is split into one axis per
    ``platform_id|proc_family`` combination, since samples processed
    differently are not on a common scale.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    geo_platform_table : str, optional
        Method assignment table. Defaults to ``"geo_platform"``.

    Returns
    -------
    list of MethodAxis
        Every usable axis, logged with its row and model counts.

    Raises
    ------
    RuntimeError
        If no axis could be built — nothing downstream can proceed.

    Notes
    -----
    Three degradations are warned about rather than raised, in keeping with
    flag-don't-drop:

    * a missing ``geo_platform`` table collapses all GEO samples into one
      method, which is wrong but visible;
    * GEO samples with no metadata record get their own
      ``UNKNOWN|unknown`` method rather than being discarded;
    * a method covering fewer than 20 cell lines is flagged, since its
      lineage strata will mostly fall back to the method-global scale.

    A source whose query fails is skipped with a warning, so the run
    continues on whichever sources are available.
    """
    axes: list[MethodAxis] = []

    # ---- hpa_rna (long) --------------------------------------------------
    try:
        m = con.execute("""
            SELECT DISTINCT model_id FROM hpa_rna WHERE model_id IS NOT NULL
        """).fetchdf()["model_id"].tolist()
        if m:
            axes.append(make_axis("hpa_rna_seq", "hpa_rna", "hpa_rna",
                                  "model_id", m, m, "tpm"))
    except Exception as e:
        log.warning("hpa_rna axis skipped: %s", e)

    # ---- depmap_expr (wide, profile_id -> model_id) -----------------------
    try:
        dm = con.execute("""
            SELECT e.profile_id, p.model_id
            FROM depmap_expr e
            JOIN depmap_profiles p USING (profile_id)
            WHERE p.model_id IS NOT NULL
        """).fetchdf()
        if len(dm):
            axes.append(make_axis("depmap_rna_seq", "depmap", "depmap_expr",
                                  "profile_id", dm.profile_id, dm.model_id,
                                  "log2tpm1"))
    except Exception as e:
        log.warning("depmap axis skipped: %s", e)

    # ---- geo_expr (wide); method = platform_id x proc_family --------------
    has_gp = geo_platform_table in set(con.execute("SHOW TABLES").fetchdf().iloc[:, 0])
    if not has_gp:
        log.warning("%s missing -- run --build-geo-platform first; all GEO "
                    "samples will collapse into one method", geo_platform_table)
    join_gp = (f'LEFT JOIN "{geo_platform_table}" g '
               'ON upper(trim(g.gsm_id)) = upper(trim(e.gsm_id))') if has_gp else ""
    sel_gp = ("coalesce(g.method,'UNKNOWN|unknown') AS method, "
              "coalesce(g.scale_hint,'unknown') AS scale_hint"
              ) if has_gp else ("'UNKNOWN|unknown' AS method, "
                                "'unknown' AS scale_hint")
    try:
        geo = con.execute(f"""
            SELECT e.gsm_id AS gsm_id, i.model_id AS model_id, {sel_gp}
            FROM geo_expr e
            JOIN geo_info i
              ON upper(trim(i.geo_accession)) = upper(trim(e.gsm_id))
            {join_gp}
            WHERE i.model_id IS NOT NULL
        """).fetchdf()
    except Exception as e:
        log.warning("geo axis skipped: %s", e)
        geo = pd.DataFrame()

    if len(geo):
        n_nometa = int((geo.method == "UNKNOWN|unknown").sum())
        if n_nometa:
            log.warning("%d GEO samples without a geo_meta record -> own method, "
                        "flagged not dropped", n_nometa)
        for meth, g in geo.groupby("method", sort=True):
            n_models = g.model_id.nunique()
            if n_models < 20:
                log.warning("method %s: only %d cell lines -- narrow panel; its "
                            "lineage strata will fall back toward the "
                            "method-global scale", meth, n_models)
            axes.append(make_axis(f"geo_{meth}", "geo", "geo_expr", "gsm_id",
                                  g.gsm_id, g.model_id, g.scale_hint.iloc[0]))

    for a in axes:
        log.info("  axis %-30s src=%-8s raw=%6d models=%5d units=%s",
                 a.name, a.source, len(a.raw_ids), len(a.models), a.units)
    if not axes:
        raise RuntimeError("no usable expression axes found")
    return axes


def load_lineage_map(con) -> pd.Series:
    """
    Load ``model_id`` -> lineage from ``sample_info``.

    Lineage is the reference population the z is computed against, not an
    output key: one cell line has exactly one lineage.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.

    Returns
    -------
    pandas.Series
        Lineage indexed by ``model_id``, NULLs filled as ``"unknown"`` and
        the count of those warned about — an ``unknown`` stratum is still a
        stratum and will be scored as one.
    """
    s = (con.execute("""
            SELECT model_id, lineage FROM sample_info WHERE model_id IS NOT NULL
         """).fetchdf()
         .drop_duplicates("model_id").set_index("model_id")["lineage"])
    miss = int(s.isna().sum())
    if miss:
        log.warning("%d models in sample_info have NULL lineage -> 'unknown'", miss)
    return s.fillna("unknown")


def normalise_ensg(value) -> str | None:
    """Return the canonical ENSG prefix, stripping version suffixes such as .7.

    Parameters
    ----------
    value : any
        Candidate identifier; None and empty strings are tolerated.

    Returns
    -------
    str or None
        Uppercase ENSG ID without version, or None if it does not match.
    """
    if value is None:
        return None
    s = str(value).strip().upper()
    if not s:
        return None
    m = _ENSG.fullmatch(s)
    return m.group(1).upper() if m else None


def filter_valid_ensg_ids(genes, valid_ids=None) -> list[str]:
    """Keep only canonical, roster-approved ENSG IDs in a consistent order.

    Normalises, optionally restricts to a roster, and de-duplicates while
    preserving input order — so a caller-supplied gene subset keeps the
    order it was given.

    Parameters
    ----------
    genes : iterable
        Candidate identifiers, in any case and with or without versions.
    valid_ids : iterable, optional
        Allowed IDs. When omitted or empty, no restriction is applied.

    Returns
    -------
    list of str
        Canonical IDs, de-duplicated, in first-seen order.

    Notes
    -----
    Silently drops anything unrecognised, so a ``--genes`` argument with a
    typo yields a shorter list rather than an error. Compare lengths if
    that matters.
    """
    valid = {str(v).strip().upper() for v in (valid_ids or []) if str(v).strip()}
    seen = set()
    out = []
    for raw in genes or []:
        norm = normalise_ensg(raw)
        if not norm:
            continue
        if valid and norm not in valid:
            continue
        if norm in seen:
            continue
        seen.add(norm)
        out.append(norm)
    return out


def gene_roster_ensg_ids(con) -> set[str]:
    """Gene roster IDs, normalised to the canonical ENSG form used in scoring.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.

    Returns
    -------
    set of str
        Canonical ENSG IDs. Empty if ``gene_roster`` is unavailable, which
        :func:`filter_valid_ensg_ids` treats as no restriction.
    """
    try:
        ids = con.execute("""
            SELECT DISTINCT upper(trim(gene_id)) AS gene_id
            FROM gene_roster
            WHERE gene_id IS NOT NULL
        """).fetchdf().gene_id.dropna().tolist()
    except Exception:
        return set()
    return {g for g in map(normalise_ensg, ids) if g is not None}


def gene_universe(con, axes: list[MethodAxis],
                  colmaps: dict[str, dict]) -> list[str]:
    """Union of gene axes restricted to canonical ENSG IDs from the gene roster.

    A gene measured by any one method is scoreable, so the union is taken
    rather than the intersection — coverage is recorded per row as
    ``k_src`` instead of being enforced up front.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    axes : list of MethodAxis
        Axes to draw genes from.
    colmaps : dict
        Axis name to its ENSG column map; unused for ``hpa_rna``, which is
        long-format and queried directly.

    Returns
    -------
    list of str
        Sorted canonical ENSG IDs. Empty if no axis contributed any.
    """
    valid = gene_roster_ensg_ids(con)
    sets = []
    for ax in axes:
        if ax.source == "hpa_rna":
            cand = con.execute("""
                SELECT DISTINCT regexp_extract(upper(gene_id),'(ENSG[0-9]+)',1) e
                FROM hpa_rna
            """).fetchdf().e.dropna().tolist()
        else:
            cand = list((colmaps.get(ax.name) or {}).keys())
        norm = filter_valid_ensg_ids(cand, valid)
        if norm:
            sets.append(set(norm))
    return sorted(set().union(*sets)) if sets else []


# ======================================================================
# 4. Chunked fetch
# ======================================================================

def fetch_wide_chunk(con, table: str, id_col: str, genes: list[str],
                     colmap: dict[str, str]) -> pd.DataFrame:
    """One read of a wide table for a gene chunk. Indexed by upper(id).

    Selects only the requested gene columns, aliased to their canonical
    ENSG names, so a thousand-column matrix is never materialised whole.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table, id_col : str
        Wide table and its sample ID column.
    genes : list of str
        Canonical ENSG IDs for this chunk.
    colmap : dict
        Canonical ID to actual column name, from :func:`ensg_columns`.

    Returns
    -------
    pandas.DataFrame
        Indexed by uppercased, trimmed sample ID, columns named by
        canonical ENSG. Empty if no requested gene exists in this table.
    """
    present = [(g, colmap[g]) for g in genes if g in colmap]
    if not present:
        return pd.DataFrame()
    sel = ", ".join(f'"{c}"::DOUBLE AS "{g}"' for g, c in present)
    df = con.execute(
        f'SELECT upper(trim("{id_col}")) AS _id, {sel} FROM "{table}"').fetchdf()
    return df.set_index("_id")


def fetch_hpa_chunk(con, genes: list[str]) -> pd.DataFrame:
    """Long -> wide for a gene chunk. median() collapses transcript duplicates.

    HPA is stored long, with several transcript rows per gene; the median
    is taken rather than the sum or max, so one high-expressing transcript
    cannot stand in for the gene.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    genes : list of str
        Canonical ENSG IDs for this chunk, bound as query parameters.

    Returns
    -------
    pandas.DataFrame
        Indexed by uppercased ``model_id``, one column per gene. Empty if
        nothing matched.
    """
    if not genes:
        return pd.DataFrame()
    q = """
        SELECT upper(trim(model_id)) AS _id,
               regexp_extract(upper(gene_id),'(ENSG[0-9]+)',1) AS ensg,
               median(tpm) AS v
        FROM hpa_rna
        WHERE model_id IS NOT NULL
          AND regexp_extract(upper(gene_id),'(ENSG[0-9]+)',1) IN ({})
        GROUP BY 1,2
    """.format(",".join("?" * len(genes)))
    df = con.execute(q, genes).fetchdf()
    if df.empty:
        return pd.DataFrame()
    return df.pivot(index="_id", columns="ensg", values="v")


def normalize_ids(ids) -> np.ndarray:
    """Fast uppercase/strip conversion for large arrays without pandas apply.

    Parameters
    ----------
    ids : iterable
        Identifiers to normalise.

    Returns
    -------
    numpy.ndarray
        Object array of uppercased, stripped strings, same length as input.
    """
    arr = np.asarray(list(ids), dtype=object)
    out = np.empty(arr.shape, dtype=object)
    for i, v in enumerate(arr):
        s = str(v).strip().upper()
        out[i] = s
    return out


def align_axis_rows(ax: MethodAxis, wide: pd.DataFrame,
                    genes: list[str]) -> np.ndarray:
    """Slice a shared table chunk down to this axis's raw rows.

    Several GEO axes read the same ``geo_expr`` chunk, so each one
    reindexes the shared frame onto its own sample list rather than
    re-querying.

    Parameters
    ----------
    ax : MethodAxis
        Axis whose raw rows are wanted.
    wide : pandas.DataFrame or None
        Shared chunk indexed by normalised sample ID.
    genes : list of str
        Gene column order for the output.

    Returns
    -------
    numpy.ndarray
        ``(n_raw_ids, n_genes)`` float32, NaN where a sample or gene is
        absent — missing stays missing rather than becoming zero.
    """
    X = np.full((len(ax.raw_ids), len(genes)), np.nan, np.float32)
    if wide is None or wide.empty:
        return X
    keys = normalize_ids(ax.raw_ids)
    sub = wide.reindex(keys)
    for j, g in enumerate(genes):
        if g in sub.columns:
            arr = sub[g].to_numpy(dtype=np.float32, copy=False)
            X[:, j] = arr
    return X


def fetch_chunk_by_table(con, axes: list[MethodAxis], genes: list[str],
                         colmaps: dict[str, dict]) -> dict[str, pd.DataFrame]:
    """Read each physical table once per gene chunk, not once per method.

    All GEO methods share geo_expr, so without this the table is scanned once
    per platform x processing combination.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    axes : list of MethodAxis
        Axes needing data this chunk.
    genes : list of str
        Canonical ENSG IDs for this chunk.
    colmaps : dict
        Axis name to its ENSG column map.

    Returns
    -------
    dict
        Table name to its chunk frame, to be sliced per axis by
        :func:`align_axis_rows`.
    """
    cache = {}
    for ax in axes:
        if ax.table in cache:
            continue
        if ax.source == "hpa_rna":
            cache[ax.table] = fetch_hpa_chunk(con, genes)
        else:
            cache[ax.table] = fetch_wide_chunk(con, ax.table, ax.id_col, genes,
                                               colmaps.get(ax.name) or {})
    return cache


def collapse_replicates(ax: MethodAxis, X: np.ndarray) -> np.ndarray:
    """Median over replicate raw rows -> (n_models, n_genes).

    Not optional: some lines carry 9-14 GSMs, and without this they would
    dominate their own lineage's median.

    Parameters
    ----------
    ax : MethodAxis
        Axis supplying the raw-row-to-model index.
    X : numpy.ndarray
        ``(n_raw_ids, n_genes)`` values.

    Returns
    -------
    numpy.ndarray
        ``(n_models, n_genes)`` float32. Returned unchanged when the axis
        has no replicates.

    Notes
    -----
    Rows are sorted once by model index and the per-model spans located
    with ``searchsorted``, so the collapse is a single pass rather than a
    groupby per gene. Single-row models bypass ``nanmedian`` entirely.
    """
    if len(ax.models) == len(ax.raw_ids):
        return X
    out = np.full((len(ax.models), X.shape[1]), np.nan, np.float32)
    order = np.argsort(ax.inv, kind="stable")
    inv_sorted = ax.inv[order]
    starts = np.searchsorted(inv_sorted, np.arange(len(ax.models)), side="left")
    ends = np.searchsorted(inv_sorted, np.arange(len(ax.models)), side="right")
    for k, (s, e) in enumerate(zip(starts, ends)):
        if e <= s:
            continue
        blk = X[order[s:e]]
        out[k] = blk[0] if e - s == 1 else np.nanmedian(blk, axis=0)
    return out


# ======================================================================
# 5. Scale detection and detection floors
# ======================================================================

def detect_scale_from_values(v: np.ndarray) -> str:
    """Linear vs log2 from the data alone. GEO metadata is not trusted here.

    Three signatures decide it: negative values can only be log; a very high
    99th percentile, absolutely or relative to the median, indicates the
    long right tail of linear expression; otherwise a low ceiling implies
    log.

    Parameters
    ----------
    v : numpy.ndarray
        Finite-filtered value sample.

    Returns
    -------
    str
        ``"log"``, ``"linear"``, or ``"unknown"`` when fewer than 50 finite
        values were available to judge on.
    """
    v = v[np.isfinite(v)]
    if v.size < 50:
        return "unknown"
    lo, hi, med = float(v.min()), float(np.percentile(v, 99)), float(np.median(v))
    if lo < -0.5:
        return "log"
    if hi > 100 or (med > 0 and hi / max(med, EPS) > 50):
        return "linear"
    return "log" if hi < 25 else "linear"


def profile_method(con, ax: MethodAxis, colmap, probe_genes: list[str],
                   floor_q: float) -> None:
    """Set ax.scale_detected and ax.floor in place. Warns on SOFT/data conflict.

    Scale is taken from the declared units where those are authoritative
    (the two RNA-seq sources) and from the data otherwise. Where SOFT
    metadata and the data disagree, the data wins and the conflict is
    logged — a mislabelled processing string is more likely than a
    misbehaving distribution.

    Floors differ by source for the same reason: RNA-seq uses HPA's own
    ``>=1 TPM`` criterion, a published detection threshold, while arrays
    have no such criterion and get an empirical quantile of their own
    values.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    ax : MethodAxis
        Axis to profile. Mutated in place.
    colmap : dict or None
        ENSG column map; None for ``hpa_rna``.
    probe_genes : list of str
        Random gene sample, large enough to characterise the distribution
        without reading the full matrix.
    floor_q : float
        Quantile for the array detection floor.

    Returns
    -------
    None
        ``ax.scale_detected`` and ``ax.floor`` are set as side effects.

    Notes
    -----
    An axis with fewer than 50 finite probe values gets scale
    ``"unknown"`` and a floor of ``-inf``, which makes every observed value
    count as detected rather than silently censoring the method.
    """
    if ax.source == "hpa_rna":
        wide = fetch_hpa_chunk(con, probe_genes)
    else:
        wide = fetch_wide_chunk(con, ax.table, ax.id_col, probe_genes, colmap or {})
    X = collapse_replicates(ax, align_axis_rows(ax, wide, probe_genes))
    v = X[np.isfinite(X)]
    if v.size < 50:
        ax.scale_detected, ax.floor = "unknown", -np.inf
        log.warning("%s: too few finite values to profile", ax.name)
        return

    empirical = detect_scale_from_values(v)
    if ax.units == "log2tpm1":
        scale = "log"
    elif ax.units == "tpm":
        scale = "linear"
    elif ax.units in ("log", "linear"):
        scale = empirical
        if empirical != ax.units:
            log.warning("%s: SOFT metadata says %s, data looks %s "
                        "(min=%.2f p99=%.1f) -> using data",
                        ax.name, ax.units, empirical, v.min(),
                        np.percentile(v, 99))
    else:
        scale = empirical

    ax.scale_detected = scale
    vv = np.log2(np.clip(v, 0, None) + 1.0) if scale == "linear" else v
    ax.floor = (float(np.log2(TPM_FLOOR + 1.0))
                if ax.source in ("hpa_rna", "depmap")
                else float(np.quantile(vv, floor_q)))
    log.info("  %-30s units=%-9s scale=%-7s floor=%.3f (log2)",
             ax.name, ax.units, scale, ax.floor)


def to_log2(ax: MethodAxis, X: np.ndarray) -> np.ndarray:
    """
    Put an axis's values on the common log2 scale.

    Linear values become ``log2(x + 1)``, with negatives clipped to zero
    first; values already on a log scale pass through. This is what makes
    RNA-seq and array methods comparable before z-scoring.

    Parameters
    ----------
    ax : MethodAxis
        Axis, already profiled so ``scale_detected`` is set.
    X : numpy.ndarray
        Values in the axis's native scale.

    Returns
    -------
    numpy.ndarray
        Values in log2 space.

    Notes
    -----
    An axis whose scale came out ``"unknown"`` is left untransformed.
    """
    if ax.scale_detected == "linear":
        return np.log2(np.clip(X, 0, None) + 1.0)
    return X


# ======================================================================
# 6. Z-score core
# ======================================================================

def nan_mad(X: np.ndarray, axis: int = 0):
    """(median, 1.4826*MAD) ignoring NaN. Rousseeuw & Croux (1993).

    All-NaN slices are legitimate -- a gene absent from a method's platform,
    or a lineage where nothing was measured. They return NaN, which the
    branching below treats as no_data. The warning is suppressed; 

    Parameters
    ----------
    X : numpy.ndarray
        Values, typically ``(n_models, n_genes)``.
    axis : int, optional
        Axis to reduce over. Default 0, i.e. per gene.

    Returns
    -------
    med : numpy.ndarray
        Median ignoring NaN.
    mad : numpy.ndarray
        MAD scaled by 1.4826, making it a consistent estimator of sigma
        under normality.
    """
    with np.errstate(invalid="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        med = np.nanmedian(X, axis=axis)
        mad = 1.4826 * np.nanmedian(np.abs(X - med), axis=axis)
    return med, mad


def method_z(V: np.ndarray, lineage: np.ndarray, floor: float, p: Params,
             stratify: bool):
    """(n_models, n_genes) log2 values -> (Z, status, det_frac, scale_src).

    Three regimes per (lineage, gene), decided by detection fraction:
      det >= upper_cut : robust z, lineage scale shrunk toward method-global
      lower<=det<upper : rank inverse-normal; floor ties share a midrank, so
                         censored lines get a lower bound not a point estimate
      det <  lower_cut : NaN + below_detection (must not vote as z=0)

    Parameters
    ----------
    V : numpy.ndarray
        ``(n_models, n_genes)`` values already on the log2 scale.
    lineage : numpy.ndarray
        Lineage label per model row.
    floor : float
        The axis's detection floor in log2 space.
    p : Params
        Cutoffs and the shrinkage constant.
    stratify : bool
        Score within lineage when True, over all models when False. The
        caller runs both to produce ``z_lineage`` and ``z_global``.

    Returns
    -------
    Z : numpy.ndarray
        float32 z-scores, NaN where not scoreable.
    S : numpy.ndarray
        int8 status codes, keyed by :data:`STATUS`.
    D : numpy.ndarray
        float32 detection fraction per cell.
    C : numpy.ndarray
        int8 scale-source codes, keyed by :data:`SCALE_SRC`.

    Notes
    -----
    Two guards keep the parametric branch honest on thin data. The scale is
    floored at a tenth of the method-global MAD, so a lineage with almost
    no spread cannot manufacture huge z-scores; and the lineage median is
    only used as the centre when the stratum has at least five members,
    falling back to the global median otherwise. Both conditions are
    recorded in ``C`` rather than hidden.
    """
    n, G = V.shape
    Z = np.full((n, G), np.nan, np.float32)
    S = np.full((n, G), 3, np.int8)
    D = np.full((n, G), np.nan, np.float32)
    C = np.full((n, G), 0, np.int8)

    med_g, mad_g = nan_mad(V)
    mad_g = np.where(np.isfinite(mad_g) & (mad_g > EPS), mad_g, np.nan)
    gl = np.where(np.isfinite(mad_g), mad_g, 0.0)

    groups = ([(l, np.flatnonzero(lineage == l)) for l in np.unique(lineage)]
              if stratify else [("__all__", np.arange(n))])

    for _, idx in groups:
        if idx.size == 0:
            continue
        Xb = V[idx]
        nb = idx.size
        fin = np.isfinite(Xb)
        nfin = fin.sum(axis=0)
        with np.errstate(invalid="ignore"):
            det = np.where(nfin > 0,
                           np.nansum(Xb >= floor, axis=0) / np.maximum(nfin, 1),
                           np.nan)
        D[idx] = det

        # branch A: parametric robust z, scale shrunk toward method-global
        med_L, mad_L = nan_mad(Xb)
        w = nb / (nb + p.n0)
        base = np.where(np.isfinite(mad_L), mad_L, 0.0)
        scale = np.sqrt(w * base ** 2 + (1.0 - w) * gl ** 2)
        scale = np.maximum(scale, np.maximum(0.1 * gl, EPS))    # MAD floor
        use_lineage_centre = np.isfinite(med_L) & (nb >= 5)
        ctr = np.where(use_lineage_centre, med_L, med_g)
        Za = (Xb - ctr) / scale

        # branch B: rank inverse-normal (Blom / Van der Waerden form)
        with np.errstate(invalid="ignore"):
            R = rankdata(Xb, axis=0, nan_policy="omit")
        Zb = norm.ppf((R - 0.5) / np.maximum(nfin, 1))
        Zb = np.where(np.isfinite(Zb), Zb, np.nan)

        okA = det >= p.upper_cut
        okB = (det >= p.lower_cut) & ~okA
        dead = det < p.lower_cut

        blk = np.where(okA[None, :], Za, np.where(okB[None, :], Zb, np.nan))
        st = np.where(okA, 0, np.where(okB, 1, np.where(dead, 2, 3))).astype(np.int8)
        if nb < 5:
            src = np.full(G, 3, np.int8)
        elif w >= 0.9:
            src = np.where(use_lineage_centre, 1, 3).astype(np.int8)
        else:
            src = np.where(use_lineage_centre, 2, 3).astype(np.int8)

        Z[idx] = np.where(fin, blk, np.nan).astype(np.float32)
        S[idx] = np.where(fin, st[None, :], 3)
        C[idx] = np.where(fin, src[None, :], 0)

    return Z, S, D, C


def stouffer(Z: np.ndarray, W: np.ndarray | None = None):
    """Z (n, G, k) -> (combined (n,G), k_used (n,G)).

    Denominator is sqrt(sum w^2), not sum(w) -- that is what keeps the output
    unit-variance. Missing entries reduce k; they never vote as z=0.

    Parameters
    ----------
    Z : numpy.ndarray
        Z-scores stacked on a third axis, one slice per contributor.
    W : numpy.ndarray or None, optional
        Matching weights. None means equal weight — one contributor, one
        vote.

    Returns
    -------
    combined : numpy.ndarray
        float32 combined z, NaN where no contributor supplied a value.
    k_used : numpy.ndarray
        int16 count of contributors that actually voted per cell.
    """
    m = np.isfinite(Z)
    Zf = np.where(m, Z, 0.0)
    Wf = np.where(m, 1.0 if W is None else W, 0.0)
    num = (Wf * Zf).sum(axis=2)
    den = np.sqrt((Wf ** 2).sum(axis=2))
    out = np.where(den > EPS, num / np.maximum(den, EPS), np.nan)
    return out.astype(np.float32), m.sum(axis=2).astype(np.int16)


def process_chunk(con, axes, colmaps, genes, models, midx, lineage_map_s,
                  p: Params, stratify: bool = True) -> dict:
    """Per-source stacks of method-level z, aligned to the global model axis.

    Each axis is fetched, replicate-collapsed, put on log2, scored, then
    scattered onto the full model axis so every method's output shares one
    row order and the stacks can be combined directly.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    axes : list of MethodAxis
        Axes to score, already profiled.
    colmaps : dict
        Axis name to ENSG column map.
    genes : list of str
        Canonical ENSG IDs for this chunk.
    models : list of str
        Global model axis.
    midx : dict
        ``model_id`` to its row index in ``models``.
    lineage_map_s : pandas.Series
        Lineage indexed by ``model_id``.
    p : Params
        Scoring parameters.
    stratify : bool, optional
        Passed through to :func:`method_z`. Default True.

    Returns
    -------
    dict
        Source name to a list of per-method dicts, each with ``name``,
        ``z``, ``status``, ``det``, ``scale_src`` and ``n``, all on the
        global ``(n_models, n_genes)`` shape.

    Notes
    -----
    Models on an axis but not in ``midx`` are dropped by the ``rows >= 0``
    mask; models in ``midx`` but absent from an axis keep their NaN and
    status 3, so the two directions of mismatch are handled distinctly.
    """
    nM, G = len(models), len(genes)
    cache = fetch_chunk_by_table(con, axes, genes, colmaps)
    per_src = {}
    for ax in axes:
        X = collapse_replicates(ax, align_axis_rows(ax, cache.get(ax.table), genes))
        V = to_log2(ax, X)
        lin = lineage_map_s.reindex(ax.models).fillna("unknown").to_numpy()
        Z, S, D, C = method_z(V, lin, ax.floor, p, stratify)

        rows = np.array([midx.get(m, -1) for m in ax.models])
        keep = rows >= 0
        Zg = np.full((nM, G), np.nan, np.float32)
        Sg = np.full((nM, G), 3, np.int8)
        Dg = np.full((nM, G), np.nan, np.float32)
        Cg = np.zeros((nM, G), np.int8)
        Ng = np.zeros((nM, G), np.float32)
        Zg[rows[keep]] = Z[keep]
        Sg[rows[keep]] = S[keep]
        Dg[rows[keep]] = D[keep]
        Cg[rows[keep]] = C[keep]
        Ng[rows[keep]] = ax.n_samples[keep][:, None]
        per_src.setdefault(ax.source, []).append(
            dict(name=ax.name, z=Zg, status=Sg, det=Dg, scale_src=Cg, n=Ng))
    return per_src


def combine_sources(per_src: dict):
    """Stage 1: GEO methods -> z_geo (sqrt-n weighted).
       Stage 2: the 3 sources -> z_final, equal weight (one source, one vote).

    Deliberately less powerful than flat Stouffer over all methods if they were
    equally trustworthy, and more robust to any single source's artefacts.
    With GEO split on platform x processing there may be 10+ GEO methods; flat
    Stouffer would hand GEO ~80% of the weight.

    Parameters
    ----------
    per_src : dict
        Source name to per-method dicts, from :func:`process_chunk`.

    Returns
    -------
    z_final : numpy.ndarray
        ``(n_models, n_genes)`` combined z across sources.
    k_src : numpy.ndarray
        Sources that voted per cell — 1 to 3.
    src_z : dict
        Source name to its stage-1 combined z, written out as the
        ``z_hpa_rna`` / ``z_depmap`` / ``z_geo`` columns.
    extra : dict
        Per-source diagnostics: ``__k``, ``__status``, ``__scale_src``
        stacks and the mean ``__det``.
    order : list of str
        Sources present, in fixed ``hpa_rna``, ``depmap``, ``geo`` order —
        the column order the output frame relies on.

    Notes
    -----
    Only GEO is weighted within a source, by ``sqrt(n_samples)``, since its
    methods differ in how many replicates back each cell line. The RNA-seq
    sources have one method each, so weighting would be a no-op.
    """
    src_z, extra = {}, {}
    for src, items in per_src.items():
        Z = np.stack([d["z"] for d in items], axis=2)
        W = np.stack([np.sqrt(d["n"]) for d in items], axis=2)
        z, k = stouffer(Z, W if src == "geo" else None)
        src_z[src] = z
        extra[f"{src}__k"] = k
        extra[f"{src}__status"] = np.stack([d["status"] for d in items], axis=2)
        extra[f"{src}__scale_src"] = np.stack([d["scale_src"] for d in items], axis=2)
        with np.errstate(invalid="ignore"):
            extra[f"{src}__det"] = np.nanmean(
                np.stack([d["det"] for d in items], axis=2), axis=2)

    order = [s for s in ("hpa_rna", "depmap", "geo") if s in src_z]
    Zs = np.stack([src_z[s] for s in order], axis=2)
    z_final, k_src = stouffer(Zs)
    return z_final, k_src, src_z, extra, order


# ======================================================================
# 7. Sensitivity analysis
# ======================================================================

def grid_combinations(grid: dict = None) -> list[tuple]:
    """
    Expand the parameter grid into valid combinations.

    Parameters
    ----------
    grid : dict, optional
        Values per parameter. Defaults to :data:`PARAM_GRID`.

    Returns
    -------
    list of tuple
        ``(floor_q, n0, upper_cut, lower_cut)`` for every combination where
        ``lower_cut < upper_cut`` — the others describe no valid censoring
        band.
    """
    grid = grid or PARAM_GRID
    keys = ["floor_q", "n0", "upper_cut", "lower_cut"]
    return [c for c in itertools.product(*[grid[k] for k in keys]) if c[3] < c[2]]


def sample_calibration_genes(all_genes, n=N_CALIB_GENES, rng=None) -> list[str]:
    """
    Draw a random gene sample for calibration.

    Parameters
    ----------
    all_genes : sequence
        Gene universe to sample from.
    n : int, optional
        Sample size, capped at the universe size. Defaults to
        :data:`N_CALIB_GENES`.
    rng : numpy.random.Generator, optional
        Seeded generator, so calibration is reproducible. Defaults to seed
        0.

    Returns
    -------
    list of str
        Genes sampled without replacement.
    """
    rng = rng or np.random.default_rng(0)
    n = min(n, len(all_genes))
    return list(rng.choice(np.asarray(all_genes, dtype=object), size=n, replace=False))


def _score_params(con, axes, colmaps, genes, models, midx, lineage_s,
                  p: Params, holdouts, ref_cache) -> float | None:
    """Held-out cross-source rank concordance for one parameter combination.

    For each held-out RNA-seq source, the remaining axes are scored under
    ``p`` and correlated against that source's reference scores. The
    parameters never see the source they are judged against.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    axes : list of MethodAxis
        All axes; the held-out source is excluded per fold.
    colmaps, genes, models, midx, lineage_s
        Universe and axis plumbing, as elsewhere.
    p : Params
        Combination under test.
    holdouts : list of str
        Source names to hold out in turn.
    ref_cache : dict
        Source name to its reference z, scored under the fixed
        :data:`REF_PARAMS_RAW`.

    Returns
    -------
    float or None
        Mean Spearman rho across folds, or None if no fold had enough
        comparable cells.

    Notes
    -----
    Scoring is restricted to the ambiguous detection band (0.01 to 0.90).
    Outside it every parameter setting agrees, so including those cells
    would dilute the signal the search is trying to resolve. Folds with
    fewer than :data:`MIN_CALIB_CELLS` comparable cells are skipped rather
    than contributing a noisy correlation.
    """
    rhos = []
    for h in holdouts:
        fit = [a for a in axes if a.source != h]
        if not fit:
            continue
        pf = process_chunk(con, fit, colmaps, genes, models, midx, lineage_s, p)
        zp, _, _, ex, order = combine_sources(pf)
        zr = ref_cache[h]
        det = None
        for s in order:
            d = ex.get(f"{s}__det")
            det = d if det is None else np.fmax(det, d)
        if det is None:
            det = np.full_like(zp, 0.5)
        # score only in the ambiguous band -- elsewhere no setting disagrees
        m = np.isfinite(zp) & np.isfinite(zr) & (det > 0.01) & (det < 0.90)
        if m.sum() < MIN_CALIB_CELLS:
            continue
        rhos.append(float(spearmanr(zp[m], zr[m]).statistic))
    return float(np.mean(rhos)) if rhos else None


def calibrate_params(con, axes, colmaps, all_genes, models, midx, lineage_s,
                     rng=None, grid: dict = None):
    """Grid search + 1-SE rule. Returns (Params, grid_frame, genes_used).

    Each combination is scored by held-out concordance against an RNA-seq
    source that did not contribute to the fit. The reference is itself
    scored under fixed parameters (:data:`REF_PARAMS_RAW`), so the target
    does not move with the setting being tested — otherwise the search
    would be circular.

    Selection is by the 1-SE rule rather than the argmax: every combination
    within one standard error of the best is treated as tied, and the tie
    is broken toward the conservative end — higher cuts, more rank-based
    and less parametric scoring.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    axes : list of MethodAxis
        All axes. Mutated by profiling as the search proceeds.
    colmaps : dict
        Axis name to ENSG column map.
    all_genes : sequence
        Gene universe to sample calibration genes from.
    models, midx, lineage_s
        Model axis plumbing.
    rng : numpy.random.Generator, optional
        Seeded generator. Defaults to seed 0.
    grid : dict, optional
        Parameter grid. Defaults to :data:`PARAM_GRID`.

    Returns
    -------
    p : Params
        Chosen parameters.
    gdf : pandas.DataFrame
        Every scorable combination with its rho, descending, with the
        standard error in ``.attrs["se"]``. Empty on total failure.
    genes : list of str
        The calibration gene sample, returned for reproducibility.

    Raises
    ------
    RuntimeError
        If no RNA-seq source is available to hold out — with only array
        methods there is no trustworthy reference to calibrate against.

    Notes
    -----
    If no combination scores, the reference parameters are returned with a
    warning rather than raising, so a run can still proceed on documented
    defaults. Check the log: a run on reference parameters is not a
    calibrated run and should not be described as one.

    The argmax of a noisy surface is largely an artefact of the noise,
    which is why the 1-SE rule is applied. The flatness of the ``n0`` sweep
    (see :func:`plot_n0_sweep`) is the visual check on whether the data
    identify that parameter at all.
    """
    rng = rng or np.random.default_rng(0)
    genes = sample_calibration_genes(all_genes, rng=rng)
    holdouts = [s for s in ("depmap", "hpa_rna")
                if any(a.source == s for a in axes)]
    if not holdouts:
        raise RuntimeError("need at least one RNA-seq source to hold out")
    log.info("calibrating on %d genes; holdouts=%s", len(genes), holdouts)

    ref_p = Params(**REF_PARAMS_RAW)
    ref_cache = {}
    for h in holdouts:
        hax = [a for a in axes if a.source == h]
        for a in hax:
            profile_method(con, a, colmaps.get(a.name), genes[:400], ref_p.floor_q)
        ph = process_chunk(con, hax, colmaps, genes, models, midx, lineage_s, ref_p)
        ref_cache[h], _, _, _, _ = combine_sources(ph)

    combos = grid_combinations(grid)
    log.info("grid: %d combinations", len(combos))
    rows = []
    for i, (fq, n0, uc, lc) in enumerate(combos, 1):
        p = Params(fq, int(n0), uc, lc)
        for a in axes:
            profile_method(con, a, colmaps.get(a.name), genes[:400], fq)
        rho = _score_params(con, axes, colmaps, genes, models, midx, lineage_s,
                            p, holdouts, ref_cache)
        if rho is not None:
            rows.append(dict(**asdict(p), rho=rho))
        if i % 20 == 0:
            log.info("  grid %d/%d", i, len(combos))

    if not rows:
        log.warning("grid search produced no scorable combination -> reference "
                    "parameters retained")
        return ref_p, pd.DataFrame(), genes

    gdf = pd.DataFrame(rows).sort_values("rho", ascending=False).reset_index(drop=True)
    best = gdf.iloc[0]
    se = float(gdf.rho.std(ddof=1) / np.sqrt(len(gdf))) if len(gdf) > 1 else 0.0
    within = gdf[gdf.rho >= best.rho - se]
    # the argmax of a noisy surface is an artefact; tie-break toward the
    # conservative end (more rank-based, less parametric)
    chosen = within.sort_values(["upper_cut", "lower_cut", "n0"],
                                ascending=[False, False, False]).iloc[0]
    log.info("best rho=%.4f  1SE=%.4f  %d/%d combos within 1 SE",
             best.rho, se, len(within), len(gdf))
    p = Params(float(chosen.floor_q), int(chosen.n0),
               float(chosen.upper_cut), float(chosen.lower_cut))
    log.info("CHOSEN (1-SE rule): %s", asdict(p))
    gdf.attrs["se"] = se
    return p, gdf, genes


# ======================================================================
# 8. Assembling and writing output
# ======================================================================

def chunk_to_frame(models, lineage_vec, genes, z_lin, z_glb, k_src, extra,
                   order, src_z: dict, min_k: int = 1) -> pd.DataFrame:
    """Long output frame for one gene chunk: ensg | model_id | lineage | ...

    Gene-major traversal (all cell lines for gene 1, then gene 2, ...) so
    DuckDB can skip row groups when filtering on ensg.

    The source-level z columns are built from the SAME index as the rest of
    the frame. Splitting this across two functions once let the two traversal
    orders diverge, which misaligns z_depmap/z_geo silently.

    Parameters
    ----------
    models, lineage_vec : sequence
        Global model axis and its lineage labels, same order.
    genes : sequence
        Gene axis for this chunk.
    z_lin, z_glb : numpy.ndarray
        Lineage-stratified and global combined z.
    k_src : numpy.ndarray
        Sources voting per cell.
    extra : dict
        Per-source diagnostic stacks from :func:`combine_sources`.
    order : list of str
        Source order; also fixes the ``z_<source>`` column order.
    src_z : dict
        Per-source combined z.
    min_k : int, optional
        Drop rows backed by fewer than this many sources. Default 1.

    Returns
    -------
    pandas.DataFrame
        One row per surviving (gene, cell line), with ``status`` and
        ``scale_source`` mapped to their labels. Empty if nothing survives.

    Notes
    -----
    ``status`` and ``scale_source`` are collapsed across methods by
    priority, in opposite directions. Status takes the *best* available
    (any ``ok`` makes the cell ok), since one method scoring properly is
    enough. Scale source takes the *worst* (any ``fallback_global`` marks
    the cell as fallback), so the reported provenance is never rosier than
    the weakest contributor.
    """
    st_all = np.concatenate([extra[f"{s}__status"] for s in order], axis=2)
    st = np.where((st_all == 0).any(axis=2), 0,
         np.where((st_all == 1).any(axis=2), 1,
         np.where((st_all == 2).any(axis=2), 2, 3))).astype(np.int8)

    sc_all = np.concatenate([extra[f"{s}__scale_src"] for s in order], axis=2)
    sc = np.where((sc_all == 3).any(axis=2), 3,
         np.where((sc_all == 2).any(axis=2), 2,
         np.where((sc_all == 1).any(axis=2), 1, 0))).astype(np.int8)

    with np.errstate(invalid="ignore"):
        det = np.nanmean(np.stack([extra[f"{s}__det"] for s in order], axis=2),
                         axis=2)
    k_meth = np.zeros_like(k_src)
    for s in order:
        k_meth = k_meth + extra[f"{s}__k"]

    keep = (k_src >= min_k) & np.isfinite(z_lin)
    c, r = np.nonzero(keep.T)          # gene-major; c indexes genes, r models
    if r.size == 0:
        return pd.DataFrame()

    df = pd.DataFrame({
        "ensg":         np.asarray(genes, dtype=object)[c],
        "model_id":     np.asarray(models, dtype=object)[r],
        "lineage":      np.asarray(lineage_vec, dtype=object)[r],
        "z_lineage":    z_lin[r, c],
        "z_global":     z_glb[r, c],
        "k_src":        k_src[r, c].astype(np.int8),
        "k_methods":    k_meth[r, c].astype(np.int16),
        "det_frac":     det[r, c],
        "status":       pd.Series(st[r, c]).map(STATUS).to_numpy(),
        "scale_source": pd.Series(sc[r, c]).map(SCALE_SRC).to_numpy(),
    })
    for s in order:
        df[f"z_{s}"] = src_z[s][r, c]
    return df


def write_chunk(con, df: pd.DataFrame, table: str, first: bool) -> None:
    """Append a gene chunk to the DuckDB output table.

    CREATE OR REPLACE on the first chunk gives idempotency across reruns;
    later chunks INSERT. DuckDB is the only sink -- at ~30M rows CSV would be
    ~4.5 GB and unusable.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open, writable connection.
    df : pandas.DataFrame
        Chunk frame from :func:`chunk_to_frame`. Empty frames are a no-op.
    table : str
        Output table name.
    first : bool
        True for the first written chunk, which creates or replaces the
        table.

    Returns
    -------
    None

    Notes
    -----
    The temporary registration is dropped in a ``finally`` block, so a
    failed insert does not leave a stale view bound to the connection.
    """
    if df.empty:
        return
    con.register("_chunk_tmp", df)
    try:
        if first:
            con.execute(
                f'CREATE OR REPLACE TABLE "{table}" AS SELECT * FROM _chunk_tmp')
        else:
            con.execute(f'INSERT INTO "{table}" SELECT * FROM _chunk_tmp')
    finally:
        con.unregister("_chunk_tmp")


def summarise_chunk(df: pd.DataFrame) -> pd.DataFrame:
    """
    Reduce one chunk to diagnostic counts.

    Counting per chunk and summing afterwards keeps the diagnostics small
    regardless of how many rows were written.

    Parameters
    ----------
    df : pandas.DataFrame
        Chunk frame from :func:`chunk_to_frame`.

    Returns
    -------
    pandas.DataFrame
        Row counts by ``lineage``, ``status``, ``scale_source`` and
        ``k_src``. Empty for an empty chunk.
    """
    if df.empty:
        return pd.DataFrame()
    return (df.groupby(["lineage", "status", "scale_source", "k_src"],
                       observed=True).size().rename("n").reset_index())


def verify_output(con, table: str) -> pd.DataFrame:
    """Post-run sanity: row/gene/line counts and the k_src breakdown.

    The mean and SD of ``z_lineage`` are the check that matters: a mean
    near 0 and an SD near 1 indicate the combination preserved the z-scale,
    while a markedly smaller SD suggests the Stouffer denominator or the
    shrinkage is over-damping.

    Parameters
    ----------
    con : duckdb.DuckDBPyConnection
        Open connection.
    table : str
        Output table to inspect.

    Returns
    -------
    pandas.DataFrame
        Row counts and percentages by ``k_src``. The headline summary is
        logged rather than returned.
    """
    head = con.execute(f"""
        SELECT count(*) AS rows,
               count(DISTINCT ensg) AS genes,
               count(DISTINCT model_id) AS cell_lines,
               round(avg(k_src), 2) AS mean_k_src,
               round(avg(z_lineage), 3) AS mean_z,
               round(stddev_samp(z_lineage), 3) AS sd_z
        FROM "{table}"
    """).fetchdf()
    log.info("output summary:\n%s", head.to_string(index=False))
    return con.execute(f"""
        SELECT k_src, count(*) AS n,
               round(100.0 * count(*) / sum(count(*)) OVER (), 1) AS pct
        FROM "{table}" GROUP BY 1 ORDER BY 1
    """).fetchdf()


# ======================================================================
# 9. Plots
# ======================================================================

def _plt():
    """
    Import pyplot with a headless backend.

    Imported lazily and forced to ``Agg`` so the library never requires a
    display, and importing this module does not pull in matplotlib at all.

    Returns
    -------
    module
        ``matplotlib.pyplot``.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_grid_rho(grid: pd.DataFrame, p: Params, out: Path) -> None:
    """
    Heatmap of held-out concordance over the two censoring cuts.

    Sliced at the chosen ``floor_q`` and ``n0``, with the selected
    combination starred. The 1-SE value in the title is the scale against
    which apparent differences on the map should be judged.

    Parameters
    ----------
    grid : pandas.DataFrame
        Grid results from :func:`calibrate_params`.
    p : Params
        Chosen parameters.
    out : Path
        Output directory. Writes ``01_grid_rho.png``.

    Returns
    -------
    None
        Nothing is plotted for an empty grid or a slice of fewer than four
        points.
    """
    if grid is None or grid.empty:
        return
    plt = _plt()
    sub = grid[(grid.floor_q == p.floor_q) & (grid.n0 == p.n0)]
    if len(sub) < 4:
        return
    piv = sub.pivot_table(index="upper_cut", columns="lower_cut", values="rho")
    fig, ax = plt.subplots(figsize=(5.4, 4))
    im = ax.imshow(piv.values, cmap="magma", origin="lower", aspect="auto")
    ax.set_xticks(range(len(piv.columns)), [f"{c:g}" for c in piv.columns])
    ax.set_yticks(range(len(piv.index)), [f"{i:g}" for i in piv.index])
    ax.set_xlabel("lower_cut"); ax.set_ylabel("upper_cut")
    ax.set_title(f"held-out Spearman rho  (1 SE = {grid.attrs.get('se', 0):.3f})",
                 fontsize=10)
    if p.lower_cut in list(piv.columns) and p.upper_cut in list(piv.index):
        ax.scatter([list(piv.columns).index(p.lower_cut)],
                   [list(piv.index).index(p.upper_cut)],
                   marker="*", s=240, c="cyan", edgecolor="k")
    plt.colorbar(im, ax=ax); fig.tight_layout()
    fig.savefig(out / "01_grid_rho.png", dpi=150); plt.close(fig)


def plot_n0_sweep(grid: pd.DataFrame, p: Params, out: Path) -> None:
    """
    Concordance against the shrinkage constant, one line per upper cut.

    The honest reading is in the subtitle: a flat sweep means the data do
    not identify ``n0``, and whichever value the search picked is
    arbitrary within that range. Worth knowing before defending the number.

    Parameters
    ----------
    grid : pandas.DataFrame
        Grid results from :func:`calibrate_params`.
    p : Params
        Chosen parameters; ``n0`` is marked.
    out : Path
        Output directory. Writes ``02_n0_sweep.png``.

    Returns
    -------
    None
    """
    if grid is None or grid.empty:
        return
    plt = _plt()
    fig, ax = plt.subplots(figsize=(5.4, 3.6))
    for uc, g in grid[grid.floor_q == p.floor_q].groupby("upper_cut"):
        s = g.groupby("n0").rho.max()
        ax.plot(s.index, s.values, "o-", lw=1, label=f"upper={uc:g}")
    ax.axvline(p.n0, color="crimson", ls="--")
    ax.set_xlabel("n0 (lineage shrinkage constant)")
    ax.set_ylabel("held-out rho")
    ax.set_title("flat here means the data do not identify n0", fontsize=9)
    ax.legend(fontsize=7, ncol=2); fig.tight_layout()
    fig.savefig(out / "02_n0_sweep.png", dpi=150); plt.close(fig)


def plot_floors(axes: list[MethodAxis], out: Path) -> None:
    """
    Detection floor per method, annotated with its detected scale.

    Makes it visible when one method's floor sits far from its peers, which
    usually means its scale was misdetected rather than that its detection
    genuinely differs.

    Parameters
    ----------
    axes : list of MethodAxis
        Profiled axes.
    out : Path
        Output directory. Writes ``03_floors.png``.

    Returns
    -------
    None

    Notes
    -----
    An unset or infinite floor is drawn as 0 with its scale label shown as
    is, so an unprofiled method is visible rather than breaking the axis.
    """
    plt = _plt()
    fig, ax = plt.subplots(figsize=(max(6, 0.45 * len(axes)), 3.8))
    names = [a.name for a in axes]
    vals = [a.floor if (a.floor is not None and np.isfinite(a.floor)) else 0.0
            for a in axes]
    ax.bar(names, vals, color="0.55")
    for i, a in enumerate(axes):
        ax.text(i, vals[i], a.scale_detected or "?", ha="center", va="bottom",
                fontsize=6)
    ax.set_ylabel("detection floor (log2)")
    ax.tick_params(axis="x", rotation=70, labelsize=6)
    fig.tight_layout(); fig.savefig(out / "03_floors.png", dpi=150); plt.close(fig)


def plot_status_by_lineage(diag: pd.DataFrame, out: Path) -> None:
    """
    Status composition per lineage, as proportions.

    Shows which lineages are carried by properly scored cells and which
    lean on censored ranks — the lineages whose z-scores deserve the least
    weight downstream.

    Parameters
    ----------
    diag : pandas.DataFrame
        Aggregated chunk diagnostics.
    out : Path
        Output directory. Writes ``04_status_by_lineage.png``.

    Returns
    -------
    None
        Limited to the 18 largest lineages for legibility.
    """
    if diag is None or diag.empty:
        return
    plt = _plt()
    ct = diag.pivot_table(index="lineage", columns="status", values="n",
                          aggfunc="sum").fillna(0)
    ct = ct.loc[ct.sum(1).sort_values(ascending=False).index[:18]]
    fig, ax = plt.subplots(figsize=(7.5, 5))
    ct.div(ct.sum(1), axis=0).plot(kind="barh", stacked=True, ax=ax,
                                   colormap="Set2")
    ax.set_xlabel("proportion of (cell line, gene) cells"); ax.set_ylabel("")
    ax.legend(fontsize=7, loc="lower right"); fig.tight_layout()
    fig.savefig(out / "04_status_by_lineage.png", dpi=150); plt.close(fig)


def plot_coverage(diag: pd.DataFrame, out: Path) -> None:
    """
    Two bars: cells by contributing sources, and by scale provenance.

    Together these answer how much of the output rests on a single source
    and how much on a global rather than lineage scale — the two coverage
    caveats that qualify the whole layer.

    Parameters
    ----------
    diag : pandas.DataFrame
        Aggregated chunk diagnostics.
    out : Path
        Output directory. Writes ``05_coverage.png``.

    Returns
    -------
    None
    """
    if diag is None or diag.empty:
        return
    plt = _plt()
    fig, axs = plt.subplots(1, 2, figsize=(10, 3.6))
    kk = diag.groupby("k_src").n.sum()
    axs[0].bar(kk.index.astype(str), kk.values, color="#3182bd")
    axs[0].set_xlabel("sources contributing (k_src)"); axs[0].set_ylabel("cells")
    ss = diag.groupby("scale_source").n.sum()
    axs[1].bar(ss.index.astype(str), ss.values, color="#756bb1")
    axs[1].set_xlabel("scale source"); axs[1].tick_params(axis="x", labelsize=7)
    fig.tight_layout(); fig.savefig(out / "05_coverage.png", dpi=150); plt.close(fig)


def make_plots(grid, p, axes, diag, out: Path) -> None:
    """
    Render the full diagnostic figure set.

    Parameters
    ----------
    grid : pandas.DataFrame
        Grid results; may be empty when parameters were supplied rather
        than calibrated, in which case the two grid plots are skipped.
    p : Params
        Chosen parameters.
    axes : list of MethodAxis
        Profiled axes.
    diag : pandas.DataFrame
        Aggregated diagnostics; may be empty in calibrate-only mode, in
        which case the two coverage plots are skipped.
    out : Path
        Output directory, created if absent.

    Returns
    -------
    None
        Writes ``01_grid_rho.png`` through ``05_coverage.png``.
    """
    out.mkdir(parents=True, exist_ok=True)
    plot_grid_rho(grid, p, out)
    plot_n0_sweep(grid, p, out)
    plot_floors(axes, out)
    plot_status_by_lineage(diag, out)
    plot_coverage(diag, out)