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

STATUS = {0: "ok", 1: "censored_ranked", 2: "below_detection", 3: "no_data"}
SCALE_SRC = {0: "none", 1: "lineage", 2: "shrunk", 3: "fallback_global"}

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
    floor_q: float
    n0: int
    upper_cut: float
    lower_cut: float


@dataclass
class MethodAxis:
    """Row axis for one method: raw sample ids, their model_id, replicate map."""
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
    """Parse one GSMxxxx.txt SOFT sample record into a flat dict."""
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
    """!Sample_data_processing -> (proc_family, scale_hint)."""
    t = (text or "").lower()
    for name, pat, scale in PROC_RULES:
        if re.search(pat, t):
            return name, scale
    return "unknown", "unknown"


def build_geo_platform_frame(meta_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Parse geo_meta/*.txt -> (records, parse_failures)."""
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
    con.register("_geo_platform_tmp", df)
    try:
        con.execute(
            f'CREATE OR REPLACE TABLE "{table}" AS SELECT * FROM _geo_platform_tmp')
    finally:
        con.unregister("_geo_platform_tmp")
    log.info("wrote %s (%d rows)", table, len(df))


def audit_geo_platform(con, df: pd.DataFrame) -> pd.DataFrame:
    """Before/after coverage audit. Returns the method summary table."""
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
    s = (con.execute("""
            SELECT model_id, lineage FROM sample_info WHERE model_id IS NOT NULL
         """).fetchdf()
         .drop_duplicates("model_id").set_index("model_id")["lineage"])
    miss = int(s.isna().sum())
    if miss:
        log.warning("%d models in sample_info have NULL lineage -> 'unknown'", miss)
    return s.fillna("unknown")


def normalise_ensg(value) -> str | None:
    """Return the canonical ENSG prefix, stripping version suffixes such as .7."""
    if value is None:
        return None
    s = str(value).strip().upper()
    if not s:
        return None
    m = _ENSG.fullmatch(s)
    return m.group(1).upper() if m else None


def filter_valid_ensg_ids(genes, valid_ids=None) -> list[str]:
    """Keep only canonical, roster-approved ENSG IDs in a consistent order."""
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
    """Gene roster IDs, normalised to the canonical ENSG form used in scoring."""
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
    """Union of gene axes restricted to canonical ENSG IDs from the gene roster."""
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
    """One read of a wide table for a gene chunk. Indexed by upper(id)."""
    present = [(g, colmap[g]) for g in genes if g in colmap]
    if not present:
        return pd.DataFrame()
    sel = ", ".join(f'"{c}"::DOUBLE AS "{g}"' for g, c in present)
    df = con.execute(
        f'SELECT upper(trim("{id_col}")) AS _id, {sel} FROM "{table}"').fetchdf()
    return df.set_index("_id")


def fetch_hpa_chunk(con, genes: list[str]) -> pd.DataFrame:
    """Long -> wide for a gene chunk. median() collapses transcript duplicates."""
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
    """Fast uppercase/strip conversion for large arrays without pandas apply."""
    arr = np.asarray(list(ids), dtype=object)
    out = np.empty(arr.shape, dtype=object)
    for i, v in enumerate(arr):
        s = str(v).strip().upper()
        out[i] = s
    return out


def align_axis_rows(ax: MethodAxis, wide: pd.DataFrame,
                    genes: list[str]) -> np.ndarray:
    """Slice a shared table chunk down to this axis's raw rows."""
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
    """Linear vs log2 from the data alone. GEO metadata is not trusted here."""
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
    """Set ax.scale_detected and ax.floor in place. Warns on SOFT/data conflict."""
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
    """Per-source stacks of method-level z, aligned to the global model axis."""
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
    grid = grid or PARAM_GRID
    keys = ["floor_q", "n0", "upper_cut", "lower_cut"]
    return [c for c in itertools.product(*[grid[k] for k in keys]) if c[3] < c[2]]


def sample_calibration_genes(all_genes, n=N_CALIB_GENES, rng=None) -> list[str]:
    rng = rng or np.random.default_rng(0)
    n = min(n, len(all_genes))
    return list(rng.choice(np.asarray(all_genes, dtype=object), size=n, replace=False))


def _score_params(con, axes, colmaps, genes, models, midx, lineage_s,
                  p: Params, holdouts, ref_cache) -> float | None:
    """Held-out cross-source rank concordance for one parameter combination."""
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
    """Grid search + 1-SE rule. Returns (Params, grid_frame, genes_used)."""
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
    if df.empty:
        return pd.DataFrame()
    return (df.groupby(["lineage", "status", "scale_source", "k_src"],
                       observed=True).size().rename("n").reset_index())


def verify_output(con, table: str) -> pd.DataFrame:
    """Post-run sanity: row/gene/line counts and the k_src breakdown."""
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
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    return plt


def plot_grid_rho(grid: pd.DataFrame, p: Params, out: Path) -> None:
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
    out.mkdir(parents=True, exist_ok=True)
    plot_grid_rho(grid, p, out)
    plot_n0_sweep(grid, p, out)
    plot_floors(axes, out)
    plot_status_by_lineage(diag, out)
    plot_coverage(diag, out)