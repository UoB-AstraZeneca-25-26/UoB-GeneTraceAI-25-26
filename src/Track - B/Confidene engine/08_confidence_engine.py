"""
08_confidence_engine.py  —  CELL 1: CONFIDENCE-INPUT CHARACTERISATION (measures only)
CellLineSelector - Confidence Engine. Companion to the Confidence Engine Whiteboard.

WHAT THIS IS
------------
Before any weight, threshold or band edge is chosen, this cell turns the confidence
inputs already built into features/expression_features.parquet into concrete numbers,
in the same "characterise before deciding" discipline as 03_expression_diagnostics.py.
It answers the questions that decide the v1 scoring rules:

  1. SCHEMA        What columns / rows / lines / genes are actually present?
  2. COVERAGE      Distribution of n_sources, corroboration_available, single_source_gene,
                   and per-source presence -> the evidence-depth dimension.
  3. AGREEMENT     THE key check. The stored detection_agreement/n_sources_detected fold in
                   detected_geo. Per decision E3 (GEO = corroboration-only, no present/absent
                   call), this RECOMPUTES agreement over DepMap+HPA only, measures how far the
                   stored column drifts, and - crucially - shows how often a second detection
                   call (HPA) exists at all. Where it does not, agreement is trivially 1.0 and
                   the engine cannot lean on it.
  4. PROVENANCE    biotype_coarse / in_default_universe -> the quality dimension.
  5. RELIABILITY   expr_percentile_depmap vs detected_depmap -> can percentile proxy the
                   near-threshold "distance from detection" signal (z-score is not stored)?
  6. MODALITY      Per-modality line-level coverage from the frozen coverage contract
                   (four-modality scope): expression, mutation, fusion, signatures.
  7. MOD INPUTS    (optional) If mutation/fusion/signature feature tables are passed, report
                   their confidence-relevant columns - incl. whether Arriba confidence/
                   split_reads survived fusion aggregation (an open item in the whiteboard).

WHAT THIS IS *NOT*
------------------
It does not score, weight, threshold, band, gate, or combine anything. Every function
returns measured numbers; the scoring rules stay for the next cell, decided against these
numbers, not in the abstract. GEO is treated as corroboration-only throughout (E3): it
counts toward coverage (n_sources) but is excluded from detection agreement.

USAGE (in your live Jupyter session)
------------------------------------
    import importlib
    ce = importlib.import_module("08_confidence_engine")

    # Preferred: pass tables already in memory
    summary = ce.run_characterisation(
        expr       = expression_features,     # or leave None to load the parquet
        coverage   = coverage_contract,       # or None -> loads verification_out/coverage_contract.parquet
        mutations  = mut_features,            # optional
        fusions    = fus_features,            # optional (checks for Arriba confidence/split_reads)
        signatures = sig_features,            # optional
    )

    # Or let it load from disk:
    summary = ce.run_characterisation()       # reads features/expression_features.parquet

It prints a digest and writes features/confidence/confidence_input_characterisation.json.
Paste that JSON (or the printed digest) back and we fix the v1 rules against it.

Author: (CellLineSelector team)   |   Depends on: pandas, numpy
"""

from __future__ import annotations

import json
import os
from typing import Iterable

import numpy as np
import pandas as pd


# =============================================================================
# CONFIG  (paths used only when a table is not passed in memory)
# =============================================================================

EXPRESSION_FEATURES_PATH = "features/expression_features.parquet"
COVERAGE_CONTRACT_PATH   = "verification_out/coverage_contract.parquet"
OUTDIR                   = "features/confidence"

# The four built modalities we care about for the coverage dimension (contract columns).
_BUILT_MODALITY_COLS = ("gene_expression", "mutation", "fusion", "signatures")

# Candidate column names for optional modality-input introspection (char 7).
_MUT_CONF_COLS = ("is_mutated", "is_damaging", "has_hotspot", "n_altering", "n_variants")
_FUS_CONF_COLS = ("confidence", "split_reads1", "split_reads2", "reading_frame",
                  "retained_protein_domains", "tags")
_SIG_CONF_COLS = ("wgd", "ploidy", "msi", "aneuploidy", "cin", "loh")


# =============================================================================
# SHARED HELPERS  (mirroring 03_expression_diagnostics.py)
# =============================================================================

def _find_col(df: pd.DataFrame, candidates: Iterable[str], required: bool = True,
              what: str = "column") -> str | None:
    """First candidate present in df (case-insensitive), else None (or raise if required)."""
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    if required:
        raise KeyError(f"Could not resolve {what}: tried {list(candidates)}; "
                       f"available = {list(df.columns)[:40]}")
    return None


def _as_bool_num(s: pd.Series) -> pd.Series:
    """Coerce a bool-ish column (True/False/1/0/'True'/NaN) to float 1.0/0.0/NaN.

    Present-but-not-measured stays NaN so it is excluded from call counts.
    """
    if s.dtype == bool:
        return s.astype("float64")
    def _m(v):
        if pd.isna(v):
            return np.nan
        if isinstance(v, (bool, np.bool_)):
            return 1.0 if v else 0.0
        if isinstance(v, (int, float, np.integer, np.floating)):
            return float(v != 0)
        t = str(v).strip().lower()
        if t in ("true", "t", "1", "yes", "detected", "present"):
            return 1.0
        if t in ("false", "f", "0", "no", "absent", "not detected"):
            return 0.0
        return np.nan
    return s.map(_m).astype("float64")


def _summ(values) -> dict:
    """Distribution summary mirroring the verification run (NaNs excluded)."""
    v = pd.to_numeric(pd.Series(values), errors="coerce").to_numpy(dtype="float64")
    v = v[~np.isnan(v)]
    if v.size == 0:
        return {"n": 0}
    s = pd.Series(v)
    return {
        "n": int(v.size),
        "min": float(v.min()), "max": float(v.max()),
        "mean": float(v.mean()), "median": float(np.median(v)),
        "q01": float(np.quantile(v, 0.01)), "q25": float(np.quantile(v, 0.25)),
        "q75": float(np.quantile(v, 0.75)), "q99": float(np.quantile(v, 0.99)),
        "pct_zero": float((v == 0).mean() * 100.0),
        "skew": float(s.skew()) if v.size > 2 else float("nan"),
    }


def _vc_pct(s: pd.Series, dropna: bool = False) -> dict:
    """Value counts as {value: (count, pct)} - small integer/boolean distributions."""
    vc = s.value_counts(dropna=dropna)
    n = int(vc.sum())
    return {str(k): [int(v), round(v / n * 100, 2) if n else 0.0] for k, v in vc.items()}


def _pct_true(s: pd.Series) -> float:
    """Percent of non-null entries that are truthy (bool-ish tolerant)."""
    b = _as_bool_num(s)
    b = b[~b.isna()]
    return float(b.mean() * 100.0) if len(b) else float("nan")


def _outdir(path: str | None) -> str | None:
    if not path:
        return None
    os.makedirs(path, exist_ok=True)
    return path


def _load(df, path: str, what: str):
    """Return df if given, else read parquet at path, else None with a clear message."""
    if df is not None:
        return df
    if path and os.path.exists(path):
        try:
            return pd.read_parquet(path)
        except Exception as e:
            print(f"[load] could not read {what} at {path}: {type(e).__name__}: {e}")
            return None
    print(f"[load] {what} not provided and not found at {path} - skipped")
    return None


# =============================================================================
# CHAR 1 - SCHEMA
# =============================================================================

def char_1_schema(expr: pd.DataFrame) -> dict:
    r = {"n_rows": int(len(expr)), "columns": list(expr.columns)}
    mcol = _find_col(expr, ("model_id",), required=False)
    gcol = _find_col(expr, ("gene", "gene_id", "ensembl_gene_id"), required=False)
    if mcol:
        r["n_lines"] = int(expr[mcol].nunique())
    if gcol:
        r["n_genes"] = int(expr[gcol].nunique())
    print(f"[schema] rows={r['n_rows']:,}  lines={r.get('n_lines','?')}  genes={r.get('n_genes','?')}")
    print(f"[schema] columns: {', '.join(expr.columns)}")
    return r


# =============================================================================
# CHAR 2 - COVERAGE  (evidence-depth dimension)
# =============================================================================

def char_2_coverage(expr: pd.DataFrame) -> dict:
    r = {}
    ns = _find_col(expr, ("n_sources",), required=False)
    if ns:
        r["n_sources"] = _vc_pct(expr[ns])
        print(f"[coverage] n_sources: {r['n_sources']}")
    for col, key in (("corroboration_available", "corroboration_available_pct"),
                     ("single_source_gene", "single_source_gene_pct")):
        c = _find_col(expr, (col,), required=False)
        if c:
            r[key] = round(_pct_true(expr[c]), 2)
            print(f"[coverage] {col}: {r[key]}% of rows")
    # per-source presence (value present, not detection)
    pres = {}
    for col, name in (("depmap_log2tpm1", "depmap"), ("geo_intensity", "geo"), ("hpa_ntpm", "hpa")):
        c = _find_col(expr, (col,), required=False)
        if c:
            pres[name] = round(float(expr[c].notna().mean()) * 100, 2)
    if pres:
        r["per_source_present_pct"] = pres
        print(f"[coverage] per-source present: {pres}")
    print("[coverage] NOTE the table is DepMap-anchored: every row has DepMap, so n_sources "
          "counts DepMap + any GEO/HPA corroboration. single_source == DepMap-only.")
    return r


# =============================================================================
# CHAR 3 - DETECTION AGREEMENT  (the decisive check; GEO excluded per E3)
# =============================================================================

def char_3_agreement(expr: pd.DataFrame) -> dict:
    r = {}
    c_dm  = _find_col(expr, ("detected_depmap",), required=False)
    c_hpa = _find_col(expr, ("detected_hpa",), required=False)
    c_geo = _find_col(expr, ("detected_geo",), required=False)
    c_sag = _find_col(expr, ("detection_agreement",), required=False)
    c_snd = _find_col(expr, ("n_sources_detected",), required=False)

    # (a) stored values, as they currently sit in the parquet (GEO folded in)
    if c_sag:
        r["stored_agreement_summary"] = _summ(expr[c_sag])
        r["stored_agreement_valuecounts"] = _vc_pct(expr[c_sag].round(4))
        print(f"[agree] STORED detection_agreement (GEO-included): "
              f"median={r['stored_agreement_summary'].get('median')}, "
              f"values={r['stored_agreement_valuecounts']}")
    if c_geo is not None:
        print("[agree] NOTE detected_geo IS present in the table and IS folded into the stored "
              "columns; per decision E3 it is excluded below.")

    if c_dm is None:
        print("[agree] detected_depmap missing - cannot recompute; skipped")
        return r

    # (b) recompute agreement over DepMap + HPA only (E3: GEO excluded)
    dm  = _as_bool_num(expr[c_dm])
    hpa = _as_bool_num(expr[c_hpa]) if c_hpa else pd.Series(np.nan, index=expr.index)
    calls = pd.concat({"dm": dm, "hpa": hpa}, axis=1)
    n_calls = calls.notna().sum(axis=1).astype("float64")          # present present/absent calls
    n_det   = calls.sum(axis=1, skipna=True)                       # of those, how many "detected"
    with np.errstate(invalid="ignore", divide="ignore"):
        agree = np.where(n_calls > 0,
                         np.maximum(n_det, n_calls - n_det) / n_calls,
                         np.nan)
    agree = pd.Series(agree, index=expr.index)

    r["n_detection_calls_dist"] = _vc_pct(n_calls)                 # how many rows have 1 vs 2 calls
    r["detection_corroborated_pct"] = round(float((n_calls >= 2).mean() * 100), 2)
    r["single_detection_call_pct"] = round(float((n_calls == 1).mean() * 100), 2)
    r["agreement_corrected_summary"] = _summ(agree)
    print(f"[agree] n_detection_calls (DepMap+HPA only): {r['n_detection_calls_dist']}")
    print(f"[agree] >>> detection-corroborated (>=2 calls): {r['detection_corroborated_pct']}% of rows")
    print(f"[agree] >>> single detection call (agreement trivially 1.0): "
          f"{r['single_detection_call_pct']}% of rows")

    # (c) the informative subset: agreement WHERE a second call exists
    real = agree[n_calls >= 2]
    if len(real):
        r["agreement_where_corroborated_valuecounts"] = _vc_pct(real.round(4))
        r["disagreement_pct_of_corroborated"] = round(float((real < 1.0).mean() * 100), 2)
        print(f"[agree] where corroborated: {r['agreement_where_corroborated_valuecounts']} "
              f"({r['disagreement_pct_of_corroborated']}% show DepMap/HPA disagreement)")

    # (d) drift of the corrected agreement vs the stored (GEO-included) column
    if c_sag:
        stored = pd.to_numeric(expr[c_sag], errors="coerce")
        changed = (stored.round(6) != agree.round(6)) & ~(stored.isna() & agree.isna())
        r["rows_changed_vs_stored"] = int(changed.sum())
        r["rows_changed_vs_stored_pct"] = round(float(changed.mean() * 100), 2)
        print(f"[agree] excluding GEO changes detection_agreement on "
              f"{r['rows_changed_vs_stored']:,} rows ({r['rows_changed_vs_stored_pct']}%) "
              f"-> confirms the stored column must be recomputed (whiteboard L5).")

    # (e) per-source detection rates
    rates = {}
    if c_dm:
        rates["detected_depmap_pct"] = round(_pct_true(expr[c_dm]), 2)
    if c_hpa:
        rates["detected_hpa_pct_where_measured"] = round(_pct_true(expr[c_hpa]), 2)
    if rates:
        r["per_source_detection_pct"] = rates
        print(f"[agree] per-source detection: {rates}")

    # (f) the gap between level-corroboration (GEO counts) and detection-corroboration
    ca = _find_col(expr, ("corroboration_available",), required=False)
    if ca:
        level_corr = _as_bool_num(expr[ca]) == 1.0
        det_corr = (n_calls >= 2)
        gap = level_corr & ~det_corr
        r["level_corr_but_not_detection_corr_pct"] = round(float(gap.mean() * 100), 2)
        print(f"[agree] rows with corroboration_available but NO 2nd detection call "
              f"(GEO-only corroboration): {r['level_corr_but_not_detection_corr_pct']}% "
              f"-> gate detection-agreement on 2 calls, not on corroboration_available.")
    return r


# =============================================================================
# CHAR 4 - PROVENANCE / QUALITY
# =============================================================================

def char_4_provenance(expr: pd.DataFrame) -> dict:
    r = {}
    bc = _find_col(expr, ("biotype_coarse", "gene_biotype"), required=False)
    if bc:
        r["biotype_coarse"] = _vc_pct(expr[bc].astype(str))
        print(f"[prov] {bc}: {r['biotype_coarse']}")
    idu = _find_col(expr, ("in_default_universe",), required=False)
    if idu:
        r["in_default_universe_pct"] = round(_pct_true(expr[idu]), 2)
        r["outside_default_universe_pct"] = round(100 - r["in_default_universe_pct"], 2)
        print(f"[prov] in_default_universe: {r['in_default_universe_pct']}% "
              f"(outside: {r['outside_default_universe_pct']}%)")
    return r


# =============================================================================
# CHAR 5 - MEASUREMENT-RELIABILITY PROXY  (percentile vs detection; z not stored)
# =============================================================================

def char_5_reliability(expr: pd.DataFrame) -> dict:
    r = {}
    pc = _find_col(expr, ("expr_percentile_depmap",), required=False)
    dm = _find_col(expr, ("detected_depmap",), required=False)
    if pc is None:
        print("[reliab] expr_percentile_depmap missing - skipped")
        return r
    r["percentile_summary"] = _summ(expr[pc])
    print(f"[reliab] expr_percentile_depmap: {r['percentile_summary']}")
    if dm:
        b = _as_bool_num(expr[dm])
        r["percentile_where_detected"] = _summ(expr.loc[b == 1.0, pc])
        r["percentile_where_not_detected"] = _summ(expr.loc[b == 0.0, pc])
        md = r["percentile_where_detected"].get("median")
        mn = r["percentile_where_not_detected"].get("median")
        print(f"[reliab] percentile median | detected={md} vs not-detected={mn}")
        print("[reliab] NOTE the zFPKM z-score is not stored; if percentile separates "
              "detected/not cleanly it can proxy 'distance from threshold' for C2, else z "
              "must be recomputed from depmap_log2tpm1.")
    return r


# =============================================================================
# CHAR 6 - MODALITY COVERAGE  (four-modality, from the frozen coverage contract)
# =============================================================================

def char_6_modality_coverage(coverage: pd.DataFrame) -> dict:
    r = {}
    n = int(len(coverage))
    r["n_lines"] = n
    per_mod = {}
    for col in coverage.columns:
        if col in ("model_id", "n_modalities", "complete_cycle"):
            continue
        present = int(_as_bool_num(coverage[col]).fillna(0).sum()) if coverage[col].dtype != object \
            else int(coverage[col].astype(bool).sum())
        per_mod[col] = {"present": present, "missing": n - present,
                        "present_pct": round(present / n * 100, 2) if n else 0.0}
    r["per_modality"] = per_mod
    print(f"[modcov] lines={n:,}")
    for col in _BUILT_MODALITY_COLS:
        if col in per_mod:
            m = per_mod[col]
            print(f"[modcov]   {col:16} present={m['present']:>5} "
                  f"missing={m['missing']:>4} ({m['present_pct']}%)")
    nm = _find_col(coverage, ("n_modalities",), required=False)
    if nm:
        r["n_modalities_dist"] = _vc_pct(coverage[nm])
        print(f"[modcov] n_modalities: {r['n_modalities_dist']}")
    cc = _find_col(coverage, ("complete_cycle",), required=False)
    if cc:
        r["complete_cycle_pct"] = round(_pct_true(coverage[cc]), 2)
        print(f"[modcov] complete_cycle: {r['complete_cycle_pct']}%")
    return r


# =============================================================================
# CHAR 7 - MODALITY CONFIDENCE INPUTS  (optional; incl. Arriba open item)
# =============================================================================

def _describe_modality(df: pd.DataFrame, cols: Iterable[str], tag: str) -> dict:
    r = {"n_rows": int(len(df))}
    mcol = _find_col(df, ("model_id",), required=False)
    if mcol:
        r["n_lines"] = int(df[mcol].nunique())
    found = {}
    for c in cols:
        real = _find_col(df, (c,), required=False)
        if real is None:
            continue
        s = df[real]
        non_null_pct = round(float(s.notna().mean() * 100), 2)
        entry = {"present_column": True, "non_null_pct": non_null_pct}
        # small-cardinality -> value counts; else numeric summary
        if s.dropna().nunique() <= 12:
            entry["valuecounts"] = _vc_pct(s.astype(str))
        else:
            entry["summary"] = _summ(s)
        found[c] = entry
    r["confidence_columns"] = found
    print(f"[{tag}] rows={r['n_rows']:,} lines={r.get('n_lines','?')}; "
          f"confidence cols found: {list(found.keys()) or 'NONE'}")
    return r


def char_7_modality_inputs(mutations, fusions, signatures) -> dict:
    r = {}
    if mutations is not None:
        r["mutations"] = _describe_modality(mutations, _MUT_CONF_COLS, "mut")
    if fusions is not None:
        r["fusions"] = _describe_modality(fusions, _FUS_CONF_COLS, "fus")
        # explicit answer to the whiteboard open item
        fc = r["fusions"]["confidence_columns"]
        arriba = {k: (k in fc) for k in ("confidence", "split_reads1", "split_reads2", "reading_frame")}
        r["fusions"]["arriba_strength_fields_retained"] = arriba
        print(f"[fus] Arriba strength fields retained through aggregation? {arriba}")
    if signatures is not None:
        r["signatures"] = _describe_modality(signatures, _SIG_CONF_COLS, "sig")
    if not r:
        print("[modinputs] no modality tables passed - skipped (pass mutations=/fusions=/signatures=)")
    return r


# =============================================================================
# ORCHESTRATOR
# =============================================================================

def run_characterisation(expr=None, coverage=None,
                         mutations=None, fusions=None, signatures=None,
                         expr_path: str = EXPRESSION_FEATURES_PATH,
                         coverage_path: str = COVERAGE_CONTRACT_PATH,
                         outdir: str = OUTDIR, write: bool = True) -> dict:
    """Run every characterisation, print a digest, write a paste-friendly JSON. Chooses nothing."""
    print("=" * 82)
    print("CONFIDENCE-INPUT CHARACTERISATION  -  measures only, chooses no weight/threshold")
    print("GEO is corroboration-only (E3): counts toward coverage, excluded from detection agreement")
    print("=" * 82)

    expr = _load(expr, expr_path, "expression_features")
    coverage = _load(coverage, coverage_path, "coverage_contract")

    summary: dict = {}
    if expr is not None:
        for key, fn in (("schema", char_1_schema),
                        ("coverage", char_2_coverage),
                        ("agreement", char_3_agreement),
                        ("provenance", char_4_provenance),
                        ("reliability", char_5_reliability)):
            print(f"\n--- {key} ---")
            try:
                summary[key] = fn(expr)
            except Exception as e:
                summary[key] = {"error": f"{type(e).__name__}: {e}"}
                print(f"[{key}] ERROR: {type(e).__name__}: {e}")
    else:
        summary["expression"] = {"error": "expression_features not available"}

    if coverage is not None:
        print("\n--- modality_coverage ---")
        try:
            summary["modality_coverage"] = char_6_modality_coverage(coverage)
        except Exception as e:
            summary["modality_coverage"] = {"error": f"{type(e).__name__}: {e}"}
            print(f"[modcov] ERROR: {type(e).__name__}: {e}")

    if any(x is not None for x in (mutations, fusions, signatures)):
        print("\n--- modality_inputs ---")
        try:
            summary["modality_inputs"] = char_7_modality_inputs(mutations, fusions, signatures)
        except Exception as e:
            summary["modality_inputs"] = {"error": f"{type(e).__name__}: {e}"}
            print(f"[modinputs] ERROR: {type(e).__name__}: {e}")

    if write:
        od = _outdir(outdir)
        if od:
            path = os.path.join(od, "confidence_input_characterisation.json")
            with open(path, "w") as f:
                json.dump(summary, f, indent=2, default=str)
            print(f"\n[run] wrote {path}")

    print("\n" + "=" * 82)
    print("READ THESE FIRST, THEN WE SET RULES (not before):")
    print("  agreement.detection_corroborated_pct  -> if low, detection agreement is mostly")
    print("     trivial; weight shifts to coverage + provenance + reliability.")
    print("  agreement.rows_changed_vs_stored_pct  -> size of the GEO-exclusion correction (L5).")
    print("  coverage.n_sources                    -> how much corroboration exists at all.")
    print("  reliability percentile split          -> whether percentile can proxy near-threshold.")
    print("  modality_coverage.per_modality        -> the four-modality coverage dimension.")
    print("Paste confidence_input_characterisation.json back and we fix the v1 rules to it.")
    print("=" * 82)
    return summary


if __name__ == "__main__":
    print(__doc__)
    print("Import this module and call run_characterisation(...) from your live Jupyter session.")
