"""
01_verification_checks.py
CellLineSelector - Phase 1 -> Feature Engineering: pre-flight verification & scope freeze.

PURPOSE
-------
Before a single feature is engineered, the Phase 1 planning document (Section L)
requires four things to be settled:

    L.1  Freeze the identity / coverage contract.
    L.2  Resolve the expression-scale question (Q1) - the blocking dependency.
    L.3  Verify the scales / fields marked "Needs verification" (Section K).
    L.4  Formally decide which modalities are in scope and descope the rest.

This module implements those as small, independently runnable, documented checks.
Each check consumes an already-harmonised DataFrame (the output of
`00_data_harmonisation_pipeline.ipynb`), returns a structured `dict` result, and
prints a human-readable verdict. Nothing here mutates your data - it only inspects.

WHY A MODULE (not a notebook)
-----------------------------
The checks are reusable and should produce a machine-readable audit log that can be
cited in the dissertation. Import them into a notebook and call cell-by-cell, or run
`python 01_verification_checks.py` after filling in the loader at the bottom.

RESULT CONTRACT
---------------
Every check returns a dict with at least:
    {
      "check":        str,   # function name
      "planning_ref": str,   # which planning-doc item this resolves
      "status":       str,   # "PASS" | "FLAG" | "INFO" | "ERROR"
      "verdict":      str,   # one-line plain-English conclusion
      "details":      dict,  # the supporting numbers
    }

`status` semantics:
    PASS  - expectation confirmed, no action needed.
    FLAG  - a real issue confirmed; a design decision is now required.
    INFO  - characterised successfully; feeds a later decision, not pass/fail.
    ERROR - the check could not run (usually a missing/renamed column).

IMPORTANT
---------
The scale verdicts (log vs linear vs log-ratio) are *heuristics* to be cross-checked
against each source's official documentation, not proofs. They are designed to catch
the confirmed DepMap-log / GEO-linear / HPA-linear mismatch and to characterise the
modalities whose scale is currently unknown - not to replace reading the source methods.

Author: (CellLineSelector team)   |   Depends on: pandas, numpy   |   Optional: scipy
"""

from __future__ import annotations

import json
import hashlib
from datetime import datetime, timezone
from typing import Iterable, Sequence

import numpy as np
import pandas as pd

try:
    from scipy.stats import skew as _scipy_skew
    _HAVE_SCIPY = True
except Exception:  # scipy is optional
    _HAVE_SCIPY = False


PLANNING_DOC = "CellLineSelector_FeatureEngineering_Phase1_Planning"


# =============================================================================
# SECTION 0 - MODALITY SCOPE REGISTRY  (implements L.4 - the descope decision)
# =============================================================================

# Status vocabulary (keep these exact strings - downstream code branches on them):
#   IN_SCOPE            - present, usable, build features from it.
#   IN_SCOPE_PENDING    - present, but a named verification must PASS first.
#   BLOCKED             - present but unusable until a specific blocker is cleared;
#                         if the blocker cannot be cleared, flip to DESCOPED.
#   PENDING_ACQUISITION - NOT in the harmonised data yet, but a standard, low-cost
#                         source exists (same DepMap portal + ACH key). Acquire it,
#                         run it through the existing cleaning, then promote to
#                         IN_SCOPE. NO feature code may depend on it until acquired.
#   DESCOPED            - not in the harmonised data and NOT planned for acquisition;
#                         NO feature code may depend on it. Re-scope only by acquiring.
#   EXTERNAL_RESOURCE   - not a cell-line measurement; an external knowledge base to be
#                         acquired at the specific feature stage that needs it.
#
# This dict is the single source of truth for what the feature layer is allowed to
# build. It is deliberately explicit so the descoping is traceable in the dissertation.

MODALITY_SCOPE: dict[str, dict] = {
    # ---- present and usable -------------------------------------------------
    "rna_depmap": {
        "status": "IN_SCOPE_PENDING",
        "coverage_note": "~1428 lines; log2(TPM+1)",
        "blocker_or_gate": "compare_expression_scales must PASS (Q1 scale reconciliation)",
        "rationale": "Primary expression source; highest gene coverage (53,961 ENSG).",
    },
    "rna_geo": {
        "status": "IN_SCOPE_PENDING",
        "coverage_note": "~2308 GSM samples; linear scale",
        "blocker_or_gate": "compare_expression_scales + batch audit (Q1, Q14)",
        "rationale": "Adds samples/studies but on a different scale from DepMap.",
    },
    "rna_hpa": {
        "status": "IN_SCOPE_PENDING",
        "coverage_note": "~1105 lines; linear TPM/pTPM/nTPM",
        "blocker_or_gate": "verify_hpa_variants must choose tpm/ptpm/ntpm",
        "rationale": "Independent RNA source useful for cross-source agreement.",
    },
    "protein": {
        "status": "IN_SCOPE_PENDING",
        "coverage_note": "375 lines (20.4%); TMT log-ratio",
        "blocker_or_gate": "verify_proteomics_scale + UniProt->gene map (protein_map)",
        "rationale": "Enables RNA->protein confirmation (D3) despite low coverage.",
    },
    "mutation": {
        "status": "IN_SCOPE",
        "coverage_note": "1697 lines (92.2%)",
        "blocker_or_gate": None,
        "rationale": "High-coverage genomic suitability signal (D5).",
    },
    "fusion": {
        "status": "IN_SCOPE",
        "coverage_note": "1428 lines (77.6%); gene-PAIR grain",
        "blocker_or_gate": "grain: emit per-partner or per-pair (design decision)",
        "rationale": "Structural-variant suitability signal (D5).",
    },
    "signatures": {
        "status": "IN_SCOPE",
        "coverage_note": "1704 lines (92.6%); mixed value types",
        "blocker_or_gate": "per-column handling (binary WGD, count aneuploidy, cont. ploidy)",
        "rationale": "Genome-instability context (D11); highest coverage after mutation.",
    },
    "lineage": {
        "status": "IN_SCOPE_PENDING",
        "coverage_note": "1840 lines (from sample_info)",
        "blocker_or_gate": "audit_lineage_vocabulary (controlled-vocab check)",
        "rationale": "Cell-line context (D6); needs vocabulary consistency.",
    },
    "metabolomics": {
        "status": "IN_SCOPE_PENDING",
        "coverage_note": "927 lines (50.4%); 225 metabolites",
        "blocker_or_gate": "verify_metabolomics_scale (units/normalisation unknown)",
        "rationale": "Per-line metabolite abundance; scale must be characterised first.",
    },

    # ---- present but blocked ------------------------------------------------
    "mirna": {
        "status": "BLOCKED",
        "coverage_note": "914 lines (49.7%); ~953 features",
        "blocker_or_gate": "verify_mirna_identifiers: opaque 'nmir*' IDs must map to miRBase",
        "rationale": "Optional modality; unusable until identifiers are resolvable. "
                     "If mapping fails, flip to DESCOPED.",
    },

    # ---- NOT in the harmonised data YET -> acquire (standard DepMap files) ---
    # These are high-priority acquisitions, NOT descopes: both are standard DepMap
    # release files keyed on the same model_id (ACH) already resolved by this
    # pipeline, so integration cost is near-zero (download -> existing cleaning ->
    # join on model_id). Verify the exact filenames/coverage against your DepMap
    # release, as they can change between releases.
    "copy_number": {
        "status": "PENDING_ACQUISITION",
        "coverage_note": "ABSENT now; DepMap covers ~nearly all lines once acquired",
        "acquire": "DepMap OmicsCNGene.csv (gene-level relative copy number); "
                   "segment-level file optional. Same ACH key.",
        "blocker_or_gate": "acquire + clean + join on model_id, then promote to IN_SCOPE",
        "rationale": "Directly required by D5 (the CNA in 'mutation, CNA, or fusion'). "
                     "Amplifications/deep deletions (e.g. ERBB2, MYC amp; CDKN2A, PTEN "
                     "loss) are first-class selection signals that mutation status alone "
                     "misses. HIGH impact, near-zero acquisition cost.",
    },
    "crispr_dependency": {
        "status": "PENDING_ACQUISITION",
        "coverage_note": "ABSENT now; DepMap CRISPR typically ~1000+ lines once acquired",
        "acquire": "DepMap CRISPRGeneEffect.csv (Chronos gene-effect scores); "
                   "CRISPRGeneDependency.csv (dependency probability) optional. Same ACH key.",
        "blocker_or_gate": "acquire + clean + join on model_id, then promote to IN_SCOPE",
        "rationale": "The functional-evidence axis (functional-dependency features) and, "
                     "for AZ target-selection, plausibly the single most decision-relevant "
                     "signal: it answers whether a target is functionally IMPORTANT in a "
                     "line, not merely expressed. HIGHEST impact of the acquisitions.",
    },
    "pathway_genesets": {
        "status": "EXTERNAL_RESOURCE",
        "coverage_note": "Not a cell-line modality; external gene-set knowledge",
        "blocker_or_gate": "acquire MSigDB/Reactome at pathway-feature stage (D4/Q5)",
        "rationale": "Needed to infer pathway activity; sourced when that feature is built.",
    },
    "marker_panels": {
        "status": "EXTERNAL_RESOURCE",
        "coverage_note": "Not a cell-line modality; curated marker lists",
        "blocker_or_gate": "acquire/curate marker panels at marker-feature stage (D7)",
        "rationale": "Receptor/marker features require curated gene panels.",
    },
}

# Derived, so they stay in sync with MODALITY_SCOPE if you edit it:
#   ACQUISITION_TARGETS - not present yet, but planned for acquisition. Finding one of
#                         these PRESENT is GOOD NEWS (promote to IN_SCOPE), not a problem.
#   EXTERNAL_MODALITIES  - knowledge bases, never harmonised tables; not presence-checked.
ACQUISITION_TARGETS = tuple(
    k for k, v in MODALITY_SCOPE.items() if v["status"] == "PENDING_ACQUISITION"
)
EXTERNAL_MODALITIES = tuple(
    k for k, v in MODALITY_SCOPE.items() if v["status"] == "EXTERNAL_RESOURCE"
)


def print_scope_report(scope: dict[str, dict] = MODALITY_SCOPE) -> None:
    """Print the modality scope registry grouped by status.

    Use this at the top of any feature-engineering notebook so the current scope
    (and every descoping decision) is visible and auditable before work begins.

    Args:
        scope: The MODALITY_SCOPE registry (or a copy you have edited).

    Returns:
        None. Prints to stdout.

    Resolves:
        Planning-doc L.4 / Section K - the explicit, traceable descope decision.
    """
    order = ["IN_SCOPE", "IN_SCOPE_PENDING", "PENDING_ACQUISITION", "BLOCKED",
             "EXTERNAL_RESOURCE", "DESCOPED"]
    print(f"\n=== MODALITY SCOPE ({PLANNING_DOC}) ===")
    for status in order:
        rows = {k: v for k, v in scope.items() if v["status"] == status}
        if not rows:
            continue
        print(f"\n[{status}]")
        for name, meta in rows.items():
            print(f"  - {name:18s} {meta['coverage_note']}")
            if meta.get("acquire"):
                print(f"       acquire: {meta['acquire']}")
            gate = meta.get("blocker_or_gate")
            if gate:
                print(f"       gate: {gate}")


def save_scope_registry(path: str = "verification_out/modality_scope.json",
                        scope: dict[str, dict] = MODALITY_SCOPE) -> str:
    """Persist the scope registry as JSON for the audit trail / dissertation appendix.

    Args:
        path: Output JSON path. Parent directory is created if needed.
        scope: The registry to save.

    Returns:
        The path written.
    """
    import os
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    payload = {
        "planning_doc": PLANNING_DOC,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": scope,
    }
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2)
    print(f"[saved] scope registry -> {path}")
    return path


# =============================================================================
# SECTION 1 - SHARED HELPERS
# =============================================================================

def _numeric_values(df: pd.DataFrame,
                    exclude_cols: Sequence[str] = (),
                    max_cells: int = 5_000_000,
                    random_state: int = 0) -> np.ndarray:
    """Flatten a DataFrame's numeric columns into a 1-D array for scale analysis.

    Large matrices (e.g. HPA long ~24M rows, DepMap ~77M cells) are down-sampled to
    `max_cells` values so the checks stay memory-safe. Sampling is only for computing
    summary statistics; it never touches your stored data.

    Args:
        df: Any DataFrame.
        exclude_cols: Column names to drop before selecting numeric columns
            (typically identifier columns such as 'model_id', 'ccle_id').
        max_cells: Maximum number of numeric values to analyse; excess is sampled.
        random_state: Seed for reproducible sampling.

    Returns:
        1-D numpy array of finite numeric values (NaNs removed).

    Raises:
        ValueError: If no numeric columns remain after exclusion.
    """
    cols = [c for c in df.columns if c not in set(exclude_cols)]
    num = df[cols].select_dtypes(include=[np.number])
    if num.shape[1] == 0:
        raise ValueError("No numeric columns found after exclusion - check column names.")
    vals = num.to_numpy(dtype="float64", copy=False).ravel()
    vals = vals[np.isfinite(vals)]
    if vals.size > max_cells:
        rng = np.random.default_rng(random_state)
        vals = rng.choice(vals, size=max_cells, replace=False)
    return vals


def _skew(vals: np.ndarray) -> float:
    """Return skewness, using scipy if available and a numpy fallback otherwise."""
    if vals.size < 3:
        return float("nan")
    if _HAVE_SCIPY:
        return float(_scipy_skew(vals, bias=False))
    m = vals.mean()
    s = vals.std(ddof=0)
    return float(np.nan) if s == 0 else float(((vals - m) ** 3).mean() / s ** 3)


def _result(check: str, planning_ref: str, status: str, verdict: str, **details) -> dict:
    """Assemble the standard result dict and print a one-line verdict."""
    icon = {"PASS": "PASS ", "FLAG": "FLAG ", "INFO": "INFO ", "ERROR": "ERROR"}.get(status, status)
    print(f"[{icon}] {check}: {verdict}")
    return {"check": check, "planning_ref": planning_ref, "status": status,
            "verdict": verdict, "details": details}


def _find_col(df: pd.DataFrame, candidates: Iterable[str]) -> str | None:
    """Return the first candidate column present in `df` (case-insensitive), else None."""
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    return None


def _classify_scale(vals: np.ndarray) -> tuple[str, dict]:
    """Heuristically classify a numeric distribution as log / linear / log-ratio.

    Heuristic (NOT a proof - confirm against source docs):
        - >1% negative values           -> "LOG_RATIO_OR_CENTERED"  (e.g. TMT proteomics)
        - non-negative and max <= ~30    -> "LIKELY_LOG"             (e.g. log2(TPM+1))
        - max > ~100                     -> "LIKELY_LINEAR"          (e.g. raw TPM / counts)
        - otherwise                      -> "AMBIGUOUS_INSPECT"

    Returns:
        (label, stats) where stats holds the numbers behind the label.
    """
    stats = {
        "n": int(vals.size),
        "min": float(np.min(vals)) if vals.size else float("nan"),
        "max": float(np.max(vals)) if vals.size else float("nan"),
        "mean": float(np.mean(vals)) if vals.size else float("nan"),
        "median": float(np.median(vals)) if vals.size else float("nan"),
        "q99": float(np.quantile(vals, 0.99)) if vals.size else float("nan"),
        "pct_zero": float(np.mean(vals == 0) * 100) if vals.size else float("nan"),
        "pct_negative": float(np.mean(vals < 0) * 100) if vals.size else float("nan"),
        "skew": _skew(vals),
    }
    if stats["pct_negative"] > 1.0:
        label = "LOG_RATIO_OR_CENTERED"
    elif stats["min"] >= 0 and stats["max"] <= 30:
        label = "LIKELY_LOG"
    elif stats["max"] > 100:
        label = "LIKELY_LINEAR"
    else:
        label = "AMBIGUOUS_INSPECT"
    return label, stats


# =============================================================================
# SECTION 2 - EXPRESSION-SCALE CHECKS   (implements L.2 / Q1 - the blocking gate)
# =============================================================================

def verify_expression_scale(df: pd.DataFrame,
                            name: str,
                            id_cols: Sequence[str] = ("model_id", "profileid",
                                                      "profile_id", "sample", "index"),
                            expected: str | None = None) -> dict:
    """Characterise a single expression matrix and label its measurement scale.

    Runs the scale heuristic over all numeric (gene) columns and reports whether the
    matrix looks log-transformed, linear, or log-ratio. This is the atom used by
    `compare_expression_scales` to confirm the cross-source mismatch that currently
    blocks every expression feature.

    Args:
        df: A wide expression matrix (lines/samples x genes) OR a long table whose
            numeric column is the expression value.
        name: Human label for reporting, e.g. "depmap_expr".
        id_cols: Identifier columns to exclude from the numeric scan.
        expected: Optional expected label ("LIKELY_LOG" / "LIKELY_LINEAR" /
            "LOG_RATIO_OR_CENTERED"). If given, status is PASS when it matches.

    Returns:
        Result dict; details include min/max/mean/median/q99/pct_zero/pct_negative/skew
        and the classified `scale_label`.

    Resolves:
        Planning-doc L.2, Q1, Appendix scale-reference row for this source.
    """
    try:
        vals = _numeric_values(df, exclude_cols=id_cols)
    except ValueError as exc:
        return _result("verify_expression_scale", "Q1", "ERROR",
                       f"{name}: {exc}", source=name)
    label, stats = _classify_scale(vals)
    stats["scale_label"] = label
    stats["source"] = name
    if expected is None:
        status, verdict = "INFO", f"{name} looks {label} (max={stats['max']:.2f})"
    elif label == expected:
        status, verdict = "PASS", f"{name} confirmed {label} (as expected)"
    else:
        status, verdict = "FLAG", f"{name} looks {label} but expected {expected}"
    return _result("verify_expression_scale", "Q1", status, verdict, **stats)


def compare_expression_scales(depmap_expr: pd.DataFrame | None = None,
                              geo_expr: pd.DataFrame | None = None,
                              hpa_rna: pd.DataFrame | None = None,
                              hpa_value_col: str = "ntpm") -> dict:
    """Confirm the DepMap-log / GEO-linear / HPA-linear scale mismatch (THE blocker).

    This is the single most important pre-feature check. `expr_long` concatenates
    DepMap (log2) and GEO (linear) values; if their scales differ, no expression
    threshold, similarity, or cross-source agreement feature is valid until the
    scales are reconciled. This function characterises each source and reports
    whether they are mutually comparable.

    Args:
        depmap_expr: DepMap expression matrix (expected LIKELY_LOG).
        geo_expr: GEO expression matrix (expected LIKELY_LINEAR).
        hpa_rna: HPA long RNA table; the numeric column named by `hpa_value_col`
            is analysed (expected LIKELY_LINEAR).
        hpa_value_col: Which HPA column to test ("tpm", "ptpm", or "ntpm").

    Returns:
        Result dict; details["per_source"] maps source -> scale_label, and
        details["comparable"] is True only if all present sources share a label.

    Resolves:
        Planning-doc C.1 (confirmed problem), L.2, Q1. A FLAG here is the expected,
        correct outcome given the notebook - it authorises the normalisation work.
    """
    per_source: dict[str, str] = {}
    sub: dict[str, dict] = {}

    if depmap_expr is not None:
        r = verify_expression_scale(depmap_expr, "depmap_expr", expected="LIKELY_LOG")
        per_source["depmap_expr"] = r["details"].get("scale_label", "ERROR")
        sub["depmap_expr"] = r["details"]
    if geo_expr is not None:
        r = verify_expression_scale(geo_expr, "geo_expr", expected="LIKELY_LINEAR")
        per_source["geo_expr"] = r["details"].get("scale_label", "ERROR")
        sub["geo_expr"] = r["details"]
    if hpa_rna is not None:
        col = _find_col(hpa_rna, [hpa_value_col, "ntpm", "ptpm", "tpm"])
        if col is None:
            per_source["hpa_rna"] = "ERROR"
        else:
            r = verify_expression_scale(hpa_rna[[col]], f"hpa_rna[{col}]",
                                        id_cols=(), expected="LIKELY_LINEAR")
            per_source["hpa_rna"] = r["details"].get("scale_label", "ERROR")
            sub["hpa_rna"] = r["details"]

    labels = {v for v in per_source.values() if v != "ERROR"}
    comparable = len(labels) <= 1
    if comparable and labels:
        status, verdict = "PASS", "expression sources share one scale - directly comparable"
    else:
        status = "FLAG"
        verdict = ("expression sources are on DIFFERENT scales "
                   f"({per_source}) - reconciliation (Q1) required before any expression feature")
    return _result("compare_expression_scales", "Q1/C.1", status, verdict,
                   per_source=per_source, comparable=comparable, stats=sub)


# =============================================================================
# SECTION 3 - OTHER SCALE / FIELD VERIFICATIONS   (implements L.3 / Section K)
# =============================================================================

def verify_proteomics_scale(proteomics: pd.DataFrame,
                            id_cols: Sequence[str] = ("model_id", "depmap_id")) -> dict:
    """Confirm proteomics is TMT log-ratio (centred, signed) rather than abundance.

    Log-ratio data (Nusinow et al. 2020) is expected to be roughly symmetric and
    centred near zero, with a meaningful fraction of negative values. This check
    verifies that signature so the transform/threshold choice for protein features
    is made on the correct assumption.

    Args:
        proteomics: Proteomics matrix (lines x UniProt proteins).
        id_cols: Identifier columns to exclude.

    Returns:
        Result dict; PASS if the log-ratio signature (negatives present, mean ~0)
        holds, FLAG otherwise.

    Resolves:
        Planning-doc K ("proteomics exact scale unconfirmed"), Q4, Appendix row.
    """
    try:
        vals = _numeric_values(proteomics, exclude_cols=id_cols)
    except ValueError as exc:
        return _result("verify_proteomics_scale", "Q4/K", "ERROR", str(exc))
    label, stats = _classify_scale(vals)
    centred = abs(stats["mean"]) < 1.0
    has_neg = stats["pct_negative"] > 1.0
    if has_neg and centred:
        status, verdict = "PASS", (f"proteomics confirmed log-ratio-like "
                                   f"(mean={stats['mean']:.3f}, {stats['pct_negative']:.1f}% negative)")
    else:
        status, verdict = "FLAG", (f"proteomics does NOT match log-ratio signature "
                                   f"(label={label}, mean={stats['mean']:.3f}) - verify source")
    return _result("verify_proteomics_scale", "Q4/K", status, verdict, **stats)


def verify_metabolomics_scale(metabolomics: pd.DataFrame,
                              id_cols: Sequence[str] = ("model_id", "ccle_id", "depmap_id")) -> dict:
    """Characterise the metabolomics scale (units/normalisation currently unknown).

    Metabolomics scale was marked Unknown in the plan. This check reports the label
    (log-like vs linear vs centred), zero/negative fractions, and skew so the
    normalisation decision for metabolite features rests on evidence.

    Args:
        metabolomics: Metabolomics matrix (lines x metabolites).
        id_cols: Identifier columns to exclude.

    Returns:
        Result dict (status INFO - this characterises rather than pass/fails).

    Resolves:
        Planning-doc K ("metabolomics scale unconfirmed"), Appendix row.
    """
    try:
        vals = _numeric_values(metabolomics, exclude_cols=id_cols)
    except ValueError as exc:
        return _result("verify_metabolomics_scale", "K", "ERROR", str(exc))
    label, stats = _classify_scale(vals)
    stats["scale_label"] = label
    verdict = (f"metabolomics looks {label} "
               f"(min={stats['min']:.2f}, max={stats['max']:.2f}, "
               f"{stats['pct_negative']:.1f}% neg, skew={stats['skew']:.2f}) - confirm vs CCLE docs")
    return _result("verify_metabolomics_scale", "K", "INFO", verdict, **stats)


def verify_hpa_variants(hpa_rna: pd.DataFrame,
                       variant_cols: Sequence[str] = ("tpm", "ptpm", "ntpm")) -> dict:
    """Compare HPA tpm / pTPM / nTPM columns to inform which to use downstream.

    HPA publishes three RNA columns; nTPM is HPA's within-dataset normalised value
    and is usually the right choice for cross-sample comparison. This check confirms
    which columns exist and summarises each, so the choice is recorded with evidence
    rather than assumed.

    Args:
        hpa_rna: HPA long RNA table.
        variant_cols: Candidate value-column names to look for.

    Returns:
        Result dict; details["variants"] maps each present column to its stats.

    Resolves:
        Planning-doc K ("HPA variant choice"), rna_hpa gate in MODALITY_SCOPE.
    """
    present = {c: _find_col(hpa_rna, [c]) for c in variant_cols}
    present = {k: v for k, v in present.items() if v is not None}
    if not present:
        return _result("verify_hpa_variants", "K", "ERROR",
                       f"none of {variant_cols} found in hpa_rna columns")
    variants = {}
    for key, col in present.items():
        vals = hpa_rna[col].to_numpy(dtype="float64")
        vals = vals[np.isfinite(vals)]
        _, stats = _classify_scale(vals)
        variants[key] = stats
    verdict = (f"found HPA variants {list(present)}; "
               "nTPM is HPA's cross-sample-normalised value (confirm on HPA methods page)")
    return _result("verify_hpa_variants", "K", "INFO", verdict, variants=variants)


def audit_geo_batch_fields(geo_info: pd.DataFrame,
                           candidate_fields: Sequence[str] = (
                               "platform", "platform_id", "gpl", "series", "series_id",
                               "gse", "study", "instrument", "library_strategy")) -> dict:
    """Check whether GEO metadata carries fields needed to diagnose batch effects.

    Cross-study/platform batch effects (Q14) can only be diagnosed and handled if the
    sample metadata records platform/series. This check reports which such fields are
    present and how many distinct values each has.

    Args:
        geo_info: GEO sample-metadata table.
        candidate_fields: Field names that could encode batch/platform/study.

    Returns:
        Result dict; PASS if at least one usable batch field is found, FLAG otherwise.

    Resolves:
        Planning-doc C ("GEO batch effects - suspected"), Q14.
    """
    found = {}
    for f in candidate_fields:
        col = _find_col(geo_info, [f])
        if col is not None:
            found[col] = int(geo_info[col].nunique(dropna=True))
    if found:
        status, verdict = "PASS", f"batch-relevant fields present: {found}"
    else:
        status, verdict = "FLAG", ("no platform/series/study field found in geo_info - "
                                   "batch diagnosis (Q14) needs an external join to GEO series metadata")
    return _result("audit_geo_batch_fields", "Q14", status, verdict, found=found,
                   available_columns=list(geo_info.columns))


def audit_lineage_vocabulary(sample_info: pd.DataFrame,
                             lineage_cols: Sequence[str] = (
                                 "lineage", "lineage_subtype", "subtype",
                                 "oncotree_lineage", "primary_disease")) -> dict:
    """Assess whether lineage/subtype fields need a controlled vocabulary.

    Free-text lineage labels often contain case/whitespace variants of the same value,
    which would fragment a categorical feature. This check reports cardinality and
    detects near-duplicates (values that collapse to the same normalised string).

    Args:
        sample_info: The authoritative cell-line roster.
        lineage_cols: Candidate lineage/subtype column names.

    Returns:
        Result dict; FLAG if near-duplicate labels are detected, INFO otherwise.

    Resolves:
        Planning-doc C ("lineage vocabulary - suspected"), D6, lineage gate.
    """
    report = {}
    any_dupe = False
    for c in lineage_cols:
        col = _find_col(sample_info, [c])
        if col is None:
            continue
        raw = sample_info[col].dropna().astype(str)
        norm = raw.str.strip().str.lower()
        n_raw, n_norm = raw.nunique(), norm.nunique()
        collapsed = n_raw - n_norm
        any_dupe = any_dupe or collapsed > 0
        report[col] = {"distinct_raw": int(n_raw), "distinct_normalised": int(n_norm),
                       "collapsible_variants": int(collapsed)}
    if not report:
        return _result("audit_lineage_vocabulary", "D6", "ERROR",
                       f"none of {lineage_cols} found in sample_info")
    status = "FLAG" if any_dupe else "INFO"
    verdict = ("case/whitespace variants detected - a controlled vocabulary is needed"
               if any_dupe else "lineage labels look internally consistent (still confirm mapping)")
    return _result("audit_lineage_vocabulary", "D6", status, verdict, fields=report)


def verify_mirna_identifiers(mirna: pd.DataFrame,
                             id_cols: Sequence[str] = ("model_id", "ccle_name")) -> dict:
    """Confirm miRNA feature IDs are opaque and unmappable to standard miRNA names.

    CCLE miRNA columns arrive as opaque 'nmir00001.1'-style IDs. Features cannot be
    tied to biology without a map to miRBase (hsa-miR-* / MIMAT*). This check measures
    how many columns already match a standard miRNA pattern; if essentially none do,
    miRNA stays BLOCKED and should be flipped to DESCOPED unless a map is obtained.

    Args:
        mirna: miRNA matrix (lines x miRNA features).
        id_cols: Non-feature identifier columns to exclude from the column scan.

    Returns:
        Result dict; FLAG if <5% of feature columns look like standard miRNA names.

    Resolves:
        Planning-doc K / Q13, mirna BLOCKED status in MODALITY_SCOPE.
    """
    feat_cols = [c for c in mirna.columns if c not in set(id_cols)]
    if not feat_cols:
        return _result("verify_mirna_identifiers", "Q13", "ERROR", "no feature columns found")
    import re
    std = re.compile(r"(hsa[-_]?mir|hsa[-_]?let|MIMAT\d)", re.IGNORECASE)
    n_std = sum(bool(std.search(str(c))) for c in feat_cols)
    frac_std = n_std / len(feat_cols)
    examples = [str(c) for c in feat_cols[:5]]
    if frac_std < 0.05:
        status, verdict = "FLAG", (f"{n_std}/{len(feat_cols)} miRNA IDs are standard "
                                   f"(examples {examples}) - mapping to miRBase required; keep DESCOPED/BLOCKED")
    else:
        status, verdict = "PASS", f"{frac_std:.0%} of miRNA IDs are standard - mapping feasible"
    return _result("verify_mirna_identifiers", "Q13", status, verdict,
                   n_features=len(feat_cols), n_standard=n_std, examples=examples)


def quantify_zero_inflation(df: pd.DataFrame, name: str,
                            id_cols: Sequence[str] = ("model_id", "ccle_id", "profileid")) -> dict:
    """Measure zero-fraction and skew of a modality (affects distance/threshold features).

    Heavy zero-inflation and skew distort any distance- or threshold-based feature.
    This check reports the numbers so feature design accounts for them.

    Args:
        df: A numeric modality matrix.
        name: Label for reporting.
        id_cols: Identifier columns to exclude.

    Returns:
        Result dict (INFO); details include pct_zero and skew.

    Resolves:
        Planning-doc C ("outliers / zero-inflation - suspected").
    """
    try:
        vals = _numeric_values(df, exclude_cols=id_cols)
    except ValueError as exc:
        return _result("quantify_zero_inflation", "C", "ERROR", f"{name}: {exc}")
    _, stats = _classify_scale(vals)
    verdict = f"{name}: {stats['pct_zero']:.1f}% zeros, skew={stats['skew']:.2f}"
    return _result("quantify_zero_inflation", "C", "INFO", verdict, source=name, **stats)


def quantify_identity_ambiguity(connection: pd.DataFrame,
                                model_col: str = "model_id",
                                rrid_col: str = "rrids",
                                measurement_tables: dict[str, pd.DataFrame] | None = None) -> dict:
    """Size the 'credit both' identity-ambiguity population (a confidence signal).

    The harmonisation deliberately credits both candidate lines when a source ID maps
    to more than one ACH. Any per-line count feature can therefore be inflated for those
    lines, so the ambiguous set must be flagged for the confidence layer. This check
    counts (a) roster lines sharing an RRID/CVCL, and (b) rows in supplied measurement
    tables whose model_id is list-valued or delimited.

    Args:
        connection: The `cell_line_connection` bridge table.
        model_col: model_id column name.
        rrid_col: Column holding the set/list of RRIDs per line.
        measurement_tables: Optional {name: df} of re-keyed tables to scan for
            list-valued / delimited model_id entries.

    Returns:
        Result dict; FLAG if any ambiguity is found, PASS if none.

    Resolves:
        Planning-doc C ("identity ambiguity - confirmed"), G.3 confidence flag.
    """
    details: dict = {}

    # (a) roster lines that share an RRID/CVCL with another line
    shared = None
    if rrid_col in connection.columns:
        def _explode(v):
            if isinstance(v, (list, tuple, set)):
                return list(v)
            if pd.isna(v):
                return []
            return [s.strip() for s in str(v).replace(";", ",").split(",") if s.strip()]
        long = connection[[model_col, rrid_col]].copy()
        long["_rrid"] = long[rrid_col].map(_explode)
        long = long.explode("_rrid").dropna(subset=["_rrid"])
        counts = long.groupby("_rrid")[model_col].nunique()
        shared = int((counts > 1).sum())
        details["rrids_shared_by_multiple_models"] = shared

    # (b) measurement rows with list-valued / delimited model_id
    ambiguous_rows = {}
    if measurement_tables:
        for name, tbl in measurement_tables.items():
            col = _find_col(tbl, [model_col])
            if col is None:
                continue
            s = tbl[col]
            is_listlike = s.map(lambda v: isinstance(v, (list, tuple, set))).sum()
            is_delim = s.map(lambda v: isinstance(v, str) and ("," in v or ";" in v)).sum()
            ambiguous_rows[name] = int(is_listlike + is_delim)
        details["ambiguous_measurement_rows"] = ambiguous_rows

    total_ambig = (shared or 0) + sum(ambiguous_rows.values())
    if total_ambig > 0:
        status, verdict = "FLAG", (f"{total_ambig} ambiguity signals found - "
                                   "flag these lines/rows for the confidence layer")
    else:
        status, verdict = "PASS", "no identity ambiguity detected in supplied tables"
    return _result("quantify_identity_ambiguity", "C/G.3", status, verdict, **details)


# =============================================================================
# SECTION 4 - SCOPE CONFIRMATION & CONTRACT FREEZE   (implements L.1 / L.4)
# =============================================================================

def check_acquisition_status(available_table_names: Iterable[str],
                             targets: Sequence[str] = ACQUISITION_TARGETS,
                             aliases: dict[str, Sequence[str]] | None = None) -> dict:
    """Report whether the PENDING_ACQUISITION modalities have been acquired yet.

    Copy-number and CRISPR dependency are not descoped - they are high-priority
    acquisitions (standard DepMap files, same ACH key). This check looks for them in
    the tables you have loaded and reports acquisition progress. Unlike a descope
    audit, finding one PRESENT is the desired outcome: it means you can promote that
    modality to IN_SCOPE. Anything still absent is simply not-yet-acquired.

    Args:
        available_table_names: Names of the tables you actually loaded (e.g. the keys
            of your `tables` dict from the harmonisation pipeline).
        targets: Modalities with status PENDING_ACQUISITION (defaults to the registry).
        aliases: Optional {modality: [substrings]} to recognise a modality by table
            name (defaults cover common DepMap naming, e.g. copy_number -> OmicsCNGene).

    Returns:
        Result dict; status INFO. details["acquired"] lists targets now present (promote
        them to IN_SCOPE), details["still_needed"] lists targets to acquire before their
        dependent features (D5 for copy_number, F9 for crispr_dependency) can be built.

    Resolves:
        Planning-doc K, L.4 - tracks acquisition rather than validating a descope.
    """
    default_aliases = {
        "copy_number": ["cnv", "copy_number", "copynumber", "_cn", "gene_cn", "omicscn", "cngene"],
        "crispr_dependency": ["crispr", "chronos", "ceres", "dependency", "gene_effect", "geneeffect"],
    }
    aliases = {**default_aliases, **(aliases or {})}
    names_lower = [str(n).lower() for n in available_table_names]

    acquired, still_needed = {}, []
    for mod in targets:
        subs = aliases.get(mod, [mod])
        hits = [n for n in names_lower if any(s in n for s in subs)]
        if hits:
            acquired[mod] = hits
        else:
            still_needed.append(mod)

    if acquired and not still_needed:
        verdict = f"all acquisition targets present: {list(acquired)} - promote to IN_SCOPE"
    elif acquired:
        verdict = (f"acquired: {list(acquired)} (promote to IN_SCOPE); "
                   f"still to acquire: {still_needed}")
    else:
        verdict = (f"not yet acquired: {still_needed} - acquire before their dependent "
                   "features (D5 needs copy_number; functional-dependency features "
                   "(category 16) need crispr_dependency)")
    return _result("check_acquisition_status", "K/L.4", "INFO", verdict,
                   acquired=acquired, still_needed=still_needed,
                   external_not_checked=list(EXTERNAL_MODALITIES))


def freeze_coverage_contract(coverage: pd.DataFrame,
                             model_col: str = "model_id",
                             out_path: str = "verification_out/coverage_contract.parquet",
                             ambiguity_flags: pd.DataFrame | None = None) -> dict:
    """Persist the per-line modality-coverage matrix as the canonical confidence substrate.

    Freezes the boolean/tallied coverage per model_id (plus optional identity-ambiguity
    flags) to a versioned parquet file, and records a content hash so the exact substrate
    used to build confidence features is reproducible and citable in the dissertation.

    Args:
        coverage: One row per model_id with per-modality presence (bool or count) columns.
        model_col: The model_id column.
        out_path: Destination parquet path (parent dir created if needed).
        ambiguity_flags: Optional per-line flags (from quantify_identity_ambiguity) to
            join in as confidence inputs.

    Returns:
        Result dict; details include row/column counts, the SHA-256 content hash, and
        the written path. Status PASS on success, ERROR on failure.

    Resolves:
        Planning-doc L.1 (freeze the identity/coverage contract), G.3 (confidence inputs).
    """
    import os
    if model_col not in coverage.columns:
        return _result("freeze_coverage_contract", "L.1", "ERROR",
                       f"'{model_col}' not in coverage columns")
    frozen = coverage.copy()
    if ambiguity_flags is not None and model_col in ambiguity_flags.columns:
        frozen = frozen.merge(ambiguity_flags, on=model_col, how="left")

    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    try:
        frozen.to_parquet(out_path, index=False)
    except Exception as exc:  # pragma: no cover - environment dependent
        return _result("freeze_coverage_contract", "L.1", "ERROR",
                       f"failed to write parquet: {exc}")

    content_hash = hashlib.sha256(
        pd.util.hash_pandas_object(frozen, index=True).values.tobytes()
    ).hexdigest()[:16]
    meta = {
        "planning_doc": PLANNING_DOC,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "rows": int(frozen.shape[0]),
        "cols": int(frozen.shape[1]),
        "columns": list(frozen.columns),
        "content_sha256_16": content_hash,
        "path": out_path,
    }
    with open(out_path.replace(".parquet", ".meta.json"), "w") as fh:
        json.dump(meta, fh, indent=2)
    return _result("freeze_coverage_contract", "L.1", "PASS",
                   f"coverage contract frozen: {meta['rows']} lines x {meta['cols']} cols "
                   f"(hash {content_hash})", **meta)


# =============================================================================
# SECTION 5 - RUNNER
# =============================================================================

def run_all(tables: dict[str, pd.DataFrame]) -> list[dict]:
    """Run every applicable verification in dependency order and return the audit log.

    Only runs checks for tables that are present in `tables`, so it degrades gracefully
    if you have not loaded everything. Prints verdicts as it goes and returns the list of
    result dicts for saving to the audit log.

    Args:
        tables: Dict mapping canonical names to harmonised DataFrames. Recognised keys:
            'depmap_expr', 'geo_expr', 'hpa_rna', 'proteomics', 'metabolomics', 'mirna',
            'sample_info', 'geo_info', 'signatures', 'cell_line_connection', 'coverage'.

    Returns:
        List of result dicts (the machine-readable audit log).

    Resolves:
        Orchestrates Planning-doc L.1-L.4 in the recommended order.
    """
    log: list[dict] = []

    # Step 1 - scope: track acquisitions, then print/save the scope decisions
    print("\n--- Step 1: scope confirmation, acquisitions & descope ---")
    log.append(check_acquisition_status(available_table_names=tables.keys()))
    print_scope_report()

    # Step 2 - THE blocking gate: expression scales
    print("\n--- Step 2: expression-scale reconciliation (Q1, blocking) ---")
    log.append(compare_expression_scales(
        depmap_expr=tables.get("depmap_expr"),
        geo_expr=tables.get("geo_expr"),
        hpa_rna=tables.get("hpa_rna"),
    ))

    # Step 3 - other uncertain scales
    print("\n--- Step 3: verify remaining uncertain scales (Section K) ---")
    if "proteomics" in tables:
        log.append(verify_proteomics_scale(tables["proteomics"]))
    if "metabolomics" in tables:
        log.append(verify_metabolomics_scale(tables["metabolomics"]))
    if "hpa_rna" in tables:
        log.append(verify_hpa_variants(tables["hpa_rna"]))

    # Step 4 - suspected quality issues
    print("\n--- Step 4: confirm suspected quality issues (Section C) ---")
    if "geo_info" in tables:
        log.append(audit_geo_batch_fields(tables["geo_info"]))
    if "sample_info" in tables:
        log.append(audit_lineage_vocabulary(tables["sample_info"]))
    if "mirna" in tables:
        log.append(verify_mirna_identifiers(tables["mirna"]))
    for nm in ("depmap_expr", "geo_expr", "mirna", "metabolomics"):
        if nm in tables:
            log.append(quantify_zero_inflation(tables[nm], nm))
    if "cell_line_connection" in tables:
        log.append(quantify_identity_ambiguity(
            tables["cell_line_connection"],
            measurement_tables={k: v for k, v in tables.items()
                                if k in ("geo_expr", "hpa_rna")},
        ))

    # Step 5 - freeze the contract
    print("\n--- Step 5: freeze identity/coverage contract (L.1) ---")
    if "coverage" in tables:
        log.append(freeze_coverage_contract(tables["coverage"]))
    else:
        print("[INFO ] freeze_coverage_contract skipped: no 'coverage' table supplied")

    return log


def save_audit_log(log: list[dict], path: str = "verification_out/verification_log.json") -> str:
    """Write the run_all() results to a JSON audit log for the dissertation trail."""
    import os
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    payload = {
        "planning_doc": PLANNING_DOC,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_checks": len(log),
        "summary": {s: sum(1 for r in log if r["status"] == s)
                    for s in ("PASS", "FLAG", "INFO", "ERROR")},
        "results": log,
    }
    with open(path, "w") as fh:
        json.dump(payload, fh, indent=2, default=str)
    print(f"\n[saved] audit log -> {path}  ({payload['summary']})")
    return path


if __name__ == "__main__":
    # ------------------------------------------------------------------ #
    # FILL IN: load your harmonised tables here (paths are local to you). #
    # The keys must match those recognised by run_all() above.            #
    # ------------------------------------------------------------------ #
    #
    # import pandas as pd
    # BASE = "path/to/harmonised/"
    # tables = {
    #     "depmap_expr":          pd.read_parquet(BASE + "depmap_expr.parquet"),
    #     "geo_expr":             pd.read_parquet(BASE + "geo_expr.parquet"),
    #     "hpa_rna":              pd.read_parquet(BASE + "hpa_rna.parquet"),
    #     "proteomics":           pd.read_parquet(BASE + "proteomics.parquet"),
    #     "metabolomics":         pd.read_parquet(BASE + "metabolomics.parquet"),
    #     "mirna":                pd.read_parquet(BASE + "mirna.parquet"),
    #     "sample_info":          pd.read_parquet(BASE + "sample_info.parquet"),
    #     "geo_info":             pd.read_parquet(BASE + "geo_info.parquet"),
    #     "signatures":           pd.read_parquet(BASE + "signatures.parquet"),
    #     "cell_line_connection": pd.read_parquet(BASE + "cell_line_connection.parquet"),
    #     "coverage":             pd.read_parquet(BASE + "coverage_matrix.parquet"),
    # }
    #
    # log = run_all(tables)
    # save_audit_log(log)
    # save_scope_registry()
    #
    print(__doc__)
    print("Fill in the loader in the __main__ block, then run again to execute the checks.")