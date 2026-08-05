"""
08_confidence_engine.py  —  CELL 2: SCORING (rule-based, GRADE-shaped, no weighted sum)
CellLineSelector - Confidence Engine. Paste below CELL 1 (characterisation), or run standalone.

WHAT THIS IS
------------
Turns the confidence inputs in features/expression_features.parquet into a per-(line, gene)
ordinal confidence judgement, using a transparent rule table grounded in the CELL 1
characterisation of the PROTEIN-CODING scope (E6 default universe). The design is GRADE-shaped:
a base tier set by corroboration, downgraded a level for measurement uncertainty, with a
distinct Conflicting state and a hard Low floor. Evidence is TALLIED, never averaged - there
is no weighted sum, consistent with the evidence-assessment framing (dimensions measure
different things and are not interchangeable units).

WHY THESE RULES (from the protein-coding characterisation)
----------------------------------------------------------
  - Corroboration is the majority case (75.7%); detection-corroborated (HPA call present) 72.7%.
    -> agreement is a real, co-primary signal here (unlike the full set, where it was trivial).
  - Of corroborated rows, 86.4% concordant, 13.6% discordant -> ~9.9% of all rows are
    detection-discordant -> a genuine Conflicting band, flagged not buried.
  - Biotype is 100% constant under this scope -> NOT a scoring dimension; it is a precondition.
  - Percentile separates detected (median 0.64) from not-detected (median 0.20; q99 0.33)
    -> percentile proxies 'distance from detection threshold'; overlap zone ~0.15-0.35.
  - GEO makes no present/absent call (E3): it counts toward coverage/rank-corroboration only.

WHAT THIS IS *NOT*
------------------
No weighted sum; no ML; no calibration claim. Confidence is documented-not-validated: it is
checked by the built-in face-validity battery (known archetypes land in the expected bands),
not against ground truth. Bands are the reported unit; the evidence tally is an ordering aid
for within-band ranking, not a probability.

USAGE (live Jupyter session)
----------------------------
    import importlib
    ce2 = importlib.import_module("08_confidence_engine")   # if pasted into the same file
    conf = ce2.run_scoring(expr=expression_features)        # or omit expr to load the parquet
    # -> features/confidence/expression_confidence.parquet  (+ summary JSON)
    print(ce2.format_rationale(conf.iloc[0]))               # human-readable explanation

Author: (CellLineSelector team)   |   Depends on: pandas, numpy
"""

from __future__ import annotations

from ast import expr
import json
import os
from typing import Iterable
from altair import expr

import numpy as np
import pandas as pd


# =============================================================================
# CONFIG  (every constant is a documented, tunable design choice - not literature)
# =============================================================================

EXPRESSION_FEATURES_PATH = "features/expression_features.parquet"
OUTDIR                   = "features/confidence"

# Reliability cutpoints on expr_percentile_depmap (protein-coding scope, from CELL 1):
#   not-detected q99 = 0.328 ; detected q01 = 0.255, q25 = 0.454.
# The [LO, HI] band is the overlap/uncertain zone where a call is near the threshold.
NEAR_THRESHOLD_LO = 0.15
NEAR_THRESHOLD_HI = 0.35

# THE policy decision (see CELL 1: ~63% of protein-coding rows stay High-eligible):
#   True  -> High requires an independent, concordant HPA detection call (honest/strict).
#   False -> a confident DepMap-only call (deep in distribution) may also reach High.
HIGH_REQUIRES_HPA_CONCORDANCE = True

# Ordinal ladder. Conflicting is a separate state, not a rung (ClinGen precedent).
_TIER_TO_BAND = {3: "High", 2: "Moderate", 1: "Low"}
_BANDS = ("High", "Moderate", "Low", "Conflicting")


# =============================================================================
# HELPERS (self-contained so this cell runs alone)
# =============================================================================

def _find_col(df: pd.DataFrame, candidates: Iterable[str], required: bool = True,
              what: str = "column") -> str | None:
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    if required:
        raise KeyError(f"Could not resolve {what}: tried {list(candidates)}; "
                       f"available = {list(df.columns)[:40]}")
    return None


def _as_bool_num(s: pd.Series) -> pd.Series:
    """Coerce bool-ish (True/False, 1.0/0.0, '1'/'true', NaN) to float 1.0/0.0/NaN.

    NaN = not measured. Numeric floats (e.g. 1.0/0.0 after a parquet round-trip) MUST be
    handled here - if they are dropped, the whole agreement dimension silently vanishes.
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


def _load(df, path: str, what: str):
    if df is not None:
        return df
    if path and os.path.exists(path):
        return pd.read_parquet(path)
    raise FileNotFoundError(f"{what} not provided and not found at {path}")


def _outdir(path: str | None) -> str | None:
    if not path:
        return None
    os.makedirs(path, exist_ok=True)
    return path


# =============================================================================
# THE RUBRIC  (vectorised; base tier -> downgrade -> band; Conflicting routed out)
# =============================================================================

def score(expr: pd.DataFrame) -> pd.DataFrame:
    """Return a per-(line, gene) confidence frame with band, tally, assessments.

    Rule table (GRADE-shaped, evidence tallied not averaged):
      Conflicting : DepMap and HPA both call, and disagree.
      Base tier   : concordant DepMap+HPA -> 3 ; any other valid call -> 2.
      Downgrade -1: near-threshold measurement; identity-ambiguous line; method-fallback line.
      High gate   : if HIGH_REQUIRES_HPA_CONCORDANCE, tier 3 needs concordant HPA (default).
      Floor       : Low (tier cannot fall below 1).
    """
    n = len(expr)
    mcol = _find_col(expr, ("model_id",))
    gcol = _find_col(expr, ("gene", "gene_id"))
    c_dm  = _find_col(expr, ("detected_depmap",))
    c_hpa = _find_col(expr, ("detected_hpa",), required=False)
    # c_geo = _find_col(expr, ("detected_geo", "geo_intensity"), required=False)
    # was: c_geo = _find_col(expr, ("geo_intensity", "detected_geo"), required=False)
    c_geo = _find_col(expr, ("geo_intensity",), required=False)   # presence = value present; detected_geo is empty
    c_pct = _find_col(expr, ("expr_percentile_depmap",))

    dm  = _as_bool_num(expr[c_dm])
    hpa = _as_bool_num(expr[c_hpa]) if c_hpa else pd.Series(np.nan, index=expr.index)
    hpa_present = hpa.notna()
    geo_present = (expr[c_geo].notna() if c_geo else pd.Series(False, index=expr.index))
    pct = pd.to_numeric(expr[c_pct], errors="coerce")

    # --- evidence descriptors (independent signals) ---
    concordant = hpa_present & (dm == hpa)                 # both call, agree
    discordant = hpa_present & (dm != hpa)                 # both call, disagree -> Conflicting
    called_present = (dm == 1.0)                           # DepMap present/absent call
    # reliability: near-threshold when percentile sits in the DepMap/HPA overlap zone
    near_threshold = (pct >= NEAR_THRESHOLD_LO) & (pct <= NEAR_THRESHOLD_HI)

    # --- coverage assessment (categorical, for the rationale) ---
    cov = np.select(
        [hpa_present & geo_present, hpa_present & ~geo_present,
         ~hpa_present & geo_present],
        ["depmap_geo_hpa", "depmap_hpa", "depmap_geo"],
        default="depmap_only")

    agr = np.select([discordant, concordant],
                    ["discordant", "concordant"],
                    default="single_call")   # ~hpa_present -> no second call
    rel = np.where(near_threshold, "near_threshold", "confident")

    # --- residual provenance downgrades (rare; resolved defensively if columns absent) ---
    amb_col = _find_col(expr, ("identity_ambiguous", "model_id_ambiguous"), required=False)
    ambiguous = _as_bool_num(expr[amb_col]).fillna(0).astype(bool) if amb_col \
        else pd.Series(False, index=expr.index)
    fb_col = _find_col(expr, ("zfpkm_fallback", "tpm_fallback", "method_fallback"), required=False)
    fallback = _as_bool_num(expr[fb_col]).fillna(0).astype(bool) if fb_col \
        else pd.Series(False, index=expr.index)

    # --- base tier ---
    base = np.where(concordant, 3, 2).astype("int64")
    if not HIGH_REQUIRES_HPA_CONCORDANCE:
        # allow a confident DepMap-only call, deep in distribution, to also start at 3
        deep_single = (~hpa_present) & called_present & (pct >= NEAR_THRESHOLD_HI)
        base = np.where(deep_single, np.maximum(base, 3), base)

    # --- downgrades (one level each; GRADE-style, then floor at 1) ---
    downgrade = near_threshold.astype(int) + ambiguous.astype(int) + fallback.astype(int)
    tier = np.clip(base - downgrade, 1, 3)

    band = np.where(discordant, "Conflicting",
                    np.vectorize(_TIER_TO_BAND.get)(tier))

    # --- evidence tally (a transparent COUNT of positive signals, not a weighted sum) ---
    tally = (
        1                                                   # DepMap present (baseline)
        + concordant.astype(int)                            # +1 concordant HPA
        + geo_present.astype(int)                            # +1 GEO rank-corroboration
        + (~near_threshold).astype(int)                     # +1 confident (deep)
        - discordant.astype(int)                            # -1 discordant
    ).astype("int64")

    out = pd.DataFrame({
        mcol: expr[mcol].values,
        gcol: expr[gcol].values,
        "confidence_band": pd.Categorical(band, categories=_BANDS),
        "evidence_tally": tally,
        "cov_assessment": cov,
        "agr_assessment": agr,
        "rel_assessment": rel,
        "expr_percentile_depmap": pct.round(4).values,
    })
    return out


# =============================================================================
# EXPLAINABILITY  (rationale generated on demand, not materialised over 28M rows)
# =============================================================================

_COV_TXT = {"depmap_geo_hpa": "3 sources (DepMap, GEO, HPA)",
            "depmap_hpa": "2 sources (DepMap, HPA)",
            "depmap_geo": "2 sources (DepMap, GEO; GEO corroborates level only)",
            "depmap_only": "DepMap only"}
_AGR_TXT = {"concordant": "DepMap and HPA agree on the call",
            "discordant": "DepMap and HPA disagree on the call",
            "single_call": "no second present/absent call"}
_REL_TXT = {"confident": "deep in the expression distribution",
            "near_threshold": "near the detection threshold"}


def format_rationale(row: pd.Series) -> str:
    """Human-readable, fully traceable explanation for one (line, gene) row."""
    band = row["confidence_band"]
    pct = row.get("expr_percentile_depmap")
    parts = [_COV_TXT.get(row["cov_assessment"], row["cov_assessment"]),
             _AGR_TXT.get(row["agr_assessment"], row["agr_assessment"]),
             f"{_REL_TXT.get(row['rel_assessment'], row['rel_assessment'])} (percentile {pct})"]
    return f"{band} - " + "; ".join(parts) + "."


# =============================================================================
# FACE-VALIDITY BATTERY  (documented-not-validated: check archetypes, not accuracy)
# =============================================================================

def _face_validity(conf: pd.DataFrame) -> dict:
    checks, ok = [], True
    def _band_of(cov, agr, rel):
        sub = conf[(conf.cov_assessment == cov) & (conf.agr_assessment == agr)
                   & (conf.rel_assessment == rel)]
        if not len(sub):
            return None
        return sub["confidence_band"].value_counts().idxmax()

    expectations = [
        ("depmap_geo_hpa", "concordant", "confident", "High"),
        ("depmap_hpa", "concordant", "near_threshold", "Moderate"),
        ("depmap_only", "single_call", "confident", "Moderate"),
        ("depmap_only", "single_call", "near_threshold", "Low"),
    ]
    for cov, agr, rel, want in expectations:
        got = _band_of(cov, agr, rel)
        passed = (got is None) or (got == want)
        ok = ok and passed
        checks.append({"case": f"{cov}/{agr}/{rel}", "expected": want,
                       "got": str(got), "pass": bool(passed)})
    # discordant must always be Conflicting
    disc = conf[conf.agr_assessment == "discordant"]
    disc_ok = (not len(disc)) or bool((disc["confidence_band"] == "Conflicting").all())
    ok = ok and disc_ok
    checks.append({"case": "any discordant", "expected": "Conflicting",
                   "got": "Conflicting" if disc_ok else "MIXED", "pass": disc_ok})
    return {"all_pass": bool(ok), "checks": checks}


# =============================================================================
# ORCHESTRATOR
# =============================================================================

def run_scoring(expr=None, expr_path: str = EXPRESSION_FEATURES_PATH,
                outdir: str = OUTDIR, write: bool = True) -> pd.DataFrame:
    print("=" * 82)
    print("CONFIDENCE SCORING (rule-based, GRADE-shaped) - protein-coding scope (E6)")
    print(f"HIGH_REQUIRES_HPA_CONCORDANCE = {HIGH_REQUIRES_HPA_CONCORDANCE}   "
          f"near-threshold zone = [{NEAR_THRESHOLD_LO}, {NEAR_THRESHOLD_HI}]")
    print("=" * 82)

    expr = _load(expr, expr_path, "expression_features")
    conf = score(expr)

    # band distribution
    dist = conf["confidence_band"].value_counts(dropna=False)
    n = len(conf)
    print("\n[bands] distribution:")
    for b in _BANDS:
        c = int(dist.get(b, 0))
        print(f"[bands]   {b:12} {c:>12,}  ({c / n * 100:5.2f}%)")

    # band x evidence structure (a sanity cross-tab)
    xt = pd.crosstab(conf["cov_assessment"], conf["confidence_band"])
    print("\n[bands] coverage x band:\n" + xt.to_string())

    # face validity
    fv = _face_validity(conf)
    print(f"\n[face-validity] all_pass = {fv['all_pass']}")
    for c in fv["checks"]:
        print(f"[face-validity]   {c['case']:>42}  expect {c['expected']:<11} "
              f"got {c['got']:<11} {'OK' if c['pass'] else 'FAIL'}")

    # example rationales (sampled, one per band where available)
    print("\n[examples] one rationale per band:")
    for b in _BANDS:
        sub = conf[conf["confidence_band"] == b]
        if len(sub):
            print(f"[examples]   {format_rationale(sub.iloc[0])}")

    summary = {
        "scope": "protein_coding_default_universe_E6",
        "n_rows": int(n),
        "policy": {"HIGH_REQUIRES_HPA_CONCORDANCE": HIGH_REQUIRES_HPA_CONCORDANCE,
                   "near_threshold_zone": [NEAR_THRESHOLD_LO, NEAR_THRESHOLD_HI]},
        "band_distribution": {b: int(dist.get(b, 0)) for b in _BANDS},
        "band_distribution_pct": {b: round(int(dist.get(b, 0)) / n * 100, 3) for b in _BANDS},
        "coverage_x_band": xt.to_dict(),
        "face_validity": fv,
    }

    if write:
        od = _outdir(outdir)
        if od:
            conf.to_parquet(os.path.join(od, "expression_confidence.parquet"), index=False)
            with open(os.path.join(od, "confidence_scoring_summary.json"), "w") as f:
                json.dump(summary, f, indent=2, default=str)
            print(f"\n[run] wrote {od}/expression_confidence.parquet (+ summary JSON)")

    print("\n" + "=" * 82)
    print("Bands are the reported unit; evidence_tally orders within a band (not a probability).")
    print("Confidence is documented-not-validated: face validity only, no accuracy claim.")
    print("=" * 82)
    return conf


if __name__ == "__main__":
    print(__doc__)
    print("Import and call run_scoring(expr=expression_features) from your live session.")
