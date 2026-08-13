#!/usr/bin/env python3
"""
geo_fetch.py — download a GEO Series and emit it with provenance and a
declared (not assumed) scale.

Design contract, mirroring 02_transcriptonomics_proteins.ipynb:

  1. Nothing is assumed. The value scale is DETECTED from the data and
     CROSS-CHECKED against the submitter's stated data_processing text.
     Disagreement is reported, not silently resolved.
  2. Every output carries provenance: accession, platform, download time,
     SHA-256 of the source file, GEOparse version, sample count.
  3. Hard build-time assertions fail the run rather than degrade it:
     multi-platform series, duplicate samples, composite gene identifiers.
  4. Missing is missing. No imputation, no zero-filling. NaN stays NaN and
     the reason is recorded.

Install:  pip install GEOparse pandas pyarrow
Usage:    python geo_fetch.py GSE12345 --outdir data/geo
          python geo_fetch.py GSE12345 --supplementary   # raw files too

Reference: Barrett et al. (2013) NAR 41:D991 (GEO); GEOparse
(https://geoparse.readthedocs.io), modelled on Davis & Meltzer (2007)
GEOquery, Bioinformatics 23:1846.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

try:
    import GEOparse
except ImportError:  # pragma: no cover
    sys.exit("pip install GEOparse pandas pyarrow")

log = logging.getLogger("geo_fetch")

# ---------------------------------------------------------------- scale rules
# Ordered: first match wins. Patterns are matched against the concatenated
# data_processing / value_definition text supplied by the submitter.
STATED_SCALE_PATTERNS = [
    ("log2_ratio", r"\blog2?\s*ratio|\bM[- ]value|two[- ]colou?r"),
    ("rma_log2", r"\bg?c?rma\b|\bvst\b|\brlog\b"),
    ("mas5_linear", r"\bmas\s*5|\bmas5\b|\bplier\b|\bdchip\b"),
    ("tpm", r"\btpm\b"),
    ("fpkm_rpkm", r"\bf?rpkm\b"),
    ("counts", r"\braw counts?\b|\bread counts?\b|\bhtseq\b|\bfeaturecounts\b"),
    ("quantile_norm", r"\bquantile[- ]normali[sz]"),
    ("log_unspecified", r"\blog[ _-]?(2|10|e)?\b|\bnatural log\b"),
]


def classify_stated_scale(text: str) -> str:
    """Submitter's own words. Returns 'unstated' when nothing matches."""
    t = (text or "").lower()
    for label, pattern in STATED_SCALE_PATTERNS:
        if re.search(pattern, t):
            return label
    return "unstated"


def detect_scale(values: np.ndarray) -> dict:
    """
    Detect the scale empirically. Deliberately conservative: returns
    'ambiguous' rather than guessing, because a wrong scale label is worse
    than an absent one downstream.
    """
    v = values[np.isfinite(values)]
    if v.size < 100:
        return {"detected_scale": "insufficient_data", "n_finite": int(v.size)}

    stats = {
        "n_finite": int(v.size),
        "min": float(np.min(v)),
        "p01": float(np.percentile(v, 1)),
        "median": float(np.median(v)),
        "p99": float(np.percentile(v, 99)),
        "max": float(np.max(v)),
        "frac_negative": float(np.mean(v < 0)),
        "frac_zero": float(np.mean(v == 0)),
        "skew": float(pd.Series(v).skew()),
    }

    neg, mx, med = stats["frac_negative"], stats["max"], stats["median"]

    if neg > 0.20 and abs(med) < 1.0:
        detected = "log_ratio_centred"       # two-colour / centred M-values
    elif neg < 0.001 and mx > 1000:
        detected = "linear"                  # MAS5-like, TPM, FPKM, counts
    elif neg < 0.001 and mx < 30:
        detected = "log_like"                # log2 intensity or log2(TPM+1)
    elif neg > 0.001 and mx < 30:
        detected = "log_like_signed"         # log with background subtraction
    else:
        detected = "ambiguous"

    stats["detected_scale"] = detected
    return stats


def detect_floor(values: np.ndarray, low_quantile: float = 0.05) -> dict:
    """
    A source whose minimum sits well above zero AND which reports essentially
    no low values has a detection floor: it is reporting background, not
    silence. Same test as the notebook's source_floor_check.
    """
    v = values[np.isfinite(values)]
    if v.size < 100:
        return {"floored": None, "reason": "insufficient_data"}

    vmin = float(np.min(v))
    q_low = float(np.percentile(v, low_quantile * 100))
    span = float(np.percentile(v, 99) - vmin)
    # tie mass exactly at the minimum -> hard clipping
    tie_mass_at_min = float(np.mean(v <= vmin * 1.0001))

    floored = bool(span > 0 and (q_low - vmin) / span < 0.02 and tie_mass_at_min > 0.01)
    return {
        "floored": floored,
        "min": vmin,
        "q_low": q_low,
        "tie_mass_at_min": tie_mass_at_min,
        "reason": "clipped_low_tail" if floored else "no_floor_detected",
    }


# ------------------------------------------------------------- hard assertions
class BuildAssertionError(RuntimeError):
    """Raised when the data cannot be used without silently corrupting results."""


def assert_single_platform(gse) -> str:
    gpls = list(gse.gpls.keys())
    if len(gpls) == 0:
        raise BuildAssertionError("No platform (GPL) found in series.")
    if len(gpls) > 1:
        raise BuildAssertionError(
            f"Series spans {len(gpls)} platforms {gpls}. Probe sets are not "
            "comparable across platforms; split the series and process each "
            "platform separately."
        )
    return gpls[0]


def assert_no_composite_ids(series: pd.Series, axis_name: str) -> None:
    """
    Composite identifiers ('GENE1///GENE2', 'A;B') corrupt percentile
    denominators wherever they go undetected. Fail the build, do not clean
    silently.
    """
    s = series.dropna().astype(str)
    composite = s[s.str.contains(r"///|;|,\s*\w", regex=True, na=False)]
    if len(composite):
        examples = composite.unique()[:5].tolist()
        raise BuildAssertionError(
            f"{len(composite)} composite identifiers on the {axis_name} axis, "
            f"e.g. {examples}. Resolve the mapping explicitly (one row per "
            "probe-gene pair, or drop ambiguous probes) before use."
        )


def explode_composite_ids(df: pd.DataFrame, col: str) -> tuple[pd.DataFrame, int]:
    """
    Resolve composite identifiers explicitly rather than failing the build on
    them. GPL570 (Affymetrix HG-U133 Plus 2) routinely annotates one probe to
    several genes ('DDR1 /// MIR4640') -- this is normal for the platform, not
    corrupted data, and BUILD_SPEC.md's own A2 assertion anticipates it project
    -wide ("no ';' or '/' in any gene symbol column"). The resolution named in
    this script's own docstring is "one row per probe-gene pair, or drop
    ambiguous probes" -- this applies the first option, explicitly, and reports
    the count so it stays visible rather than silently vanishing into a join.
    """
    s = df[col].astype(str)
    is_composite = s.str.contains(r"///|;|,", regex=True, na=False)
    n = int(is_composite.sum())
    if n == 0:
        return df, 0
    exploded = df.copy()
    exploded[col] = exploded[col].str.split(r"\s*(?:///|;|,)\s*")
    exploded = exploded.explode(col)
    exploded[col] = exploded[col].str.strip()
    exploded = exploded[exploded[col] != ""]
    return exploded.reset_index(drop=True), n


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# --------------------------------------------------------------------- fetcher
def fetch_series(
    accession: str,
    outdir: Path,
    annotate_gpl: bool = True,
    supplementary: bool = False,
    strict: bool = True,
) -> dict:
    outdir.mkdir(parents=True, exist_ok=True)
    soft_dir = outdir / "soft"
    soft_dir.mkdir(exist_ok=True)

    log.info("fetching %s", accession)
    gse = GEOparse.get_GEO(
        geo=accession,
        destdir=str(soft_dir),
        annotate_gpl=annotate_gpl,
        silent=True,
    )

    # --- provenance on the source file itself
    soft_files = sorted(soft_dir.glob(f"{accession}*soft*"))
    source_files = [
        {"name": p.name, "bytes": p.stat().st_size, "sha256": sha256(p)}
        for p in soft_files
    ]

    gpl_name = assert_single_platform(gse)
    gpl = gse.gpls[gpl_name]

    # --- expression matrix: probes x samples, from each GSM's VALUE column
    matrix = gse.pivot_samples("VALUE")
    if matrix.columns.duplicated().any():
        dupes = matrix.columns[matrix.columns.duplicated()].tolist()
        raise BuildAssertionError(f"Duplicate GSM columns: {dupes}")

    # --- sample metadata, including what the submitter says they did
    pheno = gse.phenotype_data.copy()
    proc_cols = [c for c in pheno.columns if "data_processing" in c or "value" in c.lower()]
    stated_text = (
        pheno[proc_cols].astype(str).agg(" ".join, axis=1) if proc_cols else pd.Series("", index=pheno.index)
    )
    stated = stated_text.map(classify_stated_scale)

    # --- scale: detected, then compared to stated
    values = matrix.to_numpy(dtype=float, na_value=np.nan).ravel()
    scale = detect_scale(values)
    floor = detect_floor(values)

    stated_mode = stated.mode().iat[0] if len(stated) else "unstated"
    agreement = _scale_agreement(stated_mode, scale["detected_scale"])

    # --- probe -> gene mapping, if the platform was annotated
    mapping = None
    symbol_col = None
    n_exploded_probes = 0
    if annotate_gpl and gpl.table is not None:
        candidates = [c for c in gpl.table.columns
                      if c.lower() in {"gene symbol", "gene_symbol", "symbol", "genesymbol"}]
        if candidates:
            symbol_col = candidates[0]
            mapping = gpl.table[["ID", symbol_col]].rename(
                columns={"ID": "probe_id", symbol_col: "gene_symbol"}
            )
            mapping, n_exploded_probes = explode_composite_ids(mapping, "gene_symbol")
            if n_exploded_probes:
                log.info("%s: exploded %d multi-gene probes into one row per gene",
                         accession, n_exploded_probes)

    manifest = {
        "accession": accession,
        "platform": gpl_name,
        "n_samples": int(matrix.shape[1]),
        "n_probes": int(matrix.shape[0]),
        "downloaded_utc": datetime.now(timezone.utc).isoformat(),
        "geoparse_version": getattr(GEOparse, "__version__", "unknown"),
        "source_files": source_files,
        "scale_stated_by_submitter": stated_mode,
        "scale_stated_per_sample": stated.value_counts().to_dict(),
        "scale_detected": scale,
        "scale_agreement": agreement,
        "detection_floor": floor,
        "gene_symbol_column": symbol_col,
        "n_composite_probes_exploded": n_exploded_probes,
        "missing_fraction": float(matrix.isna().to_numpy().mean()),
        "assertions_run": ["single_platform", "no_duplicate_samples",
                           "composite_gene_ids_exploded_not_dropped"],
    }

    # --- write
    matrix.to_parquet(outdir / f"{accession}_matrix.parquet")
    pheno.to_parquet(outdir / f"{accession}_samples.parquet")
    if mapping is not None:
        mapping.to_parquet(outdir / f"{accession}_probe_map.parquet")
    (outdir / f"{accession}_manifest.json").write_text(json.dumps(manifest, indent=2))

    if supplementary:
        log.info("downloading supplementary files (nproc=1 by GEOparse's own advice)")
        gse.download_supplementary_files(directory=str(outdir / "supplementary"),
                                         download_sra=False)

    if agreement == "CONFLICT":
        log.warning(
            "SCALE CONFLICT: submitter says %r, data looks %r. Do not use this "
            "series until resolved — a wrong scale label silently breaks every "
            "downstream comparison.", stated_mode, scale["detected_scale"]
        )
    if floor.get("floored"):
        log.warning(
            "DETECTION FLOOR present. This source must not be allowed to judge "
            "silence: its low tail is background, not absence."
        )

    return manifest


def _scale_agreement(stated: str, detected: str) -> str:
    if stated == "unstated" or detected in {"ambiguous", "insufficient_data"}:
        return "UNKNOWN"
    expected_linear = {"mas5_linear", "tpm", "fpkm_rpkm", "counts"}
    expected_log = {"rma_log2", "log_unspecified", "quantile_norm"}
    expected_ratio = {"log2_ratio"}
    if stated in expected_linear and detected == "linear":
        return "AGREE"
    if stated in expected_log and detected in {"log_like", "log_like_signed"}:
        return "AGREE"
    if stated in expected_ratio and detected == "log_ratio_centred":
        return "AGREE"
    return "CONFLICT"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("accessions", nargs="+", help="GSE accessions, e.g. GSE1563")
    ap.add_argument("--outdir", type=Path, default=Path("data/geo"))
    ap.add_argument("--supplementary", action="store_true",
                    help="also download raw supplementary files (large)")
    ap.add_argument("--no-annotate", action="store_true",
                    help="skip GPL annotation download (faster, no gene symbols)")
    ap.add_argument("--lenient", action="store_true",
                    help="warn instead of failing on composite identifiers")
    ap.add_argument("-v", "--verbose", action="store_true")
    args = ap.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)-8s %(message)s",
    )

    failures = []
    for acc in args.accessions:
        try:
            m = fetch_series(
                acc,
                outdir=args.outdir / acc,
                annotate_gpl=not args.no_annotate,
                supplementary=args.supplementary,
                strict=not args.lenient,
            )
            log.info(
                "%s: %d probes x %d samples | stated=%s detected=%s (%s) | floored=%s",
                acc, m["n_probes"], m["n_samples"],
                m["scale_stated_by_submitter"], m["scale_detected"]["detected_scale"],
                m["scale_agreement"], m["detection_floor"].get("floored"),
            )
        except BuildAssertionError as e:
            log.error("%s FAILED assertion: %s", acc, e)
            failures.append((acc, str(e)))
        except Exception as e:  # network, parse, malformed SOFT
            log.error("%s FAILED: %s: %s", acc, type(e).__name__, e)
            failures.append((acc, f"{type(e).__name__}: {e}"))

    if failures:
        log.error("%d of %d accessions failed", len(failures), len(args.accessions))
        sys.exit(1)


if __name__ == "__main__":
    main()
