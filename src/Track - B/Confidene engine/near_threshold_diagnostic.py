"""
near_threshold_diagnostic.py
CellLineSelector - Confidence Engine: where do NEAR_THRESHOLD_LO/HI come from?

Splits expr_percentile_depmap by the DepMap detection call and shows why [0.15, 0.35]
is the "near threshold" (measurement-uncertain) zone:

  * detected and not-detected genes OVERLAP in percentile - a call landing in the overlap
    could plausibly have gone either way, so it is less reliable.
  * P(detected | percentile) rises from ~0 to ~1 across that zone; where it crosses 0.5 is
    the empirical decision boundary; where it sits between p_lo and p_hi is the uncertain band.

This replaces "I rounded the quantiles" with a measured band you can cite. Measures only.

USAGE (protein-coding-scoped frame, as in Cell 1)
-------------------------------------------------
    import importlib
    nt = importlib.import_module("near_threshold_diagnostic")
    res = nt.visualize_near_threshold(expr=expression_features)   # or omit to load parquet
    print(res["suggested_band"])   # data-driven [lo, hi]

Author: (CellLineSelector team)   |   Depends on: pandas, numpy, matplotlib
"""

from __future__ import annotations

import json
import os
from typing import Iterable

import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    _HAVE_MPL = True
except Exception:
    _HAVE_MPL = False


EXPRESSION_FEATURES_PATH = "features/expression_features.parquet"
OUTDIR = "features/confidence"

# the hand-set band currently used by the scoring cell (drawn for comparison)
CURRENT_LO, CURRENT_HI = 0.15, 0.35
# a call is "uncertain" while P(detected|pct) is between these (data-driven band edges)
P_LO, P_HI = 0.10, 0.90


def _find_col(df: pd.DataFrame, candidates: Iterable[str], what: str = "column") -> str:
    lower = {c.lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in lower:
            return lower[cand.lower()]
    raise KeyError(f"Could not resolve {what}: tried {list(candidates)}; "
                   f"available = {list(df.columns)[:40]}")


def _as_bool_num(s: pd.Series) -> pd.Series:
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
        if t in ("true", "t", "1", "yes"):
            return 1.0
        if t in ("false", "f", "0", "no"):
            return 0.0
        return np.nan
    return s.map(_m).astype("float64")


def _q(a: np.ndarray, q: float) -> float:
    return float(np.quantile(a, q)) if a.size else float("nan")


def _cross(x: np.ndarray, y: np.ndarray, level: float) -> float:
    """First percentile x where the monotone-ish curve y crosses `level` (linear interp)."""
    d = y - level
    sign = np.sign(d)
    idx = np.where(np.diff(sign) != 0)[0]
    if not len(idx):
        return float("nan")
    i = idx[0]
    y0, y1, x0, x1 = y[i], y[i + 1], x[i], x[i + 1]
    if y1 == y0:
        return float(x0)
    return float(x0 + (level - y0) * (x1 - x0) / (y1 - y0))


def visualize_near_threshold(expr=None, expr_path: str = EXPRESSION_FEATURES_PATH,
                             bins: int = 100, lo: float = CURRENT_LO, hi: float = CURRENT_HI,
                             p_lo: float = P_LO, p_hi: float = P_HI,
                             outdir: str = OUTDIR, show: bool = True) -> dict:
    if expr is None:
        if not os.path.exists(expr_path):
            raise FileNotFoundError(f"expression_features not found at {expr_path}")
        expr = pd.read_parquet(expr_path, columns=None)

    c_pct = _find_col(expr, ("expr_percentile_depmap",), "percentile")
    c_dm  = _find_col(expr, ("detected_depmap",), "detected_depmap")

    pct = pd.to_numeric(expr[c_pct], errors="coerce").to_numpy()
    det = _as_bool_num(expr[c_dm]).to_numpy()
    m = ~np.isnan(pct) & ~np.isnan(det)
    pct, det = pct[m], det[m]
    p_det, p_not = pct[det == 1.0], pct[det == 0.0]

    # --- quantiles that motivate the band ---
    qs = {
        "not_detected": {"n": int(p_not.size), "median": _q(p_not, .5), "q99": _q(p_not, .99)},
        "detected":     {"n": int(p_det.size), "median": _q(p_det, .5),
                         "q01": _q(p_det, .01), "q25": _q(p_det, .25)},
    }
    print(f"[nt] not-detected: n={qs['not_detected']['n']:,} "
          f"median={qs['not_detected']['median']:.3f} q99={qs['not_detected']['q99']:.3f}")
    print(f"[nt] detected:     n={qs['detected']['n']:,} "
          f"median={qs['detected']['median']:.3f} q01={qs['detected']['q01']:.3f} "
          f"q25={qs['detected']['q25']:.3f}")

    # --- P(detected | percentile) over common bins (fast, exact; no KDE on 28M) ---
    edges = np.linspace(0, 1, bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2
    h_det, _ = np.histogram(p_det, bins=edges)
    h_not, _ = np.histogram(p_not, bins=edges)
    total = h_det + h_not
    with np.errstate(invalid="ignore", divide="ignore"):
        p_curve = np.where(total > 0, h_det / total, np.nan)

    valid = ~np.isnan(p_curve)
    xc, yc = centers[valid], p_curve[valid]
    crossover = _cross(xc, yc, 0.50)                    # empirical decision boundary
    band_lo = _cross(xc, yc, p_lo)                       # data-driven uncertain-band edges
    band_hi = _cross(xc, yc, p_hi)
    print(f"[nt] P(detected|pct) crosses 0.5 at percentile = {crossover:.3f}  (decision boundary)")
    print(f"[nt] data-driven uncertain band [P in {p_lo}-{p_hi}] = "
          f"[{band_lo:.3f}, {band_hi:.3f}]   (current hand-set: [{lo}, {hi}])")

    result = {
        "quantiles": qs,
        "decision_boundary_pct": crossover,
        "suggested_band": [band_lo, band_hi],
        "current_band": [lo, hi],
        "p_curve_levels": {"p_lo": p_lo, "p_hi": p_hi},
    }

    # --- figure ---
    if _HAVE_MPL:
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 7), sharex=True)

        ax1.hist(p_not, bins=edges, density=True, alpha=0.55, color="#c0392b",
                 label=f"not detected (n={p_not.size:,})")
        ax1.hist(p_det, bins=edges, density=True, alpha=0.55, color="#27ae60",
                 label=f"detected (n={p_det.size:,})")
        ax1.axvspan(lo, hi, color="#888", alpha=0.15, label=f"current band [{lo}, {hi}]")
        for x, lab, c in ((qs["not_detected"]["q99"], "not-det q99", "#c0392b"),
                          (qs["detected"]["q01"], "det q01", "#27ae60"),
                          (qs["detected"]["q25"], "det q25", "#1e8449")):
            ax1.axvline(x, color=c, ls=":", lw=1)
            ax1.text(x, ax1.get_ylim()[1] * 0.92, f" {lab}={x:.2f}", fontsize=7, color=c, rotation=90,
                     va="top")
        ax1.set_ylabel("density (each group normalised)")
        ax1.set_title("Percentile by DepMap detection call — the overlap is the uncertain zone")
        ax1.legend(fontsize=8, loc="upper center")

        ax2.plot(xc, yc, color="#2c3e50", lw=1.6, label="P(detected | percentile)")
        ax2.axhline(0.5, color="#999", lw=0.8, ls="--")
        ax2.axvspan(band_lo, band_hi, color="#f39c12", alpha=0.20,
                    label=f"data-driven band [{band_lo:.2f}, {band_hi:.2f}]")
        ax2.axvspan(lo, hi, color="#888", alpha=0.12)
        ax2.axvline(crossover, color="#8e44ad", lw=1.2,
                    label=f"decision boundary = {crossover:.2f}")
        ax2.set_ylim(0, 1)
        ax2.set_xlim(0, 1)
        ax2.set_xlabel("expr_percentile_depmap")
        ax2.set_ylabel("P(detected)")
        ax2.set_title("P(detected | percentile): where it is far from 0/1, the call is uncertain")
        ax2.legend(fontsize=8, loc="center right")

        fig.tight_layout()
        if outdir:
            os.makedirs(outdir, exist_ok=True)
            path = os.path.join(outdir, "near_threshold_band.png")
            fig.savefig(path, dpi=130, bbox_inches="tight")
            print(f"[nt] wrote {path}")
            with open(os.path.join(outdir, "near_threshold_band.json"), "w") as f:
                json.dump(result, f, indent=2, default=str)
        if show:
            plt.show()
        else:
            plt.close(fig)
    else:
        print("[nt] matplotlib not available - numbers computed, plot skipped")

    return result


if __name__ == "__main__":
    print(__doc__)
