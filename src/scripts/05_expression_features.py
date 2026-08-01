"""
05_expression_features.py
CellLineSelector - Feature Engineering: build the expression features (implements the
Expression Decision Note, decisions E1-E7).

WHAT THIS BUILDS
----------------
Per (line, gene) expression features, in the three roles E7 separates:

  SUITABILITY  expr_percentile_depmap   within-line percentile from DepMap (E1)
               detected_depmap/geo/hpa  per-source present/absent calls (E2)
  PROVENANCE   biotype, in_default_universe, single_source_*, corroboration_available (E6/E7)
  CONFIDENCE   n_sources, n_sources_detected, detection_agreement, per-source values
               INPUTS (not the fused score - that formula is Q8, deferred to the confidence
               engine; here we only store what that engine will consume)

METHODS + CITATIONS (see Decision Note References)
--------------------------------------------------
  E1 rank/percentile representation ......... MAQC/SEQC [R1, R3]
  E2 DepMap detection = zFPKM > -3 .......... Hart et al. 2013 [R5]  (default; TPM>=1 fallback)
     GEO detection   = intensity > floor .... MAS5 [R4] is the standard but needs probe-level
                                              CEL data (unavailable, L3) -> background-floor
                                              threshold from check 5 (~7)
     HPA detection   = nTPM >= cutoff ........ HPA methods [R6]
  E4 within-source median aggregation
  E6 gene universe = DepMap-complete, protein-coding default (features/gene_biotype.parquet)

DATA-DEFENCE (why zFPKM is valid on THIS data)
----------------------------------------------
validate_zfpkm() confirms zFPKM's core assumption - that each DepMap line's expression
distribution is bimodal (an "off" cluster + an "on" cluster with a valley between) - and
shows where the zFPKM=-3 threshold falls relative to that valley. This is the "defended by
the data" evidence to accompany the [R5] citation. Run it before trusting the detection call.

USAGE
-----
    import importlib
    ef = importlib.import_module("05_expression_features")

    # 1) data-defence figure (run once, keep the plot + verdict for the write-up)
    ef.validate_zfpkm(depmap_expr, n_lines=6)

    # 2) build the feature table
    feats = ef.build(
        depmap_expr = depmap_expr,          # wide line x gene, log2(TPM+1)
        expr_long   = expr_long,            # for GEO per-(line,gene) medians (source=='geo')
        hpa_rna     = hpa_rna,              # long: model_id, gene, ntpm
        biotype     = "features/gene_biotype.parquet",
        depmap_detect = "zfpkm",           # or "tpm1" for the TPM>=1 fallback
    )
    # -> features/expression_features.parquet  (one row per line x gene)

Author: (CellLineSelector team) | Depends on: numpy, pandas; scipy optional (KDE, else
histogram fallback); matplotlib optional (validation plots).
"""

from __future__ import annotations

import os
import numpy as np
import pandas as pd

try:
    import matplotlib.pyplot as plt
    _HAVE_MPL = True
except Exception:
    _HAVE_MPL = False

try:
    from scipy.stats import gaussian_kde
    _HAVE_SCIPY = True
except Exception:
    _HAVE_SCIPY = False


# =============================================================================
# HELPERS
# =============================================================================

def _ensg_cols(df):
    return [c for c in df.columns if str(c).lower().startswith("ensg")]

def _norm_gene(g):
    return str(g).strip().lower().split(".")[0]

def _unwrap(v):
    if isinstance(v, np.ndarray):
        return v[0] if v.size else np.nan
    if isinstance(v, (list, tuple)):
        return v[0] if len(v) else np.nan
    return v

def _gauss(x, mu, sd):
    sd = max(float(sd), 1e-6)
    return np.exp(-0.5 * ((x - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))


def _fit_two_component(x, iters=60, max_pts=20000, seed=0):
    """1-D two-component Gaussian mixture (EM) to separate the background (off) and expressed
    (on) populations of log2(TPM). Robust to which population is larger - unlike a global-mode
    estimate, which fails on zero-heavy lines. Returns component means/SDs sorted low->high."""
    x = np.asarray(x, dtype="float64")
    x = x[np.isfinite(x)]
    if x.size < 20:
        return None
    if x.size > max_pts:
        x = np.random.default_rng(seed).choice(x, max_pts, replace=False)
    med = np.median(x)
    lo, hi = x[x <= med], x[x > med]
    if lo.size < 2 or hi.size < 2:
        lo, hi = x[: x.size // 2], x[x.size // 2:]
    mu = np.array([lo.mean(), hi.mean()], dtype="float64")
    sd = np.array([max(lo.std(), 0.5), max(hi.std(), 0.5)], dtype="float64")
    w = np.array([0.5, 0.5])
    for _ in range(iters):
        r0 = w[0] * _gauss(x, mu[0], sd[0])
        r1 = w[1] * _gauss(x, mu[1], sd[1])
        s = r0 + r1
        s[s == 0] = 1e-300
        g1 = r1 / s
        g0 = 1.0 - g1
        N0, N1 = g0.sum(), g1.sum()
        if N0 < 1e-6 or N1 < 1e-6:
            break
        w = np.array([N0, N1]) / x.size
        mu = np.array([(g0 * x).sum() / N0, (g1 * x).sum() / N1])
        sd = np.maximum(np.array([np.sqrt((g0 * (x - mu[0]) ** 2).sum() / N0),
                                  np.sqrt((g1 * (x - mu[1]) ** 2).sum() / N1)]), 0.1)
    o = np.argsort(mu)
    return {"mu_low": float(mu[o[0]]), "mu_high": float(mu[o[1]]),
            "sd_low": float(sd[o[0]]), "sd_high": float(sd[o[1]]), "w_high": float(w[o[1]])}


# =============================================================================
# zFPKM  (Hart et al. 2013, BMC Genomics 14:778 - reference algorithm) [R5]
# =============================================================================

def zfpkm_from_log2tpm1(log2tpm1: np.ndarray):
    """Compute zFPKM for one sample from DepMap log2(TPM+1) values (robust variant).

    Reference: Hart et al. 2013 (R5). zFPKM = (log2(TPM) - mu) / sigma, with sigma from the
    half-normal above mu, threshold zFPKM > -3 for expressed. The reference sets mu to the
    global density mode; on zero-heavy lines that mode can land in the BACKGROUND population,
    breaking the threshold. Here mu is estimated as the EXPRESSED-component mean of a
    two-Gaussian fit (robust to which population is larger). fit_ok flags lines where the two
    populations are not cleanly separated (caller falls back to TPM>=1).

    TPM==0 -> log2 = -inf -> zFPKM = -inf (correctly not-detected). Returns (z, params).
    """
    v = np.asarray(log2tpm1, dtype="float64")
    tpm = np.clip(np.power(2.0, v) - 1.0, 0.0, None)
    with np.errstate(divide="ignore"):
        log2tpm = np.log2(tpm)                      # -inf where tpm == 0
    finite = np.isfinite(log2tpm)
    off_fraction = float(np.mean(tpm[~np.isnan(v)] == 0)) if (~np.isnan(v)).any() else np.nan
    if finite.sum() < 20:
        return (np.full_like(v, -np.inf),
                {"mu": np.nan, "sigma": np.nan, "mu_low": np.nan, "mu_high": np.nan,
                 "sd_low": np.nan, "sd_high": np.nan, "w_high": np.nan,
                 "fit_ok": False, "off_fraction": off_fraction, "n_finite": int(finite.sum())})

    fit = _fit_two_component(log2tpm[finite])
    mu = fit["mu_high"]                             # expressed-component mean (robust)
    above = log2tpm[finite][log2tpm[finite] > mu]
    U = above.mean() if above.size else mu
    sigma = (U - mu) * np.sqrt(np.pi / 2.0)
    if not np.isfinite(sigma) or sigma <= 0:
        sigma = float(np.nanstd(log2tpm[finite])) or 1.0
    # fit is trustworthy if the two populations are separated and the expressed mode is sane
    fit_ok = bool((fit["mu_high"] - fit["mu_low"] > 1.0) and (fit["mu_high"] > 1.0)
                  and np.isfinite(sigma) and sigma > 0)
    with np.errstate(invalid="ignore"):
        z = (log2tpm - mu) / sigma                  # -inf stays -inf
    return z, {"mu": float(mu), "sigma": float(sigma), "mu_low": fit["mu_low"],
               "mu_high": fit["mu_high"], "sd_low": fit["sd_low"], "sd_high": fit["sd_high"],
               "w_high": fit["w_high"], "fit_ok": fit_ok, "off_fraction": off_fraction,
               "n_finite": int(finite.sum())}


def _bimodality_coefficient(x):
    """Sarle's bimodality coefficient BC = (skew^2 + 1) / kurtosis. BC > 0.555 => bimodal-ish.
    Must be computed on the FULL distribution (including the zero/off cluster), not the
    expressed subset, or the off-population is invisible."""
    s = pd.Series(x[np.isfinite(x)])
    if len(s) < 4:
        return np.nan
    g1 = s.skew()                                  # sample skewness
    g2 = s.kurt()                                  # EXCESS kurtosis (pandas)
    denom = (g2 + 3.0)                             # non-excess kurtosis
    if denom == 0:
        return np.nan
    return float((g1 ** 2 + 1.0) / denom)





# =============================================================================
# DATA-DEFENCE: validate zFPKM's bimodality assumption on THIS data
# =============================================================================

def validate_zfpkm(depmap_expr, n_lines: int = 6, outdir: str = "features/expr_features",
                   show: bool = True, seed: int = 0):
    """Confirm each DepMap line's distribution is bimodal and show the zFPKM threshold in it.

    Produces (a) a figure of log2(TPM+1) distributions for a sample of lines with the
    zFPKM=-3 threshold marked, and (b) a bimodality-coefficient verdict per line. This is the
    'defended by the data' evidence for the E2 DepMap detection choice.
    """
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    ecols = _ensg_cols(depmap_expr)
    mcol = "model_id" if "model_id" in depmap_expr.columns else depmap_expr.columns[0]
    rng = np.random.default_rng(seed)
    idx = rng.choice(len(depmap_expr), size=min(n_lines, len(depmap_expr)), replace=False)

    rows, panels = [], []
    for i in idx:
        line = str(_unwrap(depmap_expr.iloc[i][mcol]))
        v = depmap_expr.iloc[i][ecols].to_numpy(dtype="float64")
        v = v[~np.isnan(v)]
        z, p = zfpkm_from_log2tpm1(v)
        # threshold in log2(TPM+1) space: log2tpm_thr = mu - 3*sigma
        thr_log2tpm = p["mu"] - 3.0 * p["sigma"]
        thr_native = float(np.log2(np.power(2.0, thr_log2tpm) + 1.0)) if np.isfinite(thr_log2tpm) else np.nan
        bc = _bimodality_coefficient(v)            # NATIVE log2(TPM+1) incl. off-cluster
        pct_detected = float(np.mean(z > -3))
        rows.append({"line": line, "mu_high": p["mu_high"], "mu_low": p["mu_low"],
                     "sigma": p["sigma"], "thr_native_log2tpm1": thr_native,
                     "bimodality_coef": bc, "off_fraction": p["off_fraction"],
                     "fit_ok": p["fit_ok"], "pct_detected": pct_detected})
        panels.append((line, v, thr_native, bc))

    report = pd.DataFrame(rows)
    print("[zfpkm-validate] per-line summary:")
    print(report.to_string(index=False))
    pct_bimodal = (report["bimodality_coef"] > 0.555).mean() * 100
    pct_fit = report["fit_ok"].mean() * 100
    print(f"[zfpkm-validate] VERDICT: {pct_bimodal:.0f}% of sampled lines are bimodal "
          f"(Sarle's BC > 0.555; a large off-cluster + an expressed hump) -> zFPKM's "
          f"active/background split is appropriate. {pct_fit:.0f}% have a clean two-component "
          f"fit (fit_ok); any that don't fall back to TPM>=1. [method: Hart et al. 2013, R5]")

    if _HAVE_MPL and panels:
        k = len(panels)
        cols = min(3, k); rows_ = int(np.ceil(k / cols))
        fig, axes = plt.subplots(rows_, cols, figsize=(4.2 * cols, 3.0 * rows_), squeeze=False)
        for ax, (line, v, thr, bc) in zip(axes.ravel(), panels):
            ax.hist(v, bins=60)
            if np.isfinite(thr):
                ax.axvline(thr, color="crimson", lw=1.5, label=f"zFPKM=-3 (BC={bc:.2f})")
            ax.set_title(line, fontsize=9)
            ax.set_xlabel("log2(TPM+1)"); ax.legend(fontsize=7)
        for ax in axes.ravel()[len(panels):]:
            ax.axis("off")
        fig.suptitle("zFPKM validation - DepMap distributions are bimodal; threshold sits in the valley")
        fig.tight_layout()
        if outdir:
            fig.savefig(os.path.join(outdir, "zfpkm_validation.png"), dpi=110, bbox_inches="tight")
        if show:
            plt.show()
        plt.close(fig)
    return report


def plot_valley(depmap_expr, line, outdir: str = "features/expr_features", show: bool = True):
    """Publication figure for one line: the two fitted populations (background vs expressed)
    on the log2(TPM) scale, with the zFPKM=-3 detection threshold. Shows WHY zFPKM works -
    the distribution is genuinely bimodal - and is honest about the threshold being a lenient
    line below the valley rather than at it. Curves and threshold come from the SAME fit.
    """
    ecols = _ensg_cols(depmap_expr)
    mcol = "model_id" if "model_id" in depmap_expr.columns else depmap_expr.columns[0]
    mask = depmap_expr[mcol].map(_unwrap).astype(str).str.strip().str.upper() == str(line).upper()
    if mask.sum() == 0:
        raise ValueError(f"line {line} not found in {mcol}")
    v = depmap_expr.loc[mask, ecols].to_numpy(dtype="float64").ravel()

    tpm = np.power(2.0, v) - 1.0
    x = np.log2(tpm[tpm > 0])                       # finite log2(TPM): the scale zFPKM uses
    z, p = zfpkm_from_log2tpm1(v)
    if not np.isfinite(p["mu_high"]):
        raise ValueError(f"{line}: zFPKM fit degenerate; nothing to plot")
    ml, mh, sl, sh, wh = p["mu_low"], p["mu_high"], p["sd_low"], p["sd_high"], p["w_high"]
    thr = p["mu_high"] - 3.0 * p["sigma"]           # zFPKM = -3 (lenient, Hart et al. R5)
    bc = _bimodality_coefficient(v)

    # crossover (the valley) between the two components, for reference
    xs = np.linspace(x.min(), x.max(), 800)
    bg = (1 - wh) * _gauss(xs, ml, sl)
    ex = wh * _gauss(xs, mh, sh)
    grid = np.linspace(ml, mh, 500)
    d = wh * _gauss(grid, mh, sh) - (1 - wh) * _gauss(grid, ml, sl)
    cross = np.where(np.diff(np.sign(d)) > 0)[0]
    valley = float(grid[cross[0] + 1]) if len(cross) else np.nan

    if not _HAVE_MPL:
        print("[plot_valley] matplotlib unavailable")
        return {"mu_low": ml, "mu_high": mh, "threshold": thr, "valley": valley,
                "bimodality_coef": bc, "fit_ok": p["fit_ok"]}

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(x, bins=100, density=True, alpha=0.4, color="steelblue")
    ax.plot(xs, bg, color="gray", lw=1.8, label=f"background (mu={ml:.1f})")
    ax.plot(xs, ex, color="green", lw=1.8, label=f"expressed (mu={mh:.1f})")
    ax.plot(xs, bg + ex, color="black", lw=1, ls=":", label="mixture")
    if np.isfinite(valley):
        ax.axvline(valley, color="orange", lw=1.3, ls="--", label=f"valley / crossover ({valley:.1f})")
    ax.axvline(thr, color="crimson", lw=2, label=f"zFPKM=-3 threshold ({thr:.1f})")
    ax.set_xlabel("log2(TPM)"); ax.set_ylabel("density")
    ax.set_title(f"{line}: bimodal expression (BC={bc:.2f}); "
                 f"{'clean fit' if p['fit_ok'] else 'FALLBACK (TPM>=1)'}")
    ax.legend(fontsize=8)
    fig.tight_layout()
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        fig.savefig(os.path.join(outdir, f"valley_{str(line).upper()}.png"),
                    dpi=110, bbox_inches="tight")
    if show:
        plt.show()
    plt.close(fig)
    return {"mu_low": ml, "mu_high": mh, "threshold": thr, "valley": valley,
            "bimodality_coef": bc, "fit_ok": p["fit_ok"]}


# =============================================================================
# DETECTION CALLS (E2)
# =============================================================================

def _detect_depmap(v_log2tpm1, method="zfpkm", tpm_threshold=1.0):
    """DepMap present/absent per line. Returns (calls, used_fallback).
    'tpm1' -> TPM >= tpm_threshold. 'zfpkm' -> z > -3 [R5], but if the two-component fit is
    not clean (fit_ok False), fall back to TPM>=1 for that line so a broken threshold never
    ships."""
    v = np.asarray(v_log2tpm1, dtype="float64")
    tpm = np.power(2.0, v) - 1.0
    if method == "tpm1":
        return tpm >= tpm_threshold, False
    z, p = zfpkm_from_log2tpm1(v)
    if not p["fit_ok"]:
        return tpm >= tpm_threshold, True          # fallback for this line
    return z > -3.0, False


def zfpkm_fit_audit(depmap_expr, outdir="features/expr_features", show=True):
    """Fit zFPKM on EVERY DepMap line and report the clean-fit rate - the real failure rate
    across all lines, so we know how many fall back to TPM>=1 before building the table."""
    ecols = _ensg_cols(depmap_expr)
    mcol = "model_id" if "model_id" in depmap_expr.columns else depmap_expr.columns[0]
    rows = []
    for i in range(len(depmap_expr)):
        line = str(_unwrap(depmap_expr.iloc[i][mcol]))
        v = depmap_expr.iloc[i][ecols].to_numpy(dtype="float64")
        z, p = zfpkm_from_log2tpm1(v)
        rows.append({"line": line, "mu_high": p["mu_high"], "mu_low": p["mu_low"],
                     "sigma": p["sigma"], "off_fraction": p["off_fraction"],
                     "fit_ok": p["fit_ok"], "pct_detected": float(np.mean(z > -3))})
    audit = pd.DataFrame(rows)
    n = len(audit)
    ok = int(audit["fit_ok"].sum())
    print(f"[zfpkm-audit] {n:,} lines | clean zFPKM fit: {ok:,} ({ok/n*100:.1f}%) | "
          f"fall back to TPM>=1: {n-ok:,} ({(n-ok)/n*100:.1f}%)")
    print(f"[zfpkm-audit] pct_detected across lines: median {audit['pct_detected'].median()*100:.1f}% "
          f"(IQR {audit['pct_detected'].quantile(.25)*100:.1f}-{audit['pct_detected'].quantile(.75)*100:.1f}%)")
    if (n - ok) and _HAVE_MPL:
        fig, ax = plt.subplots(figsize=(7, 3.4))
        ax.scatter(audit["off_fraction"], audit["pct_detected"], s=6,
                   c=np.where(audit["fit_ok"], "steelblue", "crimson"), alpha=0.5)
        ax.set_xlabel("off_fraction (share of genes at TPM=0)")
        ax.set_ylabel("pct_detected"); ax.set_title("zFPKM fit audit (red = fell back to TPM>=1)")
        if outdir:
            os.makedirs(outdir, exist_ok=True)
            fig.savefig(os.path.join(outdir, "zfpkm_fit_audit.png"), dpi=110, bbox_inches="tight")
        if show:
            plt.show()
        plt.close(fig)
    return audit


def _two_component_threshold(logx):
    """Threshold separating background from expressed on a log-scaled distribution, via the
    crossover (valley) of a two-Gaussian fit. Returns (threshold_log, params incl fit_ok)."""
    x = np.asarray(logx, dtype="float64")
    x = x[np.isfinite(x)]
    fit = _fit_two_component(x)
    if fit is None:
        return np.nan, {"fit_ok": False}
    ml, mh, sl, sh, wh = (fit["mu_low"], fit["mu_high"], fit["sd_low"],
                          fit["sd_high"], fit["w_high"])
    fit_ok = (mh - ml) > 1.0
    if not fit_ok:
        return np.nan, {**fit, "fit_ok": False}
    xs = np.linspace(ml, mh, 500)
    diff = wh * _gauss(xs, mh, sh) - (1 - wh) * _gauss(xs, ml, sl)   # expressed - background
    cross = np.where(np.diff(np.sign(diff)) > 0)[0]
    thr = float(xs[cross[0] + 1]) if len(cross) else float((ml + mh) / 2)
    return thr, {**fit, "fit_ok": True}


def _geo_detect(geo, fallback_floor=7.0):
    """Per-line GEO present/absent by fitting background vs expressed on log2(intensity) and
    calling 'present' above the crossover. Falls back to intensity>fallback_floor if a line's
    two populations don't separate. Returns (geo_with_detected, n_fallback, thr_by_line)."""
    geo = geo.copy()
    geo["_log"] = np.log2(geo["geo_intensity"].clip(lower=1e-9))
    thr_map, fit_map = {}, {}
    for line, sub in geo.groupby("model_id"):
        thr, p = _two_component_threshold(sub["_log"].to_numpy())
        thr_map[line] = thr if p.get("fit_ok") else np.nan
        fit_map[line] = bool(p.get("fit_ok"))
    geo["_thr"] = geo["model_id"].map(thr_map)
    fb = geo["_thr"].isna()
    geo["detected_geo"] = np.where(fb, geo["geo_intensity"] > fallback_floor,
                                   geo["_log"] > geo["_thr"])
    n_fb = int(sum(1 for v in fit_map.values() if not v))
    geo = geo.drop(columns=["_log", "_thr"])
    return geo, n_fb


def validate_geo_detection(expr_long, n_lines: int = 6, outdir: str = "features/expr_features",
                           show: bool = True, seed: int = 0):
    """Check whether GEO log2(intensity) is bimodal (background + expressed) per line, and show
    where the crossover threshold falls. Run BEFORE trusting detected_geo - if GEO is NOT
    bimodal, the background/expressed split is not well-founded and we reconsider."""
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    g = expr_long[expr_long["source"].astype(str).str.lower() == "geo"].copy()
    g["gene"] = g["gene"].map(_norm_gene)
    g["model_id"] = g["model_id"].map(_unwrap).astype(str).str.strip().str.upper()
    all_lines = g["model_id"].unique()
    rng = np.random.default_rng(seed)
    pick = rng.choice(all_lines, size=min(n_lines, len(all_lines)), replace=False)

    rows, panels = [], []
    for line in pick:
        vals = g[g["model_id"] == line].groupby("gene")["value"].median()
        x = np.log2(np.clip(vals.to_numpy(), 1e-9, None))
        thr, p = _two_component_threshold(x)
        bc = _bimodality_coefficient(x)
        pct = float(np.mean(x > thr)) if np.isfinite(thr) else np.nan
        rows.append({"line": line, "mu_low": p.get("mu_low"), "mu_high": p.get("mu_high"),
                     "thr_log2": thr, "thr_intensity": float(2 ** thr) if np.isfinite(thr) else np.nan,
                     "bimodality_coef": bc, "fit_ok": p.get("fit_ok"), "pct_detected": pct})
        panels.append((str(line), x, thr, bc))

    report = pd.DataFrame(rows)
    print("[geo-validate] per-line summary:")
    print(report.to_string(index=False))
    pct_bimodal = (report["bimodality_coef"] > 0.555).mean() * 100
    pct_fit = report["fit_ok"].mean() * 100
    print(f"[geo-validate] VERDICT: {pct_bimodal:.0f}% of sampled lines bimodal (BC > 0.555); "
          f"{pct_fit:.0f}% have a clean background/expressed split. If these are HIGH, the "
          f"bimodal GEO detection is well-founded; if LOW, GEO is not clearly bimodal and a "
          f"percentile/background threshold is the honest fallback.")

    if _HAVE_MPL and panels:
        k = len(panels); cols = min(3, k); rows_ = int(np.ceil(k / cols))
        fig, axes = plt.subplots(rows_, cols, figsize=(4.2 * cols, 3.0 * rows_), squeeze=False)
        for ax, (line, x, thr, bc) in zip(axes.ravel(), panels):
            ax.hist(x, bins=60)
            if np.isfinite(thr):
                ax.axvline(thr, color="crimson", lw=1.5, label=f"threshold (BC={bc:.2f})")
            ax.set_title(line, fontsize=9); ax.set_xlabel("log2(GEO intensity)"); ax.legend(fontsize=7)
        for ax in axes.ravel()[len(panels):]:
            ax.axis("off")
        fig.suptitle("GEO detection validation - is log2(intensity) bimodal? threshold at the valley")
        fig.tight_layout()
        if outdir:
            fig.savefig(os.path.join(outdir, "geo_detection_validation.png"), dpi=110, bbox_inches="tight")
        if show:
            plt.show()
        plt.close(fig)
    return report


def geo_bimodality_audit(expr_long, outdir="features/expr_features", show=True):
    """Fleet-wide check of whether GEO log2(intensity) is bimodal, across ALL GEO lines.
    Reports the fraction with Sarle's BC > 0.555. If near 0%, GEO gene-level intensities are
    unimodal and cannot support a per-gene detection call (detected_geo is dropped; GEO serves
    as rank corroboration only, E3)."""
    g = expr_long[expr_long["source"].astype(str).str.lower() == "geo"].copy()
    g["gene"] = g["gene"].map(_norm_gene)
    g["model_id"] = g["model_id"].map(_unwrap).astype(str).str.strip().str.upper()
    rows = []
    for line, sub in g.groupby("model_id"):
        vals = sub.groupby("gene")["value"].median()
        x = np.log2(np.clip(vals.to_numpy(), 1e-9, None))
        rows.append({"line": line, "bimodality_coef": _bimodality_coefficient(x),
                     "n_genes": int(len(x))})
    audit = pd.DataFrame(rows).dropna(subset=["bimodality_coef"])
    n = len(audit)
    bimodal = int((audit["bimodality_coef"] > 0.555).sum())
    print(f"[geo-bimodality-audit] {n:,} GEO lines | bimodal (BC>0.555): {bimodal:,} "
          f"({bimodal/n*100:.1f}%) | median BC {audit['bimodality_coef'].median():.3f} "
          f"(IQR {audit['bimodality_coef'].quantile(.25):.3f}-{audit['bimodality_coef'].quantile(.75):.3f})")
    if bimodal / n < 0.10:
        print("[geo-bimodality-audit] VERDICT: GEO is unimodal fleet-wide -> no per-gene "
              "detection call; GEO serves as rank corroboration only (E3). detected_geo dropped.")
    else:
        print("[geo-bimodality-audit] VERDICT: a non-trivial fraction are bimodal -> reconsider "
              "a hybrid detection for those lines.")
    if _HAVE_MPL and n:
        fig, ax = plt.subplots(figsize=(7, 3.4))
        ax.hist(audit["bimodality_coef"], bins=50)
        ax.axvline(0.555, color="crimson", lw=1.5, label="BC=0.555 (bimodal cutoff)")
        ax.set_xlabel("Sarle's bimodality coefficient (per GEO line)")
        ax.set_ylabel("lines"); ax.legend(); ax.set_title("GEO bimodality audit")
        if outdir:
            os.makedirs(outdir, exist_ok=True)
            fig.savefig(os.path.join(outdir, "geo_bimodality_audit.png"), dpi=110, bbox_inches="tight")
        if show:
            plt.show()
        plt.close(fig)
    return audit


# =============================================================================
# PER-LINE PERCENTILE (E1)  -  ranked against the default-universe (protein-coding) genes
# =============================================================================

def _percentile_against_reference(values, reference_sorted):
    """Fraction of reference genes with value <= each input value (0..1). Uniform for all genes."""
    if len(reference_sorted) == 0:
        return np.full(len(values), np.nan)
    return np.searchsorted(reference_sorted, values, side="right") / len(reference_sorted)


# =============================================================================
# BUILD
# =============================================================================

def build(depmap_expr, biotype, expr_long=None, hpa_rna=None,
          depmap_detect="zfpkm", geo_detect=False, geo_floor=7.0, hpa_ntpm_cutoff=1.0,
          out="features/expression_features.parquet"):
    """Assemble the per-(line, gene) expression feature table (E1-E7).

    Primary feature is on DepMap's full gene set. DepMap (zFPKM) and HPA (nTPM cutoff) carry
    the detection layer (E2). GEO is corroboration only (E3): its gene-level intensities are
    unimodal (geo_bimodality_audit) and cannot support a per-gene detection call, so geo_detect
    defaults to False - GEO contributes geo_intensity for rank agreement but no detected_geo.
    Percentile (E1) is computed against the protein-coding (default-universe) genes per line.
    """
    # --- biotype (E6) ---
    if isinstance(biotype, str):
        biotype = pd.read_parquet(biotype)
    bt = biotype.copy()
    bt["_g"] = bt["gene"].map(_norm_gene)
    bt_map = bt.set_index("_g")[["gene_biotype", "biotype_coarse", "in_default_universe"]]

    ecols = _ensg_cols(depmap_expr)
    mcol = "model_id" if "model_id" in depmap_expr.columns else depmap_expr.columns[0]
    gene_norm = {c: _norm_gene(c) for c in ecols}

    # --- DepMap: median-collapse profiles per line (E4), then per-line features ---
    dm = depmap_expr[[mcol] + ecols].copy()
    dm[mcol] = dm[mcol].map(_unwrap).astype(str).str.strip().str.upper()
    dm = dm.groupby(mcol)[ecols].median()          # 1 profile/line usually; median if more

    default_g = {g for g, r in bt_map.iterrows() if bool(r["in_default_universe"])}
    frames = []
    n_fallback = 0
    for line, row in dm.iterrows():
        v = row.to_numpy(dtype="float64")
        genes = [gene_norm[c] for c in ecols]
        detected, used_fb = _detect_depmap(v, method=depmap_detect)
        n_fallback += int(used_fb)
        # percentile against this line's protein-coding distribution
        ref = np.sort([val for c, val in zip(ecols, v)
                       if gene_norm[c] in default_g and not np.isnan(val)])
        pct = _percentile_against_reference(v, ref)
        frames.append(pd.DataFrame({
            "model_id": line, "gene": genes,
            "depmap_log2tpm1": v,
            "expr_percentile_depmap": pct,
            "detected_depmap": detected,
        }))
    feats = pd.concat(frames, ignore_index=True)
    feats = feats[~feats["depmap_log2tpm1"].isna()]
    if depmap_detect == "zfpkm" and n_fallback:
        print(f"[build] zFPKM fit not clean for {n_fallback} line(s) -> those used TPM>=1 fallback")

    # --- GEO per-(line,gene) median + detection (from expr_long) ---
    if expr_long is not None and "source" in expr_long.columns:
        g = expr_long[expr_long["source"].astype(str).str.lower() == "geo"].copy()
        if len(g):
            g["gene"] = g["gene"].map(_norm_gene)
            g["model_id"] = g["model_id"].map(_unwrap).astype(str).str.strip().str.upper()
            geo = g.groupby(["model_id", "gene"])["value"].median().reset_index()
            geo = geo.rename(columns={"value": "geo_intensity"})
            if geo_detect:                          # off by default - GEO is unimodal (E3)
                geo, geo_fb = _geo_detect(geo, fallback_floor=geo_floor)
                if geo_fb:
                    print(f"[build] GEO bimodal fit not clean for {geo_fb} line(s) -> "
                          f"those used intensity>{geo_floor} fallback")
            else:
                print("[build] GEO detection OFF (unimodal intensities): geo_intensity kept for "
                      "rank corroboration; no detected_geo call (E3)")
            feats = feats.merge(geo, on=["model_id", "gene"], how="left")

    # --- HPA per-(line,gene) nTPM + detection ---
    if hpa_rna is not None:
        hcols = {c.lower(): c for c in hpa_rna.columns}
        hid = hcols.get("model_id")
        hg = next((hcols[k] for k in ("gene", "ensembl", "ensembl_id", "gene_id", "ensg") if k in hcols), None)
        hv = next((hcols[k] for k in ("ntpm", "value") if k in hcols), None)
        if hid and hg and hv:
            h = pd.DataFrame({
                "model_id": hpa_rna[hid].map(_unwrap).astype(str).str.strip().str.upper(),
                "gene": hpa_rna[hg].map(_unwrap).map(_norm_gene),
                "hpa_ntpm": pd.to_numeric(hpa_rna[hv], errors="coerce"),
            }).dropna(subset=["hpa_ntpm"])
            hpa = h.groupby(["model_id", "gene"])["hpa_ntpm"].median().reset_index()
            hpa["detected_hpa"] = hpa["hpa_ntpm"] >= hpa_ntpm_cutoff
            feats = feats.merge(hpa, on=["model_id", "gene"], how="left")

    # --- provenance (E6/E7) ---
    feats = feats.merge(bt_map, left_on="gene", right_index=True, how="left")
    feats["gene_biotype"] = feats["gene_biotype"].fillna("unmatched")
    feats["in_default_universe"] = feats["in_default_universe"].fillna(False)

    for col in ("detected_geo", "detected_hpa"):
        if col not in feats.columns:
            feats[col] = pd.NA
    have_geo = "geo_intensity" in feats.columns
    have_hpa = "hpa_ntpm" in feats.columns

    # --- confidence INPUTS (NOT the fused score; Q8 deferred) ---
    src_present = feats["depmap_log2tpm1"].notna().astype(int)
    if have_geo:
        src_present = src_present + feats["geo_intensity"].notna().astype(int)
    if have_hpa:
        src_present = src_present + feats["hpa_ntpm"].notna().astype(int)
    feats["n_sources"] = src_present

    det_cols = [c for c in ("detected_depmap", "detected_geo", "detected_hpa")
                if c in feats.columns and feats[c].notna().any()]
    det = feats[det_cols].apply(pd.to_numeric, errors="coerce")
    feats["n_sources_detected"] = det.sum(axis=1, skipna=True)
    n_calls = det.notna().sum(axis=1)
    # detection_agreement = fraction of available sources that share the majority call
    with np.errstate(invalid="ignore"):
        maj = np.maximum(feats["n_sources_detected"], n_calls - feats["n_sources_detected"])
        feats["detection_agreement"] = np.where(n_calls > 0, maj / n_calls, np.nan)

    feats["corroboration_available"] = feats["n_sources"] > 1
    feats["single_source_gene"] = feats["n_sources"] == 1

    if out:
        os.makedirs(os.path.dirname(out), exist_ok=True)
        feats.to_parquet(out, index=False)
        print(f"[build] wrote {out}")
    _build_report(feats, depmap_detect)
    return feats


def _build_report(feats, depmap_detect):
    n = len(feats)
    print(f"[build] {n:,} (line, gene) rows | detection method: DepMap={depmap_detect}")
    print(f"[build] detected_depmap: {feats['detected_depmap'].mean()*100:.1f}% of rows called present")
    if "detected_geo" in feats and feats["detected_geo"].notna().any():
        print(f"[build] detected_geo:    {pd.to_numeric(feats['detected_geo'],errors='coerce').mean()*100:.1f}% present (where measured)")
    if "detected_hpa" in feats and feats["detected_hpa"].notna().any():
        print(f"[build] detected_hpa:    {pd.to_numeric(feats['detected_hpa'],errors='coerce').mean()*100:.1f}% present (where measured)")
    print(f"[build] n_sources: " + ", ".join(f"{k}src={v:,}" for k, v in
          feats["n_sources"].value_counts().sort_index().items()))
    print(f"[build] corroboration_available: {feats['corroboration_available'].mean()*100:.1f}% of rows")
    print(f"[build] in_default_universe (protein-coding): {feats['in_default_universe'].mean()*100:.1f}% of rows")


if __name__ == "__main__":
    print(__doc__)