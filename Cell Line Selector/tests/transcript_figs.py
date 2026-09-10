#!/usr/bin/env python3
"""
generate_transcriptomic_figures.py

Figures 2-4 for the Transcriptomic Evidence section, plus the numerical
audit. Figure 1 (execution pipeline) is drawn natively in LaTeX/TikZ and
has no data dependency, so it is not produced here.

  Figure 2  fig2_scale_harmonisation   per-method distributions, 3 stages
  Figure 3  fig3_detection_regimes     regime allocation
  Figure 4  fig4_calibration_example   (a) calibration surface
                                       (b) worked biological example

This script REUSES the project implementation: it imports
src.scripts.transcriptomics and calls the same functions the scoring run
calls. It reimplements no estimator. Every number quoted in the manuscript
is emitted here into numerical_audit.csv, and the provenance block above
each figure records exactly how it was derived.

Run:
    python generate_transcriptomic_figures.py \
        --db db/celllineselector.duckdb \
        --results results/transcriptomics \
        --table transcriptomics_z \
        --out figs

Prerequisites:
    python 05_transcriptomics_stat_layer.py --build-geo-platform
    python 05_transcriptomics_stat_layer.py --calibrate-only
    python 05_transcriptomics_stat_layer.py \
        --params results/transcriptomics/chosen_params.json \
        --write-table transcriptomics_z
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# --- project root discovery (same idiom as the pipeline scripts) --------
PROJECT_ROOT = Path(__file__).resolve().parent
for _cand in [PROJECT_ROOT, *PROJECT_ROOT.parents]:
    if (_cand / "src").exists():
        PROJECT_ROOT = _cand
        break
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from src.scripts import transcriptomics as tx   # canonical implementation

SEED = 0
N_PROBE_GENES = 800     # gene sample for Fig 2, and Fig 3 fallback only

plt.rcParams.update({
    "figure.dpi": 150,
    "font.family": "serif",
    "font.serif": ["Times New Roman", "DejaVu Serif", "serif"],
    "font.size": 13,
    "axes.labelsize": 13,
    "axes.titlesize": 14,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "legend.fontsize": 11,
    "figure.titlesize": 16,
    "savefig.bbox": "tight",
    "axes.spines.top": False,
    "axes.spines.right": False,
})

REGIME_ORDER = ["ok", "censored_ranked", "below_detection", "no_data"]
REGIME_COLOUR = {"ok": "#5f9e5f", "censored_ranked": "#e6b34d",
                 "below_detection": "#c9605f", "no_data": "#b8b8b8"}
SOURCE_COLOUR = {"hpa_rna": "#4878a8", "depmap": "#5f9e5f", "geo": "#b07aa1"}

# Haematological markers, fixed BEFORE any score is inspected, so the
# worked example is not selected on its outcome. The candidate with the
# widest coverage in transcriptomics_z is the one used.
WORKED_EXAMPLE_CANDIDATES = {
    "ENSG00000081237": "PTPRC (CD45)",
    "ENSG00000156738": "MS4A1 (CD20)",
    "ENSG00000177455": "CD19",
    "ENSG00000198851": "CD3E",
    "ENSG00000010610": "CD4",
}

HAEM_PATTERN = "lymph|myeloid|leuk|blood|haem|hema"

AUDIT: dict[str, dict] = {}


def note(key, value, source):
    """Record a number and where it came from. Nothing is reported unaudited."""
    AUDIT[key] = {"value": value, "source": source}
    print(f"  [audit] {key} = {value}   <- {source}")


# ======================================================================
# Shared loading
# ======================================================================

def load_params(results_dir: Path) -> tx.Params:
    """Calibrated parameters. Never silently defaulted."""
    p = results_dir / "chosen_params.json"
    if not p.exists():
        raise FileNotFoundError(
            f"{p} not found. Run 05_transcriptomics_stat_layer.py "
            f"--calibrate-only first. Parameters must not be assumed.")
    params = tx.Params(**json.loads(p.read_text()))
    for k, v in asdict(params).items():
        note(f"param.{k}", v, str(p))
    return params


def build_context(con, params: tx.Params, rng):
    """Axes, colmaps, lineage map and a seeded probe sample, via project code."""
    axes = tx.build_method_axes(con)
    lineage_s = tx.load_lineage_map(con)
    colmaps = {ax.name: (None if ax.source == "hpa_rna"
                         else tx.ensg_columns(con, ax.table)) for ax in axes}
    genes_all = tx.gene_universe(con, axes, colmaps)
    if not genes_all:
        raise RuntimeError("gene universe is empty; check gene_roster join")

    probe = list(rng.choice(np.asarray(genes_all, dtype=object),
                            size=min(N_PROBE_GENES, len(genes_all)),
                            replace=False))
    for ax in axes:
        tx.profile_method(con, ax, colmaps.get(ax.name), probe, params.floor_q)

    models = sorted(set().union(*[set(ax.models) for ax in axes]))
    midx = {m: i for i, m in enumerate(models)}

    note("n_sources", len({a.source for a in axes}),
         "tx.build_method_axes -> distinct .source")
    note("n_methods", len(axes), "tx.build_method_axes -> len(axes)")
    note("n_methods_geo", sum(a.source == "geo" for a in axes),
         "tx.build_method_axes")
    note("n_samples_total", int(sum(len(a.raw_ids) for a in axes)),
         "sum(len(axis.raw_ids))")
    note("n_cell_lines_axis_union", len(models), "union of axis.models")
    note("n_genes_universe", len(genes_all), "tx.gene_universe")

    per_source_samples = {}
    for ax in axes:
        per_source_samples[ax.source] = \
            per_source_samples.get(ax.source, 0) + len(ax.raw_ids)
    note("n_samples_by_source", per_source_samples, "sum len(raw_ids) by source")

    return axes, colmaps, lineage_s, genes_all, probe, models, midx


# ======================================================================
# FIGURE 2 - method / scale harmonisation
# ----------------------------------------------------------------------
# Input      : depmap_expr, geo_expr, hpa_rna, geo_platform, sample_info
# Columns    : ENSG value columns; axis id columns (profile_id/gsm_id/model_id)
# Filtering  : seeded probe sample of N_PROBE_GENES genes from
#              tx.gene_universe. No value filtering; NaN preserved.
# Aggregation: tx.collapse_replicates (median per cell line x method)
# Parameters : q_floor from chosen_params.json (sets array detection floors)
# Output     : figs/fig2_scale_harmonisation.{png,pdf}
# ======================================================================

def figure2(con, axes, colmaps, probe, out: Path):
    import matplotlib.patches as mpatches
    import matplotlib.colors as mc

    print("\n[Fig 2] method / scale harmonisation")
    cache = tx.fetch_chunk_by_table(con, axes, probe, colmaps)

    rows = []
    for ax in axes:
        X = tx.collapse_replicates(
            ax, tx.align_axis_rows(ax, cache.get(ax.table), probe))
        raw = X[np.isfinite(X)]
        if raw.size < 50:
            print(f"  skip {ax.name}: only {raw.size} finite probe values")
            continue

        V = tx.to_log2(ax, X)                     # linear -> log2(x+1)
        logv = V[np.isfinite(V)]

        # within-method standardisation, unstratified (global) reference
        med, mad = tx.nan_mad(V)
        mad = np.where(np.isfinite(mad) & (mad > tx.EPS), mad, np.nan)
        Z = (V - med) / mad
        zz = Z[np.isfinite(Z)]

        rows.append(dict(name=ax.name, source=ax.source,
                         scale=ax.scale_detected, units=ax.units,
                         floor=ax.floor, raw=raw, log=logv, z=zz,
                         n_models=len(ax.models)))

    if not rows:
        raise RuntimeError("no method produced >=50 finite probe values")

    rows.sort(key=lambda r: (r["source"], r["name"]))
    n = len(rows)
    xs = np.arange(1, n + 1)

    # ------------------------------------------------------------------
    # Palette: colorblind-safe (blue / amber / green triad, distinct hue
    # AND luminance so it still separates in grayscale print). Hatching
    # is layered on top as a second, redundant encoding of source so the
    # figure survives photocopying / B&W supplementary PDFs.
    # ------------------------------------------------------------------
    _SRC_COL = {
        "depmap":  "#2B6C8F",   # steel blue
        "geo":     "#C97B2E",   # warm amber
        "hpa_rna": "#4B7A57",   # muted forest green
    }
    _SRC_HATCH = {"depmap": "", "geo": "////", "hpa_rna": "...."}
    _REF_GOLD  = "#8A5A00"   # detection-floor reference line
    _REF_GREY  = "#3F3F3F"   # zero-median reference line

    def _col(src):
        return _SRC_COL.get(src, "#888888")

    def _edge(src):
        r_v, g_v, b_v = mc.to_rgb(_col(src))
        return mc.to_hex((r_v * 0.55, g_v * 0.55, b_v * 0.55))

    labels_wrapped = [r["source"].replace("_", " ").title() for r in rows]

    # contiguous source groups (rows are already sorted by source, name)
    # -- used for faint background banding + separators between sources
    groups, start, cur = [], 0, rows[0]["source"]
    for i in range(1, n):
        if rows[i]["source"] != cur:
            groups.append((start, i - 1, cur))
            start, cur = i, rows[i]["source"]
    groups.append((start, n - 1, cur))

    panels = [
        ("raw", "Native scale",                          "Expression value"),
        ("log", r"Scale detected + log$_2$(x+1)",        "Expression value"),
        ("z",   "Robust within-method standardisation",  "Robust z-score"),
    ]
    panel_letters = ["a", "b", "c"]

    _rc = {
        "font.family":       "serif",
        "font.serif":        ["Times New Roman", "Times", "DejaVu Serif"],
        "font.size":         13,
        "axes.labelsize":    14,
        "axes.titlesize":    15,
        "axes.titleweight":  "medium",
        "xtick.labelsize":   12,
        "ytick.labelsize":   12,
        "legend.fontsize":   12,
        "axes.linewidth":    0.6,
        "axes.edgecolor":    "#4a4a4a",
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "xtick.major.size":  2.5,
        "ytick.major.size":  2.5,
        "xtick.color":       "#4a4a4a",
        "ytick.color":       "#4a4a4a",
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "text.color":        "#1a1a1a",
        "axes.labelcolor":   "#1a1a1a",
    }

    with plt.rc_context(_rc):
        fig, axs = plt.subplots(1, 3, figsize=(14.5, 5.8),
                                gridspec_kw={"wspace": 0.34})

        for k, (key, title, ylabel) in enumerate(panels):
            axi = axs[k]
            data = [r[key] for r in rows]

            # faint per-source background bands + separators, drawn first
            for (s, e, src) in groups:
                axi.axvspan(xs[s] - 0.5, xs[e] + 0.5, color=_col(src),
                           alpha=0.07, zorder=0, lw=0)
            for (s, e, src) in groups[:-1]:
                axi.axvline(xs[e] + 0.5, color="#d5d5d5", lw=0.6, zorder=0)

            bp = axi.boxplot(
                data, positions=xs, widths=0.38, showfliers=False,
                patch_artist=True,
                medianprops=dict(color="#141414", lw=1.7),
                whiskerprops=dict(color="#555555", lw=0.9),
                capprops=dict(color="#555555", lw=0.9),
                boxprops=dict(lw=0.9),
            )
            for patch, r in zip(bp["boxes"], rows):
                patch.set_facecolor(_col(r["source"]))
                patch.set_alpha(0.50)
                patch.set_edgecolor(_edge(r["source"]))
                patch.set_hatch(_SRC_HATCH.get(r["source"], ""))
                patch.set_linewidth(0.9)

            # Panel (b): per-method detection floor segment + label once
            if key == "log":
                first_x, first_y = None, None
                for i, r in enumerate(rows):
                    f = r["floor"]
                    if f is not None and np.isfinite(f):
                        axi.plot([xs[i] - 0.29, xs[i] + 0.29], [f, f],
                                 color=_REF_GOLD, lw=1.1, zorder=4,
                                 solid_capstyle="round")
                        if first_x is None:
                            first_x, first_y = xs[i], f
                if first_y is not None:
                    axi.annotate("Detection floor",
                                xy=(first_x, first_y), xytext=(4, 4),
                                textcoords="offset points",
                                fontsize=9, color=_REF_GOLD,
                                style="italic", ha="left", va="bottom")

            # Panel (c): dashed zero line, unobtrusive corner label
            if key == "z":
                axi.axhline(0, color=_REF_GREY, lw=0.9,
                            ls=(0, (4, 3)), alpha=0.7, zorder=0)
                axi.text(0.02, 0.96, "Median = 0",
                         va="top", ha="left", fontsize=9,
                         color=_REF_GREY, style="italic",
                         transform=axi.transAxes)

            # Compact per-box annotation placed at whisker top
            for i, (r, d) in enumerate(zip(rows, data)):
                if len(d) < 4:
                    continue
                wtop = float(bp["whiskers"][2 * i + 1].get_ydata()[1])
                data_span = float(
                    np.nanpercentile(d, 95) - np.nanpercentile(d, 5))
                offset = max(data_span * 0.02, 1e-6)
                if key == "raw":
                    ann = f"n={r['n_models']}"
                elif key == "log":
                    ann = (f"IQR="
                           f"{float(np.subtract(*np.percentile(d, [75, 25]))):.1f}")
                else:
                    ann = f"med={float(np.median(d)):.2f}"
                axi.text(xs[i], wtop + offset, ann,
                         ha="center", va="bottom",
                         fontsize=8.5, color="#555555")

            # journal-style panel letter, set apart from the descriptive title
            axi.text(-0.14, 1.14, panel_letters[k], transform=axi.transAxes,
                     fontsize=16, fontweight="bold", va="top", ha="left")
            axi.set_title(title, loc="left", fontsize=15,
                          fontweight="medium", pad=4)
            axi.set_ylabel(ylabel, fontsize=14, labelpad=6)
            axi.set_xlim(0.35, n + 0.65)
            axi.set_xticks(xs)
            axi.set_xticklabels(labels_wrapped, rotation=0, ha="center",
                               linespacing=0.92, fontsize=12)
            axi.tick_params(axis="y", labelsize=12)
            axi.grid(axis="y", lw=0.4, color="#e6e6e6", zorder=0)
            axi.set_axisbelow(True)

        # Single source-colour legend below all three panels
        shown_sources = [s for s in ("depmap", "geo", "hpa_rna")
                         if any(r["source"] == s for r in rows)]
        legend_handles = [
            mpatches.Patch(facecolor=_col(s), alpha=0.55,
                           edgecolor=_edge(s), lw=0.9,
                           hatch=_SRC_HATCH.get(s, ""), label=s)
            for s in shown_sources
        ]
        fig.legend(handles=legend_handles,
                   loc="lower center", ncol=len(legend_handles),
                   fontsize=12, frameon=False,
                   bbox_to_anchor=(0.5, 0.035),
                   handlelength=1.6, handleheight=1.1,
                   title="Data source", title_fontsize=12)

        fig.suptitle(
            "Cross-Method Comparison of Expression Distributions",
            fontsize=21, fontweight="bold", y=0.975, x=0.5, ha="center")

        fig.subplots_adjust(top=0.84, bottom=0.20, left=0.075, right=0.97,
                            wspace=0.34)
        for ext in ("png", "pdf"):
            fig.savefig(out / f"fig2_scale_harmonisation.{ext}",
                        dpi=300 if ext == "png" else None,
                        facecolor="white", bbox_inches="tight")
        plt.close(fig)

    # the quantitative test of harmonisation
    spread_log = float(np.std([np.median(r["log"]) for r in rows]))
    spread_z = float(np.std([np.median(r["z"]) for r in rows]))
    note("fig2.sd_of_method_medians_log", round(spread_log, 4),
         "sd over methods of median(log2 value), panel (b)")
    note("fig2.sd_of_method_medians_z", round(spread_z, 4),
         "sd over methods of median(z), panel (c)")
    note("fig2.n_methods_plotted", len(rows), ">=50 finite probe values")
    note("fig2.scale_detected_counts",
         pd.Series([r["scale"] for r in rows]).value_counts().to_dict(),
         "MethodAxis.scale_detected")
    conflicts = [r["name"] for r in rows
                 if r["units"] in ("log", "linear") and r["scale"] != r["units"]]
    note("fig2.metadata_vs_data_scale_conflicts", conflicts or "none",
         "SOFT scale_hint != empirically detected scale")
    note("fig2.floors_by_method",
         {r["name"]: (round(float(r["floor"]), 3)
                      if r["floor"] is not None and np.isfinite(r["floor"])
                      else None) for r in rows},
         "MethodAxis.floor (log2 units)")


# ======================================================================
# FIGURE 3 - detection regimes
# ----------------------------------------------------------------------
# Input (preferred) : results/transcriptomics/regimes_prefilter.csv
#                     -- exact corpus-wide counts written by the scoring
#                     run before chunk_to_frame applies its keep filter.
# Input (fallback)  : expression tables, recomputed on the probe sample.
#
# WHY THE FALLBACK EXISTS: chunk_to_frame ends with
#     keep = (k_src >= min_k) & np.isfinite(z_lin)
# so below_detection and no_data cells never reach transcriptomics_z, and
# diagnostics.csv is built from the already-filtered frame. Without the
# prefilter file the withheld fraction is not recoverable and must be
# estimated. The basis actually used is recorded in the audit and printed
# in the figure title.
#
# Output : figs/fig3_detection_regimes.{png,pdf}
# ======================================================================

def figure3(con, axes, colmaps, lineage_s, probe, params, out: Path,
            results_dir: Path):
    import textwrap
    import matplotlib.patches as mpatches
    import matplotlib.transforms as mtransforms
    from matplotlib.ticker import FuncFormatter

    print("\n[Fig 3] detection regimes")
    rpath = results_dir / "regimes_prefilter.csv"

    if rpath.exists():
        basis = "corpus-wide (exact)"
        note("fig3.basis", basis, str(rpath))
        df = pd.read_csv(rpath)
        need = {"source", "method", "lineage", "status", "n"}
        if not need.issubset(df.columns):
            raise ValueError(f"{rpath} missing columns: {need - set(df.columns)}")
        dfm = df[["method", "source", "status", "n"]].copy()
        dfl = df[["lineage", "status", "n"]].copy()
    else:
        basis = f"probe sample, {len(probe)} genes"
        print(f"  {rpath} absent -- recomputing on {basis}.")
        note("fig3.basis", basis, "recomputed via tx.method_z, pre-filter")
        cache = tx.fetch_chunk_by_table(con, axes, probe, colmaps)
        m_parts, l_parts = [], []
        for ax in axes:
            X = tx.collapse_replicates(
                ax, tx.align_axis_rows(ax, cache.get(ax.table), probe))
            V = tx.to_log2(ax, X)
            lin = lineage_s.reindex(ax.models).fillna("unknown").to_numpy()
            _, S, _, _ = tx.method_z(V, lin, ax.floor, params, stratify=True)
            st_flat = pd.Series(S.ravel()).map(tx.STATUS).to_numpy()
            vc = pd.Series(st_flat).value_counts()
            m_parts.append(pd.DataFrame(dict(
                method=ax.name, source=ax.source,
                status=vc.index, n=vc.values)))
            lin_rep = np.repeat(lin[:, None], S.shape[1], axis=1).ravel()
            l_parts.append(pd.DataFrame(dict(
                lineage=lin_rep, status=st_flat, n=1)))
        dfm = pd.concat(m_parts, ignore_index=True)
        dfl = pd.concat(l_parts, ignore_index=True)

    overall = dfm.groupby("status").n.sum().reindex(REGIME_ORDER).fillna(0)
    total = float(overall.sum())
    if total <= 0:
        raise RuntimeError("no cells counted for Figure 3")
    pct = (100 * overall / total).round(2)

    # ------------------------------------------------------------------
    # Statuses that count as "usable" for downstream analysis vs. ones
    # that represent a withheld/uncertain cell. Adjust this set if your
    # status vocabulary differs -- everything in REGIME_ORDER not listed
    # here is treated as usable.
    # ------------------------------------------------------------------
    NOT_USABLE = {"below_detection", "no_data"}
    USABLE_ORDER = [s for s in REGIME_ORDER if s not in NOT_USABLE]
    NOTUSABLE_ORDER = [s for s in REGIME_ORDER if s in NOT_USABLE]

    # Blue (usable) / orange (not usable) -- the standard colorblind-safe
    # diverging pair, and it reads as "good vs. bad" without relying on
    # red-green, which a meaningful share of readers can't distinguish.
    _SEMANTIC = {
        "ok":              "#1B5E86",   # deep blue  - usable, high confidence
        "censored_ranked": "#7FB2D9",   # light blue - usable, lower confidence
        "below_detection": "#C1571F",   # burnt orange - not usable, low signal
        "no_data":         "#9A9A9A",   # neutral grey - not usable, missing
        "saturated":       "#8B3A3A",
        "ambiguous":       "#6B5B95",
    }
    _FALLBACK = ["#1B5E86", "#C1571F", "#4B7A57", "#8B5E83", "#6B6B6B"]
    REGIME_COL = {
        s: _SEMANTIC.get(s, _FALLBACK[i % len(_FALLBACK)])
        for i, s in enumerate(REGIME_ORDER)
    }
    _USABLE_TXT, _NOTUSABLE_TXT = "#123f5c", "#8a3f10"

    def _pretty(status):
        return status.replace("_", " ").title()

    def _wrap(label, width=16, max_lines=2):
        pretty = label.replace("_", " ").strip()
        lines = textwrap.wrap(pretty, width=width, break_long_words=False)
        if not lines:
            return pretty
        if len(lines) > max_lines:
            lines = lines[:max_lines]
            lines[-1] = lines[-1].rstrip() + "\u2026"
        return "\n".join(lines)

    def _diverging_barh(ax, prop_df, n_series, row_labels, show_n_labels=True):
        """Zero-centred diverging stacked bar: usable statuses extend
        right of 0, not-usable statuses extend left of 0."""
        n = len(prop_df)
        y = np.arange(n)
        right_cum = np.zeros(n)
        for s in USABLE_ORDER:
            vals = prop_df[s].to_numpy() if s in prop_df.columns else np.zeros(n)
            ax.barh(y, vals, left=right_cum, height=0.66,
                   color=REGIME_COL[s], edgecolor="#333333",
                   linewidth=0.5, zorder=3)
            right_cum += vals
        left_cum = np.zeros(n)
        for s in NOTUSABLE_ORDER:
            vals = prop_df[s].to_numpy() if s in prop_df.columns else np.zeros(n)
            ax.barh(y, -vals, left=-left_cum, height=0.66,
                   color=REGIME_COL[s], edgecolor="#333333",
                   linewidth=0.5, zorder=3)
            left_cum += vals

        ax.axvline(0, color="#2a2a2a", lw=0.9, zorder=4)
        ax.set_xlim(-1.02, 1.02)
        ax.set_ylim(-0.6, n - 0.4)
        ax.invert_yaxis()  # highest usable fraction at the top
        ax.set_yticks(y)
        ax.set_yticklabels([_wrap(lbl) for lbl in row_labels], fontsize=9.5)
        ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{abs(v) * 100:.0f}%"))
        ax.grid(axis="x", lw=0.3, color="#e6e6e6", zorder=0)
        ax.set_axisbelow(True)

        if show_n_labels:
            trans = mtransforms.blended_transform_factory(ax.transAxes, ax.transData)
            for i, n_val in enumerate(n_series):
                ax.text(1.03, y[i], f"n={int(n_val):,}", transform=trans,
                        fontsize=6.2, va="center", ha="left", color="#555555")

        ax.text(0.5, 1.0, "\u2190 Not usable", transform=ax.transAxes,
                fontsize=6.8, ha="left", va="bottom", color=_NOTUSABLE_TXT,
                style="italic")
        ax.text(0.5, 1.0, "Usable \u2192", transform=ax.transAxes,
                fontsize=6.8, ha="right", va="bottom", color=_USABLE_TXT,
                style="italic")

    _rc = {
        "font.family":       "sans-serif",
        "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":         11,
        "axes.labelsize":    11,
        "axes.titlesize":    11,
        "xtick.labelsize":   10,
        "ytick.labelsize":   10,
        "legend.fontsize":   10,
        "axes.linewidth":    0.6,
        "axes.edgecolor":    "#4a4a4a",
        "xtick.color":       "#4a4a4a",
        "ytick.color":       "#4a4a4a",
        "xtick.major.width": 0.5,
        "ytick.major.width": 0.5,
        "axes.spines.top":   False,
        "axes.spines.right": False,
        "text.color":        "#1a1a1a",
        "axes.labelcolor":   "#1a1a1a",
    }

    with plt.rc_context(_rc):
        fig, axs = plt.subplots(1, 3, figsize=(14.5, 5.4),
                                gridspec_kw={"width_ratios": [1, 1.35, 1.35],
                                            "wspace": 0.55})

        # ---------------- (a) overall allocation ----------------------
        order_a = USABLE_ORDER + NOTUSABLE_ORDER
        bar_colors_a = [REGIME_COL[s] for s in order_a]
        vals_a = [overall.get(s, 0) for s in order_a]
        axs[0].axvspan(-0.5, len(USABLE_ORDER) - 0.5, color="#1B5E86", alpha=0.05, zorder=0)
        axs[0].axvspan(len(USABLE_ORDER) - 0.5, len(order_a) - 0.5, color="#C1571F", alpha=0.05, zorder=0)
        axs[0].bar(range(len(order_a)), vals_a, color=bar_colors_a,
                  edgecolor="#333333", linewidth=0.8, width=0.62, zorder=3)
        axs[0].grid(axis="y", lw=0.3, color="#e6e6e6", zorder=0)
        axs[0].set_axisbelow(True)
        for i, s in enumerate(order_a):
            axs[0].text(i, vals_a[i], f"{int(vals_a[i]):,}\n{pct.get(s, 0):.1f}%",
                       ha="center", va="bottom", fontsize=7.3, linespacing=1.15)
        axs[0].set_ylim(0, max(vals_a) * 1.22)
        axs[0].set_xticks(range(len(order_a)))
        axs[0].set_xticklabels([_wrap(s, width=10) for s in order_a], fontsize=9.5)
        axs[0].yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{int(v):,}"))
        axs[0].set_ylabel("(Cell line, gene) pairs")
        axs[0].text(-0.20, 1.14, "a", transform=axs[0].transAxes,
                   fontsize=12, fontweight="bold", va="top", ha="left")
        axs[0].set_title("Overall allocation", loc="left", fontsize=11,
                        fontweight="medium", pad=12)
        axs[0].text(0.5, 1.0, f"{pct.get('ok', 0) + pct.get('censored_ranked', 0):.1f}% usable overall",
                   transform=axs[0].transAxes, ha="center", va="bottom",
                   fontsize=7.2, style="italic", color=_USABLE_TXT)

        # ---------------- (b) by profiling method ----------------------
        bym = (dfm.pivot_table(index="method", columns="status", values="n",
                               aggfunc="sum")
                  .reindex(columns=REGIME_ORDER).fillna(0))
        bym_n = bym.sum(1)
        bym_prop = bym.div(bym_n.replace(0, np.nan), axis=0)
        usable_frac_m = bym_prop[USABLE_ORDER].sum(1)
        order_m = usable_frac_m.sort_values(ascending=False).index
        bym_prop, bym_n = bym_prop.loc[order_m], bym_n.loc[order_m]

        # short codes on the axis instead of full method names -- the
        # names go into a small legend inside the panel instead
        codes_m = [f"M{i + 1}" for i in range(len(bym_prop))]
        _diverging_barh(axs[1], bym_prop, bym_n.values, codes_m,
                       show_n_labels=False)
        legend_lines = [f"{c}   {name}   (n={int(bym_n[name]):,})"
                        for c, name in zip(codes_m, bym_prop.index)]
        axs[1].text(0.98, 0.03, "\n".join(legend_lines),
                   transform=axs[1].transAxes, fontsize=6.5,
                   ha="right", va="bottom", color="#333333", linespacing=1.6,
                   bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                            edgecolor="#cfcfcf", linewidth=0.7, alpha=0.95))
        axs[1].set_xlabel("Share of pairs")
        axs[1].text(-0.14, 1.14, "b", transform=axs[1].transAxes,
                   fontsize=12, fontweight="bold", va="top", ha="left")
        axs[1].set_title("By profiling method", loc="left", fontsize=11,
                        fontweight="medium", pad=14)

        # ---------------- (c) by cancer lineage -------------------------
        byl_all = (dfl.pivot_table(index="lineage", columns="status", values="n",
                                   aggfunc="sum")
                      .reindex(columns=REGIME_ORDER).fillna(0))
        byl_n_all = byl_all.sum(1)
        top_idx = byl_n_all.sort_values(ascending=False).index[:18]
        byl = byl_all.loc[top_idx]
        byl_n = byl_n_all.loc[top_idx]
        byl_prop = byl.div(byl_n.replace(0, np.nan), axis=0)
        usable_frac_l = byl_prop[USABLE_ORDER].sum(1)
        order_l = usable_frac_l.sort_values(ascending=False).index
        byl_prop, byl_n = byl_prop.loc[order_l], byl_n.loc[order_l]

        _diverging_barh(axs[2], byl_prop, byl_n.values, list(byl_prop.index))
        axs[2].set_xlabel("Share of pairs")
        axs[2].text(-0.22, 1.14, "c", transform=axs[2].transAxes,
                   fontsize=12, fontweight="bold", va="top", ha="left")
        axs[2].set_title("By cancer lineage (18 largest)", loc="left",
                        fontsize=11, fontweight="medium", pad=14)

        # ---------------- match row thickness between (b) and (c) -------
        fig.canvas.draw()
        pos_b, pos_c = axs[1].get_position(), axs[2].get_position()
        row_h = pos_c.height / len(byl_prop)
        new_h_b = row_h * len(bym_prop)
        axs[1].set_position([pos_b.x0, pos_b.y1 - new_h_b, pos_b.width, new_h_b])

        # ---------------- shared status legend + title -------------------
        legend_handles = [
            mpatches.Patch(facecolor=REGIME_COL[s], edgecolor="#333333",
                           linewidth=0.6, label=_pretty(s))
            for s in REGIME_ORDER
        ]
        fig.legend(handles=legend_handles, loc="lower center",
                  ncol=len(legend_handles), fontsize=7.8, frameon=False,
                  bbox_to_anchor=(0.5, -0.02), title="Detection status",
                  title_fontsize=7.8)

        fig.suptitle(
            "Figure 3. Detection-Regime Allocation Across (Cell Line, Gene) Pairs",
            fontsize=12.5, fontweight="semibold", x=0.02, ha="left", y=1.05)
        fig.text(0.02, 1.005,
                 f"Basis: {basis}   "
                 rf"($q_{{upper}}$={params.upper_cut:g}, "
                 rf"$q_{{lower}}$={params.lower_cut:g})   "
                 f"N = {int(total):,} pairs evaluated",
                 fontsize=8, color="#555555", ha="left")

        fig.subplots_adjust(top=0.83, bottom=0.16, wspace=0.55)
        for ext in ("png", "pdf"):
            fig.savefig(out / f"fig3_detection_regimes.{ext}",
                       dpi=300 if ext == "png" else None,
                       facecolor="white", bbox_inches="tight")
        plt.close(fig)

    note("fig3.n_cells_evaluated", int(total), "sum of regime counts")
    for s in REGIME_ORDER:
        note(f"fig3.count.{s}", int(overall[s]), "tx.method_z status codes")
        note(f"fig3.pct.{s}", float(pct[s]), "share of evaluated cells")


# ======================================================================
# FIGURE 4 - calibration surface and worked biological example
# ----------------------------------------------------------------------
# Panel (a)
#   Input     : results/transcriptomics/grid_search.csv, chosen_params.json
#   Columns   : floor_q, n0, upper_cut, lower_cut, rho
#   Filtering : none; the full grid is used for the 1-SE computation, and
#               the heatmap slices at the selected q_floor and n0
# Panel (b)
#   Input     : transcriptomics_z, restricted to the worked-example gene
#   Columns   : ensg, model_id, lineage, z_lineage, z_global, k_src,
#               k_methods, det_frac, status, scale_source, per-source z
#   Filtering : gene chosen from WORKED_EXAMPLE_CANDIDATES, a list fixed
#               before any score was inspected; widest-coverage candidate
#               wins. No selection on effect size.
# Output      : figs/fig4_calibration_example.{png,pdf}
# ======================================================================

def figure4(con, table: str, results_dir: Path, params: tx.Params, out: Path):
    import matplotlib.patches as mpatches
    import matplotlib.colors as mcolors
    from matplotlib.ticker import FuncFormatter

    print("\n[Fig 4] calibration surface and worked example")
    gpath = results_dir / "grid_search.csv"
    grid = pd.read_csv(gpath) if gpath.exists() else None
    if grid is None:
        print(f"  MISSING {gpath} -- panel (a) will be blank.")
        note("fig4a", "UNAVAILABLE", f"{gpath} absent; rerun --calibrate-only")

    # ---- pick the worked-example gene, by coverage, not by outcome ----
    q = ",".join("?" * len(WORKED_EXAMPLE_CANDIDATES))
    present = con.execute(f"""
        SELECT ensg, count(*) AS n, count(DISTINCT lineage) AS n_lin
        FROM "{table}" WHERE ensg IN ({q})
        GROUP BY 1 ORDER BY n DESC
    """, list(WORKED_EXAMPLE_CANDIDATES)).fetchdf()

    df = None
    if present.empty:
        print("  none of the pre-declared candidate genes are present.")
        note("fig4b", "UNAVAILABLE", "no candidate marker gene in table")
    else:
        ensg = present.ensg.iloc[0]
        label = WORKED_EXAMPLE_CANDIDATES[ensg]
        print(f"  selected {ensg} = {label} "
              f"({int(present.n.iloc[0])} rows, {int(present.n_lin.iloc[0])} lineages)")
        df = con.execute(f'SELECT * FROM "{table}" WHERE ensg = ?',
                         [ensg]).fetchdf()

    # ------------------------------------------------------------------
    # Shared palette: continues the blue (good/primary) / amber (accent)
    # language from Figs 2-3 instead of introducing a new hue family.
    # ------------------------------------------------------------------
    _BLUE, _AMBER, _INK = "#1B5E86", "#C1571F", "#232323"
    _CMAP = mcolors.LinearSegmentedColormap.from_list(
        "calib_cmap", ["#C1571F", "#F4EFE6", "#1B5E86"])

    def _stat_box(ax, text, loc="lower left"):
        anchors = {"lower left": (0.02, 0.02, "bottom", "left"),
                   "lower right": (0.98, 0.02, "bottom", "right"),
                   "upper left": (0.02, 0.98, "top", "left")}
        x, y, va, ha = anchors[loc]
        ax.text(x, y, text, transform=ax.transAxes, fontsize=6.8,
                va=va, ha=ha, color="#333333", linespacing=1.6,
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                         edgecolor="#cfcfcf", linewidth=0.7, alpha=0.94))

    def _placeholder(ax, msg, letter, title):
        ax.set_facecolor("#fafafa")
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linestyle((0, (4, 3)))
            spine.set_color("#c9c9c9")
        ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=9.5,
                color="#9a9a9a", style="italic", transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        ax.text(-0.05, 1.1, letter, transform=ax.transAxes, fontsize=12,
                fontweight="bold", va="top", ha="left")
        ax.set_title(title, loc="left", fontsize=11, fontweight="medium", pad=12)

    _rc = {
        "font.family":       "sans-serif",
        "font.sans-serif":   ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size":         11,
        "axes.labelsize":    11,
        "axes.titlesize":    11,
        "xtick.labelsize":   10,
        "ytick.labelsize":   10,
        "legend.fontsize":   10,
        "axes.linewidth":    0.6,
        "axes.edgecolor":    "#4a4a4a",
        "xtick.color":       "#4a4a4a",
        "ytick.color":       "#4a4a4a",
        "text.color":        "#1a1a1a",
        "axes.labelcolor":   "#1a1a1a",
    }

    with plt.rc_context(_rc):
        fig, axs = plt.subplots(1, 2, figsize=(13.0, 5.2),
                                gridspec_kw={"wspace": 0.32})

        # ---------------- panel (a): calibration ----------------
        if grid is not None and len(grid):
            best = float(grid.rho.max())
            se = float(grid.rho.std(ddof=1) / np.sqrt(len(grid))) if len(grid) > 1 else 0.0
            within = grid[grid.rho >= best - se]

            sub = grid[(grid.floor_q == params.floor_q) & (grid.n0 == params.n0)]
            if len(sub) >= 4:
                piv = sub.pivot_table(index="upper_cut", columns="lower_cut",
                                      values="rho")
                im = axs[0].imshow(piv.values, cmap=_CMAP, origin="lower",
                                   aspect="auto",
                                   vmin=float(np.nanmin(piv.values)),
                                   vmax=float(np.nanmax(piv.values)))
                axs[0].set_xticks(range(len(piv.columns)),
                                  [f"{c:g}" for c in piv.columns])
                axs[0].set_yticks(range(len(piv.index)),
                                  [f"{i:g}" for i in piv.index])

                # highlight cells within 1 SE of the global max -- adds
                # information (the plateau shape) without adding colors
                mask = piv.values >= (best - se)
                for (i, j), ok in np.ndenumerate(mask):
                    if ok and np.isfinite(piv.values[i, j]):
                        axs[0].add_patch(plt.Rectangle(
                            (j - 0.5, i - 0.5), 1, 1, fill=False,
                            edgecolor="white", linewidth=1.3, zorder=4))

                cols, idx = list(piv.columns), list(piv.index)
                handles = []
                if params.lower_cut in cols and params.upper_cut in idx:
                    h1 = axs[0].scatter(
                        [cols.index(params.lower_cut)], [idx.index(params.upper_cut)],
                        marker="*", s=280, facecolor="white", edgecolor=_INK,
                        linewidth=1.1, zorder=6, label="Selected (1-SE)")
                    handles.append(h1)
                bi = grid.rho.idxmax()
                if (grid.loc[bi, "floor_q"] == params.floor_q
                        and grid.loc[bi, "n0"] == params.n0
                        and grid.loc[bi, "lower_cut"] in cols
                        and grid.loc[bi, "upper_cut"] in idx):
                    h2 = axs[0].scatter(
                        [cols.index(grid.loc[bi, "lower_cut"])],
                        [idx.index(grid.loc[bi, "upper_cut"])],
                        marker="D", s=90, facecolor="none", edgecolor=_INK,
                        linewidth=1.3, zorder=6, label="Grid maximum")
                    handles.append(h2)
                if handles:
                    axs[0].legend(handles=handles, fontsize=6.8, loc="lower left",
                                 frameon=False, labelcolor=_INK)

                cbar = fig.colorbar(im, ax=axs[0], fraction=0.046, pad=0.03)
                cbar.set_label(r"Held-out Spearman $\rho$", fontsize=10)
                cbar.ax.tick_params(labelsize=9)
                cbar.outline.set_linewidth(0.5)

                _stat_box(axs[0],
                         f"Max $\\rho$ = {best:.3f}\n"
                         f"SE = {se:.3f}\n"
                         f"{len(within)}/{len(grid)} within 1 SE",
                         loc="lower right")

            axs[0].set_xlabel(r"$q_{\mathrm{lower}}$")
            axs[0].set_ylabel(r"$q_{\mathrm{upper}}$")
            axs[0].text(-0.14, 1.12, "a", transform=axs[0].transAxes,
                       fontsize=12, fontweight="bold", va="top", ha="left")
            axs[0].set_title("Calibration surface", loc="left", fontsize=11,
                            fontweight="medium", pad=16)
            axs[0].text(0.0, 1.03,
                       rf"$q_{{\mathrm{{floor}}}}$={params.floor_q:g}, "
                       rf"$n_0$={params.n0:g}",
                       transform=axs[0].transAxes, fontsize=7, color="#666666",
                       style="italic", ha="left", va="bottom")

            note("calib.n_combinations", int(len(grid)), f"{gpath} rows")
            note("calib.rho_max", round(best, 4), f"{gpath} max rho")
            note("calib.se", round(se, 4), "sd(rho)/sqrt(n) over grid")
            note("calib.n_within_1se", int(len(within)), "rho >= max - SE")
            sel = grid[(grid.floor_q == params.floor_q) & (grid.n0 == params.n0)
                       & (grid.upper_cut == params.upper_cut)
                       & (grid.lower_cut == params.lower_cut)]
            note("calib.rho_selected",
                 round(float(sel.rho.iloc[0]), 4) if len(sel) else "NOT IN GRID",
                 "rho at the chosen_params.json setting")
            n0_best = grid.groupby("n0").rho.max()
            note("calib.n0_rho_range", round(float(n0_best.max() - n0_best.min()), 4),
                 "spread of best rho across n0; small => n0 unidentified")
        else:
            _placeholder(axs[0], "grid_search.csv not available", "a",
                        "Calibration surface")

        # ---------------- panel (b): worked example ----------------
        if df is not None and len(df):
            is_h = df.lineage.astype(str).str.contains(HAEM_PATTERN, case=False,
                                                       na=False)
            axs[1].scatter(df.z_global[~is_h], df.z_lineage[~is_h], s=7,
                          alpha=0.32, color="#A9AFB6", edgecolor="none",
                          label="Other lineages", zorder=2)
            axs[1].scatter(df.z_global[is_h], df.z_lineage[is_h], s=18,
                          alpha=0.88, color=_AMBER, edgecolor="#7a3712",
                          linewidth=0.3, label="Haematological", zorder=3)
            lo = float(np.nanmin([df.z_global.min(), df.z_lineage.min()]))
            hi = float(np.nanmax([df.z_global.max(), df.z_lineage.max()]))
            axs[1].plot([lo, hi], [lo, hi], ls=(0, (4, 3)), lw=0.9,
                       color="#888888", zorder=1, label="y = x")
            axs[1].axhline(0, color="#dddddd", lw=0.7, zorder=0)
            axs[1].axvline(0, color="#dddddd", lw=0.7, zorder=0)
            axs[1].grid(lw=0.3, color="#eeeeee", zorder=0)
            axs[1].set_axisbelow(True)
            axs[1].set_xlabel(r"$z_{\mathrm{global}}$")
            axs[1].set_ylabel(r"$z_{\mathrm{lineage}}$")
            axs[1].legend(fontsize=6.8, frameon=False, loc="upper left",
                         markerscale=1.3)

            haem = df[is_h]
            axs[1].text(-0.12, 1.12, "b", transform=axs[1].transAxes,
                       fontsize=12, fontweight="bold", va="top", ha="left")
            axs[1].set_title(f"Worked example: {label}", loc="left",
                            fontsize=11, fontweight="medium", pad=16)
            axs[1].text(0.0, 1.03,
                       f"{len(df)} cell lines, {int(is_h.sum())} haematological",
                       transform=axs[1].transAxes, fontsize=7, color="#666666",
                       style="italic", ha="left", va="bottom")
            _stat_box(axs[1],
                     f"Haem. median $z_{{global}}$ = {haem.z_global.median():.2f}\n"
                     f"Haem. median $z_{{lineage}}$ = {haem.z_lineage.median():.2f}",
                     loc="lower right")

            note("fig4b.gene", f"{label} ({ensg})",
                 "pre-declared candidate list; widest coverage")
            note("fig4b.n_cell_lines", int(len(df)), f"{table} rows for gene")
            note("fig4b.n_haematological", int(is_h.sum()),
                 f"lineage matches /{HAEM_PATTERN}/i")
            note("fig4b.median_z_lineage_haem",
                 round(float(haem.z_lineage.median()), 3) if len(haem) else "n/a",
                 "median z_lineage, haematological lines")
            note("fig4b.median_z_global_haem",
                 round(float(haem.z_global.median()), 3) if len(haem) else "n/a",
                 "median z_global, haematological lines")
            note("fig4b.mean_k_src", round(float(df.k_src.mean()), 2), table)
            note("fig4b.mean_k_methods", round(float(df.k_methods.mean()), 2), table)
            note("fig4b.status_counts", df.status.value_counts().to_dict(), table)
            note("fig4b.scale_source_counts",
                 df.scale_source.value_counts().to_dict(), table)
            for c in ("z_hpa_rna", "z_depmap", "z_geo"):
                if c in df.columns:
                    note(f"fig4b.median_{c}",
                         round(float(df[c].median()), 3)
                         if df[c].notna().any() else "n/a", table)
        else:
            _placeholder(axs[1], "Worked-example gene not available", "b",
                        "Biological validation example")

        fig.suptitle(
            "Figure 4. Parameter Calibration and Biological Validation",
            fontsize=12.5, fontweight="semibold", x=0.02, ha="left", y=1.06)

        fig.subplots_adjust(top=0.84, wspace=0.32)
        for ext in ("png", "pdf"):
            fig.savefig(out / f"fig4_calibration_example.{ext}",
                       dpi=300 if ext == "png" else None,
                       facecolor="white", bbox_inches="tight")
        plt.close(fig)

def figure4(con, table: str, results_dir: Path, params: tx.Params, out: Path):
    """Variant: one sensitivity curve per parameter instead of a 2D slice."""
    print("\n[Fig 4] calibration sensitivity and worked example")
    gpath = results_dir / "grid_search.csv"
    grid = pd.read_csv(gpath) if gpath.exists() else None
    if grid is None:
        print(f"  MISSING {gpath} -- sensitivity panels will be blank.")
        note("fig4a", "UNAVAILABLE", f"{gpath} absent; rerun --calibrate-only")

    q = ",".join("?" * len(WORKED_EXAMPLE_CANDIDATES))
    present = con.execute(f"""
        SELECT ensg, count(*) AS n, count(DISTINCT lineage) AS n_lin
        FROM "{table}" WHERE ensg IN ({q})
        GROUP BY 1 ORDER BY n DESC
    """, list(WORKED_EXAMPLE_CANDIDATES)).fetchdf()

    df = None
    if present.empty:
        print("  none of the pre-declared candidate genes are present.")
        note("fig4b", "UNAVAILABLE", "no candidate marker gene in table")
    else:
        ensg = present.ensg.iloc[0]
        label = WORKED_EXAMPLE_CANDIDATES[ensg]
        print(f"  selected {ensg} = {label} "
              f"({int(present.n.iloc[0])} rows, {int(present.n_lin.iloc[0])} lineages)")
        df = con.execute(f'SELECT * FROM "{table}" WHERE ensg = ?', [ensg]).fetchdf()

    _BLUE, _AMBER, _INK = "#1B5E86", "#C1571F", "#232323"

    def _stat_box(ax, text, loc="lower right"):
        anchors = {"lower left": (0.02, 0.02, "bottom", "left"),
                   "lower right": (0.98, 0.02, "bottom", "right"),
                   "upper left": (0.02, 0.98, "top", "left")}
        x, y, va, ha = anchors[loc]
        ax.text(x, y, text, transform=ax.transAxes, fontsize=6.8,
                va=va, ha=ha, color="#333333", linespacing=1.6,
                bbox=dict(boxstyle="round,pad=0.4", facecolor="white",
                         edgecolor="#cfcfcf", linewidth=0.7, alpha=0.94))

    def _placeholder(ax, msg, letter, title):
        ax.set_facecolor("#fafafa")
        for spine in ax.spines.values():
            spine.set_visible(True)
            spine.set_linestyle((0, (4, 3)))
            spine.set_color("#c9c9c9")
        ax.text(0.5, 0.5, msg, ha="center", va="center", fontsize=9.5,
                color="#9a9a9a", style="italic", transform=ax.transAxes)
        ax.set_xticks([]); ax.set_yticks([])
        ax.text(-0.08, 1.16, letter, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="top", ha="left")
        ax.set_title(title, loc="left", fontsize=11, fontweight="medium", pad=10)

    def _sensitivity(ax, grid, param, fixed, best, se, letter, pretty_name, show_ylabel):
        chosen_val = fixed[param]
        other = {k: v for k, v in fixed.items() if k != param}
        mask = np.ones(len(grid), dtype=bool)
        for k, v in other.items():
            mask &= np.isclose(grid[k], v)
        sl = grid[mask].sort_values(param)
        profiled = len(sl) < 3
        if profiled:
            sl = grid.groupby(param, as_index=False).rho.max().sort_values(param)

        ax.axhspan(best - se, best, color=_BLUE, alpha=0.08, zorder=0)
        ax.axhline(best - se, color=_BLUE, lw=0.6, ls=(0, (3, 2)), alpha=0.55, zorder=1)
        ax.axvline(chosen_val, color=_INK, lw=0.6, ls=(0, (1, 1.5)), alpha=0.6, zorder=1)
        ax.plot(sl[param], sl.rho, color=_BLUE, lw=1.4, marker="o", markersize=3, zorder=3)
        if chosen_val in sl[param].values:
            yv = float(sl.loc[sl[param] == chosen_val, "rho"].iloc[0])
            ax.scatter([chosen_val], [yv], marker="*", s=140, facecolor="white",
                      edgecolor=_INK, linewidth=1.0, zorder=5)
        ax.grid(lw=0.3, color="#eeeeee", zorder=0)
        ax.set_axisbelow(True)
        ax.set_xlabel(pretty_name, fontsize=11)
        if show_ylabel:
            ax.set_ylabel(r"Held-out $\rho$", fontsize=11)
        ax.text(-0.10, 1.28, letter, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="top", ha="left")
        tag = " (profiled)" if profiled else ""
        ax.set_title(f"vs. {pretty_name}{tag}", loc="left", fontsize=11, pad=14)

    _rc = {
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 11, "axes.labelsize": 11, "axes.titlesize": 12,
        "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 10,
        "axes.linewidth": 0.6, "axes.edgecolor": "#4a4a4a",
        "xtick.color": "#4a4a4a", "ytick.color": "#4a4a4a",
        "text.color": "#1a1a1a", "axes.labelcolor": "#1a1a1a",
    }

    with plt.rc_context(_rc):
        fig = plt.figure(figsize=(14.5, 7.3))
        gs = fig.add_gridspec(2, 4, height_ratios=[1, 1.35], hspace=0.65, wspace=0.28)
        top_axs = [fig.add_subplot(gs[0, i]) for i in range(4)]
        for a in top_axs[1:]:
            a.sharey(top_axs[0])
        ax_b = fig.add_subplot(gs[1, :])

        if grid is not None and len(grid):
            best = float(grid.rho.max())
            se = float(grid.rho.std(ddof=1) / np.sqrt(len(grid))) if len(grid) > 1 else 0.0
            within = grid[grid.rho >= best - se]

            fixed = {"floor_q": params.floor_q, "n0": params.n0,
                    "upper_cut": params.upper_cut, "lower_cut": params.lower_cut}
            pretty = {"floor_q": r"$q_{\mathrm{floor}}$", "n0": "$n_0$",
                     "upper_cut": r"$q_{\mathrm{upper}}$", "lower_cut": r"$q_{\mathrm{lower}}$"}
            letters = ["a", "b", "c", "d"]
            for i, (p, ax) in enumerate(zip(fixed.keys(), top_axs)):
                _sensitivity(ax, grid, p, fixed, best, se, letters[i], pretty[p], i == 0)

            pad = max((grid.rho.max() - grid.rho.min()) * 0.15, 1e-3)
            top_axs[0].set_ylim(grid.rho.min() - pad, grid.rho.max() + pad)
            _stat_box(top_axs[3],
                     f"Max $\\rho$={best:.3f}\nSE={se:.3f}\n{len(within)}/{len(grid)} within 1 SE",
                     loc="upper left")

            note("calib.n_combinations", int(len(grid)), f"{gpath} rows")
            note("calib.rho_max", round(best, 4), f"{gpath} max rho")
            note("calib.se", round(se, 4), "sd(rho)/sqrt(n) over grid")
            note("calib.n_within_1se", int(len(within)), "rho >= max - SE")
            sel = grid[(grid.floor_q == params.floor_q) & (grid.n0 == params.n0)
                       & (grid.upper_cut == params.upper_cut)
                       & (grid.lower_cut == params.lower_cut)]
            note("calib.rho_selected",
                 round(float(sel.rho.iloc[0]), 4) if len(sel) else "NOT IN GRID",
                 "rho at the chosen_params.json setting")
            n0_best = grid.groupby("n0").rho.max()
            note("calib.n0_rho_range", round(float(n0_best.max() - n0_best.min()), 4),
                 "spread of best rho across n0; small => n0 unidentified")
        else:
            for i, l in enumerate(["a", "b", "c", "d"]):
                _placeholder(top_axs[i], "n/a", l, "")

        if df is not None and len(df):
            is_h = df.lineage.astype(str).str.contains(HAEM_PATTERN, case=False, na=False)
            ax_b.scatter(df.z_global[~is_h], df.z_lineage[~is_h], s=7,
                        alpha=0.32, color="#A9AFB6", edgecolor="none",
                        label="Other lineages", zorder=2)
            ax_b.scatter(df.z_global[is_h], df.z_lineage[is_h], s=18,
                        alpha=0.88, color=_AMBER, edgecolor="#7a3712",
                        linewidth=0.3, label="Haematological", zorder=3)
            lo = float(np.nanmin([df.z_global.min(), df.z_lineage.min()]))
            hi = float(np.nanmax([df.z_global.max(), df.z_lineage.max()]))
            ax_b.plot([lo, hi], [lo, hi], ls=(0, (4, 3)), lw=0.9,
                     color="#888888", zorder=1, label="y = x")
            ax_b.axhline(0, color="#dddddd", lw=0.7, zorder=0)
            ax_b.axvline(0, color="#dddddd", lw=0.7, zorder=0)
            ax_b.grid(lw=0.3, color="#eeeeee", zorder=0)
            ax_b.set_axisbelow(True)
            ax_b.set_xlabel(r"$z_{\mathrm{global}}$")
            ax_b.set_ylabel(r"$z_{\mathrm{lineage}}$")
            ax_b.legend(fontsize=7, frameon=False, loc="upper left", markerscale=1.3)

            haem = df[is_h]
            ax_b.text(-0.055, 1.10, "e", transform=ax_b.transAxes, fontsize=12,
                     fontweight="bold", va="top", ha="left")
            ax_b.set_title(f"Worked example: {label}", loc="left", fontsize=11,
                          fontweight="medium", pad=14)
            ax_b.text(0.0, 1.02, f"{len(df)} cell lines, {int(is_h.sum())} haematological",
                     transform=ax_b.transAxes, fontsize=7.2, color="#666666",
                     style="italic", ha="left", va="bottom")
            _stat_box(ax_b,
                     f"Haem. median $z_{{global}}$ = {haem.z_global.median():.2f}\n"
                     f"Haem. median $z_{{lineage}}$ = {haem.z_lineage.median():.2f}")

            note("fig4b.gene", f"{label} ({ensg})", "pre-declared candidate list; widest coverage")
            note("fig4b.n_cell_lines", int(len(df)), f"{table} rows for gene")
            note("fig4b.n_haematological", int(is_h.sum()), f"lineage matches /{HAEM_PATTERN}/i")
            note("fig4b.median_z_lineage_haem",
                 round(float(haem.z_lineage.median()), 3) if len(haem) else "n/a", table)
            note("fig4b.median_z_global_haem",
                 round(float(haem.z_global.median()), 3) if len(haem) else "n/a", table)
            note("fig4b.mean_k_src", round(float(df.k_src.mean()), 2), table)
            note("fig4b.mean_k_methods", round(float(df.k_methods.mean()), 2), table)
            note("fig4b.status_counts", df.status.value_counts().to_dict(), table)
            note("fig4b.scale_source_counts", df.scale_source.value_counts().to_dict(), table)
            for c in ("z_hpa_rna", "z_depmap", "z_geo"):
                if c in df.columns:
                    note(f"fig4b.median_{c}",
                         round(float(df[c].median()), 3) if df[c].notna().any() else "n/a", table)
        else:
            _placeholder(ax_b, "Worked-example gene not available", "e", "Biological validation example")

        fig.suptitle("Figure 4. Parameter Calibration Sensitivity and Biological Validation",
                     fontsize=12.5, fontweight="semibold", x=0.02, ha="left", y=0.99)
        for ext in ("png", "pdf"):
            fig.savefig(out / f"fig4_calibration_example.{ext}",
                       dpi=300 if ext == "png" else None,
                       facecolor="white", bbox_inches="tight")
        plt.close(fig)


from scipy.stats import norm

PANELS = [
    ("z_global",  r"$z_{\mathrm{global}}$",  "#4878A8"),
    ("z_lineage", r"$z_{\mathrm{lineage}}$", "#D08442"),
    ("z_hpa_rna", r"$z_{\mathrm{HPA}}$",     "#6A9A6A"),
    ("z_depmap",  r"$z_{\mathrm{DepMap}}$",  "#9A7AA8"),
    ("z_geo",     r"$z_{\mathrm{GEO}}$",     "#E0A85C"),
]

TRIM_Q = (0.1, 99.9)
N_BINS = 90


def parse_args():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", default="db/celllineselector.duckdb")
    ap.add_argument("--table", default="transcriptomics_z")
    ap.add_argument("--sample", type=int, default=500_000,
                    help="reservoir sample size; 0 reads the whole table")
    ap.add_argument("--min-k", type=int, default=None,
                    help="restrict to rows with k_src >= this")
    ap.add_argument("--status", default=None,
                    help="restrict to one status, e.g. ok")
    ap.add_argument("--out", default="results/figures")
    ap.add_argument("--seed", type=int, default=0)
    return ap.parse_args()


def load(args):
    con = duckdb.connect(args.db, read_only=True)
    try:
        have = set(con.execute(
            f"PRAGMA table_info('{args.table}')").fetchdf()["name"])
        cols = [c for c, _, _ in PANELS if c in have]
        missing = [c for c, _, _ in PANELS if c not in have]
        if missing:
            print(f"[warn] absent from {args.table}, panels skipped: {missing}")
        if not cols:
            raise RuntimeError(f"no z columns found in {args.table}")

        where = []
        if args.min_k is not None:
            where.append(f"k_src >= {int(args.min_k)}")
        if args.status:
            where.append(f"status = '{args.status}'")
        clause = f"WHERE {' AND '.join(where)}" if where else ""
        smp = (f"USING SAMPLE {args.sample} ROWS (reservoir, {args.seed})"
               if args.sample else "")

        df = con.execute(
            f'SELECT {", ".join(cols)} FROM "{args.table}" {clause} {smp}'
        ).fetchdf()
        n_total = con.execute(
            f'SELECT count(*) FROM "{args.table}" {clause}').fetchone()[0]
    finally:
        con.close()
    return df, cols, n_total


def panel(ax, v_all, label, colour):
    v_all = v_all[np.isfinite(v_all)]
    if v_all.size < 100:
        ax.text(0.5, 0.5, "insufficient data", ha="center", va="center",
                transform=ax.transAxes, fontsize=8, color="0.45")
        ax.set_xlabel(label)
        return

    lo, hi = np.percentile(v_all, TRIM_Q)
    v = v_all[(v_all >= lo) & (v_all <= hi)]          # trim, do not clip
    mu, sd = float(v.mean()), float(v.std(ddof=1))
    frac0 = float(np.mean(v_all == 0.0))

    # --- background tint + cumulative sigma bands (outer → inner) ----------
    ax.set_facecolor("#f7f7f9")
    for k, alpha in [(3, 0.09), (2, 0.13), (1, 0.18)]:
        ax.axvspan(max(lo, mu - k * sd), min(hi, mu + k * sd),
                   color=colour, alpha=alpha, zorder=0, lw=0)

    # --- histogram -----------------------------------------------------------
    ax.hist(v, bins=N_BINS, density=True, color=colour,
            edgecolor="white", linewidth=0.25, alpha=0.78, zorder=2)

    # --- normal reference curve ----------------------------------------------
    x = np.linspace(lo, hi, 400)
    ax.plot(x, norm.pdf(x, mu, sd), "--", color="#1a1a1a", lw=1.5,
            alpha=0.70, label="Normal ref", zorder=5)

    # --- +/-1sigma, +/-2sigma, +/-3sigma dashed lines + top labels ----------
    xform = ax.get_xaxis_transform()   # x=data coords, y=axes [0,1]
    for k, alpha, lw in [(1, 0.85, 1.0), (2, 0.55, 0.85), (3, 0.32, 0.65)]:
        for sign in (-1, +1):
            xv = mu + sign * k * sd
            if lo < xv < hi:
                ax.axvline(xv, color=colour, lw=lw, ls=(0, (5, 3)),
                           alpha=alpha, zorder=4)
        xpos = mu + k * sd
        if lo < xpos < hi:
            ax.text(xpos, 1.01, f"$+{k}\\sigma$", transform=xform,
                    ha="center", va="bottom", fontsize=6.0,
                    color=colour, alpha=min(1.0, alpha + 0.15))

    # --- mean line -----------------------------------------------------------
    ax.axvline(mu, color="#1a1a1a", lw=1.6, ls="-", zorder=6, label="Mean")
    ax.text(mu, 1.01, r"$\mu$", transform=xform,
            ha="center", va="bottom", fontsize=7.5,
            color="#1a1a1a", fontweight="bold")

    # --- z=0 reference -------------------------------------------------------
    ax.axvline(0.0, color="#aaaaaa", lw=0.8, ls=":", zorder=3)

    # --- stat box ------------------------------------------------------------
    txt = (f"n\u2009=\u2009{v_all.size:,}\n"
           f"$\\mu$\u2009=\u2009{mu:.2f}\n"
           f"$\\sigma$\u2009=\u2009{sd:.2f}")
    if frac0 > 0.002:
        txt += f"\nP(z=0)\u2009=\u2009{100 * frac0:.1f}%"
    ax.text(0.97, 0.96, txt, transform=ax.transAxes, ha="right", va="top",
            fontsize=7.5, linespacing=1.45,
            bbox=dict(boxstyle="round,pad=0.40", fc="white", ec="#cccccc",
                      lw=0.7, alpha=0.90))

    ax.set_xlabel(label, fontsize=12, labelpad=4)
    ax.set_xlim(lo, hi)
    ax.grid(axis="y", color="#e0e0e0", lw=0.5, zorder=1)
    ax.set_axisbelow(True)
    ax.tick_params(axis="both", labelsize=10.5, colors="#444444")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#cccccc")
    ax.spines["bottom"].set_color("#cccccc")


def make_figure(df, cols, n_total, args, out: Path):
    panels = [p for p in PANELS if p[0] in cols]
    n = len(panels)
    fig, axs = plt.subplots(1, n, figsize=(4.2 * n, 6.4),
                            gridspec_kw={"wspace": 0.32})
    axs = np.atleast_1d(axs)
    fig.patch.set_facecolor("#ffffff")

    for ax, (col, label, colour) in zip(axs, panels):
        panel(ax, df[col].to_numpy(dtype=float), label, colour)

    axs[0].set_ylabel("Density", fontsize=13, labelpad=5)
    axs[0].legend(fontsize=10, frameon=True, loc="upper left",
                  framealpha=0.88, edgecolor="#cccccc",
                  handlelength=1.6, handletextpad=0.6)

    filt = []
    if args.min_k is not None:
        filt.append(f"$k_{{src}} \\geq$ {args.min_k}")
    if args.status:
        filt.append(f"status\u202f=\u202f{args.status}")
    sub = (f"Source: {args.table}\u2002|\u2002"
           f"{len(df):,} of {n_total:,} rows\u2002|\u2002"
           f"trimmed to {TRIM_Q[0]}\u2013{TRIM_Q[1]} pct\u2002|\u2002"
           f"shaded bands: \u00b11\u03c3, \u00b12\u03c3, \u00b13\u03c3")
    if filt:
        sub += "\u2002|\u2002" + ", ".join(filt)

    fig.suptitle("Distribution of RNA Expression $Z$-Scores",
                 fontsize=17, fontweight="bold", x=0.5, ha="center", y=1.04)
    fig.text(0.5, 0.995, sub, fontsize=9.5, color="#666666", ha="center",
             va="top", transform=fig.transFigure)
    fig.subplots_adjust(top=0.88, bottom=0.13, left=0.07, right=0.98)

    out.mkdir(parents=True, exist_ok=True)
    for ext in ("png", "pdf"):
        f = out / f"fig_zscore_distribution.{ext}"
        fig.savefig(f, dpi=200, facecolor="white", bbox_inches="tight")
        print(f"wrote {f}")
    plt.close(fig)


def main():
    args = parse_args()
    df, cols, n_total = load(args)
    print(f"loaded {len(df):,} rows x {len(cols)} z columns "
          f"(table has {n_total:,} matching rows)")
    for c in cols:
        v = df[c].to_numpy(dtype=float)
        v = v[np.isfinite(v)]
        print(f"  {c:<12} finite={v.size:>9,}  mean={v.mean():+.3f}  "
              f"sd={v.std(ddof=1):.3f}  exact0={100 * np.mean(v == 0):.2f}%")
    make_figure(df, cols, n_total, args, Path(args.out))


if __name__ == "__main__":
    main()