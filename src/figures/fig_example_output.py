"""
17_example_query_output.png — one answerable gene beside one the system refuses.

The point of the figure is the pair of numbers in the two headers: GAPDH's top
cell line scores HIGHER than ERBB2's, and means nothing. core_score is a
within-gene percentile, so every gene is rescaled onto the same 0-1 spread
whether or not the gene discriminates. The abstention state is what separates
them, and it is measured (top_vs_median_fold), not asserted.

Everything here is read live from the production parquet at render time; no
number in this figure is transcribed from a doc.

Run:  python src/figures/fig_example_output.py
"""
import sys
from pathlib import Path

import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "pipeline"))

from common import (ROOT, PIPE, REF, OUT, plt, BLUE, RED, GREEN, GREY, INK,
                    STAGE_FILL, STAGE_EDGE)
from evidence_state import gene_verdict, load_dispersion
from explain_pair import sort_gene_rows

TOP_N = 10
GENES = ("ERBB2", "GAPDH")

BANNER_OK = "#e6f2ea"
BANNER_BAD = "#fbe9e7"
HEAD_BG = "#eef2f7"
ROW_ALT = "#f7f9fb"

TIER_COLOUR = {"high": GREEN, "moderate": "#e0a020", "low": RED, "unknown": GREY}


# ── data ──────────────────────────────────────────────────────────────────
def short_disease(s):
    """'ncit; c4017; breast ductal carcinoma || ordo; ...' -> 'breast ductal carcinoma'."""
    if not isinstance(s, str) or not s.strip():
        return "-"
    first = s.split("||")[0]
    parts = [p.strip() for p in first.split(";")]
    name = parts[-1] if parts else first
    name = name.replace(" of no special type", "")
    return name[:36]


def collect(symbol):
    gl = pd.read_parquet(REF / "gene_lookup.parquet", columns=["ensg_id", "hgnc_symbol"])
    ensg = (gl.loc[gl["hgnc_symbol"].astype("string").str.upper() == symbol, "ensg_id"]
            .astype("string").str.lower().iloc[0])

    rows = pd.read_parquet(
        PIPE / "predictions_with_confidence.parquet",
        columns=["model_id", "ensg_id", "core_score", "stratum_rank", "n_layers",
                 "class", "regime_source", "rank_basis", "has_driver_alteration",
                 "confidence"],
        filters=[("ensg_id", "==", ensg)])
    rows["model_id"] = rows["model_id"].astype("string").str.lower()

    disp = load_dispersion(str(PIPE))
    drow = disp.loc[[ensg]] if disp is not None and ensg in disp.index else None

    verdict = gene_verdict(ensg, symbol, len(rows), dispersion_row=drow,
                           gene_class=rows["class"].iloc[0] if len(rows) else None)

    cll = pd.read_parquet(REF / "cell_line_lookup.parquet",
                          columns=["model_id", "cell_line_name", "canonical_name",
                                   "diseases"])
    cll["model_id"] = cll["model_id"].astype("string").str.lower()

    top = sort_gene_rows(rows).head(TOP_N).merge(cll, on="model_id", how="left")
    top["display"] = (top["cell_line_name"].fillna(top["canonical_name"])
                      .fillna(top["model_id"]).astype(str).str.upper())
    top["tissue"] = top["diseases"].map(short_disease)

    return {"symbol": symbol, "ensg": ensg, "rows": rows, "top": top,
            "verdict": verdict,
            "fold": float(drow["top_vs_median_fold"].iloc[0]) if drow is not None else None,
            "spread": str(drow["signal_spread"].iloc[0]) if drow is not None else "n/a",
            "median_tpm": float(drow["median_log2tpm"].iloc[0]) if drow is not None else None}


# ── drawing ───────────────────────────────────────────────────────────────
def box(fig, x, y, w, h, fc, ec, lw=1.1, r=0.006, z=1):
    fig.patches.append(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
        facecolor=fc, edgecolor=ec, linewidth=lw, transform=fig.transFigure, zorder=z))


def banner_height(d):
    """Uniform across panels so the two tables start at the same y and can be
    read against each other — the whole point of the figure."""
    return 0.058 + 0.024 * len(d["verdict"].detail[:3])


def panel(fig, d, x0, w, bh):
    sym = d["symbol"]
    v = d["verdict"]
    ok = v.state == "RANKED"
    top = d["top"]
    accent = GREEN if ok else RED

    # column x-offsets within the panel, as fractions of w
    cx = {"rank": 0.030, "line": 0.075, "tissue": 0.265, "score": 0.625,
          "layers": 0.735, "tier": 0.815}

    y = 0.885

    # ── gene header ───────────────────────────────────────────────────────
    hh = 0.072
    box(fig, x0, y - hh, w, hh, STAGE_FILL, STAGE_EDGE, lw=1.4)
    fig.text(x0 + 0.012, y - 0.026, sym, fontsize=17, fontweight="bold", color=INK,
             va="center")
    fig.text(x0 + 0.012 + 0.075, y - 0.026,
             f"{d['ensg'].upper()}   ·   {len(d['rows']):,} cell lines scored",
             fontsize=8.5, color=GREY, va="center")
    fig.text(x0 + w - 0.012, y - 0.026,
             f"top core_score  {top['core_score'].iloc[0]:.4f}",
             fontsize=12, fontweight="bold", color=accent, va="center", ha="right",
             family="monospace")
    fig.text(x0 + 0.012, y - 0.055,
             f"class: {d['rows']['class'].iloc[0]}    ·    regime_source: "
             f"{d['rows']['regime_source'].iloc[0]}    ·    rank_basis: "
             f"{top['rank_basis'].iloc[0]}",
             fontsize=8, color=GREY, va="center")
    y -= hh + 0.014

    # ── verdict banner ────────────────────────────────────────────────────
    detail = v.detail[:3]
    box(fig, x0, y - bh, w, bh, BANNER_OK if ok else BANNER_BAD, accent, lw=1.4)
    fig.patches.append(Rectangle((x0, y - bh), 0.0045, bh, facecolor=accent,
                                 edgecolor="none", transform=fig.transFigure, zorder=2))
    fig.text(x0 + 0.014, y - 0.020, v.state, fontsize=10.5, fontweight="bold",
             color=accent, va="center")
    fig.text(x0 + 0.014, y - 0.044, v.headline, fontsize=9, color=INK, va="center")
    ty = y - 0.070
    for line in detail:
        fig.text(x0 + 0.016, ty, "· " + line, fontsize=7.6, color="#4a5058", va="center")
        ty -= 0.024
    y -= bh + 0.016

    # ── table header ──────────────────────────────────────────────────────
    rh = 0.0405
    box(fig, x0, y - rh, w, rh, HEAD_BG, "#c9d2dc", lw=0.9)
    for key, label, ha in (("rank", "#", "center"), ("line", "cell line", "left"),
                           ("tissue", "disease", "left"),
                           ("score", "core_score", "right"),
                           ("layers", "layers", "center"),
                           ("tier", "confidence", "left")):
        fig.text(x0 + cx[key] * w, y - rh / 2, label, fontsize=8, fontweight="bold",
                 color="#525a63", va="center", ha=ha)
    y -= rh

    # ── rows ──────────────────────────────────────────────────────────────
    for i, r in top.iterrows():
        if i % 2 == 1:
            fig.patches.append(Rectangle((x0, y - rh), w, rh, facecolor=ROW_ALT,
                                         edgecolor="none", transform=fig.transFigure,
                                         zorder=0.5))
        yc = y - rh / 2
        fig.text(x0 + cx["rank"] * w, yc, str(int(r["sorted_position"])), fontsize=8.5,
                 color=GREY, va="center", ha="center")
        fig.text(x0 + cx["line"] * w, yc, r["display"][:16], fontsize=9,
                 fontweight="bold", color=INK, va="center")
        fig.text(x0 + cx["tissue"] * w, yc, r["tissue"], fontsize=8, color="#555c64",
                 va="center")
        fig.text(x0 + cx["score"] * w, yc, f"{r['core_score']:.4f}", fontsize=9,
                 color=INK, va="center", ha="right", family="monospace")
        fig.text(x0 + cx["layers"] * w, yc, str(int(r["n_layers"])), fontsize=8.5,
                 color=BLUE if r["n_layers"] == 2 else GREY, va="center", ha="center",
                 fontweight="bold" if r["n_layers"] == 2 else "normal")
        tier = str(r["confidence"])
        tc = TIER_COLOUR.get(tier, GREY)
        fig.patches.append(FancyBboxPatch(
            (x0 + cx["tier"] * w, y - rh + 0.009), 0.052, rh - 0.018,
            boxstyle="round,pad=0,rounding_size=0.004", facecolor=tc, edgecolor="none",
            transform=fig.transFigure, zorder=2))
        fig.text(x0 + cx["tier"] * w + 0.026, yc, tier, fontsize=7.5, color="white",
                 fontweight="bold", va="center", ha="center", zorder=3)
        y -= rh

    fig.patches.append(Rectangle((x0, y), w, top.shape[0] * rh, facecolor="none",
                                 edgecolor="#c9d2dc", linewidth=0.9,
                                 transform=fig.transFigure, zorder=2))

    # ── footer: the measured fact behind the verdict ──────────────────────
    y -= 0.020
    fold = d["fold"]
    fig.text(x0 + 0.004, y - 0.016,
             f"measured spread:  top line is {fold:.1f}× a typical line   ·   "
             f"signal_spread = {d['spread']}   ·   median log2TPM {d['median_tpm']:.1f}",
             fontsize=8.2, color=accent, fontweight="bold", va="center")
    dist = d["rows"]["confidence"].value_counts().to_dict()
    dist_s = "   ".join(f"{k} {v:,}" for k, v in sorted(dist.items()))
    fig.text(x0 + 0.004, y - 0.040,
             f"all {len(d['rows']):,} rows for this gene:   {dist_s}",
             fontsize=8.2, color=GREY, va="center")
    return y - 0.055


def main():
    data = [collect(g) for g in GENES]

    fig = plt.figure(figsize=(17.5, 10.6))

    fig.text(0.03, 0.968,
             "Example output — the same query, and the answer the system declines to give",
             fontsize=19, fontweight="bold", color=INK, va="center")
    fig.text(0.03, 0.936,
             "Left: a gene whose ranking is worth acting on. Right: a gene that is measured, ranks perfectly well arithmetically, and means nothing. "
             "Note the two top scores — the meaningless one is the higher.",
             fontsize=10, color=GREY, va="center")

    bh = max(banner_height(d) for d in data)
    ends = [panel(fig, data[0], 0.030, 0.452, bh),
            panel(fig, data[1], 0.518, 0.452, bh)]

    # ── the punchline strip ───────────────────────────────────────────────
    a, b = data[0], data[1]
    ytop = min(ends) - 0.012
    h = 0.088
    box(fig, 0.030, ytop - h, 0.940, h, "#fdf6e3", "#e0a020", lw=1.4)
    fig.text(0.044, ytop - 0.026,
             "Why a percentile alone cannot be the answer",
             fontsize=11, fontweight="bold", color="#8a6410", va="center")
    fig.text(0.044, ytop - 0.052,
             f"core_score is a WITHIN-GENE percentile across cell lines, so every gene is rescaled onto the same 0–1 range. "
             f"{b['symbol']}'s best line scores {b['top']['core_score'].iloc[0]:.4f} at {b['fold']:.1f}× a typical line; "
             f"{a['symbol']}'s best scores {a['top']['core_score'].iloc[0]:.4f} at {a['fold']:.1f}×.",
             fontsize=9.2, color="#5a4a20", va="center")
    fig.text(0.044, ytop - 0.073,
             "The scores are indistinguishable; the biology is not. The pipeline therefore measures the spread separately (top_vs_median_fold) and abstains — "
             "rather than printing a confident-looking table with a warning underneath it.",
             fontsize=9.2, color="#5a4a20", va="center")

    fig.text(0.030, 0.018,
             "Rendered live from predictions_with_confidence.parquet · gene_dispersion.parquet · cell_line_lookup.parquet. "
             "Ordering re-derived at query time by explain_pair.sort_gene_rows; verdict by evidence_state.gene_verdict.",
             fontsize=7.5, color=GREY, va="center")

    path = OUT / "17_example_query_output.png"
    fig.savefig(path, dpi=200, bbox_inches=None)
    plt.close(fig)
    print(f"wrote {path}  ({path.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
