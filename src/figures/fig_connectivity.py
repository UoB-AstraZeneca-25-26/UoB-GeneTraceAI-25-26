"""
18_connectivity_raw.png — the raw gene <-> cell line bipartite graph.

RAW, NOT PIPELINE. The tallies below are the pre-pipeline connectivity count
over the warehouse (20,163 genes x 2,052 cell lines), NOT the scored universe
(19,176 x 1,485). Nothing here has been through harmonisation, scoring or
routing -- it is the shape of the data as delivered.

The counts are transcribed from the raw connectivity tally rather than
recomputed, because the full pair-level scan is expensive and the question
being answered ("which cell line connects to more genes, or vice versa") is
fully determined by the summary. Everything DERIVED from them is computed
here and reconciles:

    mean(genes/line)  x 2,052  = 31,670,568  ~ 31,670,726   (rounding)
    mean(lines/gene)  x 20,163 = 31,676,073  ~ 31,670,726   (rounding)
    31,670,726 / (20,163 x 2,052)           = 76.55%        (stated 76.5%)

Run:  python src/figures/fig_connectivity.py
"""
import sys
from pathlib import Path

from matplotlib.patches import FancyBboxPatch

sys.path.insert(0, str(Path(__file__).resolve().parent))

from common import OUT, plt, BLUE, RED, GREEN, AMBER, GREY, INK, LIGHT

# ── the raw tally, exactly as counted ─────────────────────────────────────
RAW = {
    "pairs": 31_670_726,
    "genes": 20_163,
    "lines": 2_052,
    "density_pct": 76.5,
    "lines_per_gene": {"median": 1_582, "mean": 1_571, "max": 1_735,
                       "bins": [("0–10", 0), ("10–100", 0),
                                ("100–1,000", 0), ("1,000+", 20_163)]},
    "genes_per_line": {"median": 20_163, "mean": 15_434, "max": 20_163,
                       "bins": [("0–100", 242), ("100–1,000", 223),
                                ("1,000–15,000", 7), ("15,000+", 1_580)]},
}

# upper bound on what the sparse lines can contribute: every line in a bin
# credited with the TOP of its bin. Deliberately generous to the sparse side.
BIN_TOP = {"0–100": 100, "100–1,000": 1_000, "1,000–15,000": 15_000}


def derived(r):
    g, l, p = r["genes"], r["lines"], r["pairs"]
    sparse_lines = sum(n for lbl, n in r["genes_per_line"]["bins"] if lbl in BIN_TOP)
    sparse_max = sum(n * BIN_TOP[lbl] for lbl, n in r["genes_per_line"]["bins"]
                     if lbl in BIN_TOP)
    dense_lines = l - sparse_lines
    return {
        "possible": g * l,
        "density": 100 * p / (g * l),
        "dense_lines": dense_lines,
        "sparse_lines": sparse_lines,
        "dense_line_pct": 100 * dense_lines / l,
        "sparse_pair_pct_max": 100 * sparse_max / p,
        "dense_pair_pct_min": 100 - 100 * sparse_max / p,
        "recon_line": r["genes_per_line"]["mean"] * l,
        "recon_gene": r["lines_per_gene"]["mean"] * g,
    }


def box(fig, x, y, w, h, fc, ec, lw=1.2, r=0.006, z=1):
    fig.patches.append(FancyBboxPatch(
        (x, y), w, h, boxstyle=f"round,pad=0,rounding_size={r}",
        facecolor=fc, edgecolor=ec, linewidth=lw, transform=fig.transFigure, zorder=z))


def kpi(fig, x, y, w, value, label, colour):
    box(fig, x, y, w, 0.088, "#f4f7fa", "#c9d2dc", lw=1.0)
    fig.text(x + w / 2, y + 0.058, value, fontsize=20, fontweight="bold",
             color=colour, ha="center", va="center")
    fig.text(x + w / 2, y + 0.024, label, fontsize=9, color=GREY,
             ha="center", va="center")


def bar_panel(ax, bins, title, subtitle, colour, note):
    labels = [b[0] for b in bins][::-1]
    vals = [b[1] for b in bins][::-1]
    ypos = range(len(vals))
    top = max(vals) or 1

    cols = [colour if v else "#ffffff" for v in vals]
    ax.barh(list(ypos), vals, color=cols, edgecolor=colour, linewidth=1.3, height=0.62)

    for i, v in enumerate(vals):
        pct = 100 * v / sum(vals) if sum(vals) else 0
        if v == 0:
            ax.text(top * 0.012, i, "0", fontsize=10, color=GREY, va="center")
        else:
            inside = v > top * 0.45
            ax.text(v - top * 0.015 if inside else v + top * 0.015, i,
                    f"{v:,}   ({pct:.1f}%)", fontsize=10.5, fontweight="bold",
                    color="white" if inside else INK,
                    va="center", ha="right" if inside else "left")

    ax.set_yticks(list(ypos))
    ax.set_yticklabels(labels, fontsize=10)
    ax.set_xlim(0, top * 1.32)
    ax.set_xticks([])
    ax.grid(False)
    for s in ("top", "right", "bottom"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color("#c9d2dc")
    ax.set_title(title, fontsize=13, fontweight="bold", color=INK, loc="left", pad=16)
    ax.text(0, 1.045, subtitle, transform=ax.transAxes, fontsize=9.2, color=GREY)
    ax.text(0, -0.13, note, transform=ax.transAxes, fontsize=9, color=colour,
            fontweight="bold")


def stacked(ax, rows):
    """rows: [(title, [(label, frac, colour)])] — two comparable 100% bars."""
    for i, (title, segs) in enumerate(rows):
        left = 0.0
        y = len(rows) - 1 - i
        for label, frac, col in segs:
            ax.barh(y, frac, left=left, height=0.5, color=col, edgecolor="white",
                    linewidth=1.6)
            if frac > 6:
                ax.text(left + frac / 2, y, f"{label}\n{frac:.1f}%", fontsize=9.5,
                        color="white", fontweight="bold", ha="center", va="center")
            else:
                ax.text(left + frac / 2, y - 0.42, f"{label}  {frac:.1f}%",
                        fontsize=9, color=RED, fontweight="bold", ha="center",
                        va="center")
            left += frac
        ax.text(-1.5, y, title, fontsize=10.5, fontweight="bold", color=INK,
                ha="right", va="center")

    ax.set_xlim(0, 100)
    ax.set_ylim(-0.75, len(rows) - 0.4)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)


def main():
    r, d = RAW, derived(RAW)

    fig = plt.figure(figsize=(16.8, 9.8))

    fig.text(0.035, 0.963, "Gene ↔ cell line connectivity — raw counts, before the pipeline",
             fontsize=20, fontweight="bold", color=INK, va="center")
    fig.text(0.035, 0.930,
             "The bipartite graph as the sources deliver it: 20,163 genes × 2,052 cell lines. "
             "Not the scored universe (19,176 × 1,485) — nothing here has been harmonised, scored or routed.",
             fontsize=10, color=GREY, va="center")

    # ── answer up front ───────────────────────────────────────────────────
    box(fig, 0.035, 0.822, 0.930, 0.080, "#eef4fb", BLUE, lw=1.5)
    fig.text(0.048, 0.882, "Which side connects to more?", fontsize=10.5,
             fontweight="bold", color=BLUE, va="center")
    fig.text(0.048, 0.858,
             "Genes are all alike — every one of the 20,163 touches over 1,000 cell lines, and the spread is tiny "
             "(median 1,582, max 1,735).",
             fontsize=9.8, color=INK, va="center")
    fig.text(0.048, 0.836,
             "Cell lines are not: they split into a fully-profiled majority and a barely-measured tail, with almost "
             "nothing in between. So the question only has an answer on one side.",
             fontsize=9.8, color=INK, va="center")

    # ── KPI strip ─────────────────────────────────────────────────────────
    kpi(fig, 0.035, 0.716, 0.222, f"{r['pairs']:,}", "distinct gene–line pairs", BLUE)
    kpi(fig, 0.271, 0.716, 0.222, f"{r['genes']:,}", "genes with ≥1 cell line", GREEN)
    kpi(fig, 0.507, 0.716, 0.222, f"{r['lines']:,}", "cell lines with ≥1 gene", AMBER)
    kpi(fig, 0.743, 0.716, 0.222, f"{d['density']:.1f}%",
        f"of the {d['possible']:,} possible pairs", GREY)

    # ── distributions ─────────────────────────────────────────────────────
    axA = fig.add_axes([0.075, 0.400, 0.375, 0.220])
    bar_panel(
        axA, r["genes_per_line"]["bins"],
        "Genes per cell line  —  bimodal",
        f"median {r['genes_per_line']['median']:,}   ·   mean {r['genes_per_line']['mean']:,}"
        f"   ·   max {r['genes_per_line']['max']:,}",
        AMBER,
        "mean sits far below the median → a long left tail of barely-measured lines")

    axB = fig.add_axes([0.590, 0.400, 0.375, 0.220])
    bar_panel(
        axB, r["lines_per_gene"]["bins"],
        "Cell lines per gene  —  no structure",
        f"median {r['lines_per_gene']['median']:,}   ·   mean {r['lines_per_gene']['mean']:,}"
        f"   ·   max {r['lines_per_gene']['max']:,}",
        GREEN,
        "every gene lands in one bin: there is no interesting variation on this axis")

    # ── the asymmetry ─────────────────────────────────────────────────────
    axC = fig.add_axes([0.215, 0.155, 0.750, 0.150])
    sparse_pct = 100 - d["dense_line_pct"]
    stacked(axC, [
        ("share of cell lines",
         [(f"{d['dense_lines']:,} fully profiled", d["dense_line_pct"], BLUE),
          (f"{d['sparse_lines']:,} sparse", sparse_pct, RED)]),
        ("share of all gene–line pairs",
         [("carried by the fully-profiled lines", d["dense_pair_pct_min"], BLUE),
          ("sparse", d["sparse_pair_pct_max"], RED)]),
    ])
    axC.set_title("The asymmetry, in one line", fontsize=13, fontweight="bold",
                  color=INK, loc="left", pad=14)

    box(fig, 0.035, 0.048, 0.930, 0.078, "#fdf6e3", AMBER, lw=1.4)
    fig.text(0.048, 0.106, "So: the answer is on the cell-line side",
             fontsize=10.5, fontweight="bold", color="#8a6410", va="center")
    fig.text(0.048, 0.079,
             f"{d['dense_lines']:,} of {r['lines']:,} cell lines ({d['dense_line_pct']:.0f}%) carry at least "
             f"{d['dense_pair_pct_min']:.1f}% of every gene–line connection in the dataset. The remaining "
             f"{d['sparse_lines']:,} lines contribute at most {d['sparse_pair_pct_max']:.1f}%.",
             fontsize=9.6, color="#5a4a20", va="center")
    fig.text(0.048, 0.059,
             "That bound is generous to the sparse side — every sparse line is credited with the TOP of its bin. The real share is lower. "
             "It is why coverage, not gene choice, is what limits the pipeline.",
             fontsize=9.0, color="#7a6430", va="center")

    fig.text(0.035, 0.018,
             f"Raw connectivity tally over the warehouse; summary counts transcribed, all percentages and bounds derived here. "
             f"Reconciles: mean(genes/line)×{r['lines']:,} = {d['recon_line']:,} and mean(lines/gene)×{r['genes']:,} = "
             f"{d['recon_gene']:,}, both ≈ {r['pairs']:,} to rounding.",
             fontsize=7.6, color=GREY, va="center")

    path = OUT / "18_connectivity_raw.png"
    fig.savefig(path, dpi=200)
    plt.close(fig)
    print(f"wrote {path}  ({path.stat().st_size // 1024} KB)")


if __name__ == "__main__":
    main()
