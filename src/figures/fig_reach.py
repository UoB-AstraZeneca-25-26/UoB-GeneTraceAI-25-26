"""
Data-reach diagnostics (report section 2).

  04_gene_axis_reach     -- rows a per-gene query returns per source,
                            native ENSG axis vs after the ID bridges
  05_layer_availability  -- UpSet-style: which modality combinations the
                            cell-line panel actually has
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import BLUE, GREEN, GREY, INK, LIGHT, load_stats, save

# Which sources carry a usable ENSG column as delivered, and which only
# become reachable once the bridge tables exist. Taken from the native_key
# recorded alongside each source in _stats.json.
NATIVE = {"Expression (DepMap)", "Mutations", "Fusions", "CNA flags",
          "Chronos essentiality"}
BRIDGED = {"Proteomics (CCLE/Gygi)": "UniProt → ENSG",
           "GDSC drug targets": "symbol → ENSG"}


def fig_gene_axis_reach(S):
    src = S["gene_reach"]["sources"]
    universe = S["gene_reach"]["universe"]

    # order: bridged first (the story), then native, then unreachable
    def key(r):
        if r["source"] in BRIDGED:
            return (0, -r["genes_in_universe"])
        if r["source"] in NATIVE:
            return (1, -r["genes_in_universe"])
        return (2, -r["genes_in_universe"])

    src = sorted(src, key=key)
    labels = [r["source"] for r in src]
    after = np.array([r["genes_in_universe"] for r in src], dtype=float)
    before = np.array([r["genes_in_universe"] if r["source"] in NATIVE else 0
                       for r in src], dtype=float)

    fig, ax = plt.subplots(figsize=(12.4, 7.4))
    y = np.arange(len(labels))[::-1]
    h = 0.36

    ax.barh(y + h / 2, before, height=h, color=LIGHT, edgecolor="#9fb0c4",
            label="reachable on a native ensg_id column")
    ax.barh(y - h / 2, after, height=h, color=BLUE, edgecolor="none",
            label="reachable after the ID bridges (production)")

    for yi, b, a, r in zip(y, before, after, src):
        if a > 0:
            ax.text(a + universe * 0.008, yi - h / 2, f"{int(a):,}",
                    va="center", ha="left", fontsize=8, color=INK)
        if a == 0:
            ax.text(universe * 0.004, yi, "  0", va="center", ha="left",
                    fontsize=8, color=GREY)

    ax.axvline(universe, color=GREY, ls=":", lw=1.2)
    ax.text(universe, len(labels) - 0.2, f"  gene universe\n  {universe:,}",
            fontsize=8, color=GREY, va="top")

    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=9.5)
    ax.set_xlabel("Genes in the protein-coding universe the source can answer for")
    ax.set_xlim(0, universe * 1.16)
    ax.set_title("Genes reachable by a per-gene query, by source",
                 fontsize=11.5, loc="left")
    ax.legend(loc="lower right", fontsize=9)
    ax.grid(axis="y", visible=False)

    fig.tight_layout()
    return save(fig, "04_gene_axis_reach")


def fig_layer_availability(S):
    combos = S["modality"]["combinations"]
    sizes = S["modality"]["set_sizes"]
    order = [k for k, _ in sorted(sizes.items(), key=lambda kv: -kv[1])]
    total = S["modality"]["universe"]

    top = combos[:14]
    n_bars = len(top)

    fig = plt.figure(figsize=(13.6, 9.0))
    gs = fig.add_gridspec(2, 2, height_ratios=[2.05, 1.55],
                          width_ratios=[1.0, 3.35], hspace=0.06, wspace=0.03)
    ax_bar = fig.add_subplot(gs[0, 1])
    ax_dot = fig.add_subplot(gs[1, 1], sharex=ax_bar)
    ax_set = fig.add_subplot(gs[1, 0], sharey=ax_dot)
    fig.add_subplot(gs[0, 0]).axis("off")

    x = np.arange(n_bars)
    vals = [c["n"] for c in top]
    bars = ax_bar.bar(x, vals, color=BLUE, width=0.66)
    for b, v in zip(bars, vals):
        ax_bar.text(b.get_x() + b.get_width() / 2, v + total * 0.006,
                    f"{v}", ha="center", va="bottom", fontsize=8)
    # highlight the fully-covered combination
    for i, c in enumerate(top):
        if len(c["layers"]) == len(order):
            bars[i].set_color(GREEN)
    ax_bar.set_ylabel("cell lines with exactly\nthis combination")
    ax_bar.set_ylim(0, max(vals) * 1.16)
    ax_bar.tick_params(labelbottom=False)
    ax_bar.grid(axis="x", visible=False)
    ax_bar.set_title(
        f"Modality combinations across the {total:,} cell lines "
        f"({n_bars} most common shown)", fontsize=11.5, loc="left")

    yy = np.arange(len(order))[::-1]
    for i, c in enumerate(top):
        present = set(c["layers"])
        for j, lay in enumerate(order):
            yv = yy[j]
            on = lay in present
            ax_dot.plot(i, yv, "o", ms=8.5,
                        color=INK if on else "#dfe3e8", zorder=3)
        idx = [yy[j] for j, lay in enumerate(order) if lay in present]
        if len(idx) > 1:
            ax_dot.plot([i, i], [min(idx), max(idx)], color=INK, lw=1.6, zorder=2)

    ax_dot.set_yticks(yy)
    ax_dot.set_yticklabels([])
    ax_dot.set_xticks(x)
    ax_dot.set_xticklabels([])
    ax_dot.set_ylim(-0.7, len(order) - 0.3)
    ax_dot.grid(False)
    for s in ax_dot.spines.values():
        s.set_visible(False)
    ax_dot.tick_params(length=0)
    ax_dot.set_xlabel("modality combination (ordered by number of cell lines)")

    ax_set.barh(yy, [sizes[l] for l in order], height=0.5, color="#9fb0c4")
    # Labels live in reserved space to the right of the bars (x axis is
    # inverted, so negative x is further right).
    for j, lay in enumerate(order):
        ax_set.text(-total * 0.03, yy[j], f"{lay}  ({sizes[lay]:,})",
                    ha="left", va="center", fontsize=9)
    ax_set.set_xlim(total * 1.02, -total * 0.62)
    ax_set.set_xticks([0, 500, 1000, 1500])
    ax_set.set_yticks([])
    ax_set.set_xlabel("lines with the\nmodality at all")
    ax_set.grid(axis="y", visible=False)
    for s in ("top", "right", "left"):
        ax_set.spines[s].set_visible(False)


    return save(fig, "05_layer_availability")


def main():
    S = load_stats()
    print("Data-reach figures:")
    fig_gene_axis_reach(S)
    fig_layer_availability(S)


if __name__ == "__main__":
    main()
