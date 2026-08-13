"""
Method-justification figures (report section 3).

Each figure exists because it killed or shaped a design decision:

  06_mrna_protein        -> rejection of Noisy-OR (independence violated)
  07_normalisation       -> choice of percentile normalisation
  08_interlayer_corr     -> the correlation penalty in the weighted sum
  09_score_vs_layers     -> core_score is a within-gene ranking, not a probability
  10_tsg_deletion        -> the TSG score inversion is mechanistically real
  11_routing_statistic   -> gene-class routing, and the size of the unknown regime

Figures carry data, axis labels and legends only. The interpretation lives in
CAPTIONS.md so it can be read, edited and cited as report prose.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (AMBER, BLUE, GREEN, GREY, INK, LIGHT, REF, RED, VALID,
                    load_expr_wide, load_prot_wide, load_stats, save)


def _pct(a):
    return pd.Series(a).rank(pct=True).values


# ─────────────────────────────────────────────────────────────────────────
def fig_mrna_protein(S):
    epc = pd.read_parquet(VALID / "expr_protein_correlation_per_gene.parquet")
    # a pre-v2 validation artefact, still keyed on UPPERCASE ENSG. The gene id
    # picked out of it is used to index the expression/protein matrices, which
    # are lowercase, so normalise it here rather than at the point of use.
    epc.index = epc.index.astype(str).str.split(".").str[0].str.lower()
    rho_all = epc["rho"].dropna()
    med = float(rho_all.median())

    gene = (rho_all - med).abs().idxmin()
    gl = pd.read_parquet(REF / "gene_lookup.parquet")[["ensg_id", "hgnc_symbol"]]
    gl["ensg_id"] = gl["ensg_id"].astype(str).str.split(".").str[0].str.lower()
    sym = gl.loc[gl["ensg_id"] == gene, "hgnc_symbol"]
    sym = sym.iloc[0] if len(sym) else gene

    E = load_expr_wide(genes=[gene])
    P = load_prot_wide(genes=[gene])
    shared = E.index.intersection(P.index)
    e, p = E.loc[shared, gene], P.loc[shared, gene]
    ok = e.notna() & p.notna()
    ex, px = _pct(e[ok].values), _pct(p[ok].values)
    rho_g = spearmanr(ex, px).statistic

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8))

    ax = axes[0]
    ax.scatter(ex, px, s=20, alpha=0.40, color=BLUE, linewidths=0)
    ax.plot([0, 1], [0, 1], color=GREY, ls="--", lw=1.1, label="y = x")
    ax.set_xlabel("mRNA expression percentile")
    ax.set_ylabel("Protein abundance percentile")
    ax.set_xlim(-0.02, 1.02)
    ax.set_ylim(-0.02, 1.02)
    ax.set_title(f"{sym} —  ρ = {rho_g:.2f},  n = {int(ok.sum())}", fontsize=10.5)
    ax.legend(fontsize=8.5, loc="upper left")

    ax = axes[1]
    ax.hist(rho_all, bins=50, color=BLUE, alpha=0.85, edgecolor="none")
    ax.axvline(med, color=RED, lw=1.8, label=f"median = {med:.3f}")
    ax.axvline(0, color=GREY, lw=1.0, ls=":")
    ax.set_xlabel("Spearman ρ between mRNA and protein, per gene")
    ax.set_ylabel("genes")
    ax.set_xlim(-0.6, 1.02)
    ax.set_title(f"All genes with both layers  (n = {len(rho_all):,})", fontsize=10.5)
    ax.legend(fontsize=8.5, loc="upper left")

    fig.tight_layout()
    return save(fig, "06_mrna_protein")


# ─────────────────────────────────────────────────────────────────────────
def fig_normalisation(S):
    N = S["normalisation"]
    fig, axes = plt.subplots(2, 3, figsize=(12.0, 6.2))

    for row, key in enumerate(("expr", "prot")):
        d = N[key]
        colour = BLUE if key == "expr" else RED
        edges = np.array(d["raw_edges"])
        centres = (edges[:-1] + edges[1:]) / 2

        ax = axes[row, 0]
        ax.bar(centres, d["raw_hist"], width=np.diff(edges), color=colour,
               alpha=0.85, edgecolor="none")
        ax.set_ylabel(f"{d['label']}\n\nmeasurements", fontsize=9)
        ax.set_xlabel("native units")
        if row == 0:
            ax.set_title("raw", fontsize=10.5)

        ax = axes[row, 1]
        ax.bar(np.linspace(0.008, 0.992, 60), d["mm_hist"], width=1 / 60,
               color=colour, alpha=0.85, edgecolor="none")
        ax.set_xlabel(f"min–max value          IQR = {d['mm_iqr']:.3f}")
        ax.set_xlim(0, 1)
        if row == 0:
            ax.set_title("min–max", fontsize=10.5)

        ax = axes[row, 2]
        ax.bar(np.linspace(0.008, 0.992, 60), d["pct_hist"], width=1 / 60,
               color=colour, alpha=0.85, edgecolor="none")
        ax.set_xlabel(f"percentile within gene          IQR = {d['pct_iqr']:.3f}")
        ax.set_xlim(0, 1)
        if row == 0:
            ax.set_title("percentile rank  (production)", fontsize=10.5)

    fig.tight_layout()
    return save(fig, "07_normalisation")


# ─────────────────────────────────────────────────────────────────────────
def fig_interlayer(S):
    I = S["interlayer"]
    layers = ["Expression", "Protein", "Chronos"]
    rho = np.array([
        [1.0, I["rho_EP_median"], I["rho_EC_median"]],
        [I["rho_EP_median"], 1.0, I["rho_PC_median"]],
        [I["rho_EC_median"], I["rho_PC_median"], 1.0]])

    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.6))

    ax = axes[0]
    im = ax.imshow(rho, cmap="RdBu_r", vmin=-1, vmax=1)
    for i in range(3):
        for j in range(3):
            ax.text(j, i, f"{rho[i, j]:.3f}", ha="center", va="center",
                    fontsize=11, fontweight="bold",
                    color="white" if abs(rho[i, j]) > 0.55 else INK)
    ax.set_xticks(range(3), layers)
    ax.set_yticks(range(3), layers)
    ax.grid(False)
    ax.set_title("median per-gene Spearman ρ", fontsize=10.5)
    fig.colorbar(im, ax=ax, shrink=0.76, label="ρ")

    ax = axes[1]
    w = I["weights_three_layer"]
    names, vals = list(w), [w[n] for n in w]
    bars = ax.bar(names, vals, color=[BLUE, RED, GREEN], width=0.55, alpha=0.9)
    ax.axhline(1 / 3, color=GREY, ls="--", lw=1.3, label="equal weighting")
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.005, f"{v:.3f}",
                ha="center", fontsize=9.5, fontweight="bold")
    ax.set_ylabel("wᵢ = 1/(1+|ρ̄ᵢ|), normalised")
    ax.set_ylim(0, max(vals) * 1.20)
    ax.set_title("resulting layer weights", fontsize=10.5)
    ax.legend(fontsize=8.5)

    fig.tight_layout()
    return save(fig, "08_interlayer_correlation")


# ─────────────────────────────────────────────────────────────────────────
def fig_score_vs_layers(S):
    D = S["score_vs_layers"]
    fig, axes = plt.subplots(1, 2, figsize=(10.4, 4.8))

    for ax, col in zip(axes, ("core_score", "stratum_rank")):
        keys = sorted(D[col])
        data = [np.array(D[col][k]["sample"]) for k in keys]
        bp = ax.boxplot(data, patch_artist=True, widths=0.5,
                        medianprops=dict(color="black", lw=2),
                        flierprops=dict(marker=".", ms=2, alpha=0.22,
                                        markerfacecolor=GREY,
                                        markeredgecolor="none"))
        for b, c in zip(bp["boxes"], [LIGHT, "#9fc2e8"]):
            b.set_facecolor(c)
        ax.set_xticks(range(1, len(keys) + 1),
                      [f"n_layers = {k}\nn = {D[col][k]['n']:,}" for k in keys])
        ax.set_ylabel(col)
        ax.set_ylim(-0.02, 1.02)
        for i, k in enumerate(keys, start=1):
            ax.text(i + 0.30, D[col][k]["mean"], f"mean {D[col][k]['mean']:.3f}",
                    fontsize=8.5, va="center", color=INK)

    fig.tight_layout()
    return save(fig, "09_score_vs_layers")


# ─────────────────────────────────────────────────────────────────────────
def fig_tsg_deletion(S):
    T = S["tsg_deletion"]
    genes = [g for g, d in T.items() if d.get("deleted_n", 0) >= 5]

    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.9))

    ax = axes[0]
    d = T["PTEN"]
    dele = 1 - np.array(d["deleted_sample"])
    reta = 1 - np.array(d["retained_sample"])
    bins = np.linspace(0, 1, 34)
    ax.hist(reta, bins=bins, density=True, color=LIGHT, edgecolor="#9fb0c4",
            label=f"copy retained / no call  (n = {d['retained_n']:,})")
    ax.hist(dele, bins=bins, density=True, color=RED, alpha=0.75,
            edgecolor="none", label=f"homozygous deletion  (n = {d['deleted_n']})")
    ax.axvline(dele.mean(), color=RED, lw=2)
    ax.axvline(reta.mean(), color="#6b7683", lw=2)
    ax.set_xlabel("served score for PTEN  (1 − core_score)")
    ax.set_ylabel("density")
    ax.set_title(f"PTEN:  {reta.mean():.3f} → {dele.mean():.3f},  "
                 f"p = {d['mannwhitney_p']:.1e}", fontsize=10.5)
    ax.legend(fontsize=8.5, loc="upper left")

    ax = axes[1]
    y = np.arange(len(genes))[::-1]
    dm = [1 - T[g]["deleted_mean"] for g in genes]
    rm = [1 - T[g]["retained_mean"] for g in genes]
    ax.hlines(y, rm, dm, color="#c8cfd6", lw=2.2, zorder=1)
    ax.scatter(rm, y, s=60, color="#6b7683", zorder=3, label="copy retained / no call")
    ax.scatter(dm, y, s=60, color=RED, zorder=3, label="homozygous deletion")
    ax.set_yticks(y, [f"{g}  (n = {T[g]['deleted_n']})" for g in genes])
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("mean served score  (1 − core_score)")
    ax.legend(fontsize=8.5, loc="lower left")
    ax.grid(axis="y", visible=False)

    fig.tight_layout()
    return save(fig, "10_tsg_deletion")


# ─────────────────────────────────────────────────────────────────────────
def fig_routing(S):
    R = S["routing"]
    cur = pd.DataFrame(R["curated"])
    thr = R["threshold"]

    fig, axes = plt.subplots(1, 2, figsize=(11.8, 4.8),
                             gridspec_kw={"width_ratios": [1.45, 1.0]})

    ax = axes[0]
    ab = cur[cur["class"] == "abundance_tracking"]["stat"]
    act = cur[cur["class"] == "activation_driven"]["stat"]
    bins = np.linspace(0, cur["stat"].max() * 1.05, 34)
    ax.hist(act, bins=bins, color=RED, alpha=0.70,
            label=f"activation-driven  (n = {len(act)})")
    ax.hist(ab, bins=bins, color=BLUE, alpha=0.70,
            label=f"abundance-tracking  (n = {len(ab)})")
    ax.axvline(thr, color=AMBER, lw=2, ls="--",
               label=f"routing threshold = {thr:.4f}")
    ax.set_xlabel("flat-score hit@20 on the training genes")
    ax.set_ylabel("genes")
    ax.set_title(f"routing statistic, {len(cur)} curated genes", fontsize=10.5)
    ax.legend(fontsize=8.5)

    ax = axes[1]
    rows = R["class_rows_full"]
    order = [o for o in ["unknown", "abundance_tracking", "loss_of_function",
                         "activation_driven"] if o in rows]
    vals = [rows[o] for o in order]
    total = sum(vals)
    cols = {"unknown": "#c8cfd6", "abundance_tracking": BLUE,
            "activation_driven": RED, "loss_of_function": AMBER}
    y = np.arange(len(order))[::-1]
    ax.barh(y, vals, color=[cols[o] for o in order], height=0.6)
    ax.set_yticks(y, order)
    ax.set_xscale("log")
    ax.set_xlabel("rows in predictions_with_confidence  (log scale)")
    for yi, v in zip(y, vals):
        ax.text(v * 1.12, yi, f"{v:,}  ({100 * v / total:.2f}%)", va="center",
                fontsize=8.5)
    ax.set_xlim(1e5, total * 6)
    ax.set_title("how all 28.4M rows are routed", fontsize=10.5)
    ax.grid(axis="y", visible=False)

    fig.tight_layout()
    return save(fig, "11_routing_statistic")


def main():
    S = load_stats()
    print("Method-justification figures:")
    fig_mrna_protein(S)
    fig_normalisation(S)
    fig_interlayer(S)
    fig_score_vs_layers(S)
    fig_tsg_deletion(S)
    fig_routing(S)


if __name__ == "__main__":
    main()
