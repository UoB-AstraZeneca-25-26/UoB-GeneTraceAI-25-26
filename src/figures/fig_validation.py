"""
Validation and results figures (report section 4).

  12_recipe_forest      -- scoring recipes with bootstrap 95% CIs
  13_per_gene_hit20     -- per-gene hit@20, consistent denominator, by regime
  14_confidence_vs_hit  -- confidence tier against observed outcome
  15_ablation_waterfall -- what each component adds to held-out hit@20
  16_negative_results   -- four null results on one page

Figures carry data, axis labels and legends only. The interpretation lives in
CAPTIONS.md so it can be read, edited and cited as report prose.
"""
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

sys.path.insert(0, str(Path(__file__).resolve().parent))
from common import (AMBER, BLUE, GREEN, GREY, INK, LIGHT, RED, load_stats, save)

PRETTY = {
    "corr_penalized_weighted_sum": "percentile → correlation-penalised\nweighted sum   (production)",
    "residualized_noisy_or": "percentile → residualised Noisy-OR",
    "percentile_noisy_or_MATCHED": "percentile → plain Noisy-OR",
    "minmax_product_floor_MATCHED": "min–max → product with floor",
}


def fig_recipe_forest(S):
    R = S["recipes"]
    rows = sorted(R["variants"], key=lambda r: r["rho"])
    names = [PRETTY.get(r["variant"], r["variant"]) for r in rows]
    rho = np.array([r["rho"] for r in rows])
    lo = np.array([r["lo"] for r in rows])
    hi = np.array([r["hi"] for r in rows])
    prod_i = [i for i, r in enumerate(rows) if r["variant"] == R["production"]][0]

    fig, ax = plt.subplots(figsize=(10.4, 3.9))
    y = np.arange(len(rows))

    for i in range(len(rows)):
        col = BLUE if i == prod_i else GREY
        lw = 2.4 if i == prod_i else 1.7
        ax.plot([lo[i], hi[i]], [y[i], y[i]], color=col, lw=lw,
                solid_capstyle="round")
        for b in (lo[i], hi[i]):
            ax.plot([b, b], [y[i] - 0.10, y[i] + 0.10], color=col, lw=lw)
        ax.plot(rho[i], y[i], "o", ms=10 if i == prod_i else 7.5, color=col,
                zorder=4, markeredgecolor="white", markeredgewidth=1.1)
        ax.text(hi[i] + 0.0020, y[i],
                f"{rho[i]:+.4f}  [{lo[i]:+.4f}, {hi[i]:+.4f}]",
                va="center", fontsize=8.5, color=INK,
                fontweight="bold" if i == prod_i else "normal")

    ax.axvline(0, color=RED, ls="--", lw=1.3, label="ρ = 0")
    ax.set_yticks(y, names, fontsize=9)
    ax.set_xlabel(f"mean Spearman ρ against GDSC drug response, "
                  f"over {R['n_pairs']} drug–target pairs  "
                  f"(95% bootstrap CI, {R['n_boot']:,} resamples)")
    ax.set_xlim(-0.004, max(hi) + 0.028)
    ax.set_ylim(-0.6, len(rows) - 0.35)
    ax.legend(fontsize=8.5, loc="upper right")
    ax.grid(axis="y", visible=False)

    fig.tight_layout()
    return save(fig, "12_recipe_forest")


def fig_per_gene_hits(S):
    P = S["per_gene_hits"]
    df = pd.DataFrame(P["genes"]).sort_values("hr_driver")
    base = P["random_baseline"]

    fig, ax = plt.subplots(figsize=(9.6, 7.2))
    y = np.arange(len(df))
    cols = {"activation_driven": RED, "abundance_tracking": BLUE}

    for yi, (_, r) in zip(y, df.iterrows()):
        c = cols.get(r["class"], GREY)
        ax.plot([r["hr_flat"], r["hr_driver"]], [yi, yi], color="#c8cfd6",
                lw=2.0, zorder=1)
        ax.plot(r["hr_flat"], yi, "o", ms=6, color="white",
                markeredgecolor=c, markeredgewidth=1.7, zorder=3)
        ax.plot(r["hr_driver"], yi, "o", ms=8, color=c, zorder=4)

    ax.axvline(base, color=AMBER, ls="--", lw=1.5,
               label=f"random baseline = {base:.4f}")
    ax.axvline(P["summary"]["hr_driver"], color=GREEN, ls=":", lw=1.7,
               label=f"mean, driver-gated = {P['summary']['hr_driver']:.4f}")

    ax.set_yticks(y, [f"{s}   (n = {int(k)})" for s, k in
                      zip(df["hgnc_symbol"].fillna(df["ensg_id"]), df["known"])],
                  fontsize=8.5)
    ax.set_xlabel("hit@20   (denominator = GDSC-sensitive ∩ scored, shown as n)")
    ax.set_xlim(-0.002, max(df["hr_driver"].max(), df["hr_flat"].max()) * 1.10)
    ax.grid(axis="y", visible=False)

    handles = [
        Line2D([], [], marker="o", ls="", ms=8, color=RED, label="activation-driven"),
        Line2D([], [], marker="o", ls="", ms=8, color=BLUE, label="abundance-tracking"),
        Line2D([], [], marker="o", ls="", ms=6, mfc="white", mec=GREY, mew=1.7,
               label="flat score (same gene)"),
    ] + ax.get_legend_handles_labels()[0]
    ax.legend(handles=handles, fontsize=8.5, loc="lower right")

    fig.tight_layout()
    return save(fig, "13_per_gene_hit20")


def fig_confidence(S):
    C = S["confidence"]
    pair = pd.DataFrame(C["pair_level"])
    order = [t for t in ["high", "moderate", "low", "prior", "unknown"]
             if t in set(pair["confidence"])]
    pair = pair.set_index("confidence").loc[order].reset_index()
    gene = C["gene_level_hit20"]
    pts = pd.DataFrame(C["gene_level_points"])

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.8))
    cols = {"high": GREEN, "moderate": AMBER, "low": RED,
            "prior": GREY, "unknown": "#c8cfd6"}

    ax = axes[0]
    x = np.arange(len(pair))
    for i, r in pair.iterrows():
        ax.plot([x[i], x[i]], [r["lo"], r["hi"]], color=cols[r["confidence"]], lw=2.2)
        ax.plot(x[i], r["rate"], "o", ms=11, color=cols[r["confidence"]],
                markeredgecolor="white", markeredgewidth=1.3, zorder=4)
        ax.text(x[i] + 0.10, r["rate"], f"{r['rate']:.3f}", va="center",
                fontsize=9.5, fontweight="bold")
    ax.axhline(C["overall_rate"], color=INK, ls="--", lw=1.2,
               label=f"overall = {C['overall_rate']:.3f}")
    ax.set_xticks(x, [f"{t}\n{int(n):,} pairs\n{int(g)} genes" for t, n, g in
                      zip(pair["confidence"], pair["n_pairs"], pair["n_genes"])],
                  fontsize=9)
    ax.set_ylabel("P(GDSC-sensitive), pair level")
    ax.set_ylim(0, pair["hi"].max() * 1.15)
    ax.set_xlim(-0.5, len(pair) - 0.35)
    ax.legend(fontsize=8.5, loc="lower left")

    ax = axes[1]
    tiers = [t for t in order if t in gene]
    rng = np.random.default_rng(0)
    for i, t in enumerate(tiers):
        v = pts[pts["confidence"] == t]["hit_at_20"].values
        ax.plot(i + rng.normal(0, 0.055, len(v)), v, "o", ms=5, alpha=0.42,
                color=cols[t], markeredgewidth=0)
        ax.plot([i - 0.26, i + 0.26], [gene[t]["mean"]] * 2, color=cols[t], lw=3)
        ax.text(i + 0.30, gene[t]["mean"], f"{gene[t]['mean']:.4f}",
                fontsize=9, fontweight="bold", color=cols[t], va="center")
    ax.set_xticks(range(len(tiers)),
                  [f"{t}\n{int(gene[t]['count'])} genes" for t in tiers], fontsize=9)
    ax.set_xlim(-0.55, len(tiers) - 0.30)
    ax.set_ylabel("hit@20, gene level")

    fig.tight_layout()
    return save(fig, "14_confidence_vs_hit")


def fig_ablation(S):
    A = S["ablation"]
    steps = A["steps"]
    labels = [s["label"] for s in steps]
    vals = [s["value"] for s in steps]
    base = vals[0]

    fig, ax = plt.subplots(figsize=(10.6, 5.0))
    kind_col = {"base": GREY, "live": GREEN, "inert": LIGHT, "harmful": RED}

    for i, s in enumerate(steps):
        col = kind_col[s["kind"]]
        if i == 0:
            ax.bar(i, s["value"], color=col, width=0.56)
        else:
            prev = vals[i - 1]
            d = s["value"] - prev
            if abs(d) > 1e-9:
                ax.bar(i, abs(d), bottom=min(prev, s["value"]), color=col, width=0.56)
            else:
                ax.bar(i, base * 0.004, bottom=s["value"], color=col,
                       width=0.56, edgecolor="#9fb0c4")
            ax.plot([i - 0.72, i - 0.28], [prev, prev], color=GREY, lw=1.0, ls=":")

        if s["delta"] is not None and abs(s["delta"]) > 1e-9:
            lab = f"{s['delta']:+.5f}"
            if s["lo"] is not None:
                lab += f"\n[{s['lo']:+.5f}, {s['hi']:+.5f}]"
        elif s["delta"] is not None:
            lab = "no change"
        else:
            lab = ""
        ytop = max(s["value"], vals[i - 1] if i else s["value"])
        ax.text(i, ytop + base * 0.05, lab, ha="center", va="bottom", fontsize=8.5,
                color=GREEN if s["kind"] == "live" else
                      (RED if s["kind"] == "harmful" else GREY))
        ax.text(i, -base * 0.06, f"{s['value']:.5f}", ha="center", va="top",
                fontsize=9, fontweight="bold")

    ref = A["any_alteration_reference"]
    ax.axhline(ref["value"], color=AMBER, ls="--", lw=1.4,
               label=f"any-alteration arm (rejected) = {ref['value']:.5f}")
    ax.axhline(base, color=GREY, ls=":", lw=1.1)
    ax.set_xticks(range(len(steps)), labels, fontsize=8.5)
    ax.set_ylabel(f"held-out hit@20  ({A['n_genes']} genes)")
    ax.set_ylim(-base * 0.18, max(vals + [ref["value"]]) * 1.28)
    ax.grid(axis="x", visible=False)
    ax.legend(fontsize=8.5, loc="upper right")

    fig.tight_layout()
    return save(fig, "15_ablation_waterfall")


def fig_negatives(S):
    N = S["negatives"]
    fig, axes = plt.subplots(2, 2, figsize=(11.4, 7.4))

    # A -- mutation-signature discount
    ax = axes[0, 0]
    d = N["signature_discount"]
    ax.errorbar([d["point"]], [0],
                xerr=[[d["point"] - d["lo"]], [d["hi"] - d["point"]]],
                fmt="o", ms=10, color=RED, ecolor=RED, capsize=6, lw=2.0)
    ax.axvline(0, color=INK, ls="--", lw=1.4)
    ax.set_yticks([])
    ax.set_ylim(-1, 1)
    ax.set_xlabel("Δ hit@20 with signature discount  (paired bootstrap)")
    ax.set_title(f"A · signature discount:  {d['point']:+.6f}\n"
                 f"95% CI [{d['lo']:+.6f}, {d['hi']:+.6f}]", fontsize=10)
    ax.grid(axis="y", visible=False)

    # B -- metabolomics / miRNA lineage confound
    ax = axes[0, 1]
    m, r = N["metabolomics"], N["mirna"]
    r2 = [100 * m["anova_r2"], 100 * r["anova_r2"]]
    bars = ax.bar(["metabolomics", "miRNA"], r2, color=[RED, "#d98b86"], width=0.5)
    for b, v, s in zip(bars, r2, (m, r)):
        ax.text(b.get_x() + b.get_width() / 2, v + 0.2,
                f"{v:.1f}%\np = {s['anova_p']:.1e}", ha="center", fontsize=8.5)
    ax.set_ylabel("variance explained by lineage (%)")
    ax.set_ylim(0, max(r2) * 1.45)
    ax.set_title("B · candidate per-line scalars, one-way ANOVA", fontsize=10)
    ax.grid(axis="x", visible=False)

    # C -- ploidy-normalised CNA thresholds
    ax = axes[1, 0]
    p = N["ploidy"]
    x = np.arange(2)
    rho = [p["rho_before"], p["rho_after"]]
    wgd = [p["wgd_before"], p["wgd_after"]]
    ax.bar(x - 0.19, rho, width=0.34, color=BLUE, label="ρ(alteration rate, ploidy)")
    ax2 = ax.twinx()
    ax2.bar(x + 0.19, wgd, width=0.34, color=AMBER, label="WGD / non-WGD ratio")
    ax2.set_ylim(0, 2.2)
    ax2.set_ylabel("WGD rate ratio", color=AMBER)
    ax2.grid(False)
    for xi, v, pv in zip(x - 0.19, rho, [p["p_before"], p["p_after"]]):
        ax.text(xi, v + 0.003, f"{v:+.3f}\np = {pv:.3f}", ha="center", fontsize=8.5)
    for xi, v in zip(x + 0.19, wgd):
        ax2.text(xi, v + 0.04, f"{v:.2f}×", ha="center", fontsize=8.5, color=AMBER)
    ax.set_xticks(x, ["absolute\n(current)", "ploidy-relative\n(candidate)"])
    ax.set_ylabel("ρ(alteration rate, ploidy)", color=BLUE)
    ax.set_ylim(0, max(rho) * 1.7)
    ax.set_title("C · CNA thresholds, before and after ploidy normalisation",
                 fontsize=10)
    ax.grid(axis="x", visible=False)

    # D -- independent replication
    ax = axes[1, 1]
    ps = N["project_score"]
    vals = [ps["chronos_mean_rho"], ps["mean_rho"]]
    bars = ax.bar(["Chronos\n(DepMap)", "Project Score\n(Sanger + Broad)"], vals,
                  color=[GREY, GREEN], width=0.5)
    ax.axhline(0, color=INK, lw=1.1)
    for b, v in zip(bars, vals):
        ax.text(b.get_x() + b.get_width() / 2, v - 0.0003, f"{v:+.4f}",
                ha="center", va="top", fontsize=9.5, fontweight="bold")
    ax.set_ylabel("mean per-gene ρ, abundance vs dependency")
    ax.set_ylim(min(vals) * 1.5, abs(min(vals)) * 0.7)
    ax.set_title(f"D · replication over {ps['n_genes']:,} genes", fontsize=10)
    ax.grid(axis="x", visible=False)

    fig.tight_layout()
    return save(fig, "16_negative_results")


def main():
    S = load_stats()
    print("Validation figures:")
    fig_recipe_forest(S)
    fig_per_gene_hits(S)
    fig_confidence(S)
    fig_ablation(S)
    fig_negatives(S)


if __name__ == "__main__":
    main()
