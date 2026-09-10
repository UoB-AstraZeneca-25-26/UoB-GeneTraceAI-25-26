"""
Figure 1 - transcriptomic execution pipeline (Chapter 3).

Begins at the HARMONISED CORPUS produced by the Data Harmonisation
section of the same chapter, not at raw data.
No data dependency: this documents control flow implemented in
transcriptomics.py / 05_transcriptomics_stat_layer.py. Labels are the
actual function, column and status names in that code.
"""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

CH2, SRC, PROC, DEC, OUT = "#e4e4ef", "#dbe8f5", "#eeeeee", "#fdf0d5", "#d8ecd8"
EDGE = "#4a4a4a"

fig, ax = plt.subplots(figsize=(8.6, 11.6))
ax.set_xlim(-1.2, 12.4); ax.set_ylim(0, 27.4); ax.axis("off")

def box(x, y, w, h, title, sub, fc, ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.10",
                                fc=fc, ec=EDGE, lw=0.9, linestyle=ls))
    ax.text(x + w/2, y + h*(0.62 if sub else 0.5), title, ha="center",
            va="center", fontsize=8.6, fontweight="bold")
    if sub:
        ax.text(x + w/2, y + h*0.24, sub, ha="center", va="center",
                fontsize=6.9, style="italic", color="#333333")

def arrow(x1, y1, x2, y2, ls="-"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle="-|>",
                                 mutation_scale=11, lw=0.9, color=EDGE,
                                 linestyle=ls, shrinkA=0, shrinkB=0))

# --- section boundary ------------------------------------------------
box(1.60, 25.75, 6.80, 1.30,
    "Harmonised corpus",
    "2,127 cell lines x 20,163 genes | keyed on model id + ENSG", CH2)
ax.plot([-1.0, 11.2], [25.35, 25.35], ls=(0, (4, 3)), lw=0.8, color="#8a5a00")
ax.text(11.1, 25.5, "Data Harmonisation", fontsize=6.4, ha="right",
        va="bottom", color="#8a5a00", style="italic")
ax.text(11.1, 25.05, "Transcriptomic Evidence", fontsize=6.4, ha="right",
        va="top", color="#8a5a00", style="italic")

box(1.60, 23.55, 6.80, 1.25, "Transcriptomic records retrieved",
    "hpa_rna | depmap_expr + depmap_profiles | geo_expr + geo_info", PROC)
arrow(5.00, 25.75, 5.00, 24.80)

# --- sources -----------------------------------------------------------
box(0.15, 21.45, 2.90, 1.30, "HPA", "TPM | 1 method", SRC)
box(3.55, 21.45, 2.90, 1.30, "DepMap", "log2(TPM+1) | 1 method", SRC)
box(6.95, 21.45, 2.90, 1.30, "GEO", "arrays + seq | N methods", SRC)
for x in (1.60, 5.00, 8.40):
    arrow(5.00, 23.55, x, 22.75)

# --- method resolution -------------------------------------------------
box(6.95, 19.55, 2.90, 1.35, "Method resolution",
    "SOFT data_processing\nplatform x proc_family", PROC)
arrow(8.40, 21.45, 8.40, 20.90)
ax.text(10.05, 20.22, "unclassified ->\nown method,\nflagged;\n<20 lines ->\nnarrow panel",
        fontsize=6.0, ha="left", va="center", color="#8a5a00")

for x in (1.60, 5.00):
    arrow(x, 21.45, x, 18.60)
arrow(8.40, 19.55, 8.40, 18.60)

# --- shared trunk ------------------------------------------------------
steps = [
    (17.20, "Empirical scale detection", "data override metadata; conflict logged", DEC),
    (15.35, "Scale transformation", "linear -> log2(x + 1)", PROC),
    (13.50, "Replicate collapse", "median per cell line x method; retain n", PROC),
    (11.65, "Detection floor", "RNA-seq: 1 TPM   |   array: quantile q_floor", DEC),
    (9.80,  "Detection fraction  d(L,g)", "n at/above floor / n finite, per gene x lineage", PROC),
]
for y, t, s, c in steps:
    box(1.60, y, 6.80, 1.40, t, s, c)
for a, b in zip([18.60, 17.20, 15.35, 13.50, 11.65], [17.20, 15.35, 13.50, 11.65, 9.80]):
    arrow(5.00, a, 5.00, b + 1.40)

# --- three regimes -----------------------------------------------------
box(0.15, 7.50, 3.05, 1.55, "d >= q_upper", "robust z + shrinkage\nstatus: ok", OUT)
box(3.45, 7.50, 3.05, 1.55, "q_lower <= d < q_upper",
    "rank inverse normal\nstatus: censored_ranked", OUT)
box(6.75, 7.50, 3.10, 1.55, "d < q_lower",
    "withheld (NaN)\nbelow_detection / no_data", "#f6dcdc")
ax.text(5.00, 9.40, "regime selection", fontsize=7.4, ha="center",
        style="italic", color="#555555")
for x in (1.68, 4.95, 8.30):
    arrow(5.00, 9.80, x, 9.05)

ax.text(1.68, 7.00,
        "median / 1.4826xMAD\nlineage AND global\nshrink w=n_L/(n_L+n0)\nscale_source recorded",
        fontsize=6.2, ha="center", va="top", color="#333333")
ax.text(4.95, 7.00, "z = inv_norm((R-0.5)/n_fin)\nfloor ties share a midrank",
        fontsize=6.2, ha="center", va="top", color="#333333")
ax.text(8.30, 7.00, "absence is NOT zero\nreduces k, never votes",
        fontsize=6.2, ha="center", va="top", color="#8a2020")

# --- combination -------------------------------------------------------
box(1.60, 4.60, 6.80, 1.40, "Stage 1  -  within source",
    "Stouffer, weights sqrt(n)  ->  z_hpa_rna, z_depmap, z_geo", PROC)
for x in (1.68, 4.95, 8.30):
    arrow(x, 6.30, 5.00, 6.00)

box(1.60, 2.75, 6.80, 1.40, "Stage 2  -  across sources",
    "Stouffer, EQUAL weight (one source, one vote)", PROC)
arrow(5.00, 4.60, 5.00, 4.15)
ax.text(8.60, 3.45, "method count affects\nPRECISION, not\nINFLUENCE", fontsize=6.2,
        ha="left", va="center", color="#8a5a00")

# --- calibration feedback ---------------------------------------------
box(-1.05, 12.15, 2.35, 1.60, "Calibration",
    "held-out\ncross-source rho\n1-SE rule", DEC, ls=(0, (3, 2)))
ax.annotate("", xy=(1.60, 12.35), xytext=(1.30, 12.95),
            arrowprops=dict(arrowstyle="-|>", lw=0.9, color="#8a5a00",
                            linestyle="dashed"))
ax.annotate("", xy=(1.60, 10.50), xytext=(0.12, 12.15),
            arrowprops=dict(arrowstyle="-|>", lw=0.9, color="#8a5a00",
                            connectionstyle="angle,angleA=-90,angleB=180,rad=6",
                            linestyle="dashed"))
ax.text(-1.00, 10.9, "selects q_floor, n0,\nq_upper, q_lower", fontsize=6.0,
        ha="left", va="center", color="#8a5a00")

# --- output ------------------------------------------------------------
box(1.60, 0.75, 6.80, 1.45, "Transcriptomic evidence score",
    "z_lineage | z_global | k_src | k_methods\ndet_frac | status | scale_source | per-source z",
    OUT)
arrow(5.00, 2.75, 5.00, 2.20)

fig.savefig("figs/fig1_transcriptomic_pipeline.png", dpi=300, bbox_inches="tight")
fig.savefig("figs/fig1_transcriptomic_pipeline.pdf", bbox_inches="tight")
print("wrote figs/fig1_transcriptomic_pipeline.{png,pdf}")