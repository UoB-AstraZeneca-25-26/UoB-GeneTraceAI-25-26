import textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

# ---- palette (matches the HTML artifact) ----
BG       = "#0a0d11"
BOX_BG   = "#121820"
BOX_BRD  = "#3a4652"
LINE     = "#3c4753"
TEXT     = "#dbe2e9"
TEXT_DIM = "#9aa5b1"
RNA      = "#72d9be"
PROTEIN  = "#e8ab68"

MONO = "DejaVu Sans Mono"
SANS = "DejaVu Sans"

PAD_X, PAD_Y = 0.22, 0.16
TITLE_FS, BODY_FS, NOTE_FS = 10.5, 9.3, 8.3
TITLE_H, LINE_H, GAP_TB = 0.20, 0.185, 0.09
STEM_H = 0.30
BAR_STEM = 0.22


def wrap(text, width):
    return textwrap.wrap(text, width=width)


def box_height(title, body, wrap_w):
    lines = wrap(body, wrap_w)
    return PAD_Y * 2 + TITLE_H + GAP_TB + max(len(lines), 1) * LINE_H, lines


def draw_box(ax, x, y, w, h, title, lines, accent, io=False):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0,rounding_size=0.045",
        linewidth=1.3 if io else 0.9,
        edgecolor=accent if io else BOX_BRD,
        facecolor=BOX_BG, zorder=3))
    ax.text(x + PAD_X, y + h - PAD_Y, title.upper(), fontfamily=MONO,
            fontsize=TITLE_FS, fontweight="bold", color=accent,
            ha="left", va="top", zorder=4)
    ty = y + h - PAD_Y - TITLE_H
    for ln in lines:
        ax.text(x + PAD_X, ty, ln, fontfamily=SANS, fontsize=BODY_FS,
                 color=TEXT, ha="left", va="top", zorder=4)
        ty -= LINE_H


def vline(ax, xc, y0, y1):
    ax.plot([xc, xc], [y0, y1], color=LINE, lw=1.3, zorder=2)


class Flow:
    def __init__(self, width, x0=0.0):
        self.width = width
        self.x0 = x0
        self.xc = x0 + width / 2
        self.y = 0.0   # current cursor, grows downward (negative)
        self.ops = []  # (kind, args) deferred draw calls, y already in final coords

    def stem(self, h=STEM_H):
        y0 = self.y
        self.y -= h
        self.ops.append(("line", self.xc, y0, self.y))

    def box(self, title, body, accent, wrap_w=40, io=False, w=None):
        w = w or self.width
        h, lines = box_height(title, body, wrap_w)
        y0 = self.y - h
        self.ops.append(("box", self.x0, y0, w, h, title, lines, accent, io))
        self.y = y0

    def note(self, text, wrap_w=70):
        lines = wrap(text, wrap_w)
        h = len(lines) * (NOTE_FS / 72 * 1.55) + 0.08
        y0 = self.y - h
        ty = self.y - 0.04
        for ln in lines:
            self.ops.append(("note", self.xc, ty, ln))
            ty -= NOTE_FS / 72 * 1.55
        self.y = y0

    def branch(self, cols, accent, wrap_w, col_w, gap, io_last=False):
        n = len(cols)
        total_w = col_w * n + gap * (n - 1)
        x0 = self.xc - total_w / 2
        heights = [box_height(t, b, wrap_w)[0] for t, b in cols]
        h = max(heights)

        bar_y_top = self.y
        self.ops.append(("hbar", x0 + col_w / 2, x0 + total_w - col_w / 2, bar_y_top))
        y_stem_top = bar_y_top
        y_box_top = bar_y_top - BAR_STEM
        for i in range(n):
            cx = x0 + col_w / 2 + i * (col_w + gap)
            self.ops.append(("line", cx, y_stem_top, y_box_top))
        y_box_bot = y_box_top - h
        for i, (title, body) in enumerate(cols):
            cx0 = x0 + i * (col_w + gap)
            _, lines = box_height(title, body, wrap_w)
            self.ops.append(("box", cx0, y_box_bot, col_w, h, title, lines, accent, False))
        y_stem_bot = y_box_bot - BAR_STEM
        for i in range(n):
            cx = x0 + col_w / 2 + i * (col_w + gap)
            self.ops.append(("line", cx, y_box_bot, y_stem_bot))
        self.ops.append(("hbar", x0 + col_w / 2, x0 + total_w - col_w / 2, y_stem_bot))
        self.y = y_stem_bot

    def render(self, filename, title, accent, fig_w):
        total_h = -self.y + 0.5
        fig, ax = plt.subplots(figsize=(fig_w, total_h + 0.6))
        fig.patch.set_facecolor(BG)
        ax.set_facecolor(BG)
        ax.set_xlim(self.x0 - 0.3, self.x0 + self.width + 0.3)
        ax.set_ylim(self.y - 0.2, 0.75)
        ax.axis("off")

        ax.text(self.xc, 0.35, title, fontfamily=SANS, fontsize=15,
                 fontweight="bold", color=TEXT, ha="center", va="bottom")

        for op in self.ops:
            if op[0] == "line":
                _, xc, y0, y1 = op
                vline(ax, xc, y0, y1)
            elif op[0] == "hbar":
                _, x0, x1, y = op
                ax.plot([x0, x1], [y, y], color=LINE, lw=1.3, zorder=2)
            elif op[0] == "box":
                _, x, y, w, h, ttl, lines, acc, io = op
                draw_box(ax, x, y, w, h, ttl, lines, acc, io)
            elif op[0] == "note":
                _, xc, y, text = op
                ax.text(xc, y, text, fontfamily=MONO, fontsize=NOTE_FS,
                         color=TEXT_DIM, ha="center", va="top")

        fig.tight_layout(pad=0.4)
        fig.savefig(filename, dpi=200, facecolor=BG)
        plt.close(fig)
        print("wrote", filename, "size", fig_w, "x", total_h + 0.6)


# =====================================================================
# RNA
# =====================================================================
rna = Flow(width=6.4)
rna.box("Raw expression, per source",
        "DepMap RNA-seq · HPA RNA-seq · GEO microarray, split into one method per platform x processing pipeline",
        RNA, wrap_w=52, io=True)
rna.stem()
rna.box("Collapse replicates",
        "Median across replicate samples per cell line",
        RNA, wrap_w=52)
rna.stem()
rna.box("Scale & detection floor, per method",
        "Fixed floor at ~1 TPM for RNA-seq; empirical floor for microarray methods; all values on a common log2 scale",
        RNA, wrap_w=52)
rna.stem()
rna.box("Detection fraction, per gene within each lineage",
        "Share of a lineage's cell lines at or above the detection floor",
        RNA, wrap_w=52)
rna.stem()
rna.branch([
    ("Well detected", "Robust z-score (median/MAD); lineage scale shrunk toward the method-wide scale"),
    ("Partially detected", "Rank-based inverse-normal transform, within the lineage"),
    ("Below detection", "Excluded — not scored as zero"),
], RNA, wrap_w=21, col_w=2.05, gap=0.15)
rna.note("Lineages with too few scored cell lines are excluded outright, regardless of detection fraction.")
rna.stem()
rna.box("Combine GEO's methods",
        "Weighted by sample size",
        RNA, wrap_w=52)
rna.stem()
rna.box("Combine the three sources",
        "DepMap, HPA, GEO — equal weight, one source one vote",
        RNA, wrap_w=52)
rna.stem()
rna.box("RNA z-score",
        "Final RNA expression z-score, per gene, per cell line",
        RNA, wrap_w=52, io=True)
rna.render("rna_flow.png", "Transcript Expression Score", RNA, fig_w=7.4)

# =====================================================================
# Protein
# =====================================================================
prot = Flow(width=5.4)
prot.box("Two independent platforms",
        "ProCAN and CCLE proteomics — neither treated as the reference",
        PROTEIN, wrap_w=46, io=True)
prot.stem()
prot.box("Derive agreement thresholds from the data",
        "Matched correlations compared against a null of mismatched proteins/genes",
        PROTEIN, wrap_w=46)
prot.stem()
prot.branch([
    ("ProCAN — isoforms", "Concordant forms: max abundance. Discordant: mean"),
    ("CCLE — isoforms", "Concordant forms: max abundance. Discordant: mean"),
], PROTEIN, wrap_w=24, col_w=2.55, gap=0.3)
prot.stem()
prot.box("Platform-agreement assessment, per gene",
        "ProCAN vs CCLE correlation on shared cell lines -> per-gene confidence score",
        PROTEIN, wrap_w=46)
prot.stem()
prot.box("Per-gene z-score, within each platform",
        "Robust (median/MAD), lineage-conditioned",
        PROTEIN, wrap_w=46)
prot.stem()
prot.box("Combine platforms",
        "Coverage-weighted average of the two platforms' z-scores",
        PROTEIN, wrap_w=46)
prot.stem()
prot.branch([
    ("Both platforms measured", "Scaled by the gene's agreement confidence"),
    ("One platform only", "Passed through unweighted"),
], PROTEIN, wrap_w=22, col_w=2.55, gap=0.3)
prot.stem()
prot.box("Protein z-score",
        "Final protein abundance z-score, per gene, per cell line",
        PROTEIN, wrap_w=46, io=True)
prot.render("protein_flow.png", "Protein Abundance Score", PROTEIN, fig_w=6.4)
