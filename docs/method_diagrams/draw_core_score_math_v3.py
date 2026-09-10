import textwrap
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

BG       = "#0a0d11"
BOX_BG   = "#121820"
BOX_BRD  = "#3a4652"
LINE     = "#3c4753"
TEXT     = "#dbe2e9"
TEXT_DIM = "#9aa5b1"
CORE     = "#c792ea"
FORMULA  = "#8fd0e8"   # distinct color for equation lines vs. plain gloss text

MONO = "DejaVu Sans Mono"
SANS = "DejaVu Sans"

PAD_X, PAD_Y = 0.22, 0.16
TITLE_FS, BODY_FS, FORMULA_FS, NOTE_FS = 10.5, 8.8, 10.0, 8.3
TITLE_H, LINE_H, FORMULA_LINE_H, GAP_TB = 0.20, 0.185, 0.205, 0.09
STEM_H = 0.30
BAR_STEM = 0.22


def wrap(text, width):
    return textwrap.wrap(text, width=width)


def box_height_lines(title, lines_spec):
    """lines_spec: list of ('formula'|'gloss', text)"""
    h = PAD_Y * 2 + TITLE_H + GAP_TB
    for kind, _ in lines_spec:
        h += FORMULA_LINE_H if kind == "formula" else LINE_H
    return h


def draw_box_lines(ax, x, y, w, h, title, lines_spec, accent, io=False):
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0,rounding_size=0.045",
        linewidth=1.3 if io else 0.9,
        edgecolor=accent if io else BOX_BRD,
        facecolor=BOX_BG, zorder=3))
    ax.text(x + PAD_X, y + h - PAD_Y, title.upper(), fontfamily=MONO,
            fontsize=TITLE_FS, fontweight="bold", color=accent,
            ha="left", va="top", zorder=4)
    ty = y + h - PAD_Y - TITLE_H
    for kind, text in lines_spec:
        if kind == "formula":
            ax.text(x + PAD_X, ty, text, fontfamily=MONO, fontsize=FORMULA_FS,
                     color=FORMULA, ha="left", va="top", zorder=4)
            ty -= FORMULA_LINE_H
        else:
            ax.text(x + PAD_X, ty, text, fontfamily=SANS, fontsize=BODY_FS,
                     color=TEXT_DIM, ha="left", va="top", zorder=4, style="italic")
            ty -= LINE_H


def vline(ax, xc, y0, y1):
    ax.plot([xc, xc], [y0, y1], color=LINE, lw=1.3, zorder=2)


class Flow:
    def __init__(self, width, x0=0.0):
        self.width = width
        self.x0 = x0
        self.xc = x0 + width / 2
        self.y = 0.0
        self.ops = []

    def stem(self, h=STEM_H):
        y0 = self.y
        self.y -= h
        self.ops.append(("line", self.xc, y0, self.y))

    def mbox(self, title, lines_spec, accent, io=False, w=None):
        """lines_spec: list of ('formula'|'gloss', text)"""
        w = w or self.width
        h = box_height_lines(title, lines_spec)
        y0 = self.y - h
        self.ops.append(("mbox", self.x0, y0, w, h, title, lines_spec, accent, io))
        self.y = y0

    def note(self, text, wrap_w=70):
        lines = wrap(text, wrap_w)
        h = len(lines) * (NOTE_FS / 72 * 1.55) + 0.08
        ty = self.y - 0.04
        for ln in lines:
            self.ops.append(("note", self.xc, ty, ln))
            ty -= NOTE_FS / 72 * 1.55
        self.y -= h

    def branch(self, cols, accent, col_w, gap):
        """cols: list of (title, lines_spec)"""
        n = len(cols)
        total_w = col_w * n + gap * (n - 1)
        x0 = self.xc - total_w / 2
        heights = [box_height_lines(t, ls) for t, ls in cols]
        h = max(heights)

        bar_y_top = self.y
        self.ops.append(("hbar", x0 + col_w / 2, x0 + total_w - col_w / 2, bar_y_top))
        y_box_top = bar_y_top - BAR_STEM
        for i in range(n):
            cx = x0 + col_w / 2 + i * (col_w + gap)
            self.ops.append(("line", cx, bar_y_top, y_box_top))
        y_box_bot = y_box_top - h
        for i, (title, lines_spec) in enumerate(cols):
            cx0 = x0 + i * (col_w + gap)
            self.ops.append(("mbox", cx0, y_box_bot, col_w, h, title, lines_spec, accent, False))
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
            elif op[0] == "mbox":
                _, x, y, w, h, ttl, lines_spec, acc, io = op
                draw_box_lines(ax, x, y, w, h, ttl, lines_spec, acc, io)
            elif op[0] == "note":
                _, xc, y, text = op
                ax.text(xc, y, text, fontfamily=MONO, fontsize=NOTE_FS,
                         color=TEXT_DIM, ha="center", va="top")

        fig.tight_layout(pad=0.4)
        fig.savefig(filename, dpi=200, facecolor=BG)
        plt.close(fig)
        print("wrote", filename, "size", fig_w, "x", total_h + 0.6)


F = "formula"
G = "gloss"

# =====================================================================
core = Flow(width=6.6)

core.mbox("Inputs", [
    (F, "rna_z(gene, model, n_sources)      [Stage 1]"),
    (F, "prot_z(gene, model)                [Stage 2]"),
], CORE, io=True)
core.stem()

core.mbox("RNA standardization  —  per (gene, n_sources) stratum", [
    (F, "w = n / (n + n0),   n0 = 25"),
    (F, "σ′ = √(w·σ² + (1−w)·σg²)          if n < 5"),
    (F, "μ′ = w·μ + (1−w)·μg              if n < 5"),
    (F, "rna_std_z = clip((rna_z − μ′) / σ′,  −5, 5)"),
], CORE)
core.stem()

core.mbox("Protein residualization  —  per gene", [
    (F, "λ = K / (n + K),   K = 30"),
    (F, "ρg = λ·ρ0 + (1−λ)·ρraw,     ρ0 = 0.353"),
    (F, "βg = ρg · (σprot / σrna)"),
    (F, "resid = prot_z − βg·rna_z"),
], CORE)
core.stem()

core.mbox("Protein residual shrinkage + clip  —  per gene  (v3)", [
    (F, "w = n / (n + n0),   n0 = 25"),
    (F, "σ' = √(w·σresid² + (1−w)·σglobal²)"),
    (G, "σglobal = residual SD pooled across ALL genes"),
    (F, "prot_resid_z = clip((resid − median) / σ',  −5, 5)"),
], CORE)
core.stem()

core.mbox("Combine  —  Kish-weighted", [
    (F, "W_rna = √1.431,   W_prot = √1.45"),
    (F, "core_z = (W_rna·rna_std_z + W_prot·prot_resid_z)"),
    (F, "          ÷ √(W_rna² + W_prot²)"),
    (G, "Φ = standard normal CDF"),
], CORE)
core.stem()

core.branch([
    ("Both layers  (n≥ 2)", [(F, "core_score = Φ(core_z)")]),
    ("RNA only  (fallback)", [(F, "core_score = Φ(rna_std_z)")]),
], CORE, col_w=3.0, gap=0.3)
core.stem()

core.mbox("Rank within lineage", [
    (F, "rank = rank_desc(core_score)  within (gene, lineage)"),
    (G, "ties broken by model_id, ascending"),
], CORE)
core.stem()

core.mbox("core_score", [
    (F, "core_score ∈ [0, 1],  per (gene, cell line)"),
], CORE, io=True)

core.render("core_score_math_v3.png", "core_score — Mathematical Formulation (v3, live)", CORE, fig_w=7.6)
