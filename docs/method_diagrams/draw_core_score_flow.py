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
CORE     = "#c792ea"   # core_score accent (distinct from RNA teal / protein amber / pipeline blue)

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
        self.y = 0.0
        self.ops = []

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
        ty = self.y - 0.04
        for ln in lines:
            self.ops.append(("note", self.xc, ty, ln))
            ty -= NOTE_FS / 72 * 1.55
        self.y -= h

    def branch(self, cols, accent, wrap_w, col_w, gap):
        n = len(cols)
        total_w = col_w * n + gap * (n - 1)
        x0 = self.xc - total_w / 2
        heights = [box_height(t, b, wrap_w)[0] for t, b in cols]
        h = max(heights)

        bar_y_top = self.y
        self.ops.append(("hbar", x0 + col_w / 2, x0 + total_w - col_w / 2, bar_y_top))
        y_box_top = bar_y_top - BAR_STEM
        for i in range(n):
            cx = x0 + col_w / 2 + i * (col_w + gap)
            self.ops.append(("line", cx, bar_y_top, y_box_top))
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
# core_score
# =====================================================================
core = Flow(width=6.2)

core.box("Two independent inputs",
        "RNA z-score (Stage 1) and protein z-score (Stage 2), per gene, per cell line",
        CORE, wrap_w=50, io=True)
core.stem()

core.branch([
    ("RNA standardization", "Per (gene, n_sources) stratum; thin strata shrunk toward the gene-global scale; winsorized"),
    ("Protein residualization", "Protein regressed on RNA (per-gene shrunk slope); the residual is standardized per gene"),
], CORE, wrap_w=24, col_w=2.9, gap=0.3)
core.stem()

core.box("Combine RNA + protein",
        "Weighted by each arm's effective source count (Kish); mapped through the normal CDF onto [0,1]",
        CORE, wrap_w=50)
core.stem()

core.branch([
    ("Both layers scored", "core_score from the combined RNA + protein estimate"),
    ("Protein-free gene/line", "core_score from RNA alone (fallback, not dropped)"),
], CORE, wrap_w=22, col_w=2.9, gap=0.3)
core.stem()

core.box("Rank within lineage",
        "Cell lines ranked by core_score within each (gene, lineage) stratum",
        CORE, wrap_w=50)
core.stem()

core.box("core_score",
        "Final combined score, per gene, per cell line",
        CORE, wrap_w=50, io=True)

core.render("core_score_flow.png", "Combined Core Score", CORE, fig_w=7.2)
