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
PIPE     = "#7fb3e8"   # pipeline-level accent (distinct from RNA teal / protein amber)

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
# End-to-end pipeline
# =====================================================================
pipe = Flow(width=6.6)

pipe.box("Warehouse check — Stage 0",
        "DuckDB warehouse exists (DepMap, GEO, HPA, cell-line metadata) — harmonised once upstream, not rebuilt here.",
        PIPE, wrap_w=54, io=True)
pipe.stem()

pipe.branch([
    ("RNA score — Stage 1", "Full method: see Transcript Expression Score"),
    ("Protein score — Stage 2", "Full method: see Protein Abundance Score"),
], PIPE, wrap_w=27, col_w=3.0, gap=0.3)
pipe.note("Stage 2 also derives a fixed ProCAN/CCLE tier per protein; checked as a precondition, not used in the combination below.")
pipe.stem()

pipe.box("Combine RNA + protein — core_score",
        "RNA standardized per gene x source-count; protein residualized against RNA, then standardized; combined at fixed reliability weights into one score per gene, per cell line.",
        PIPE, wrap_w=54)
pipe.stem()

pipe.branch([
    ("Mutations", "VEP + pathogenicity + burden -> p_mutation (Noisy-OR)"),
    ("Fusions", "Confidence + FFPM + recurrence -> p_fusion (Noisy-OR)"),
    ("Copy number", "COSMIC ploidy-aware calls, gated by gene role"),
], PIPE, wrap_w=20, col_w=2.1, gap=0.15)
pipe.stem()

pipe.box("Driver-gated routing",
        "core_score merged with mutation/fusion/CNA evidence; driver flag set when any channel clears its threshold.",
        PIPE, wrap_w=54)
pipe.stem()

pipe.box("Confidence tiers",
        "HIGH / MEDIUM / CONTEXT / LOW from score percentile x driver flag, direction-aware for loss-of-function alterations.",
        PIPE, wrap_w=54)
pipe.stem()

pipe.box("Evidence ledger — Stage 5",
        "One row per gene, per cell line, for the ranking/query layer — a projection, not new modelling.",
        PIPE, wrap_w=54)
pipe.stem()

pipe.box("Held-out evaluation — Stage 6",
        "hit@20 on a 20% gene hold-out: driver-gated ranking vs. score alone, against per-drug sensitivity labels.",
        PIPE, wrap_w=54)
pipe.stem()

pipe.box("Diagnostic tests — Stage 7",
        "Read-only checks on RNA/protein source agreement and the core_score distribution.",
        PIPE, wrap_w=54, io=True)

pipe.render("pipeline_flow.png", "GeneTraceAI Scoring Pipeline", PIPE, fig_w=7.6)
