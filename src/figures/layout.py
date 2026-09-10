"""
A tiny inch-space layout helper for the schematic figures.

Matplotlib's default axes coordinates are fraction-of-axes, which makes text
overflow impossible to reason about: a box height of 0.09 means nothing until
you know the figure is 15 inches tall. Everything here works in inches with
the origin at the bottom-left of the figure, so a box can size itself to its
own text and boxes can be stacked with a known gap.
"""
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

PT = 1 / 72.0                      # points -> inches

TITLE_LEAD = 1.65                  # line height as a multiple of font size
BODY_LEAD = 1.55
PAD_TOP = 0.10
PAD_BOT = 0.11
PAD_X = 0.10


def canvas(width, height):
    """Axes whose coordinate system is inches, origin bottom-left."""
    fig = plt.figure(figsize=(width, height))
    ax = fig.add_axes([0, 0, 1, 1])
    ax.set_xlim(0, width)
    ax.set_ylim(0, height)
    ax.axis("off")
    ax.grid(False)
    ax.set_facecolor("white")
    return fig, ax


def box_height(n_lines, title_size=9.5, body_size=7.6, title=True):
    h = PAD_TOP + PAD_BOT
    if title:
        h += title_size * TITLE_LEAD * PT
    h += n_lines * body_size * BODY_LEAD * PT
    return h


def box(ax, x, y_top, w, title, lines, *, fill="#e8f0fa", edge="#2a6fb5",
        title_size=9.5, body_size=7.6, lw=1.4, ls="solid", align="center",
        title_colour=None, body_colour="#3d4550"):
    """Draw a self-sizing box with its TOP edge at y_top.
    Returns (y_bottom, y_centre, height)."""
    h = box_height(len(lines), title_size, body_size, title is not None)
    y = y_top - h
    ax.add_patch(FancyBboxPatch(
        (x, y), w, h, boxstyle="round,pad=0.02,rounding_size=0.05",
        facecolor=fill, edgecolor=edge, linewidth=lw, linestyle=ls, zorder=2))

    if align == "center":
        tx, ha = x + w / 2, "center"
    else:
        tx, ha = x + PAD_X, "left"

    cursor = y_top - PAD_TOP
    if title is not None:
        ax.text(tx, cursor, title, ha=ha, va="top", fontsize=title_size,
                fontweight="bold", color=title_colour or "#22262b", zorder=3)
        cursor -= title_size * TITLE_LEAD * PT
    for ln in lines:
        ax.text(tx, cursor, ln, ha=ha, va="top", fontsize=body_size,
                color=body_colour, zorder=3)
        cursor -= body_size * BODY_LEAD * PT
    return y, y + h / 2, h


def arrow(ax, p0, p1, *, colour="#2a6fb5", lw=1.6, style="-|>", rad=0.0,
          zorder=1):
    ax.add_patch(FancyArrowPatch(
        p0, p1, arrowstyle=style, mutation_scale=13, color=colour,
        linewidth=lw, connectionstyle=f"arc3,rad={rad}", zorder=zorder,
        shrinkA=1, shrinkB=1))


def arrow_label(ax, p0, p1, label, sub=None, *, dx=0.14, ha="left",
                fs=7.6, sub_fs=7.0, colour="#22262b", sub_colour="#7a7a7a"):
    mx = (p0[0] + p1[0]) / 2 + dx
    my = (p0[1] + p1[1]) / 2
    off = 0.075 if sub else 0.0
    ax.text(mx, my + off, label, ha=ha, va="center", fontsize=fs,
            color=colour, family="DejaVu Sans Mono", zorder=4,
            bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="none",
                      alpha=0.95))
    if sub:
        ax.text(mx, my - off - 0.02, sub, ha=ha, va="center", fontsize=sub_fs,
                color=sub_colour, style="italic", zorder=4)


def text(ax, x, y, s, **kw):
    kw.setdefault("va", "top")
    kw.setdefault("ha", "left")
    return ax.text(x, y, s, zorder=4, **kw)
