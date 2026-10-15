#!/usr/bin/env python3
"""Figure 1A: candidate-identification workflow.

Redraw of the submitted panel. Two changes of substance: the inward extension
and the deduplication step are now shown, and the decision chain is reduced from
four diamonds to two, which is what made the submitted panel dense.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch, Polygon

OUT = Path(__file__).resolve().parents[2] / "analysis" / "results" / "plots"

PROC_FC, PROC_EC = "#CFE3F5", "#2166AC"     # process step
DEC_FC, DEC_EC = "#FFE9B0", "#C8860D"       # decision
TERM_FC, TERM_EC = "#E2D5F0", "#6A51A3"     # start / stop
ARROW = "#444444"
FS = 9.2


def box(ax, x, y, w, h, text, fc=PROC_FC, ec=PROC_EC, rounding=0.035):
    ax.add_patch(FancyBboxPatch((x - w / 2, y - h / 2), w, h,
                                boxstyle=f"round,pad=0,rounding_size={rounding}",
                                fc=fc, ec=ec, lw=1.5, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=FS, zorder=4, linespacing=1.45)
    return (x, y, w, h)


def diamond(ax, x, y, w, h, text):
    ax.add_patch(Polygon([(x, y + h / 2), (x + w / 2, y), (x, y - h / 2), (x - w / 2, y)],
                         closed=True, fc=DEC_FC, ec=DEC_EC, lw=1.5, zorder=3))
    ax.text(x, y, text, ha="center", va="center", fontsize=FS - 0.5, zorder=4, linespacing=1.4)
    return (x, y, w, h)


def term(ax, x, y, text):
    return box(ax, x, y, 0.62, 0.42, text, TERM_FC, TERM_EC, rounding=0.2)


def arrow(ax, p0, p1, label=None, rad=0.0, lx=0, ly=0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle="-|>", mutation_scale=13,
                                 lw=1.4, color=ARROW, zorder=2,
                                 connectionstyle=f"arc3,rad={rad}",
                                 shrinkA=0, shrinkB=0))
    if label:
        mx, my = (p0[0] + p1[0]) / 2 + lx, (p0[1] + p1[1]) / 2 + ly
        ax.text(mx, my, label, fontsize=FS - 0.8, color=ARROW, ha="center",
                va="center", zorder=5,
                bbox=dict(fc="white", ec="none", pad=1.2))


def main() -> None:
    fig, ax = plt.subplots(figsize=(13.2, 4.3))
    ax.set_xlim(0, 13.2); ax.set_ylim(0, 4.3); ax.axis("off")

    yt, yb = 3.35, 1.30          # top and bottom row centres

    # ---- top row ----------------------------------------------------------
    s = term(ax, 0.75, yt, "Start")
    p1 = box(ax, 2.75, yt, 2.55, 0.95,
             "Precompute purity prefix sums\nand next non-ACGT arrays")
    p2 = box(ax, 5.85, yt, 2.25, 0.95,
             "Next centre $c$,\nspacer length $s$")
    d1 = diamond(ax, 8.75, yt, 2.45, 1.30,
                 "Max reachable arm\n$\\geq$ min arm, and\npurity attainable?")
    p3 = box(ax, 11.85, yt, 2.25, 1.05,
             "Extend outward:\ncompare $S[l]$, $S[r]$;\ntrack mismatches, purity")

    arrow(ax, (s[0] + s[2] / 2, yt), (p1[0] - p1[2] / 2, yt))
    arrow(ax, (p1[0] + p1[2] / 2, yt), (p2[0] - p2[2] / 2, yt))
    arrow(ax, (p2[0] + p2[2] / 2, yt), (d1[0] - d1[2] / 2, yt))
    arrow(ax, (d1[0] + d1[2] / 2, yt), (p3[0] - p3[2] / 2, yt), "yes", ly=0.22)

    # ---- top to bottom ----------------------------------------------------
    d2 = diamond(ax, 11.85, yb, 2.35, 1.25, "Valid candidate\nfound?")
    arrow(ax, (p3[0], yt - p3[3] / 2), (d2[0], yb + d2[3] / 2))

    # ---- bottom row, right to left ---------------------------------------
    p4 = box(ax, 8.60, yb, 2.75, 1.05,
             "Inward extension:\nabsorb matching pairs\nat the spacer ends")
    p5 = box(ax, 5.55, yb, 2.25, 0.95, "Add candidate\nto list $C$")

    arrow(ax, (d2[0] - d2[2] / 2, yb), (p4[0] + p4[2] / 2, yb), "yes", ly=0.22)
    arrow(ax, (p4[0] - p4[2] / 2, yb), (p5[0] + p5[2] / 2, yb))

    # ---- loop back to the centre/spacer iterator --------------------------
    arrow(ax, (p5[0], yb + p5[3] / 2), (p2[0], yt - p2[3] / 2))
    ax.text(p5[0] - 0.12, (yb + yt) / 2, "next $c$, $s$", fontsize=FS - 0.8,
            color=ARROW, ha="right", va="center",
            bbox=dict(fc="white", ec="none", pad=1.2))

    # "no" branches rejoin the same loop
    ymerge = 2.30
    ax.add_patch(FancyArrowPatch((d1[0], yt - d1[3] / 2), (d1[0], ymerge),
                                 arrowstyle="-", lw=1.4, color=ARROW, zorder=2))
    ax.text(d1[0] - 0.22, yt - d1[3] / 2 - 0.26, "no", fontsize=FS - 0.8,
            color=ARROW, ha="right", va="center")
    arrow(ax, (d2[0] + d2[2] / 2, yb), (12.98, yb), "no", ly=0.22)
    ax.add_patch(FancyArrowPatch((12.98, yb), (12.98, ymerge),
                                 arrowstyle="-", lw=1.4, color=ARROW, zorder=2))
    ax.add_patch(FancyArrowPatch((12.98, ymerge), (d1[0], ymerge),
                                 arrowstyle="-", lw=1.4, color=ARROW, zorder=2))
    ax.add_patch(FancyArrowPatch((d1[0], ymerge), (p2[0] + 0.45, ymerge),
                                 arrowstyle="-", lw=1.4, color=ARROW, zorder=2))
    arrow(ax, (p2[0] + 0.45, ymerge), (p2[0] + 0.45, yt - p2[3] / 2))

    # ---- exit path --------------------------------------------------------
    p6 = box(ax, 2.55, yb, 2.75, 0.92,
             "Deduplicate candidates\nby (start, end, arm)")
    p7 = box(ax, 2.55, 0.30, 2.75, 0.50, "Return candidate list $C$")

    arrow(ax, (p5[0] - p5[2] / 2, yb), (p6[0] + p6[2] / 2, yb),
          "all $c$, $s$ done", ly=0.24)
    arrow(ax, (p6[0], yb - p6[3] / 2), (p7[0], p7[1] + p7[3] / 2))
    t2 = term(ax, 5.55, 0.30, "Stop")
    arrow(ax, (p7[0] + p7[2] / 2, 0.30), (t2[0] - t2[2] / 2, 0.30))

    fig.text(0.008, 0.985, "A", fontsize=17, fontweight="bold", va="top")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"figure1_panelA_workflow.{ext}", dpi=300,
                    bbox_inches="tight", facecolor="white")
    print("wrote", OUT / "figure1_panelA_workflow.pdf")


if __name__ == "__main__":
    main()
