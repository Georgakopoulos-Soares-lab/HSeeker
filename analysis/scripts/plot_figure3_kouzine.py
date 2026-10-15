#!/usr/bin/env python3
"""Figure 3: HSeeker predictions vs single-stranded DNA mapped by Kouzine et al. 2017.

Replaces the former Figure 3 (E. coli insertion recovery and in-sample score threshold).
Panels:
  A  Recall: share of the authors' formed (ssDNA+) H-DNA that an HSeeker locus overlaps.
  B  Selectivity: odds that one of the authors' predicted H-DNA motifs is formed when HSeeker
     calls it vs when it does not (all analysed motifs, and mappable motifs only).
  C  Single-stranded signal at HSeeker loci relative to composition-matched controls:
     all loci with genome-wide controls (pre-planned), loci isolated from other sources of
     single-stranded DNA with local controls, the authors' formed motifs (positive control) and
     controls against controls (placebo). Log scale; 1 = no difference.
  D  Isolated loci by mirror-repeat arm length.
All values: minimum arm 8 (shipped defaults), 95 % CIs.

Inputs: analysis/results/kouzine_v1/published/{mm9,hg19}.cli_minrep8{,.p4}.json
Outputs: analysis/results/plots/figure3_kouzine.{pdf,png,svg} and figure3_kouzine_data.csv

Usage:  python3 analysis/scripts/plot_figure3_kouzine.py
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import NullFormatter  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
RES = ROOT / "analysis" / "results" / "kouzine_v1" / "published"
OUT = ROOT / "analysis" / "results" / "plots"

# Validated categorical slots (light surface): mouse = slot 1 blue, human = slot 2 orange;
# resting B cells = a lighter step of the mouse hue.
MOUSE, MOUSE_REST, HUMAN = "#2a78d6", "#86b6ef", "#eb6834"
INK, INK2, GRID, SURFACE = "#0b0b0b", "#52514e", "#e4e3df", "#fcfcfb"
KEY = "tier_B | local | C"
PC_KEY = "positive_control_tier_B | local | C"
ARM = ["8-9", "10-11", "12-14", "15-19", ">=20"]


def load(name):
    return json.loads((RES / name).read_text())


def style(ax):
    ax.set_facecolor(SURFACE)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(INK2)
        ax.spines[s].set_linewidth(0.8)
    ax.tick_params(colors=INK2, labelcolor=INK, width=0.8, labelsize=8)
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.grid(axis="x", color=GRID, linewidth=0.6)
    ax.set_axisbelow(True)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    m, h = load("mm9.cli_minrep8.json"), load("hg19.cli_minrep8.json")
    mp4, hp4 = load("mm9.cli_minrep8.p4.json"), load("hg19.cli_minrep8.p4.json")
    rows = []  # table view

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 8, "text.color": INK,
                         "axes.labelcolor": INK, "figure.facecolor": SURFACE})
    fig = plt.figure(figsize=(7.2, 6.2))
    gs = fig.add_gridspec(2, 2, width_ratios=[1, 1.35], height_ratios=[1, 1.25],
                          hspace=0.55, wspace=0.55)

    # ---- A: recall
    ax = fig.add_subplot(gs[0, 0])
    style(ax)
    vals = [(f"Mouse\nactivated B", m["P1_recall"], MOUSE), ("Human\nRaji", h["P1_recall"], HUMAN)]
    for y, (lab, p1, col) in enumerate(vals):
        frac = 100 * p1["fraction"]
        ax.barh(y, frac, height=0.5, color=col, edgecolor=SURFACE, linewidth=2)
        ax.text(frac + 1.5, y, f"{frac:.1f} %\n({p1['overlapped_by_hseeker']:,} of {p1['called_ssDNA_plus']:,})",
                va="center", fontsize=7.5, color=INK)
        rows.append(dict(panel="A recall", group=lab.replace("\n", " "), value=round(frac, 2),
                         ci_low="", ci_high="", n=p1["called_ssDNA_plus"]))
    ax.set_yticks(range(len(vals)), [v[0] for v in vals])
    ax.invert_yaxis()
    ax.set_xlim(0, 140)
    ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Formed H-DNA found by HSeeker (%)")
    ax.set_title("A  Recall", loc="left", fontsize=9, fontweight="bold", color=INK)

    # ---- B: selectivity odds ratios
    ax = fig.add_subplot(gs[0, 1])
    style(ax)
    items = [("Mouse: all motifs", m["P2_selectivity"], MOUSE),
             ("Mouse: mappable", m["P2_selectivity_mappable_sensitivity"]["mean_ge_0.5"], MOUSE),
             ("Human: all motifs", h["P2_selectivity"], HUMAN),
             ("Human: mappable", h["P2_selectivity_mappable_sensitivity"]["mean_ge_0.5"], HUMAN)]
    for y, (lab, p2, col) in enumerate(items):
        o, (lo, hi) = p2["odds_ratio"], p2["or_ci95"]
        ax.plot([lo, hi], [y, y], color=col, linewidth=2, solid_capstyle="round")
        ax.plot(o, y, "o", color=col, markersize=7, markeredgecolor=SURFACE, markeredgewidth=1.5)
        ax.text(hi * 1.06, y, f"{o:.2f}", va="center", fontsize=7.5, color=INK)
        rows.append(dict(panel="B selectivity odds ratio", group=lab, value=round(o, 3),
                         ci_low=round(lo, 3), ci_high=round(hi, 3), n=p2["analysed_motifs"]))
    ax.axvline(1, color=INK2, linestyle=(0, (3, 3)), linewidth=0.8)
    ax.set_xscale("log")
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(0.8, 8)
    ax.set_xticks([1, 2, 4, 8], ["1", "2", "4", "8"])
    ax.set_yticks(range(len(items)), [i[0] for i in items])
    ax.invert_yaxis()
    ax.set_xlabel("Odds ratio of formation, called vs not called")
    ax.set_title("B  Selectivity", loc="left", fontsize=9,
                 fontweight="bold", color=INK)

    # ---- C: signal ratios
    ax = fig.add_subplot(gs[1, 0])
    style(ax)
    pre_m, pre_h = m["P3_enrichment_authors_signal"], h["P3_enrichment_authors_signal"]
    c_items = []
    for lab_cell, col, pre, p4, key in (("Mouse activated", MOUSE, pre_m["activated"], mp4["signal"]["activated"], "activated"),
                                        ("Mouse resting", MOUSE_REST, pre_m["resting"], mp4["signal"]["resting"], "resting"),
                                        ("Human Raji", HUMAN, pre_h["raji"], hp4["signal"]["raji"], "raji")):
        c_items.append((f"{lab_cell}: all loci", col, pre["asymmetry_matched"]["ratio_of_means"],
                        pre["asymmetry_matched"]["ratio_ci95"], "o"))
        c_items.append((f"{lab_cell}: isolated loci", col, p4[KEY]["ratio_of_means"], p4[KEY]["ratio_ci95"], "D"))
        c_items.append((f"{lab_cell}: authors' formed motifs", col, p4[PC_KEY]["ratio_of_means"],
                        p4[PC_KEY]["ratio_ci95"], "^"))
        c_items.append((f"{lab_cell}: placebo", col, p4[KEY]["placebo"]["ratio_of_means"],
                        p4[KEY]["placebo"]["ratio_ci95"], "s"))
    y = 0
    yt, yl = [], []
    for i, (lab, col, r, (lo, hi), mk) in enumerate(c_items):
        if i and i % 4 == 0:
            y += 0.6
        filled = mk != "s"
        ax.plot([lo, hi], [y, y], color=col, linewidth=2, solid_capstyle="round")
        ax.plot(r, y, mk, color=col, markersize=6, markerfacecolor=col if filled else SURFACE,
                markeredgecolor=col if not filled else SURFACE, markeredgewidth=1.2)
        yt.append(y)
        yl.append(lab.split(": ")[1])
        rows.append(dict(panel="C signal ratio, loci / matched controls", group=lab, value=round(r, 4),
                         ci_low=round(lo, 4), ci_high=round(hi, 4), n=""))
        y += 1
    ax.axvline(1, color=INK2, linestyle=(0, (3, 3)), linewidth=0.8)
    ax.set_xscale("log")
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(0.5, 12)
    ax.set_xticks([0.5, 1, 2, 4, 8], ["0.5", "1", "2", "4", "8"])
    ax.set_yticks(yt, yl, fontsize=7)
    ax.invert_yaxis()
    for start, lab, col in ((0, "Mouse, activated B", MOUSE), (4.6, "Mouse, resting B", MOUSE_REST),
                            (9.2, "Human, Raji", HUMAN)):
        ax.text(11.5, start - 0.1, lab, ha="right", va="center", fontsize=7.5, color=INK, fontweight="bold")
    ax.set_xlabel("Single-stranded signal, loci / matched controls")
    ax.set_title("C  Signal at HSeeker loci", loc="left", fontsize=9, fontweight="bold", color=INK)

    # ---- D: arm length
    ax = fig.add_subplot(gs[1, 1])
    style(ax)
    ax.grid(axis="y", color=GRID, linewidth=0.6)
    ax.grid(axis="x", visible=False)
    series = (("Mouse activated B", MOUSE, mp4["signal"]["activated"][KEY], -0.18),
              ("Mouse resting B", MOUSE_REST, mp4["signal"]["resting"][KEY], 0.0),
              ("Human Raji", HUMAN, hp4["signal"]["raji"][KEY], 0.18))
    for lab, col, d, dx in series:
        xs, rs = [], []
        for j, a in enumerate(ARM):
            v = d["by_arm_length"][a]
            if not v.get("n"):
                continue
            lo, hi = v["ratio_ci95"]
            ax.plot([j + dx, j + dx], [lo, hi], color=col, linewidth=2, solid_capstyle="round")
            xs.append(j + dx)
            rs.append(v["ratio_of_means"])
            rows.append(dict(panel="D isolated loci by arm length", group=f"{lab}, arm {a} nt",
                             value=round(v["ratio_of_means"], 4), ci_low=round(lo, 4), ci_high=round(hi, 4),
                             n=v["n"]))
        ax.plot(xs, rs, "-o", color=col, linewidth=1.5, markersize=6, markeredgecolor=SURFACE,
                markeredgewidth=1.2, label=lab)
    ax.axhline(1, color=INK2, linestyle=(0, (3, 3)), linewidth=0.8)
    ax.set_yscale("log")
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.set_ylim(0.5, 4)
    ax.set_yticks([0.5, 1, 2, 4], ["0.5", "1", "2", "4"])
    ax.set_xticks(range(len(ARM)), ["8–9", "10–11", "12–14", "15–19", "≥20"])
    ax.set_xlabel("Mirror-repeat arm length (nt)")
    ax.set_ylabel("Signal, isolated loci / controls")
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.set_title("D  Dependence on arm length", loc="left", fontsize=9, fontweight="bold", color=INK)

    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"figure3_kouzine.{ext}", dpi=300, bbox_inches="tight")
    with open(OUT / "figure3_kouzine_data.csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=["panel", "group", "value", "ci_low", "ci_high", "n"])
        w.writeheader()
        w.writerows(rows)
    print("wrote", OUT / "figure3_kouzine.pdf")


if __name__ == "__main__":
    main()
