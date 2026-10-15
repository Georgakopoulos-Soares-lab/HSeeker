#!/usr/bin/env python3
"""Figure 3: recovery of benchmark sequences embedded in a genomic background.

A) Concordance of per-record calls between scanning sequences on their own and
   scanning them embedded in the E. coli K-12 genome.
B) Detection performance in the genomic background.
C) Robustness to the positional overlap criterion used to match hits to inserts.
D) Fraction of each inserted sequence covered by the matching hit, with the
   background locus rate of the unmodified genome.

Inputs: analysis/results/sensitivity_injected_2026-10-15/ and the standalone
per-record calls in analysis/results/tool_comparison/.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "analysis" / "results" / "plots"

BLUE = "#2166ac"
ORANGE = "#e08214"
GREY = "#999999"
RED = "#b2182b"


def truthy(v: str) -> int:
    return 1 if str(v).strip().lower() in ("true", "1") else 0


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path,
                    default=ROOT / "analysis" / "results" / "sensitivity_injected_2026-10-15")
    args = ap.parse_args()
    d = args.dir

    # ---- inputs ------------------------------------------------------------
    standalone = {"HSeeker": {}, "triplex": {}}
    for r in csv.DictReader(open(ROOT / "analysis/results/tool_comparison/per_record_calls.csv")):
        if r["setting"] == "defaults" and r["tool"] in standalone:
            standalone[r["tool"]][r["record_id"]] = int(r["call"])

    hs_rows = list(csv.DictReader(open(d / "hseeker_insert_matches.csv")))
    tx_rows = list(csv.DictReader(open(d / "triplex_insert_matches.csv")))
    in_genome = {
        "HSeeker": {r["sequence_id"]: truthy(r["matched"]) for r in hs_rows},
        "triplex": {r["sequence_id"]: truthy(r["matched"]) for r in tx_rows},
    }
    labels = {r["sequence_id"]: r["label"] for r in hs_rows}

    metrics = {r["method"]: r for r in csv.DictReader(open(d / "summary_metrics.csv"))}
    sweep = list(csv.DictReader(open(d / "overlap_sweep.csv")))
    bg = json.loads((d / "background.json").read_text())

    fig, axes = plt.subplots(2, 2, figsize=(11.6, 8.0))
    fig.subplots_adjust(hspace=0.42, wspace=0.30, left=0.085, right=0.975,
                        top=0.92, bottom=0.085)

    # ---- A: call concordance ----------------------------------------------
    ax = axes[0, 0]
    tools = ["HSeeker", "triplex"]
    conc, lost = [], []
    for t in tools:
        common = set(standalone[t]) & set(in_genome[t])
        same = sum(1 for k in common if standalone[t][k] == in_genome[t][k])
        conc.append(100.0 * same / len(common))
        lost.append(sum(1 for k in common
                        if standalone[t][k] == 1 and in_genome[t][k] == 0))
    y = np.arange(len(tools))[::-1]
    ax.barh(y, conc, color=[BLUE, ORANGE], height=0.5, zorder=3)
    for yi, c, l, t in zip(y, conc, lost, tools):
        n = len(set(standalone[t]) & set(in_genome[t]))
        txt = f"{c:.1f}%  ({int(round(c*n/100))}/{n})"
        if l:
            txt += f"\n{l} detections lost"
        ax.annotate(txt, (c, yi), textcoords="offset points", xytext=(6, 0),
                    va="center", fontsize=9)
    ax.set_yticks(y); ax.set_yticklabels(tools)
    ax.set_xlim(0, 136); ax.set_xticks([0, 25, 50, 75, 100])
    ax.set_xlabel("Records whose call is unchanged when embedded (%)")
    ax.set_title("A   Calls on their own vs in a genome", loc="left", fontweight="bold")
    ax.grid(alpha=0.25, axis="x", ls=":")

    # ---- B: detection in genomic background --------------------------------
    ax = axes[0, 1]
    names = [("HSeeker overlap", "HSeeker", BLUE), ("Triplex overlap", "triplex", ORANGE)]
    stats = ["sensitivity", "specificity", "MCC"]
    width = 0.36
    x = np.arange(len(stats))
    for i, (key, lab, col) in enumerate(names):
        vals = [float(metrics[key][s]) for s in stats]
        ax.bar(x + (i - 0.5) * width, vals, width, label=lab, color=col, zorder=3)
        for xi, v in zip(x + (i - 0.5) * width, vals):
            ax.annotate(f"{v:.3f}", (xi, v), textcoords="offset points",
                        xytext=(0, 3), ha="center", fontsize=8.5)
    ax.set_xticks(x)
    ax.set_xticklabels(["Sensitivity", "Specificity", "MCC"])
    ax.set_ylim(0, 1.17); ax.set_yticks([0, 0.25, 0.5, 0.75, 1.0])
    ax.set_ylabel("Value")
    ax.set_title("B   Detection in the genomic background", loc="left", fontweight="bold")
    ax.legend(fontsize=9, frameon=False, loc="upper right", ncol=2)
    ax.grid(alpha=0.25, axis="y", ls=":")

    # ---- C: overlap threshold sweep ----------------------------------------
    ax = axes[1, 0]
    for tool, col, mark in (("HSeeker", BLUE, "o"), ("triplex", ORANGE, "s")):
        rows = [r for r in sweep if r["tool"] == tool]
        xs = [float(r["overlap_fraction"]) for r in rows]
        ys = [float(r["MCC"]) for r in rows]
        ax.plot(xs, ys, marker=mark, color=col, lw=2, ms=6, label=tool, zorder=3)
    ax.axvline(0.80, ls="--", color=GREY, lw=1.3)
    ax.annotate("criterion used", (0.80, 0.06), fontsize=8.5, color="#555555",
                rotation=90, ha="right", va="bottom")
    ax.set_xlabel("Minimum overlap required to match a hit to an insert")
    ax.set_ylabel("MCC")
    ax.set_ylim(0, 1.08)
    ax.set_xticks([0.5, 0.6, 0.7, 0.8, 0.9, 0.95])
    ax.set_title("C   Robustness to the matching criterion", loc="left", fontweight="bold")
    ax.legend(fontsize=9, frameon=False, loc="lower left", ncol=2)
    ax.grid(alpha=0.25, ls=":")

    # ---- D: coverage of inserts + background rate --------------------------
    ax = axes[1, 1]
    cov = [100.0 * float(r["best_overlap_fraction_insert"])
           for r in hs_rows if truthy(r["matched"])]
    ax.hist(cov, bins=np.arange(45, 105.1, 5), color=BLUE, edgecolor="white", zorder=3)
    ax.set_xlabel("Inserted sequence covered by the matching hit (%)")
    ax.set_ylabel("Recovered sequences")
    ax.set_title("D   Extent of recovery", loc="left", fontweight="bold")
    ax.grid(alpha=0.25, axis="y", ls=":")
    med = float(np.median(cov))
    ax.axvline(med, color=RED, ls="--", lw=1.5)
    ax.annotate(f"median {med:.0f}%", (med, ax.get_ylim()[1] * 0.92),
                fontsize=9, color=RED, ha="right", va="top",
                textcoords="offset points", xytext=(-4, 0))
    ax.annotate(
        f"Unmodified genome: {bg['background_loci_unmodified_genome']} loci in "
        f"{bg['reference_length_bp']/1e6:.2f} Mb\n"
        f"({bg['loci_per_mb']:.1f} per Mb); no inserted negative was called.",
        (0.0, -0.265), xycoords="axes fraction", fontsize=8.2, color="#444444", va="top")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"figure3_ecoli_recovery.{ext}", dpi=300, bbox_inches="tight")
    print("wrote", OUT / "figure3_ecoli_recovery.pdf")

    for t, c, l in zip(tools, conc, lost):
        print(f"  {t:8s} concordance {c:5.1f}%  detections lost {l}")
    print(f"  coverage median {med:.1f}%  min {min(cov):.0f}%")
    print(f"  background {bg['background_loci_unmodified_genome']} loci "
          f"({bg['loci_per_mb']:.1f}/Mb)")


if __name__ == "__main__":
    main()
