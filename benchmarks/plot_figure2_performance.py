#!/usr/bin/env python3
"""Figure 2: computational performance of HSeeker (October 2026 measurements).

A) Wall-clock time and speedup vs worker count on hg38 chr1.
B) Peak resident set size vs worker count, with the whole-genome run marked.
C) Parallel efficiency vs worker count.

The comparator tools are deliberately not plotted: nBMST and NeSSie have no
parallel mode, so a worker sweep is undefined for them. Their timings are
printed by this script for use in the caption and Results text.

Inputs: benchmarks/results_2026-10-11/ (sweep JSON, /usr/bin/time -v output).
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
RES = ROOT / "benchmarks" / "results_2026-10-11"
OUT = ROOT / "analysis" / "results" / "plots"

BLUE = "#2166ac"
RED = "#b2182b"
GREY = "#999999"

# Peak RSS is reported as kbytes/1024, matching the published benchmark table.
TIME_RSS = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
TIME_WALL = re.compile(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)")


def parse_time_file(path: Path) -> tuple[float, float]:
    """Return (wall seconds, peak RSS in MB) from a /usr/bin/time -v file."""
    text = path.read_text()
    rss = int(TIME_RSS.search(text).group(1)) / 1024.0
    parts = [float(p) for p in TIME_WALL.search(text).group(1).split(":")]
    wall = parts[-1] + (parts[-2] * 60 if len(parts) > 1 else 0) + (parts[-3] * 3600 if len(parts) > 2 else 0)
    return wall, rss


def mean_of(pattern: str) -> tuple[float, float]:
    vals = [parse_time_file(p) for p in sorted(RES.glob(pattern))]
    return float(np.mean([v[0] for v in vals])), float(np.mean([v[1] for v in vals]))


def main() -> None:
    sweep = json.loads((RES / "sweep_greedy_minrep8.json").read_text())["results"]
    workers = np.array([r["workers"] for r in sweep])
    wall = np.array([r["wall_s"] for r in sweep])
    rss = np.array([r["peak_rss_mb"] for r in sweep])
    speedup = wall[0] / wall

    genome_wall, genome_rss = mean_of("genome_greedy_*.time")

    # Comparator timings are reported in the caption and Results, not plotted:
    # two of the four tools have no parallel mode, so a worker sweep is undefined
    # for them and a mixed plot would imply a comparison the data cannot support.
    comparators = [
        ("HSeeker", wall[-1], rss[-1], "45,927"),
        ("Triplexator", *mean_of("comparators/triplexator_*.time"), "275,476"),
        ("triplex", *mean_of("comparators/triplex_*.time"), "81,906"),
        ("nBMST", *mean_of("comparators/gfa_*.time"), "147,222"),
        ("NeSSie", *mean_of("comparators/nessie.time"), "720,811"),
    ]

    fig = plt.figure(figsize=(13.5, 4.3))
    gs = fig.add_gridspec(1, 3, wspace=0.42, left=0.055, right=0.985, top=0.86, bottom=0.17)

    # ---- A: wall time + speedup -------------------------------------------
    ax = fig.add_subplot(gs[0, 0])
    ax.plot(workers, wall, "o-", color=BLUE, lw=2, ms=7, zorder=3)
    for x, y in zip(workers, wall):
        ax.annotate(f"{y:.1f} s", (x, y), textcoords="offset points", xytext=(6, 7),
                    fontsize=8.5, color=BLUE)
    ax.set_xscale("log", base=2)
    ax.set_yscale("log")
    ax.set_xticks(workers)
    ax.set_xticklabels(workers)
    ax.set_xlabel("Worker threads")
    ax.set_ylabel("Wall-clock time (s)", color=BLUE)
    ax.tick_params(axis="y", labelcolor=BLUE)
    ax.set_title("A   Scaling on chromosome 1", loc="left", fontweight="bold")
    ax.grid(alpha=0.25, which="both", ls=":")

    ax2 = ax.twinx()
    ax2.plot(workers, workers, ls="--", color=GREY, lw=1.3, label="Ideal")
    ax2.plot(workers, speedup, "s-", color=RED, lw=2, ms=6, label="Observed")
    ax2.annotate(f"{speedup[-1]:.2f}×", (workers[-1], speedup[-1]),
                 textcoords="offset points", xytext=(-32, 4), fontsize=9,
                 color=RED, fontweight="bold")
    ax2.set_ylabel("Speedup over 1 worker", color=RED)
    ax2.tick_params(axis="y", labelcolor=RED)
    ax2.set_ylim(0, 17.5)
    ax2.legend(loc="upper left", fontsize=8.5, frameon=False)

    # ---- B: peak memory ----------------------------------------------------
    ax = fig.add_subplot(gs[0, 1])
    ax.plot(workers, rss, "o-", color=BLUE, lw=2, ms=7, label="Chromosome 1")
    ax.axhline(genome_rss, ls="--", color=RED, lw=1.6)
    ax.annotate(f"Whole genome, 16 workers: {genome_rss:,.0f} MB\n"
                f"({genome_wall:.1f} s, 559,601 loci)",
                (workers[1], genome_rss), textcoords="offset points", xytext=(0, -30),
                fontsize=8.5, color=RED)
    ax.set_xscale("log", base=2)
    ax.set_xticks(workers)
    ax.set_xticklabels(workers)
    ax.set_ylim(0, max(genome_rss, rss.max()) * 1.35)
    ax.set_xlabel("Worker threads")
    ax.set_ylabel("Peak resident set size (MB)")
    ax.set_title("B   Memory", loc="left", fontweight="bold")
    ax.grid(alpha=0.25, ls=":")
    ax.legend(loc="lower right", fontsize=8.5, frameon=False)

    # ---- C: parallel efficiency ------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    eff = 100.0 * speedup / workers
    ax.plot(workers, eff, "o-", color=BLUE, lw=2, ms=7, zorder=3)
    ax.axhline(100, ls="--", color=GREY, lw=1.3)
    ax.annotate("Ideal", (workers[-1], 100), textcoords="offset points",
                xytext=(-26, 5), fontsize=8.5, color=GREY)
    for x, y in zip(workers, eff):
        ax.annotate(f"{y:.0f}%", (x, y), textcoords="offset points", xytext=(4, -12),
                    fontsize=8.5, color=BLUE)
    ax.set_xscale("log", base=2)
    ax.set_xticks(workers)
    ax.set_xticklabels(workers)
    ax.set_ylim(0, 118)
    ax.set_xlabel("Worker threads")
    ax.set_ylabel("Parallel efficiency (%)")
    ax.set_title("C   Parallel efficiency", loc="left", fontweight="bold")
    ax.grid(alpha=0.25, ls=":")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"figure2_performance.{ext}", dpi=300, bbox_inches="tight")
    print("wrote", OUT / "figure2_performance.pdf")

    print(f"\nchr1 sweep:   {dict(zip(workers.tolist(), wall.tolist()))}")
    print(f"speedup:      {speedup[-1]:.2f}x   RSS {rss.min():.0f}-{rss.max():.0f} MB")
    print(f"whole genome: {genome_wall:.1f} s, {genome_rss:.0f} MB")
    print("\ncomparators (caption text, not plotted):")
    for n, t, m, l in comparators:
        print(f"  {n:12s} {t:6.1f} s  {m:6.0f} MB  {l}")


if __name__ == "__main__":
    main()
