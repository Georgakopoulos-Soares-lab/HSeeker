#!/usr/bin/env python3
"""Figure 2: computational performance of HSeeker (October 2026 measurements).

A) Wall-clock time and speedup vs worker count on hg38 chr1.
B) Peak resident set size vs worker count, with the whole-genome run marked.
C) Wall-clock time and peak memory at 16 cores for HSeeker and four comparators.

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

    comparators = [
        ("HSeeker", wall[-1], rss[-1], "45,927", BLUE),
        ("Triplexator", *mean_of("comparators/triplexator_*.time"), "275,476", GREY),
        ("triplex", *mean_of("comparators/triplex_*.time"), "81,906", GREY),
        ("nBMST", *mean_of("comparators/gfa_*.time"), "147,222", GREY),
        ("NeSSie", *mean_of("comparators/nessie.time"), "720,811", GREY),
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

    # ---- C: comparator bars ------------------------------------------------
    ax = fig.add_subplot(gs[0, 2])
    names = [c[0] for c in comparators]
    times = [c[1] for c in comparators]
    mems = [c[2] for c in comparators]
    loci = [c[3] for c in comparators]
    cols = [c[4] for c in comparators]
    y = np.arange(len(names))[::-1]

    ax.barh(y, times, color=cols, height=0.62, zorder=3)
    for yi, t, n, m in zip(y, times, loci, mems):
        ax.annotate(f"{t:.1f} s  ·  {m:,.0f} MB\n({n})", (t, yi),
                    textcoords="offset points", xytext=(6, 0), va="center",
                    fontsize=8.2)
    ax.set_yticks(y)
    ax.set_yticklabels(names)
    ax.set_xlim(0, max(times) * 1.72)
    ax.set_xlabel("Wall-clock time on chromosome 1 (s)")
    ax.set_title("C   Against existing tools", loc="left", fontweight="bold")
    ax.grid(alpha=0.25, axis="x", ls=":")

    note = ("16 cores, except nBMST and NeSSie, which have no parallel mode.\n"
            "Loci counts in brackets are not comparable: the tools enumerate\n"
            "different target classes.")
    ax.annotate(note, (0.0, -0.265), xycoords="axes fraction", fontsize=7.2,
                color="#444444", va="top")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"figure2_performance.{ext}", dpi=300, bbox_inches="tight")
    print("wrote", OUT / "figure2_performance.pdf")

    print(f"\nchr1 sweep:   {dict(zip(workers.tolist(), wall.tolist()))}")
    print(f"speedup:      {speedup[-1]:.2f}x   RSS {rss.min():.0f}-{rss.max():.0f} MB")
    print(f"whole genome: {genome_wall:.1f} s, {genome_rss:.0f} MB")
    for n, t, m, l, _ in comparators:
        print(f"  {n:12s} {t:6.1f} s  {m:6.0f} MB  {l}")


if __name__ == "__main__":
    main()
