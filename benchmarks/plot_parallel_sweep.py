"""
Plot HSeeker parallel worker sweep results.

Produces a 4-panel academic figure:
  (a) Wall-clock time vs workers
  (b) Speedup vs workers (actual vs ideal linear)
  (c) Parallel efficiency vs workers
  (d) Peak RSS vs workers

Usage:
    python benchmarks/plot_parallel_sweep.py [--out figure.pdf]
"""
from __future__ import annotations

import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ── Data from 2026-06-25 sweep (chr1, 252.1 MB, hseeker 0.1.0, scoring ON, purity-rmq) ──
workers    = np.array([1,     2,     4,    8,    16])
wall_s     = np.array([187.7, 105.3, 59.5, 35.9, 25.2])
rss_mb     = np.array([1612,  1614,  1599, 1719, 1965])
speedup    = wall_s[0] / wall_s          # actual speedup relative to 1 worker
efficiency = speedup / workers * 100     # %

# ── Style ────────────────────────────────────────────────────────────────────
BLUE   = "#2166ac"
RED    = "#d6604d"
GREEN  = "#4dac26"
GRAY   = "#888888"

plt.rcParams.update({
    "font.family":      "serif",
    "font.size":        13,
    "axes.titlesize":   13,
    "axes.labelsize":   13,
    "xtick.labelsize":  12,
    "ytick.labelsize":  12,
    "legend.fontsize":  12,
    "axes.linewidth":   1.0,
    "lines.linewidth":  2.2,
    "lines.markersize": 9,
    "grid.linewidth":   0.5,
    "grid.color":       "#cccccc",
    "figure.dpi":       300,
})

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle(
    "HSeeker — Parallel Scaling on chr1",
    fontsize=14, fontweight="bold", y=0.998,
)

def style_ax(ax, xlabel, ylabel, title, letter):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{letter.upper()}) {title}", loc="left", pad=4)
    ax.set_xticks(workers)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.grid(True, axis="y", linestyle="--", alpha=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# ── (a) Wall-clock time ──────────────────────────────────────────────────────
ax = axes[0, 0]
ax.plot(workers, wall_s, "o-", color=BLUE, label="Measured")
for x, y in zip(workers, wall_s):
    ax.annotate(f"{y:.0f}s", (x, y), textcoords="offset points",
                xytext=(0, 8), ha="center", fontsize=11, color=BLUE)
style_ax(ax, "Workers", "Wall-clock time (s)", "Wall-clock time", "a")
ax.set_ylim(0, wall_s.max() * 1.18)

# ── (b) Speedup ──────────────────────────────────────────────────────────────
ax = axes[0, 1]
ax.plot(workers, workers.astype(float), "--", color=GRAY, lw=1.0, label="Ideal (linear)")
ax.plot(workers, speedup, "s-", color=BLUE, label="Measured")
for x, y in zip(workers, speedup):
    ax.annotate(f"{y:.2f}×", (x, y), textcoords="offset points",
                xytext=(0, 8), ha="center", fontsize=11, color=BLUE)
ax.legend(frameon=False)
style_ax(ax, "Workers", "Speedup (×)", "Speedup over 1 worker", "b")
ax.set_ylim(0, workers[-1] * 1.15)

# ── (c) Parallel efficiency ──────────────────────────────────────────────────
ax = axes[1, 0]
ax.axhline(100, color=GRAY, lw=1.0, linestyle="--", label="Ideal (100%)")
ax.plot(workers, efficiency, "^-", color=GREEN, label="Measured")
offsets_eff = [(0, 8), (0, -18), (0, -18), (0, -18), (0, 8)]
for x, y, ofs in zip(workers, efficiency, offsets_eff):
    ax.annotate(f"{y:.1f}%", (x, y), textcoords="offset points",
                xytext=ofs, ha="center", fontsize=11, color=GREEN)
ax.legend(frameon=False)
ax.set_ylim(0, 120)
style_ax(ax, "Workers", "Efficiency (%)", "Parallel efficiency", "c")

# ── (d) Peak RSS ─────────────────────────────────────────────────────────────
ax = axes[1, 1]
ax.plot(workers, rss_mb, "D-", color=RED)
offsets_rss = [(-14, 8), (0, 18), (14, 8), (0, 8), (0, 8)]
for x, y, ofs in zip(workers, rss_mb, offsets_rss):
    ax.annotate(f"{y:,}", (x, y), textcoords="offset points",
                xytext=ofs, ha="center", fontsize=11, color=RED)
style_ax(ax, "Workers", "Peak RSS (MB)", "Peak memory usage", "d")
ax.set_ylim(1400, 2100)
ax.set_yticks([1400, 1600, 1800, 2000])

fig.tight_layout(rect=[0, 0, 1, 0.97], h_pad=4.0, w_pad=3.5)

# ── Save ─────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--out", default="benchmarks/parallel_sweep_plot.pdf")
args, _ = parser.parse_known_args()

fig.savefig(args.out, bbox_inches="tight")
print(f"Saved → {args.out}")
