"""
Plot HSeeker, triplex R, and Triplexator parallel worker sweep results.

4-panel academic figure:
  A) Wall-clock time vs workers  (log-y; all 3 tools)
  B) Speedup vs workers          (HSeeker + triplex R; Triplexator is flat)
  C) Parallel efficiency         (HSeeker + triplex R)
  D) Peak RSS vs workers         (all 3 tools)

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

# ── HSeeker sweep (chr1 hg38 248.9 Mbp, v0.1.0, purity-RMQ, 2026-06-25) ─────
hs_workers    = np.array([1,     2,     4,    8,    16])
hs_wall_s     = np.array([187.7, 105.3, 59.5, 35.9, 25.2])
hs_rss_mb     = np.array([1612,  1614,  1599, 1719, 1965])
hs_speedup    = hs_wall_s[0] / hs_wall_s
hs_efficiency = hs_speedup / hs_workers * 100

# ── triplex R v1.50.0 sweep (mclapply, 249×1 Mbp chunks, chr1) ───────────────
# 1-8 cores: live sweep; 16-core: standalone run at lower system commit load
tx_workers    = np.array([1,      2,      4,      8,      16])
tx_wall_s     = np.array([816.31, 407.70, 209.85, 112.38,  57.70])   # R internal timer
tx_rss_mb     = np.array([643.67, 605.76, 599.21, 593.48, 578.00])
tx_speedup    = tx_wall_s[0] / tx_wall_s
tx_efficiency = tx_speedup / tx_workers * 100

# ── Triplexator v1.3.2 sweep (-rm 1 -p N, chr1) ──────────────────────────────
# Native OpenMP with -rm 1 shows flat scaling on a single long chromosome
tp_workers    = np.array([1,    2,    4,    8,    16])
tp_wall_s     = np.array([68.3, 67.8, 66.5, 67.8, 67.5])
tp_rss_mb     = np.array([1709, 1709, 1709, 1710, 1709])
tp_speedup    = tp_wall_s[0] / tp_wall_s
tp_efficiency = tp_speedup / tp_workers * 100

# ── Colours ───────────────────────────────────────────────────────────────────
BLUE    = "#2166ac"   # HSeeker line
RED     = "#d6604d"   # HSeeker RSS line
ORANGE  = "#e08214"   # triplex R
PURPLE  = "#762a83"   # Triplexator
GRAY    = "#888888"   # ideal / reference

plt.rcParams.update({
    "font.family":      "serif",
    "font.size":        13,
    "axes.titlesize":   13,
    "axes.labelsize":   13,
    "xtick.labelsize":  12,
    "ytick.labelsize":  12,
    "legend.fontsize":  10,
    "axes.linewidth":   1.0,
    "lines.linewidth":  2.2,
    "lines.markersize": 8,
    "grid.linewidth":   0.5,
    "grid.color":       "#cccccc",
    "figure.dpi":       300,
})

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
fig.suptitle(
    "Parallel Scaling on chr1 hg38 — HSeeker vs. triplex R vs. Triplexator",
    fontsize=14, fontweight="bold", y=0.998,
)

WORKERS = np.array([1, 2, 4, 8, 16])

def style_ax(ax, xlabel, ylabel, title, letter):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(f"{letter.upper()}) {title}", loc="left", pad=4)
    ax.set_xticks(WORKERS)
    ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
    ax.grid(True, axis="y", linestyle="--", alpha=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)

# ── A) Wall-clock time (log y) ────────────────────────────────────────────────
ax = axes[0, 0]

ax.plot(hs_workers, hs_wall_s, "o-",  color=BLUE,   zorder=4, label="HSeeker (this work)")
ax.plot(tx_workers, tx_wall_s, "s--", color=ORANGE, zorder=3, label="triplex R v1.50.0")
ax.plot(tp_workers, tp_wall_s, "^:",  color=PURPLE, zorder=3, label="Triplexator v1.3.2")

# Annotate HSeeker only (to keep figure readable)
# x=1 and x=2 sit near the left edge on a linear axis — shift right to avoid clip
hs_wall_offsets = [(0, 8), (8, 8), (0, 8), (0, 8), (0, 8)]
hs_wall_ha      = ["center", "left", "center", "center", "center"]
for x, y, ofs, ha in zip(hs_workers, hs_wall_s, hs_wall_offsets, hs_wall_ha):
    ax.annotate(f"{y:.0f}s", (x, y), textcoords="offset points",
                xytext=ofs, ha=ha, fontsize=9, color=BLUE)

ax.set_yscale("log")
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f"{v:.0f}"))
ax.set_yticks([20, 50, 100, 200, 500, 1000])
ax.legend(frameon=False, fontsize=10, loc="upper right")
style_ax(ax, "Workers", "Wall-clock time (s)", "Wall-clock time (log scale)", "a")

# ── B) Speedup ────────────────────────────────────────────────────────────────
ax = axes[0, 1]
ax.plot(WORKERS, WORKERS.astype(float), "--", color=GRAY, lw=1.0, label="Ideal (linear)")
ax.plot(hs_workers, hs_speedup,    "o-",  color=BLUE,   label="HSeeker")
ax.plot(tx_workers, tx_speedup,    "s--", color=ORANGE, label="triplex R")
ax.plot(tp_workers, tp_speedup,    "^:",  color=PURPLE, label="Triplexator")

for x, y in zip(hs_workers, hs_speedup):
    ax.annotate(f"{y:.1f}×", (x, y), textcoords="offset points",
                xytext=(0, 8), ha="center", fontsize=9, color=BLUE)

ax.legend(frameon=False, fontsize=10)
style_ax(ax, "Workers", "Speedup (×)", "Speedup over 1 worker", "b")
ax.set_ylim(0, WORKERS[-1] * 1.15)

# ── C) Parallel efficiency ────────────────────────────────────────────────────
ax = axes[1, 0]
ax.axhline(100, color=GRAY, lw=1.0, linestyle="--", label="Ideal (100%)")
ax.plot(hs_workers, hs_efficiency, "o-",  color=BLUE,   label="HSeeker")
ax.plot(tx_workers, tx_efficiency, "s--", color=ORANGE, label="triplex R")
ax.plot(tp_workers, tp_efficiency, "^:",  color=PURPLE, label="Triplexator")

offsets_eff = [(10, 8), (0, -18), (0, -18), (0, -18), (0, 8)]
ha_eff      = ["left", "center", "center", "center", "center"]
for x, y, ofs, ha in zip(hs_workers, hs_efficiency, offsets_eff, ha_eff):
    ax.annotate(f"{y:.1f}%", (x, y), textcoords="offset points",
                xytext=ofs, ha=ha, fontsize=9, color=BLUE)

ax.legend(frameon=False, fontsize=10)
ax.set_ylim(0, 120)
style_ax(ax, "Workers", "Efficiency (%)", "Parallel efficiency", "c")

# ── D) Peak RSS ───────────────────────────────────────────────────────────────
ax = axes[1, 1]

ax.plot(hs_workers, hs_rss_mb, "D-",  color=RED,    zorder=4, label="HSeeker")
ax.plot(tx_workers, tx_rss_mb, "s--", color=ORANGE, zorder=3, label="triplex R")
ax.plot(tp_workers, tp_rss_mb, "^:",  color=PURPLE, zorder=3, label="Triplexator")

offsets_rss = [(-14, 8), (0, 18), (16, 8), (0, 8), (0, 8)]
for x, y, ofs in zip(hs_workers, hs_rss_mb, offsets_rss):
    ax.annotate(f"{y:,}", (x, y), textcoords="offset points",
                xytext=ofs, ha="center", fontsize=9, color=RED)

ax.legend(frameon=False, fontsize=10, loc="center right")
style_ax(ax, "Workers", "Peak RSS (MB)", "Peak memory usage", "d")
ax.set_ylim(0, 2300)
ax.set_yticks([0, 500, 1000, 1500, 2000])

fig.tight_layout(rect=[0, 0, 1, 0.97], h_pad=4.0, w_pad=3.5)

# ── Save ──────────────────────────────────────────────────────────────────────
parser = argparse.ArgumentParser()
parser.add_argument("--out", default="benchmarks/parallel_sweep_plot.pdf")
args, _ = parser.parse_known_args()

fig.savefig(args.out, bbox_inches="tight")
# Also save PNG for quick preview
png_out = args.out.replace(".pdf", ".png")
fig.savefig(png_out, bbox_inches="tight", dpi=200)
print(f"Saved → {args.out}")
print(f"Saved → {png_out}")
