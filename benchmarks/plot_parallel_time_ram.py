from __future__ import annotations

import json
import time
from pathlib import Path

import hseeker
import numpy as np
import matplotlib.pyplot as plt

try:
    import psutil
except Exception as exc:  # pragma: no cover
    raise RuntimeError("psutil is required for RAM measurements") from exc


ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
OUT_JSON = ROOT / "parallel_workers_time_ram.json"
OUT_PNG = ROOT / "parallel_workers_time_ram.png"
OUT_SVG = ROOT / "parallel_workers_time_ram.svg"

WORKERS = [1, 4, 7, 10, 14]
MINREP = 8


def measure_one(path: Path, workers: int) -> dict:
    proc = psutil.Process()
    rss_before = proc.memory_info().rss

    t0 = time.perf_counter()
    c0 = time.process_time()
    hits = hseeker.scan_fasta_parallel(str(path), workers=workers, minrep=MINREP)
    c1 = time.process_time()
    t1 = time.perf_counter()

    rss_after = proc.memory_info().rss
    rss_delta_mb = (rss_after - rss_before) / 1e6

    return {
        "workers": workers,
        "wall_s": t1 - t0,
        "cpu_s": c1 - c0,
        "hits": len(hits),
        "rss_delta_mb": rss_delta_mb,
    }


def main() -> None:
    fasta_files = sorted(DATA_DIR.glob("*.fa"))
    if not fasta_files:
        raise FileNotFoundError(f"No FASTA files found in {DATA_DIR}")

    payload = {
        "minrep": MINREP,
        "workers": WORKERS,
        "files": [],
    }

    for fp in fasta_files:
        print(f"Measuring {fp.name} ...")
        rows = []
        for w in WORKERS:
            print(f"  workers={w}")
            rows.append(measure_one(fp, w))
        payload["files"].append({
            "file": fp.name,
            "size_mb": fp.stat().st_size / 1e6,
            "results": rows,
        })

    OUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")

    plt.style.use("seaborn-v0_8-whitegrid")
    plt.rcParams.update(
        {
            "font.size": 12,
            "axes.titlesize": 14,
            "axes.labelsize": 13,
            "legend.fontsize": 10,
            "xtick.labelsize": 10,
            "ytick.labelsize": 10,
        }
    )

    fig, axes = plt.subplots(1, 2, figsize=(18, 7), constrained_layout=True)

    ax_t, ax_m = axes

    files = [f["file"] for f in payload["files"]]
    x = np.arange(len(files))

    # Color-blind friendly palette with a fixed mapping by worker count.
    worker_palette = {
        1: "#0072B2",
        4: "#E69F00",
        7: "#009E73",
        10: "#D55E00",
        14: "#CC79A7",
    }

    n_w = len(WORKERS)
    total_width = 0.82
    bar_w = total_width / n_w

    for i, w in enumerate(WORKERS):
        offset = (i - (n_w - 1) / 2) * bar_w
        wall_vals = []
        rss_vals = []
        for f in payload["files"]:
            by_w = {r["workers"]: r for r in f["results"]}
            wall_vals.append(by_w[w]["wall_s"])
            rss_vals.append(by_w[w]["rss_delta_mb"])

        ax_t.bar(
            x + offset,
            wall_vals,
            width=bar_w,
            color=worker_palette[w],
            edgecolor="black",
            linewidth=0.5,
            alpha=0.92,
            label=f"{w} workers",
        )
        ax_m.bar(
            x + offset,
            rss_vals,
            width=bar_w,
            color=worker_palette[w],
            edgecolor="black",
            linewidth=0.5,
            alpha=0.92,
            label=f"{w} workers",
        )

    ax_t.set_title("Execution Time by Dataset and Worker Count")
    ax_t.set_xlabel("Dataset")
    ax_t.set_ylabel("Wall time (s)")
    ax_t.set_xticks(x)
    ax_t.set_xticklabels(files, rotation=20, ha="right")
    ax_t.grid(axis="y", linestyle="--", alpha=0.35)

    ax_m.set_title("RAM Delta by Dataset and Worker Count")
    ax_m.set_xlabel("Dataset")
    ax_m.set_ylabel("RSS delta (MB)")
    ax_m.set_xticks(x)
    ax_m.set_xticklabels(files, rotation=20, ha="right")
    ax_m.grid(axis="y", linestyle="--", alpha=0.35)

    handles, labels = ax_t.get_legend_handles_labels()
    fig.legend(
        handles,
        labels,
        loc="lower center",
        ncol=5,
        frameon=False,
        title="Worker count",
    )

    fig.suptitle(
        "HSeeker Parallel Scaling: Execution Time and RAM by Dataset\n"
        "(scan_fasta_parallel only, minrep=8)",
        fontsize=16,
    )

    fig.savefig(OUT_PNG, dpi=450, bbox_inches="tight")
    fig.savefig(OUT_SVG, bbox_inches="tight")
    print(f"Saved {OUT_JSON}")
    print(f"Saved {OUT_PNG}")
    print(f"Saved {OUT_SVG}")


if __name__ == "__main__":
    main()
