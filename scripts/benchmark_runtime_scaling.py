#!/usr/bin/env python3
"""Measure HSeeker detector runtime scaling across DNA sequence lengths."""

from __future__ import annotations

import csv
import math
import random
import statistics
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
OUT = ROOT / "results" / "runtime_scaling"
PLOTS = OUT / "plots"

PARAMS = dict(minrep=10, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10)
METHODS = [
    "original",
    "safe_pruned_mismatch_only",
    "fast_gated_fast_4x",
    "fast_gated_fast_8x",
]
LENGTHS = [1_000, 2_000, 5_000, 10_000, 20_000, 50_000, 100_000, 200_000, 500_000, 1_000_000]
REPEAT_LENGTH_LIMIT = 100_000

sys.path.insert(0, str(SRC))
import hseeker  # noqa: E402


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)


def write_csv(path: Path, rows: list[dict]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def weighted_dna(rng: random.Random, n: int, weights: dict[str, int]) -> str:
    bases = list(weights)
    return "".join(rng.choices(bases, weights=[weights[b] for b in bases], k=n))


def make_sequence(sequence_class: str, n: int, rng: random.Random) -> str:
    if sequence_class == "balanced_random":
        return weighted_dna(rng, n, dict(A=1, C=1, G=1, T=1))
    if sequence_class == "ga_rich_random":
        return weighted_dna(rng, n, dict(A=45, G=45, C=5, T=5))
    if sequence_class == "pure_gaa_repeat":
        return ("GAA" * (n // 3 + 1))[:n]
    raise ValueError(sequence_class)


def method_kwargs(method: str) -> dict:
    if method == "original":
        return dict(detector="original")
    if method == "safe_pruned_mismatch_only":
        return dict(detector="safe_pruned", safe_prune_purity=False)
    if method == "fast_gated_fast_4x":
        return dict(
            detector="fast_gated",
            gate_window=PARAMS["minrep"],
            gate_search_limit=4 * PARAMS["minrep"],
            fast_mode=True,
            use_purity_prefilter=True,
        )
    if method == "fast_gated_fast_8x":
        return dict(
            detector="fast_gated",
            gate_window=PARAMS["minrep"],
            gate_search_limit=8 * PARAMS["minrep"],
            fast_mode=True,
            use_purity_prefilter=True,
        )
    raise ValueError(method)


def reps_for_length(n: int) -> int:
    if n <= 20_000:
        return 3
    if n <= 100_000:
        return 2
    return 1


def run_detector(seq: str, method: str) -> tuple[float, dict, int]:
    start = time.perf_counter()
    hits = hseeker.search(seq, **PARAMS, **method_kwargs(method), score=False, remove_overlaps=True)
    runtime = time.perf_counter() - start
    return runtime, hseeker.profiling_info(), len(hits)


def benchmark() -> list[dict]:
    rng = random.Random(20260617)
    rows = []
    classes = ["balanced_random", "ga_rich_random", "pure_gaa_repeat"]
    for sequence_class in classes:
        for n in LENGTHS:
            if sequence_class == "pure_gaa_repeat" and n > REPEAT_LENGTH_LIMIT:
                continue
            seq = make_sequence(sequence_class, n, rng)
            for method in METHODS:
                runtimes = []
                last_prof = {}
                final_hits = 0
                for rep in range(reps_for_length(n)):
                    runtime, prof, final_hits = run_detector(seq, method)
                    runtimes.append(runtime)
                    last_prof = prof
                rows.append({
                    "sequence_class": sequence_class,
                    "sequence_length": n,
                    "method": method,
                    "replicates": len(runtimes),
                    "runtime_sec_mean": statistics.mean(runtimes),
                    "runtime_sec_median": statistics.median(runtimes),
                    "runtime_sec_min": min(runtimes),
                    "runtime_sec_max": max(runtimes),
                    "ctr_sp_pairs": last_prof.get("ctr_sp_pairs", 0),
                    "inner_iterations": last_prof.get("inner_iters", 0),
                    "raw_hits": last_prof.get("raw_hits", 0),
                    "final_hits": last_prof.get("final_hits", final_hits),
                    "gate_passed": last_prof.get("gate_passed", 0),
                    "stopped_by_mismatch_impossible": last_prof.get("stopped_by_mismatch_impossible", 0),
                    "stopped_by_purity_impossible": last_prof.get("stopped_by_purity_impossible", 0),
                })
                print(f"{sequence_class:16s} {n:8d} {method:28s} {statistics.mean(runtimes):.6f}s")
    return rows


def summarize(rows: list[dict]) -> list[dict]:
    summary = []
    for sequence_class in sorted({r["sequence_class"] for r in rows}):
        original_by_len = {
            int(r["sequence_length"]): float(r["runtime_sec_mean"])
            for r in rows
            if r["sequence_class"] == sequence_class and r["method"] == "original"
        }
        for method in METHODS:
            subset = [r for r in rows if r["sequence_class"] == sequence_class and r["method"] == method]
            if len(subset) < 2:
                continue
            xs = [math.log(float(r["sequence_length"])) for r in subset if float(r["runtime_sec_mean"]) > 0]
            ys = [math.log(float(r["runtime_sec_mean"])) for r in subset if float(r["runtime_sec_mean"]) > 0]
            xbar = statistics.mean(xs)
            ybar = statistics.mean(ys)
            denom = sum((x - xbar) ** 2 for x in xs)
            slope = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys)) / denom if denom else 0.0
            largest = max(subset, key=lambda r: int(r["sequence_length"]))
            n = int(largest["sequence_length"])
            speedup = original_by_len.get(n, 0.0) / float(largest["runtime_sec_mean"]) if float(largest["runtime_sec_mean"]) else 0.0
            summary.append({
                "sequence_class": sequence_class,
                "method": method,
                "log_log_runtime_slope": slope,
                "largest_length_tested": n,
                "runtime_at_largest_length_sec": largest["runtime_sec_mean"],
                "speedup_vs_original_at_largest_length": speedup,
            })
    return summary


def make_plots(rows: list[dict], summary: list[dict]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except Exception:
        plt = None

    if plt is not None:
        for sequence_class in sorted({r["sequence_class"] for r in rows}):
            subset = [r for r in rows if r["sequence_class"] == sequence_class]
            plt.figure(figsize=(9, 5))
            for method in METHODS:
                series = sorted([r for r in subset if r["method"] == method], key=lambda r: int(r["sequence_length"]))
                if not series:
                    continue
                plt.plot(
                    [int(r["sequence_length"]) for r in series],
                    [float(r["runtime_sec_mean"]) for r in series],
                    marker="o",
                    label=method,
                )
            plt.xscale("log")
            plt.yscale("log")
            plt.xlabel("sequence length (bp)")
            plt.ylabel("mean runtime (sec)")
            plt.title(f"Runtime scaling: {sequence_class}")
            plt.legend()
            plt.tight_layout()
            plt.savefig(PLOTS / f"runtime_scaling_{sequence_class}.png", dpi=160)
            plt.close()
        return

    try:
        from PIL import Image, ImageDraw
    except Exception:
        return

    for sequence_class in sorted({r["sequence_class"] for r in rows}):
        img = Image.new("RGB", (1100, 620), "white")
        draw = ImageDraw.Draw(img)
        draw.text((40, 20), f"Runtime scaling: {sequence_class}", fill=(20, 20, 20))
        subset = [r for r in rows if r["sequence_class"] == sequence_class]
        max_rt = max(float(r["runtime_sec_mean"]) for r in subset) or 1.0
        min_rt = min(float(r["runtime_sec_mean"]) for r in subset if float(r["runtime_sec_mean"]) > 0)
        max_len = max(int(r["sequence_length"]) for r in subset)
        min_len = min(int(r["sequence_length"]) for r in subset)
        colors = {
            "original": (20, 80, 160),
            "safe_pruned_mismatch_only": (60, 140, 80),
            "fast_gated_fast_4x": (180, 90, 40),
            "fast_gated_fast_8x": (120, 70, 160),
        }
        left, top, right, bottom = 80, 70, 1040, 520
        draw.rectangle((left, top, right, bottom), outline=(180, 180, 180))
        for method in METHODS:
            series = sorted([r for r in subset if r["method"] == method], key=lambda r: int(r["sequence_length"]))
            pts = []
            for r in series:
                lx = (math.log(int(r["sequence_length"])) - math.log(min_len)) / (math.log(max_len) - math.log(min_len))
                ly = (math.log(float(r["runtime_sec_mean"])) - math.log(min_rt)) / (math.log(max_rt) - math.log(min_rt))
                pts.append((left + lx * (right - left), bottom - ly * (bottom - top)))
            for a, b in zip(pts, pts[1:]):
                draw.line((a[0], a[1], b[0], b[1]), fill=colors[method], width=3)
            for x, y in pts:
                draw.ellipse((x - 4, y - 4, x + 4, y + 4), fill=colors[method])
        y = 540
        for method in METHODS:
            draw.rectangle((80, y, 95, y + 12), fill=colors[method])
            draw.text((105, y - 2), method, fill=(20, 20, 20))
            y += 18
        img.save(PLOTS / f"runtime_scaling_{sequence_class}.png")


def write_summary(rows: list[dict], summary: list[dict]) -> None:
    lines = ["# Runtime Scaling Benchmark", ""]
    lines.append("Detector parameters: `minrep=10`, `maxrep=1000`, `maxspacer=10`, `purity=0.90`, `mismatch=0.10`; scoring disabled.")
    lines.append("")
    lines.append("## Largest-Length Speedups")
    lines.append("")
    for r in summary:
        lines.append(
            f"- {r['sequence_class']} / {r['method']}: length {r['largest_length_tested']} bp, "
            f"runtime {float(r['runtime_at_largest_length_sec']):.4f}s, "
            f"speedup vs original {float(r['speedup_vs_original_at_largest_length']):.2f}x, "
            f"log-log slope {float(r['log_log_runtime_slope']):.2f}"
        )
    lines.append("")
    lines.append("## Caveats")
    lines.append("")
    lines.append("- Pure GAA repeat was capped at 100 kb because repeat-heavy inputs are intentionally adversarial.")
    lines.append("- Log-log slopes are empirical over the tested length range, not formal asymptotic proofs.")
    lines.append("- Runtime can vary with DNA composition because pruning/gating behavior changes with match and GA/CT density.")
    (OUT / "summary.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    ensure_dirs()
    rows = benchmark()
    summary = summarize(rows)
    write_csv(OUT / "runtime_scaling.csv", rows)
    write_csv(OUT / "runtime_scaling_summary.csv", summary)
    make_plots(rows, summary)
    write_summary(rows, summary)

    print("\nRuntime scaling summary:")
    for sequence_class in sorted({r["sequence_class"] for r in summary}):
        print(f"  {sequence_class}:")
        for r in [x for x in summary if x["sequence_class"] == sequence_class]:
            print(
                f"    {r['method']}: slope={float(r['log_log_runtime_slope']):.2f}, "
                f"largest={r['largest_length_tested']}bp, "
                f"speedup={float(r['speedup_vs_original_at_largest_length']):.2f}x"
            )


if __name__ == "__main__":
    main()
