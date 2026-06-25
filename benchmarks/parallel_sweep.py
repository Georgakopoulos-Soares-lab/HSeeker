#!/usr/bin/env python3
"""
benchmarks/parallel_sweep.py — Worker-thread sweep benchmark on chr1
======================================================================

Runs the hseeker CLI against hg38 chr1 with 1, 2, 4, 8, and 16 worker
threads, recording wall time and peak RSS at each point.

Metrics reported
----------------
  Wall time   — total end-to-end duration (s)
  Peak RSS    — maximum resident set size during the run (MB)
  Throughput  — MB of FASTA processed per second of wall time
  Speedup     — wall time relative to the 1-worker baseline

Usage
-----
  python benchmarks/parallel_sweep.py
  python benchmarks/parallel_sweep.py --workers 1,2,4,8,16,32
  python benchmarks/parallel_sweep.py --runs 3          # repeat each point
  python benchmarks/parallel_sweep.py --json sweep.json
"""

from __future__ import annotations

import argparse
import datetime
import json
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

import hseeker

# ---------------------------------------------------------------------------
# Paths & constants
# ---------------------------------------------------------------------------

DATA_DIR = Path(__file__).parent / "data"
CHR1_FA  = DATA_DIR / "chr1.fa"

_TIME_RSS_RE  = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
_TIME_WALL_RE = re.compile(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)")

_W = 80  # print width


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _parse_time_str(ts: str) -> float:
    ts = ts.strip()
    parts = ts.split(":")
    if len(parts) == 2:
        return float(parts[0]) * 60 + float(parts[1])
    elif len(parts) == 3:
        return float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])
    return 0.0


def _run_once(fa: Path, workers: int, minrep: int) -> dict[str, Any]:
    """Run hseeker CLI once and return {wall_s, peak_rss_mb, throughput_mbps}."""
    fasta_mb = fa.stat().st_size / 1e6
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = str(Path(tmp) / "bench")
        cmd = [
            "/usr/bin/time", "-v",
            sys.executable, "-m", "hseeker",
            "-seq", str(fa),
            "-out", out_prefix,
            "-minrep", str(minrep),
            "-workers", str(workers),
            "-purity-rmq",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        stderr = proc.stderr

    rss_m   = _TIME_RSS_RE.search(stderr)
    wall_m  = _TIME_WALL_RE.search(stderr)
    peak_rss_mb = int(rss_m.group(1)) / 1024.0 if rss_m else 0.0
    wall_s      = _parse_time_str(wall_m.group(1)) if wall_m else 0.0

    return {
        "workers":         workers,
        "wall_s":          wall_s,
        "peak_rss_mb":     peak_rss_mb,
        "throughput_mbps": fasta_mb / wall_s if wall_s > 0 else 0.0,
        "returncode":      proc.returncode,
    }


def _separator(char: str = "─") -> None:
    print(char * _W)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--workers", default="1,2,4,8,16",
        help="Comma-separated worker counts to sweep (default: 1,2,4,8,16)",
    )
    p.add_argument(
        "--runs", type=int, default=1,
        help="Repetitions per worker count; results are averaged (default: 1)",
    )
    p.add_argument(
        "--minrep", type=int, default=10,
        help="Minimum arm length passed to hseeker (default: 10)",
    )
    p.add_argument(
        "--json", metavar="FILE",
        help="Save machine-readable results to a JSON file",
    )
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    worker_list = [int(w.strip()) for w in args.workers.split(",") if w.strip()]
    if not worker_list:
        print("ERROR: --workers must be a non-empty comma-separated list", file=sys.stderr)
        sys.exit(1)

    if not CHR1_FA.exists():
        print(
            f"ERROR: {CHR1_FA} not found.\n"
            "Run: python benchmarks/benchmark.py --chromosomes chr1",
            file=sys.stderr,
        )
        sys.exit(1)

    fasta_mb = CHR1_FA.stat().st_size / 1e6

    _separator("=")
    print("  HSeeker Parallel Worker Sweep — chr1".center(_W))
    print(datetime.datetime.now().strftime("  %Y-%m-%d %H:%M:%S").center(_W))
    _separator("=")
    print(f"  hseeker  : {hseeker.__version__}")
    print(f"  FASTA    : {CHR1_FA.name}  ({fasta_mb:.1f} MB)")
    print(f"  minrep   : {args.minrep}  |  Scoring: on  |  Runs per point: {args.runs}")
    print(f"  Workers  : {worker_list}")
    print()

    all_results: list[dict[str, Any]] = []

    for w in worker_list:
        run_times: list[float] = []
        run_rss:   list[float] = []

        for r in range(1, args.runs + 1):
            label = f"  workers={w:>2}"
            if args.runs > 1:
                label += f"  run {r}/{args.runs}"
            print(f"{label} … ", end="", flush=True)

            result = _run_once(CHR1_FA, workers=w, minrep=args.minrep)

            if result["returncode"] != 0:
                print("FAILED")
                continue

            print(
                f"{result['wall_s']:6.1f}s  "
                f"RSS={result['peak_rss_mb']:5.0f} MB  "
                f"{result['throughput_mbps']:.2f} MB/s"
            )
            run_times.append(result["wall_s"])
            run_rss.append(result["peak_rss_mb"])

        if not run_times:
            continue

        avg_wall = sum(run_times) / len(run_times)
        avg_rss  = sum(run_rss)  / len(run_rss)
        all_results.append({
            "workers":         w,
            "wall_s":          avg_wall,
            "peak_rss_mb":     avg_rss,
            "throughput_mbps": fasta_mb / avg_wall if avg_wall > 0 else 0.0,
            "fasta_mb":        fasta_mb,
            "runs":            len(run_times),
        })

    if not all_results:
        print("\n  No successful runs.", file=sys.stderr)
        sys.exit(1)

    # Baseline = smallest worker count that succeeded
    baseline_wall = all_results[0]["wall_s"]

    # ── Summary table ─────────────────────────────────────────────────────
    print()
    _separator()
    print(
        f"  {'Workers':>7}  {'Wall (s)':>9}  {'RSS (MB)':>9}  "
        f"{'MB/s':>8}  {'Speedup':>8}  {'Efficiency':>10}"
    )
    _separator("·")
    for r in all_results:
        speedup    = baseline_wall / r["wall_s"] if r["wall_s"] > 0 else 0.0
        efficiency = speedup / r["workers"] * 100
        print(
            f"  {r['workers']:>7}  "
            f"{r['wall_s']:>9.1f}  "
            f"{r['peak_rss_mb']:>9.0f}  "
            f"{r['throughput_mbps']:>8.2f}  "
            f"{speedup:>7.2f}×  "
            f"{efficiency:>9.1f}%"
        )
    _separator("=")
    print("  Done.")
    _separator("=")

    # ── JSON output ────────────────────────────────────────────────────────
    if args.json:
        payload = {
            "date":    datetime.datetime.now().isoformat(),
            "hseeker": hseeker.__version__,
            "fasta":   str(CHR1_FA),
            "fasta_mb": fasta_mb,
            "minrep":  args.minrep,
            "results": all_results,
        }
        Path(args.json).write_text(json.dumps(payload, indent=2))
        print(f"\n  JSON → {args.json}")


if __name__ == "__main__":
    main()
