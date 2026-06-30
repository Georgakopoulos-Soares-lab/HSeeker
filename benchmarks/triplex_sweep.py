#!/usr/bin/env python3
"""
benchmarks/triplex_sweep.py
============================
Worker sweep for triplex R (v1.50.0) on chr1 hg38.
Calls triplex_benchmark.R at mc.cores = 1,2,4,8,16 via Rscript
and records wall time + peak RSS from /usr/bin/time -v.

Usage:
  python benchmarks/triplex_sweep.py [--workers 1,2,4,8,16] [--json sweep.json]
"""
from __future__ import annotations

import argparse
import csv
import datetime
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

RSCRIPT     = Path("/home1/10899/kimopro/WORK/miniconda3/envs/r_triplex/bin/Rscript")
BENCH_R     = Path(__file__).parent / "triplex_benchmark.R"
CHR1_FA     = Path(__file__).parent / "data" / "chr1.fa"

_RSS_RE  = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
_WALL_RE = re.compile(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)")


def _parse_wall(ts: str) -> float:
    parts = ts.strip().split(":")
    if len(parts) == 2:
        return float(parts[0]) * 60 + float(parts[1])
    return float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])


def run_once(ncores: int, fasta: Path) -> dict:
    fasta_mb = fasta.stat().st_size / 1e6
    with tempfile.TemporaryDirectory() as tmp:
        out_csv = Path(tmp) / "hits.csv"
        cmd = [
            "/usr/bin/time", "-v",
            str(RSCRIPT),
            str(BENCH_R),
            str(fasta),
            str(out_csv),
            str(ncores),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        stderr = proc.stderr + proc.stdout

        rss_m  = _RSS_RE.search(stderr)
        wall_m = _WALL_RE.search(stderr)
        peak_rss_mb = int(rss_m.group(1)) / 1024.0 if rss_m else 0.0
        wall_s      = _parse_wall(wall_m.group(1)) if wall_m else 0.0

        # Also read R-internal wall time from the summary file
        summary_file = Path(str(out_csv) + ".summary.txt")
        r_wall_s = None
        if summary_file.exists():
            for line in summary_file.read_text().splitlines():
                if line.startswith("wall_sec="):
                    r_wall_s = float(line.split("=", 1)[1])

        hits = 0
        if out_csv.exists():
            with open(out_csv) as fh:
                reader = csv.DictReader(fh)
                hits = sum(1 for _ in reader)

    return {
        "workers":          ncores,
        "wall_s":           round(wall_s, 2),       # /usr/bin/time (incl. R startup)
        "wall_s_r":         round(r_wall_s, 2) if r_wall_s else None,  # R internal
        "peak_rss_mb":      round(peak_rss_mb, 2),
        "throughput_mbps":  round(fasta_mb / r_wall_s, 4) if r_wall_s else 0.0,
        "total_hits":       hits,
        "returncode":       proc.returncode,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", default="1,2,4,8,16")
    ap.add_argument("--json", default="benchmarks/data/triplex_sweep.json")
    args = ap.parse_args()

    worker_list = [int(w) for w in args.workers.split(",") if w.strip()]
    fasta_mb = CHR1_FA.stat().st_size / 1e6
    baseline_wall = None

    print("=" * 64)
    print("  triplex R Worker Sweep — chr1")
    print(f"  {datetime.datetime.now():%Y-%m-%d %H:%M:%S}")
    print("=" * 64)
    print(f"  FASTA   : {CHR1_FA.name}  ({fasta_mb:.1f} MB)")
    print(f"  Workers : {worker_list}")
    print()

    results = []
    for w in worker_list:
        print(f"  workers={w:>2} … ", end="", flush=True)
        r = run_once(w, CHR1_FA)
        if r["returncode"] != 0:
            print(f"FAILED (rc={r['returncode']})")
            continue
        wall = r["wall_s_r"] or r["wall_s"]
        if baseline_wall is None:
            baseline_wall = wall
        speedup    = baseline_wall / wall if wall > 0 else 0
        efficiency = speedup / w * 100
        print(f"{wall:6.1f}s  RSS={r['peak_rss_mb']:6.0f} MB  "
              f"hits={r['total_hits']:>8,}  {speedup:.2f}×  {efficiency:.1f}%")
        r["speedup"]    = round(speedup, 4)
        r["efficiency"] = round(efficiency, 2)
        results.append(r)

    print()
    print("=" * 64)

    payload = {
        "date":    datetime.datetime.now().isoformat(),
        "tool":    "triplex",
        "version": "1.50.0",
        "fasta":   str(CHR1_FA),
        "fasta_mb": round(fasta_mb, 2),
        "results": results,
    }
    Path(args.json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.json).write_text(json.dumps(payload, indent=2))
    print(f"\n  JSON → {args.json}")


if __name__ == "__main__":
    main()
