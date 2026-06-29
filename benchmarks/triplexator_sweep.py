#!/usr/bin/env python3
"""
benchmarks/triplexator_sweep.py
================================
Worker sweep for Triplexator (v1.3.2) on chr1 hg38.
Runs -p 1,2,4,8,16 with -rm 1 and records wall time + peak RSS.

Usage:
  python benchmarks/triplexator_sweep.py [--workers 1,2,4,8,16] [--json sweep.json]
"""
from __future__ import annotations

import argparse
import datetime
import json
import re
import subprocess
import sys
import tempfile
from pathlib import Path

TRIPLEXATOR = Path("/scratch/10899/kimopro/triplexator/bin/triplexator")
CHR1_FA     = Path(__file__).parent / "data" / "chr1.fa"

_RSS_RE  = re.compile(r"Maximum resident set size \(kbytes\):\s*(\d+)")
_WALL_RE = re.compile(r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)")


def _parse_wall(ts: str) -> float:
    parts = ts.strip().split(":")
    if len(parts) == 2:
        return float(parts[0]) * 60 + float(parts[1])
    return float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])


def _count_hits(path: Path) -> int:
    count = 0
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#") and line.strip():
                count += 1
    return count


def run_once(ncores: int, fasta: Path) -> dict:
    fasta_mb = fasta.stat().st_size / 1e6
    with tempfile.TemporaryDirectory() as tmp:
        cmd = [
            "/usr/bin/time", "-v",
            str(TRIPLEXATOR),
            "-ds", str(fasta),
            "-rm", "1",
            "-p",  str(ncores),
            "-of", "0",
            "-od", tmp,
            "-o",  "tts",
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True)
        stderr = proc.stderr + proc.stdout

        rss_m  = _RSS_RE.search(stderr)
        wall_m = _WALL_RE.search(stderr)
        peak_rss_mb = int(rss_m.group(1)) / 1024.0 if rss_m else 0.0
        wall_s      = _parse_wall(wall_m.group(1)) if wall_m else 0.0

        out_file = Path(tmp) / "tts"
        hits = _count_hits(out_file) if out_file.exists() else 0

    return {
        "workers":         ncores,
        "wall_s":          round(wall_s, 2),
        "peak_rss_mb":     round(peak_rss_mb, 2),
        "throughput_mbps": round(fasta_mb / wall_s, 4) if wall_s > 0 else 0.0,
        "total_hits":      hits,
        "returncode":      proc.returncode,
    }


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--workers", default="1,2,4,8,16")
    ap.add_argument("--json", default="benchmarks/data/triplexator_sweep.json")
    args = ap.parse_args()

    worker_list = [int(w) for w in args.workers.split(",") if w.strip()]
    fasta_mb = CHR1_FA.stat().st_size / 1e6
    baseline_wall = None

    print("=" * 64)
    print("  Triplexator Worker Sweep — chr1")
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
            print("FAILED")
            continue
        if baseline_wall is None:
            baseline_wall = r["wall_s"]
        speedup    = baseline_wall / r["wall_s"] if r["wall_s"] > 0 else 0
        efficiency = speedup / w * 100
        print(f"{r['wall_s']:6.1f}s  RSS={r['peak_rss_mb']:6.0f} MB  "
              f"hits={r['total_hits']:>8,}  {speedup:.2f}×  {efficiency:.1f}%")
        r["speedup"]    = round(speedup, 4)
        r["efficiency"] = round(efficiency, 2)
        results.append(r)

    print()
    print("=" * 64)

    payload = {
        "date":    datetime.datetime.now().isoformat(),
        "tool":    "triplexator",
        "version": "1.3.2",
        "fasta":   str(CHR1_FA),
        "fasta_mb": round(fasta_mb, 2),
        "results": results,
    }
    Path(args.json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.json).write_text(json.dumps(payload, indent=2))
    print(f"\n  JSON → {args.json}")


if __name__ == "__main__":
    main()
