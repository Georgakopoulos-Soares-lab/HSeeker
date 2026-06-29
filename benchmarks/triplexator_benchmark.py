#!/usr/bin/env python3
"""
benchmarks/triplexator_benchmark.py
====================================
Benchmark Triplexator (v1.3.2) on chr1 hg38.

Runs:
  /usr/bin/time -v triplexator -ds chr1.fa -rm 1 -p <ncores> -of 0 -o <out>

Records:
  - Wall-clock time   (from /usr/bin/time -v)
  - Peak RSS (MB)     (from /usr/bin/time -v)
  - Hit count         (lines in output not starting with #)
  - Throughput (Mbp/s)

Usage:
  python benchmarks/triplexator_benchmark.py [--cores 16] [--json results.json]
"""

from __future__ import annotations

import argparse
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


def run_benchmark(ncores: int, fasta: Path) -> dict:
    fasta_mb = fasta.stat().st_size / 1e6

    with tempfile.TemporaryDirectory() as tmp:
        out_dir  = tmp
        out_base = "tts"
        out_file = Path(tmp) / out_base   # triplexator writes <dir>/<base> (no extension)

        cmd = [
            "/usr/bin/time", "-v",
            str(TRIPLEXATOR),
            "-ds", str(fasta),
            "-rm", "1",          # parallelize TTSs — best for single long sequence
            "-p",  str(ncores),
            "-of", "0",          # tab-separated output
            "-od", out_dir,      # output directory
            "-o",  out_base,     # output basename (triplexator appends no ext for TTS mode)
        ]

        print(f"Command: {' '.join(cmd[3:])}", flush=True)  # skip /usr/bin/time -v
        proc = subprocess.run(cmd, capture_output=True, text=True)
        stderr = proc.stderr + proc.stdout   # /usr/bin/time writes to stderr

        rss_m  = _RSS_RE.search(stderr)
        wall_m = _WALL_RE.search(stderr)

        peak_rss_mb = int(rss_m.group(1)) / 1024.0 if rss_m else 0.0
        wall_s      = _parse_wall(wall_m.group(1))  if wall_m else 0.0

        hits = _count_hits(out_file) if out_file.exists() else 0

    return {
        "tool":            "triplexator",
        "version":         "1.3.2",
        "fasta":           str(fasta),
        "fasta_mb":        round(fasta_mb, 2),
        "ncores":          ncores,
        "wall_s":          round(wall_s, 2),
        "peak_rss_mb":     round(peak_rss_mb, 2),
        "throughput_mbps": round(fasta_mb / wall_s, 4) if wall_s > 0 else 0.0,
        "total_hits":      hits,
        "returncode":      proc.returncode,
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--cores", type=int, default=16)
    ap.add_argument("--fasta", default=str(CHR1_FA))
    ap.add_argument("--json",  default="benchmarks/data/triplexator_bench.json")
    args = ap.parse_args()

    fa = Path(args.fasta)
    if not fa.exists():
        print(f"ERROR: {fa} not found", file=sys.stderr)
        sys.exit(1)
    if not TRIPLEXATOR.exists():
        print(f"ERROR: {TRIPLEXATOR} not found", file=sys.stderr)
        sys.exit(1)

    fasta_mb = fa.stat().st_size / 1e6
    print("=" * 60)
    print("  Triplexator Benchmark — chr1")
    print("=" * 60)
    print(f"  Binary  : {TRIPLEXATOR}")
    print(f"  FASTA   : {fa.name}  ({fasta_mb:.1f} MB)")
    print(f"  Cores   : {args.cores}")
    print(f"  Mode    : -rm 1 (parallelize TTSs, best for single long seq)")
    print()

    result = run_benchmark(ncores=args.cores, fasta=fa)

    print()
    print("=" * 60)
    print(f"  Total hits   : {result['total_hits']:,}")
    print(f"  Wall time    : {result['wall_s']:.1f} s  ({result['wall_s']/60:.2f} min)")
    print(f"  Peak RSS     : {result['peak_rss_mb']:.0f} MB")
    print(f"  Throughput   : {result['throughput_mbps']:.3f} Mbp/s")
    print(f"  Return code  : {result['returncode']}")
    print("=" * 60)

    Path(args.json).parent.mkdir(parents=True, exist_ok=True)
    Path(args.json).write_text(json.dumps(result, indent=2))
    print(f"\n  JSON → {args.json}")


if __name__ == "__main__":
    main()
