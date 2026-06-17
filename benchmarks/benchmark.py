#!/usr/bin/env python3
"""
benchmarks/benchmark.py — HSeeker End-to-End CLI Benchmark
============================================================

Benchmarks the `hseeker` CLI against real hg38/GRCh38 chromosomes downloaded
from NCBI.  For each chromosome the full end-to-end pipeline is timed:
FASTA read → C-core scan → overlap removal → scoring → TSV write.

Metrics reported
----------------
  Wall time        — total end-to-end duration (s)
  Peak RSS         — maximum resident set size during the run (MB)
  TSV output size  — size of the written ``*_HDNA.tsv`` file (MB)
  FASTA/TSV ratio  — input-to-output size ratio (dimensionless)
  FASTA throughput — MB of FASTA processed per second of wall time

Usage
-----
  python benchmarks/benchmark.py                      # all 24 chromosomes
  python benchmarks/benchmark.py --chromosomes chr1   # single chromosome
  python benchmarks/benchmark.py --chromosomes chr1,chr22 --minrep 15
  python benchmarks/benchmark.py --no-cache           # force re-download
  python benchmarks/benchmark.py --json results.json  # machine-readable output
  python benchmarks/benchmark.py --report report.md   # Markdown report
"""

from __future__ import annotations

import argparse
import datetime
import gzip
import json
import os
import platform
import re
import shutil
import subprocess
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path
from typing import Any

try:
    import numpy as np
    _HAS_NUMPY = True
except ImportError:
    np = None  # type: ignore[assignment]
    _HAS_NUMPY = False

try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False

import hseeker

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATA_DIR = Path(__file__).parent / "data"

CHROMOSOMES: list[str] = [
    f"chr{i}" for i in range(1, 23)
] + ["chrX", "chrY"]

NCBI_FASTA_BASE: str = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
    "GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/"
    "Primary_Assembly/assembled_chromosomes/FASTA/"
)

_W = 100  # report width

# Regex to parse RSS (kB) from /usr/bin/time -v output
_TIME_RSS_RE = re.compile(
    r"Maximum resident set size \(kbytes\):\s*(\d+)"
)
_TIME_WALL_RE = re.compile(
    r"Elapsed \(wall clock\) time \(h:mm:ss or m:ss\):\s*([\d:.]+)"
)
_TIME_EXIT_RE = re.compile(r"Exit status:\s*(\d+)")
_TIME_FS_OUT_RE = re.compile(
    r"File system outputs:\s*(\d+)"
)


# ---------------------------------------------------------------------------
# Chromosome download
# ---------------------------------------------------------------------------

def _download_chromosome(name: str, *, force: bool = False) -> Path | None:
    """Download a single hg38 chromosome .fna.gz from NCBI and decompress.

    Returns the local ``.fa`` path on success, ``None`` on any error.
    """
    fa_path = DATA_DIR / f"{name}.fa"
    if fa_path.exists() and not force:
        size_mb = fa_path.stat().st_size / 1e6
        print(f"  {name:<6}  cached: {fa_path.name}  ({size_mb:.1f} MB)")
        return fa_path

    remote_name = f"{name}.fna.gz"
    url = NCBI_FASTA_BASE + remote_name
    gz_path = DATA_DIR / remote_name
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    print(f"  {name:<6}  downloading {remote_name} … ", end="", flush=True)
    t0 = time.perf_counter()
    try:
        req = urllib.request.Request(
            url, headers={"User-Agent": "hseeker-benchmark/1.0"},
        )
        with urllib.request.urlopen(req, timeout=600) as resp, \
                open(gz_path, "wb") as out:
            shutil.copyfileobj(resp, out)
    except Exception as exc:
        print(f"FAILED ({exc})")
        gz_path.unlink(missing_ok=True)
        return None

    dl_s = time.perf_counter() - t0
    dl_mb = gz_path.stat().st_size / 1e6
    print(f"{dl_mb:.1f} MB  {dl_s:.1f}s", end="", flush=True)

    # Decompress
    print("  →  decompressing … ", end="", flush=True)
    t1 = time.perf_counter()
    try:
        with gzip.open(gz_path, "rb") as f_in, open(fa_path, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)
    except Exception as exc:
        print(f"FAILED ({exc})")
        gz_path.unlink(missing_ok=True)
        fa_path.unlink(missing_ok=True)
        return None

    gz_path.unlink()
    fa_mb = fa_path.stat().st_size / 1e6
    print(f"{fa_mb:.1f} MB  ({time.perf_counter() - t1:.1f}s)")
    return fa_path


# ---------------------------------------------------------------------------
# CLI benchmark runner
# ---------------------------------------------------------------------------

def bench_cli(path: Path, *, minrep: int = 10) -> dict[str, Any] | None:
    """Run ``python -m hseeker`` on *path* and measure wall/RSS/disk I/O.

    Uses ``/usr/bin/time -v`` to capture peak RSS.  Returns a dict with
    keys ``wall_s``, ``peak_rss_mb``, ``output_mb``, ``compression_ratio``,
    ``throughput_mbps``, and ``returncode``, or *None* if the CLI fails.
    """
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = str(Path(tmp) / "bench")
        fasta_mb = path.stat().st_size / 1e6

        cmd = [
            "/usr/bin/time", "-v",
            sys.executable, "-m", "hseeker",
            "-seq", str(path),
            "-out", out_prefix,
            "-minrep", str(minrep),
        ]
        proc = subprocess.run(
            cmd, capture_output=True, text=True,
        )
        # /usr/bin/time -v writes to stderr
        time_output = proc.stderr
        rc = proc.returncode

        # Parse /usr/bin/time -v output
        rss_match = _TIME_RSS_RE.search(time_output)
        wall_match = _TIME_WALL_RE.search(time_output)

        peak_rss_kb = int(rss_match.group(1)) if rss_match else 0
        wall_s = _parse_time_str(wall_match.group(1)) if wall_match else 0.0

        tsv = Path(out_prefix + "_HDNA.tsv")
        out_mb = tsv.stat().st_size / 1e6 if tsv.exists() else 0.0

        return {
            "wall_s":            wall_s,
            "peak_rss_mb":       peak_rss_kb / 1024.0,
            "output_mb":         out_mb,
            "compression_ratio": fasta_mb / out_mb if out_mb > 0 else 0.0,
            "throughput_mbps":   fasta_mb / wall_s if wall_s > 0 else 0.0,
            "returncode":        rc,
        }


def _parse_time_str(ts: str) -> float:
    """Parse /usr/bin/time elapsed string like '0:15.93' or '1:23:45'."""
    ts = ts.strip()
    parts = ts.split(":")
    if len(parts) == 2:  # m:ss.ss
        return float(parts[0]) * 60 + float(parts[1])
    elif len(parts) == 3:  # h:mm:ss
        return float(parts[0]) * 3600 + float(parts[1]) * 60 + float(parts[2])
    return 0.0


# ---------------------------------------------------------------------------
# Report formatting
# ---------------------------------------------------------------------------

def _print_separator(char: str = "─") -> None:
    print(char * _W)


# ---------------------------------------------------------------------------
# System info
# ---------------------------------------------------------------------------

def _get_system_info() -> dict[str, Any]:
    info: dict[str, Any] = {
        "os":                 platform.platform(),
        "python":             sys.version.split()[0],
        "hseeker":            hseeker.__version__,
        "cpu_cores_logical":  os.cpu_count(),
    }
    if _HAS_NUMPY and np is not None:
        info["numpy"] = np.__version__

    for cmd in (
        ["sysctl", "-n", "machdep.cpu.brand_string"],
        ["grep", "-m1", "model name", "/proc/cpuinfo"],
    ):
        try:
            out = subprocess.run(cmd, capture_output=True, text=True, timeout=3)
            if out.returncode == 0 and out.stdout.strip():
                brand = out.stdout.strip()
                if "model name" in brand:
                    brand = brand.split(":", 1)[1].strip()
                info["cpu_brand"] = brand
                break
        except Exception:
            pass

    if _HAS_PSUTIL:
        vm = psutil.virtual_memory()
        info["ram_total_gb"]     = round(vm.total    / 1e9, 1)
        info["ram_available_gb"] = round(vm.available / 1e9, 1)
        info["cpu_cores_physical"] = psutil.cpu_count(logical=False)
        try:
            freq = psutil.cpu_freq()
            if freq and freq.max and freq.max > 0:
                info["cpu_freq_max_ghz"] = round(freq.max / 1000, 2)
        except Exception:
            pass
        try:
            info["disk_free_gb"] = round(
                psutil.disk_usage(str(DATA_DIR.parent)).free / 1e9, 1
            )
        except Exception:
            pass

    return info


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------

def _generate_markdown_report(
    sys_info: dict[str, Any],
    results: list[dict[str, Any]],
    minrep: int,
) -> str:
    lines: list[str] = []
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    lines += [
        "# HSeeker CLI Performance Benchmark Report",
        "",
        f"**Date**: {now}  ",
        f"**hseeker version**: `{sys_info.get('hseeker', 'unknown')}`  ",
        f"**minrep**: {minrep}  ",
        f"**Scoring**: enabled  ",
        "",
        "## System Specifications",
        "",
        "| Property | Value |",
        "|---|---|",
    ]
    if "cpu_brand" in sys_info:
        lines.append(f"| CPU | {sys_info['cpu_brand']} |")
    phys = sys_info.get("cpu_cores_physical", "—")
    logi = sys_info.get("cpu_cores_logical", "—")
    lines.append(f"| CPU cores | {phys} physical / {logi} logical |")
    if sys_info.get("cpu_freq_max_ghz"):
        lines.append(f"| CPU max clock | {sys_info['cpu_freq_max_ghz']} GHz |")
    if "ram_total_gb" in sys_info:
        lines.append(f"| RAM (total) | {sys_info['ram_total_gb']} GB |")
    if "ram_available_gb" in sys_info:
        lines.append(f"| RAM (available) | {sys_info['ram_available_gb']} GB |")
    lines.append(f"| OS | {sys_info.get('os', '—')} |")
    lines.append(f"| Python | {sys_info.get('python', '—')} |")
    if "numpy" in sys_info:
        lines.append(f"| NumPy | {sys_info['numpy']} |")
    if "disk_free_gb" in sys_info:
        lines.append(f"| Disk free | {sys_info['disk_free_gb']} GB |")
    lines += ["", "## Methodology", "",
              "The `hseeker` CLI is invoked against each hg38 chromosome "
              "FASTA file.  Every run is measured with `/usr/bin/time -v` "
              "which reports wall-clock time and peak RSS.  "
              "Scoring is always enabled (default).",
              ""]

    # Results table
    lines += [
        "## Results",
        "",
        "| Chromosome | FASTA (MB) | Wall (s) | Peak RSS (MB) | "
        "TSV (MB) | Ratio | Throughput (MB/s) |",
        "|---|---|---|---|---|---|---|",
    ]
    total_mb = 0.0
    total_wall = 0.0
    total_rss = 0.0
    total_tsv = 0.0
    count = 0

    for r in results:
        name = r["chromosome"]
        fm = r["fasta_mb"]
        w  = r["wall_s"]
        rss = r["peak_rss_mb"]
        om = r["output_mb"]
        ratio = r["compression_ratio"]
        tp  = r["throughput_mbps"]
        lines.append(
            f"| {name} | {fm:.1f} | {w:.1f} | {rss:.0f} | "
            f"{om:.1f} | {ratio:.1f}× | {tp:.2f} |"
        )
        total_mb += fm
        total_wall += w
        total_rss = max(total_rss, rss)
        total_tsv += om
        count += 1

    if count > 0:
        avg_tp = total_mb / total_wall if total_wall > 0 else 0.0
        avg_ratio = total_mb / total_tsv if total_tsv > 0 else 0.0
        lines += [
            "",
            "## Summary",
            "",
            f"| Metric | Value |",
            "|---|---|",
            f"| Chromosomes benchmarked | {count} |",
            f"| Total FASTA size | {total_mb:.1f} MB |",
            f"| Total wall time | {total_wall:.1f} s ({total_wall/60:.1f} min) |",
            f"| Peak RSS (max) | {total_rss:.0f} MB |",
            f"| Total TSV output | {total_tsv:.1f} MB |",
            f"| Aggregate FASTA/TSV ratio | {avg_ratio:.1f}× |",
            f"| Aggregate throughput | {avg_tp:.2f} MB/s |",
            "",
            f"**Estimated time for full genome (3.1 GB):** "
            f"{3100 / avg_tp / 3600:.1f} hours" if avg_tp > 0 else "",
            "",
        ]

    lines += [
        "---",
        f"*Generated by HSeeker benchmark suite — {now}*",
        "",
    ]
    return "\n".join(lines)


def _build_concat_if_needed(chrom_paths: dict[str, "Path | None"]) -> None:
    """Concatenate chr1 + chr2 + chr3 into ``data/chr1_2_3.fa`` if all three
    are available and the combined file is not already up to date.
    """
    concat_path = DATA_DIR / "chr1_2_3.fa"
    sources = [chrom_paths.get(n) for n in ("chr1", "chr2", "chr3")]
    if any(p is None or not p.exists() for p in sources):
        return  # not all three are present; skip silently

    # Re-build if any source is newer than the concat file.
    if concat_path.exists():
        concat_mtime = concat_path.stat().st_mtime
        if all(p.stat().st_mtime <= concat_mtime for p in sources):  # type: ignore[union-attr]
            size_mb = concat_path.stat().st_size / 1e6
            print(f"  chr1_2_3 cached: {concat_path.name}  ({size_mb:.1f} MB)")
            return

    print("  Building chr1_2_3.fa (chr1 + chr2 + chr3) … ", end="", flush=True)
    t0 = time.perf_counter()
    with open(concat_path, "wb") as out:
        for p in sources:
            with open(p, "rb") as src:  # type: ignore[arg-type]
                shutil.copyfileobj(src, out)
    size_mb = concat_path.stat().st_size / 1e6
    print(f"{size_mb:.1f} MB  ({time.perf_counter() - t0:.1f}s)")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--minrep", type=int, default=10,
        help="Minimum arm length (default: 10)",
    )
    p.add_argument(
        "--no-cache", action="store_true",
        help="Re-download chromosome FASTA files",
    )
    p.add_argument(
        "--json", metavar="FILE",
        help="Save machine-readable results to a JSON file",
    )
    p.add_argument(
        "--report", metavar="FILE",
        help="Write a Markdown performance report to FILE",
    )
    p.add_argument(
        "--chromosomes", type=str, default=None, metavar="chr1,chr2,...",
        help=(
            "Comma-separated list of chromosomes to benchmark "
            "(default: all 24 — chr1–chr22, chrX, chrY)"
        ),
    )
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    # ── Resolve chromosome list ──────────────────────────────────────────
    if args.chromosomes:
        chrom_list = [
            c.strip() for c in args.chromosomes.split(",") if c.strip()
        ]
        unknown = set(chrom_list) - set(CHROMOSOMES)
        if unknown:
            print(
                f"ERROR: unknown chromosome(s): {', '.join(sorted(unknown))}",
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        chrom_list = list(CHROMOSOMES)

    # ── Header ────────────────────────────────────────────────────────────
    sys_info = _get_system_info()
    _print_separator("=")
    print("  HSeeker CLI Performance Benchmark".center(_W))
    print(datetime.datetime.now().strftime("  %Y-%m-%d %H:%M:%S").center(_W))
    _print_separator("=")
    print(f"  hseeker    : {sys_info['hseeker']}")
    print(f"  Python     : {sys_info['python']}")
    if "cpu_brand" in sys_info:
        print(f"  CPU        : {sys_info['cpu_brand']}")
    phys = sys_info.get("cpu_cores_physical", "?")
    logi = sys_info.get("cpu_cores_logical", "?")
    print(f"  CPU cores  : {logi} logical / {phys} physical")
    print(f"  minrep     : {args.minrep}  |  Scoring: on")
    if "ram_total_gb" in sys_info:
        print(
            f"  RAM        : {sys_info['ram_total_gb']} GB total, "
            f"{sys_info.get('ram_available_gb', '?')} GB available"
        )
    print()

    # ── Download chromosomes ──────────────────────────────────────────────
    print("── Downloading hg38 chromosomes ──────────────────────────────────")
    print(f"  Data directory: {DATA_DIR}")
    print()

    chrom_paths: dict[str, Path | None] = {}
    for name in chrom_list:
        path = _download_chromosome(name, force=args.no_cache)
        chrom_paths[name] = path

    available = [
        (name, p) for name, p in chrom_paths.items()
        if p is not None and p.exists()
    ]
    if not available:
        print("\n  ERROR: No chromosomes available.", file=sys.stderr)
        sys.exit(1)

    print(f"\n  {len(available)} chromosome(s) ready.")

    # ── Build chr1+chr2+chr3 concat if all three are present ─────────────
    _build_concat_if_needed(chrom_paths)
    print()

    # ── Benchmark ─────────────────────────────────────────────────────────
    print("── Benchmarking (CLI: python -m hseeker) ──────────────────────────")

    all_results: list[dict[str, Any]] = []
    for name, path in available:
        size_mb = path.stat().st_size / 1e6
        print(f"\n  [{name}]  {size_mb:.1f} MB  … ", end="", flush=True)
        r = bench_cli(path, minrep=args.minrep)
        if r is None or r["returncode"] != 0:
            print("FAILED")
            continue
        r["chromosome"] = name
        r["fasta_mb"] = size_mb
        all_results.append(r)
        print(
            f"{r['wall_s']:.1f}s  "
            f"RSS={r['peak_rss_mb']:.0f}MB  "
            f"TSV={r['output_mb']:.1f}MB  "
            f"ratio={r['compression_ratio']:.1f}×  "
            f"{r['throughput_mbps']:.2f} MB/s"
        )

    if not all_results:
        print("\n  No successful benchmarks.", file=sys.stderr)
        sys.exit(1)

    # ── Summary table ─────────────────────────────────────────────────────
    total_mb = sum(r["fasta_mb"] for r in all_results)
    total_wall = sum(r["wall_s"] for r in all_results)
    max_rss = max(r["peak_rss_mb"] for r in all_results)
    total_tsv = sum(r["output_mb"] for r in all_results)
    avg_tp = total_mb / total_wall if total_wall > 0 else 0

    print("\n\n── Summary ────────────────────────────────────────────────────────")
    _print_separator("·")
    print(
        f"  {'Chromosome':<12}  {'FASTA MB':>9}  {'Wall':>8}  "
        f"{'RSS MB':>7}  {'TSV MB':>8}  {'Ratio':>7}  {'MB/s':>8}"
    )
    _print_separator("·")
    for r in all_results:
        print(
            f"  {r['chromosome']:<12}  "
            f"{r['fasta_mb']:>8.1f}  "
            f"{r['wall_s']:>7.1f}s  "
            f"{r['peak_rss_mb']:>6.0f}  "
            f"{r['output_mb']:>7.1f}  "
            f"{r['compression_ratio']:>6.1f}×  "
            f"{r['throughput_mbps']:>7.2f}"
        )
    _print_separator("·")
    print(
        f"  {'TOTAL/MAX':<12}  "
        f"{total_mb:>8.1f}  "
        f"{total_wall:>7.1f}s  "
        f"{max_rss:>6.0f}  "
        f"{total_tsv:>7.1f}  "
        f"{total_mb/total_tsv if total_tsv > 0 else 0:>6.1f}×  "
        f"{avg_tp:>7.2f}"
    )
    print(
        f"\n  Estimated full genome (3.1 GB): "
        f"{3100 / avg_tp / 3600:.1f} hours" if avg_tp > 0 else ""
    )

    # ── JSON output ───────────────────────────────────────────────────────
    if args.json:
        out = Path(args.json)
        payload = {
            "system_info": sys_info,
            "minrep":      args.minrep,
            "results":     all_results,
        }
        out.write_text(json.dumps(payload, indent=2, default=str))
        print(f"\n  JSON → {out}")

    # ── Markdown report ───────────────────────────────────────────────────
    if args.report:
        report_path = Path(args.report)
        report_text = _generate_markdown_report(
            sys_info, all_results, args.minrep,
        )
        report_path.write_text(report_text, encoding="utf-8")
        print(f"  Report → {report_path}")

    _print_separator("=")
    print("  Done.")
    _print_separator("=")


if __name__ == "__main__":
    main()
