#!/usr/bin/env python3
"""
benchmarks/benchmark.py — HSeeker Performance Benchmark Suite
==============================================================

Generates synthetic FASTA datasets at three size tiers and measures wall
time, peak RAM, CPU utilisation, parallelism speedup, and disk I/O for
each public API path.

Sequence profiles
-----------------
  uniform    — equal ACGT probability (sparse H-DNA; baseline throughput)
  ga_biased  — 45 % G + 45 % A (dense hits; stresses hit buffer & overlap
                removal code)
  realistic  — slight AT bias + GC blocks + embedded motifs (mimics human
                chr in both content and H-DNA site density)

Size tiers
----------
  small   (default) — 30  MB  total: 5 × 6 MB   and 24 × 1.25 MB records
  medium  (--medium)— 300 MB  total: 5 × 60 MB  and 24 × 12.5 MB records
  large   (--large) — 3   GB  total: 5 × 600 MB and 24 × 125  MB records
                      ⚠ requires ~9 GB free disk space and 10–30 min

Methods benchmarked
-------------------
  scan_fasta            — sequential, returns list
  scan_fasta_iter       — sequential streaming; peak-RAM compared to above
  scan_fasta_parallel   — threaded; benchmarked at workers=1 and workers=N
  CLI (python -m hseeker) — end-to-end including TSV disk write

Usage
-----
  python benchmarks/benchmark.py                   # small only (~1 min)
  python benchmarks/benchmark.py --medium          # small + medium
  python benchmarks/benchmark.py --large           # all tiers
  python benchmarks/benchmark.py --real genome.fa  # also bench a real file
  python benchmarks/benchmark.py --workers 8       # override thread count
  python benchmarks/benchmark.py --no-cli          # skip disk-write bench
  python benchmarks/benchmark.py --json out.json   # save machine-readable results
  python benchmarks/benchmark.py --no-cache        # regenerate FASTA files
"""

from __future__ import annotations

import argparse
import contextlib
import gc
import json
import os
import random
import subprocess
import sys
import tempfile
import time
import tracemalloc
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Any, Generator

# numpy is always available: pandas (required dep) pulls it in
import numpy as np

try:
    import psutil
    _HAS_PSUTIL = True
except ImportError:
    _HAS_PSUTIL = False

import hseeker

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DATA_DIR   = Path(__file__).parent / "data"
LINE_WIDTH = 60          # bases per FASTA line
GEN_CHUNK  = 10_000_000  # bases generated per numpy call (≈10 MB in RAM)

# ACGT weights: (A, C, G, T)
_PROFILES: dict[str, tuple[float, float, float, float]] = {
    "uniform":   (0.25,  0.25,  0.25,  0.25),
    "ga_biased": (0.05,  0.05,  0.45,  0.45),  # 90 % purine → very dense hits
    "realistic": (0.295, 0.205, 0.205, 0.295),  # ~41 % GC, like human chr
}

# H-DNA motifs embedded at regular intervals — chosen to produce real hits
# at minrep=6 with default purity/mismatch settings
_MOTIFS: list[bytes] = [
    b"GAGAGAGAGAGAGAGA",      # pure GA 8-mer, is_perfect candidate
    b"GGGAAAGGGTTTTCCCAAACCC",# mixed GA, arm≈10
    b"AAAAAAAAAAAAAA",         # poly-A 14-mer, is_perfect
    b"GAGAGAGTTTTCTCTCTC",     # impure GA/CT mirror
    b"CCCTTTAATTTCCC",         # CT mirror (14 bp)
]

_BASE_ASCII = np.array([ord("A"), ord("C"), ord("G"), ord("T")], dtype=np.uint8)

# ---------------------------------------------------------------------------
# FASTA generation
# ---------------------------------------------------------------------------

def _write_fasta(
    path: Path,
    n_records: int,
    record_length: int,
    profile: str,
    seed: int = 42,
) -> None:
    """Write a synthetic FASTA file using numpy for fast sequence generation."""
    weights   = np.array(_PROFILES[profile])
    rng       = np.random.default_rng(seed)
    py_rng    = random.Random(seed)
    # Embed one motif every ~motif_period bases
    motif_period = max(300, record_length // 3000)

    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", buffering=1 << 20) as fh:
        for rec in range(n_records):
            fh.write(f">chr{rec + 1}\n")

            remaining = record_length
            genome_offset = 0   # position within this record

            while remaining > 0:
                chunk = min(GEN_CHUNK, remaining)

                # Generate random bases
                idx = rng.choice(4, size=chunk, p=weights).astype(np.uint8)
                bases_arr = _BASE_ASCII[idx]   # uint8 ASCII array

                # Embed H-DNA motifs
                pos = py_rng.randint(0, motif_period)
                while pos < chunk:
                    motif = py_rng.choice(_MOTIFS)
                    end   = min(pos + len(motif), chunk)
                    motif_arr = np.frombuffer(motif, dtype=np.uint8)
                    bases_arr[pos:end] = motif_arr[: end - pos]
                    pos += py_rng.randint(
                        motif_period * 3 // 4,
                        motif_period * 5 // 4,
                    )

                # Write line-wrapped: reshape into rows of LINE_WIDTH
                n_full = chunk // LINE_WIDTH
                partial = chunk % LINE_WIDTH

                if n_full > 0:
                    grid = bases_arr[: n_full * LINE_WIDTH].reshape(
                        n_full, LINE_WIDTH
                    )
                    # Append newline column
                    nl = np.full((n_full, 1), ord("\n"), dtype=np.uint8)
                    grid_nl = np.hstack([grid, nl])
                    fh.write(grid_nl.tobytes().decode("ascii"))

                if partial:
                    fh.write(bases_arr[-partial:].tobytes().decode("ascii"))
                    fh.write("\n")

                remaining -= chunk
                genome_offset += chunk


def _ensure_dataset(
    name: str,
    n_records: int,
    record_length: int,
    profile: str,
    seed: int,
    force: bool = False,
) -> Path:
    """Return path to cached FASTA file, generating it if absent."""
    path = DATA_DIR / f"{name}.fa"
    if path.exists() and not force:
        return path

    total_mb = n_records * record_length / 1e6
    print(
        f"  generating {name}.fa  "
        f"({n_records} records × {record_length / 1e6:.1f} MB = {total_mb:.0f} MB) …",
        flush=True,
    )
    t0 = time.perf_counter()
    _write_fasta(path, n_records, record_length, profile, seed)
    elapsed = time.perf_counter() - t0
    size_mb = path.stat().st_size / 1e6
    print(f"    → {size_mb:.1f} MB written in {elapsed:.1f} s  "
          f"({size_mb / elapsed:.0f} MB/s)")
    return path


# ---------------------------------------------------------------------------
# Measurement utilities
# ---------------------------------------------------------------------------

def _rss_bytes() -> int | None:
    if _HAS_PSUTIL:
        return psutil.Process().memory_info().rss
    return None


@contextlib.contextmanager
def _measure() -> Generator[dict[str, Any], None, None]:
    """Context manager; yields a mutable result dict populated on exit."""
    result: dict[str, Any] = {}
    gc.collect()
    tracemalloc.start()
    rss_before  = _rss_bytes()
    cpu_before  = time.process_time()
    wall_before = time.perf_counter()
    try:
        yield result
    finally:
        result["wall_s"]  = time.perf_counter() - wall_before
        result["cpu_s"]   = time.process_time() - cpu_before
        _, peak = tracemalloc.get_traced_memory()
        tracemalloc.stop()
        result["peak_ram_mb"] = peak / 1e6
        rss_after = _rss_bytes()
        if rss_before is not None and rss_after is not None:
            result["rss_delta_mb"] = (rss_after - rss_before) / 1e6


# ---------------------------------------------------------------------------
# Benchmark runners
# ---------------------------------------------------------------------------

def bench_scan_fasta(path: Path, minrep: int = 6) -> dict[str, Any]:
    with _measure() as r:
        hits = hseeker.scan_fasta(str(path), minrep=minrep)
    r["hits"] = len(hits)
    return r


def bench_scan_fasta_iter(path: Path, minrep: int = 6) -> dict[str, Any]:
    with _measure() as r:
        hits = list(hseeker.scan_fasta_iter(str(path), minrep=minrep))
    r["hits"] = len(hits)
    return r


def bench_scan_fasta_parallel(
    path: Path,
    workers: int,
    minrep: int = 6,
) -> dict[str, Any]:
    with _measure() as r:
        hits = hseeker.scan_fasta_parallel(
            str(path), minrep=minrep, workers=workers
        )
    r["hits"] = len(hits)
    return r


def bench_cli(path: Path, minrep: int = 6) -> dict[str, Any]:
    """End-to-end CLI benchmark including TSV disk write."""
    with tempfile.TemporaryDirectory() as tmp:
        out_prefix = str(Path(tmp) / "bench")
        cmd = [
            sys.executable, "-m", "hseeker",
            "-seq", str(path),
            "-out", out_prefix,
            "-minrep", str(minrep),
        ]
        t0 = time.perf_counter()
        proc = subprocess.run(cmd, capture_output=True, text=True)
        wall_s = time.perf_counter() - t0

        tsv = Path(out_prefix + "_HDNA.tsv")
        out_mb = tsv.stat().st_size / 1e6 if tsv.exists() else 0.0
        fasta_mb = path.stat().st_size / 1e6

        return {
            "wall_s":          wall_s,
            "returncode":      proc.returncode,
            "output_mb":       out_mb,
            "compression_ratio": fasta_mb / out_mb if out_mb > 0 else 0.0,
            "throughput_mbps": fasta_mb / wall_s if wall_s > 0 else 0.0,
        }


# ---------------------------------------------------------------------------
# Dataset specs
# ---------------------------------------------------------------------------

@dataclass
class DatasetSpec:
    name:          str
    n_records:     int
    record_length: int    # bases
    profile:       str
    seed:          int = 42

    @property
    def total_mb(self) -> float:
        return self.n_records * self.record_length / 1e6


# fmt: off
SMALL_DATASETS: list[DatasetSpec] = [
    DatasetSpec("small_uniform_5rec",    5,  6_000_000, "uniform",   seed=1),
    DatasetSpec("small_ga_biased_5rec",  5,  6_000_000, "ga_biased", seed=2),
    DatasetSpec("small_realistic_24rec", 24, 1_250_000, "realistic", seed=3),
]

MEDIUM_DATASETS: list[DatasetSpec] = [
    DatasetSpec("medium_uniform_5rec",    5,  60_000_000, "uniform",   seed=4),
    DatasetSpec("medium_ga_biased_5rec",  5,  60_000_000, "ga_biased", seed=5),
    DatasetSpec("medium_realistic_24rec", 24, 12_500_000, "realistic", seed=6),
]

LARGE_DATASETS: list[DatasetSpec] = [
    DatasetSpec("large_uniform_5rec",    5,  600_000_000, "uniform",   seed=7),
    DatasetSpec("large_ga_biased_5rec",  5,  600_000_000, "ga_biased", seed=8),
    DatasetSpec("large_realistic_24rec", 24, 125_000_000, "realistic", seed=9),
]
# fmt: on


# ---------------------------------------------------------------------------
# Report formatting
# ---------------------------------------------------------------------------

_W = 100  # report width


def _print_separator(char: str = "─") -> None:
    print(char * _W)


def _fmt_ram(r: dict[str, Any]) -> str:
    s = f"{r['peak_ram_mb']:6.1f} MB (tracemalloc)"
    if "rss_delta_mb" in r:
        s += f"  RSS Δ {r['rss_delta_mb']:+.0f} MB"
    return s


def _print_results_table(
    results: dict[str, dict[str, Any]],
    baseline_key: str = "scan_fasta",
) -> None:
    baseline_wall = results.get(baseline_key, {}).get("wall_s")
    header = (
        f"  {'Method':<38}  {'Wall':>7}  {'CPU':>7}  "
        f"{'PeakRAM':>8}  {'Speedup':>7}  {'Hits':>11}"
    )
    print(header)
    print("  " + "·" * (_W - 2))

    for method, r in results.items():
        if "wall_s" not in r:   # CLI result format is different
            continue
        speedup_str = "  —   "
        if baseline_wall and r["wall_s"] > 0 and method != baseline_key:
            speedup_str = f"{baseline_wall / r['wall_s']:5.2f}×"

        hits_str = f"{r.get('hits', '—'):>11,}" if isinstance(r.get("hits"), int) else f"{'—':>11}"

        rss = f"  RSS Δ{r['rss_delta_mb']:+.0f}MB" if "rss_delta_mb" in r else ""
        print(
            f"  {method:<38}  "
            f"{r['wall_s']:>6.2f}s  "
            f"{r['cpu_s']:>6.2f}s  "
            f"{r['peak_ram_mb']:>7.1f}MB  "
            f"{speedup_str:>7}  "
            f"{hits_str}"
            f"{rss}"
        )

    # CLI row (different keys)
    if "cli" in results:
        c = results["cli"]
        print(
            f"  {'cli (incl. TSV write)':<38}  "
            f"{c['wall_s']:>6.2f}s  "
            f"{'—':>7}   "
            f"{'—':>7}   "
            f"{'—':>7}  "
            f"{'—':>11}   "
            f"→ {c['output_mb']:.1f} MB TSV  "
            f"{c['throughput_mbps']:.1f} MB/s FASTA throughput"
        )


def _print_dataset_block(
    spec: DatasetSpec,
    path: Path,
    results: dict[str, dict[str, Any]],
    workers: int,
) -> None:
    _print_separator()
    print(
        f"  Dataset  : {spec.name}"
        f"   profile={spec.profile}"
        f"   {spec.n_records} records × {spec.record_length / 1e6:.1f} MB"
        f"   = {path.stat().st_size / 1e6:.1f} MB on disk"
    )
    _print_separator("·")

    _print_results_table(results)

    # Parallelism summary
    seq  = results.get("scan_fasta")
    par1 = results.get("parallel_1w")
    parN = results.get(f"parallel_{workers}w")
    if seq and par1 and parN and seq["wall_s"] > 0 and parN["wall_s"] > 0:
        overhead   = par1["wall_s"] - seq["wall_s"]  # may be negative (noise)
        speedup    = seq["wall_s"] / parN["wall_s"]
        efficiency = speedup / workers * 100
        print(
            f"\n  Parallelism  workers=1 overhead: {overhead:+.3f}s | "
            f"workers={workers}: {speedup:.2f}× speedup = {efficiency:.0f}% efficiency"
        )


def _print_scaling_table(
    path: Path,
    workers_list: list[int],
    minrep: int,
    label: str,
) -> None:
    """Show how wall time and speedup scale with worker count."""
    _print_separator()
    print(f"  Parallelism scaling — {label}  ({path.stat().st_size / 1e6:.1f} MB)")
    _print_separator("·")
    print(f"  {'workers':>8}  {'wall_s':>8}  {'cpu_s':>8}  {'speedup':>8}  {'efficiency':>10}")

    baseline = None
    for w in workers_list:
        print(f"  {w:>8}  ", end="", flush=True)
        r = bench_scan_fasta_parallel(path, workers=w, minrep=minrep)
        if baseline is None and w == 1:
            baseline = r["wall_s"]
        speedup    = baseline / r["wall_s"] if baseline and r["wall_s"] > 0 else 0
        efficiency = speedup / w * 100
        print(
            f"{r['wall_s']:>8.2f}  "
            f"{r['cpu_s']:>8.2f}  "
            f"{speedup:>8.2f}×  "
            f"{efficiency:>9.0f}%"
        )


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--medium", action="store_true",
        help="Include 300 MB datasets (in addition to 30 MB)",
    )
    p.add_argument(
        "--large", action="store_true",
        help="Include 3 GB datasets (requires ~9 GB free disk, 10–30 min)",
    )
    p.add_argument(
        "--real", metavar="FASTA",
        help="Also benchmark against this existing FASTA file",
    )
    p.add_argument(
        "--workers", type=int, default=None,
        help="Thread count for parallel scan (default: os.cpu_count())",
    )
    p.add_argument(
        "--minrep", type=int, default=6,
        help="Minimum arm length passed to hseeker (default 6)",
    )
    p.add_argument(
        "--no-cli", action="store_true",
        help="Skip the CLI / disk-write benchmark",
    )
    p.add_argument(
        "--no-cache", action="store_true",
        help="Regenerate FASTA files even if they already exist",
    )
    p.add_argument(
        "--json", metavar="FILE",
        help="Save machine-readable results to a JSON file",
    )
    p.add_argument(
        "--scaling", action="store_true",
        help="Run parallelism scaling table on the largest dataset of each tier",
    )
    return p.parse_args()


def _run_dataset(
    spec: DatasetSpec,
    workers: int,
    minrep: int,
    force: bool,
    cli: bool,
) -> tuple[Path, dict[str, dict[str, Any]]]:
    path    = _ensure_dataset(
        spec.name, spec.n_records, spec.record_length,
        spec.profile, spec.seed, force=force,
    )
    results: dict[str, dict[str, Any]] = {}

    for label, fn in [
        ("scan_fasta",      lambda: bench_scan_fasta(path, minrep)),
        ("scan_fasta_iter", lambda: bench_scan_fasta_iter(path, minrep)),
        ("parallel_1w",     lambda: bench_scan_fasta_parallel(path, 1, minrep)),
    ]:
        print(f"    {label} … ", end="", flush=True)
        r = fn()
        results[label] = r
        print(f"{r['wall_s']:.2f}s  ({r.get('hits', '?'):,} hits)")

    if workers > 1:
        label = f"parallel_{workers}w"
        print(f"    {label} … ", end="", flush=True)
        r = bench_scan_fasta_parallel(path, workers, minrep)
        results[label] = r
        print(f"{r['wall_s']:.2f}s")

    if cli:
        print(f"    cli … ", end="", flush=True)
        r = bench_cli(path, minrep)
        results["cli"] = r
        rc = r["returncode"]
        if rc != 0:
            print(f"FAILED (exit {rc})")
        else:
            print(
                f"{r['wall_s']:.2f}s  "
                f"→ {r['output_mb']:.1f} MB TSV  "
                f"{r['throughput_mbps']:.1f} MB/s"
            )

    return path, results


def main() -> None:
    args    = _parse_args()
    workers = args.workers or os.cpu_count() or 1

    # ── header ────────────────────────────────────────────────────────────
    import datetime
    _print_separator("=")
    print("  HSeeker Performance Benchmark".center(_W))
    print(datetime.datetime.now().strftime("  %Y-%m-%d %H:%M:%S").center(_W))
    _print_separator("=")
    print(f"  hseeker    : {hseeker.__version__}")
    print(f"  Python     : {sys.version.split()[0]}")
    print(f"  CPU cores  : {os.cpu_count()}")
    print(f"  Workers    : {workers}")
    print(f"  minrep     : {args.minrep}")
    if _HAS_PSUTIL:
        vm = psutil.virtual_memory()
        print(
            f"  RAM        : {vm.total / 1e9:.1f} GB total, "
            f"{vm.available / 1e9:.1f} GB available"
        )
    else:
        print("  RAM        : install psutil for RSS delta measurements")
    if not args.no_cli:
        print("  CLI bench  : yes (measures TSV disk-write throughput)")
    print()

    # ── dataset selection ─────────────────────────────────────────────────
    specs: list[DatasetSpec] = list(SMALL_DATASETS)
    if args.medium or args.large:
        specs += MEDIUM_DATASETS
    if args.large:
        specs += LARGE_DATASETS

    # ── run all datasets ──────────────────────────────────────────────────
    all_output: list[dict[str, Any]] = []
    completed: list[tuple[DatasetSpec, Path, dict]] = []

    print("\n── Generating datasets ──────────────────────────────────────────")
    print(
        f"  Data directory: {DATA_DIR}\n"
        f"  (Delete this directory to force regeneration)\n"
    )

    for spec in specs:
        print(f"\n[{spec.name}]  profile={spec.profile}  "
              f"{spec.total_mb:.0f} MB total")
        path, results = _run_dataset(
            spec, workers, args.minrep, args.no_cache, not args.no_cli
        )
        completed.append((spec, path, results))
        all_output.append({
            "dataset":   spec.name,
            "profile":   spec.profile,
            "total_mb":  spec.total_mb,
            "n_records": spec.n_records,
            "workers":   workers,
            "results":   results,
        })

    # ── optional: real FASTA ──────────────────────────────────────────────
    if args.real:
        real_path = Path(args.real)
        if not real_path.exists():
            print(f"\n[real] ERROR: {real_path} not found", file=sys.stderr)
        else:
            print(
                f"\n[real: {real_path.name}]  "
                f"{real_path.stat().st_size / 1e6:.1f} MB"
            )
            fake_spec = DatasetSpec(
                name=f"real:{real_path.name}",
                n_records=0,   # not used
                record_length=0,
                profile="(real)",
            )
            for label, fn in [
                ("scan_fasta",      lambda: bench_scan_fasta(real_path, args.minrep)),
                ("scan_fasta_iter", lambda: bench_scan_fasta_iter(real_path, args.minrep)),
                ("parallel_1w",     lambda: bench_scan_fasta_parallel(real_path, 1, args.minrep)),
                (f"parallel_{workers}w",
                 lambda: bench_scan_fasta_parallel(real_path, workers, args.minrep)),
            ]:
                print(f"    {label} … ", end="", flush=True)
                r = fn()
                print(f"{r['wall_s']:.2f}s  ({r.get('hits', '?'):,} hits)")

    # ── detailed report ───────────────────────────────────────────────────
    print("\n\n── Detailed Results ─────────────────────────────────────────────")
    for spec, path, results in completed:
        _print_dataset_block(spec, path, results, workers)

    # ── parallelism scaling table ─────────────────────────────────────────
    if args.scaling and workers > 1:
        print("\n\n── Parallelism Scaling ──────────────────────────────────────────")
        # Use the realistic multi-record dataset from each tier (best for
        # record-level parallelism since it has the most records)
        scale_targets = [
            (s, p) for s, p, _ in completed if "realistic" in s.name
        ]
        # Candidate worker counts up to cpu_count
        cpu = os.cpu_count() or 1
        w_list = sorted({1, 2, min(4, cpu), min(8, cpu), cpu})
        for spec, path in scale_targets:
            _print_scaling_table(
                path, w_list, args.minrep, spec.name
            )

    # ── summary: parallelism benefit table ───────────────────────────────
    print("\n\n── Parallelism Benefit Summary ──────────────────────────────────")
    _print_separator("·")
    print(
        f"  {'Dataset':<38}  {'seq wall':>9}  "
        f"{'par wall':>9}  {'speedup':>8}  {'efficiency':>10}"
    )
    _print_separator("·")
    for spec, path, results in completed:
        seq = results.get("scan_fasta")
        par = results.get(f"parallel_{workers}w")
        if not seq or not par:
            continue
        speedup    = seq["wall_s"] / par["wall_s"] if par["wall_s"] > 0 else 0
        efficiency = speedup / workers * 100
        print(
            f"  {spec.name:<38}  "
            f"{seq['wall_s']:>8.2f}s  "
            f"{par['wall_s']:>8.2f}s  "
            f"{speedup:>8.2f}×  "
            f"{efficiency:>9.0f}%"
        )

    # ── summary: RAM comparison ───────────────────────────────────────────
    print("\n\n── Peak RAM Comparison (tracemalloc) ────────────────────────────")
    _print_separator("·")
    print(
        f"  {'Dataset':<38}  {'scan_fasta':>12}  "
        f"{'iter':>12}  {'iter saving':>12}"
    )
    _print_separator("·")
    for spec, path, results in completed:
        seq  = results.get("scan_fasta")
        it   = results.get("scan_fasta_iter")
        if not seq or not it:
            continue
        saving = seq["peak_ram_mb"] - it["peak_ram_mb"]
        print(
            f"  {spec.name:<38}  "
            f"{seq['peak_ram_mb']:>10.1f}MB  "
            f"{it['peak_ram_mb']:>10.1f}MB  "
            f"{saving:>+10.1f}MB"
        )

    # ── summary: disk I/O ─────────────────────────────────────────────────
    if not args.no_cli:
        print("\n\n── CLI Disk I/O ─────────────────────────────────────────────────")
        _print_separator("·")
        print(
            f"  {'Dataset':<38}  {'FASTA MB':>9}  "
            f"{'TSV MB':>9}  {'ratio':>7}  {'FASTA MB/s':>11}"
        )
        _print_separator("·")
        for spec, path, results in completed:
            c = results.get("cli")
            if not c or c.get("returncode", 1) != 0:
                continue
            print(
                f"  {spec.name:<38}  "
                f"{path.stat().st_size / 1e6:>8.1f}  "
                f"{c['output_mb']:>8.1f}  "
                f"{c['compression_ratio']:>6.1f}×  "
                f"{c['throughput_mbps']:>10.1f}"
            )

    # ── JSON output ───────────────────────────────────────────────────────
    if args.json:
        out = Path(args.json)
        out.write_text(json.dumps(all_output, indent=2, default=str))
        print(f"\n  Results saved → {out}")

    _print_separator("=")
    print("  Done.")
    _print_separator("=")


if __name__ == "__main__":
    main()
