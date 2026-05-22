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
import gzip
import hashlib
import json
import datetime
import os
import platform
import random
import shutil
import subprocess
import sys
import tempfile
import time
import tracemalloc
import urllib.error
import urllib.request
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

DATA_DIR          = Path(__file__).parent / "data"
ZENODO_RECORD_FILE = Path(__file__).parent / "zenodo_record.json"

# hg38 chromosome 16 from NCBI FTP (GRCh38 / GCA_000001405.15, ~90 MB uncompressed)
CHR16_FASTA_URL = (
    "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/"
    "GCA_000001405.15_GRCh38/GCA_000001405.15_GRCh38_assembly_structure/"
    "Primary_Assembly/assembled_chromosomes/FASTA/chr16.fna.gz"
)
CHR16_LOCAL_NAME = "hg38_chr16.fa"

LINE_WIDTH = 60          # bases per FASTA line
GEN_CHUNK  = 10_000_000  # bases generated per numpy call (≈10 MB in RAM)

# ACGT weights: (A, C, G, T)
_PROFILES: dict[str, tuple[float, float, float, float]] = {
    "uniform":   (0.25,  0.25,  0.25,  0.25),
    "ga_biased": (0.05,  0.05,  0.45,  0.45),  # 90 % purine    → dense GA-mirror hits
    "ct_biased": (0.05,  0.45,  0.05,  0.45),  # 90 % pyrimidine → dense CT-mirror hits
    "realistic": (0.295, 0.205, 0.205, 0.295),  # ~41 % GC, like human chr
}

# minrep values used in the optional parameter-sensitivity sweep
MINREP_SWEEP_VALUES: list[int] = [6, 8, 10, 15, 20]

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


# ---------------------------------------------------------------------------
# Seed manifest — tracks which seed was used per dataset so that a seed
# change triggers automatic regeneration of the affected file.
# ---------------------------------------------------------------------------

_SEED_MANIFEST_PATH = DATA_DIR / "SEED_MANIFEST.json"


def _load_seed_manifest() -> dict[str, int]:
    """Load {dataset_name: seed} mapping from disk; return {} on any error."""
    if _SEED_MANIFEST_PATH.exists():
        try:
            return json.loads(_SEED_MANIFEST_PATH.read_text())
        except Exception:
            return {}
    return {}


def _save_seed_manifest(manifest: dict[str, int]) -> None:
    DATA_DIR.mkdir(parents=True, exist_ok=True)
    _SEED_MANIFEST_PATH.write_text(json.dumps(manifest, indent=2, sort_keys=True))


# ---------------------------------------------------------------------------
# Zenodo download helpers
# ---------------------------------------------------------------------------


def _read_zenodo_record() -> dict | None:
    """Load benchmarks/zenodo_record.json; return None if not populated."""
    if not ZENODO_RECORD_FILE.exists():
        return None
    try:
        data = json.loads(ZENODO_RECORD_FILE.read_text())
    except Exception:
        return None
    if not data.get("record_id"):
        return None
    return data


def _verify_md5(path: Path, expected: str) -> bool:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest() == expected


def _download_from_zenodo(name: str, record: dict) -> Path | None:
    """
    Download <name>.fa.gz from a published Zenodo record and decompress to
    DATA_DIR/<name>.fa.  Returns the local .fa path on success, None if the
    file is not listed in the record manifest.
    """
    remote_name = f"{name}.fa.gz"
    entry = next((f for f in record.get("files", []) if f["name"] == remote_name), None)
    if entry is None:
        return None

    base_url  = record.get("base_url", "https://zenodo.org").rstrip("/")
    record_id = record["record_id"]
    dl_url    = f"{base_url}/records/{record_id}/files/{remote_name}?download=1"

    gz_path  = DATA_DIR / remote_name
    fa_path  = DATA_DIR / f"{name}.fa"
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    size_mb = entry.get("size_bytes", 0) / 1e6
    print(
        f"  downloading {remote_name} ({size_mb:.1f} MB) from Zenodo … ",
        end="", flush=True,
    )
    t0 = time.perf_counter()
    try:
        req = urllib.request.Request(
            dl_url,
            headers={"User-Agent": f"hseeker-benchmark/1.0"},
        )
        with urllib.request.urlopen(req) as resp, open(gz_path, "wb") as out:
            shutil.copyfileobj(resp, out)
    except Exception as exc:
        print(f"FAILED ({exc})")
        gz_path.unlink(missing_ok=True)
        return None

    elapsed = time.perf_counter() - t0
    print(f"{elapsed:.1f}s")

    # Verify MD5 if available
    expected_md5 = entry.get("md5")
    if expected_md5 and not _verify_md5(gz_path, expected_md5):
        print(f"  ERROR: MD5 mismatch for {remote_name} — file may be corrupt.")
        gz_path.unlink(missing_ok=True)
        return None

    # Decompress
    print(f"  decompressing {remote_name} … ", end="", flush=True)
    t1 = time.perf_counter()
    with gzip.open(gz_path, "rb") as f_in, open(fa_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    gz_path.unlink()   # remove .gz; keep only .fa
    size_fa = fa_path.stat().st_size / 1e6
    print(f"{size_fa:.1f} MB ({time.perf_counter() - t1:.1f}s)")
    return fa_path


def _download_chr16(force: bool = False) -> "Path | None":
    """Download hg38 chr16 FASTA from UCSC and decompress to DATA_DIR/hg38_chr16.fa.

    Returns the local .fa path on success, None on failure.  The .gz file is
    removed after decompression; only the .fa is kept as the local cache.
    """
    fa_path = DATA_DIR / CHR16_LOCAL_NAME
    if fa_path.exists() and not force:
        print(f"  cached: {fa_path.name}  ({fa_path.stat().st_size / 1e6:.1f} MB)")
        return fa_path

    gz_path = DATA_DIR / (CHR16_LOCAL_NAME + ".gz")
    DATA_DIR.mkdir(parents=True, exist_ok=True)

    print("  downloading hg38 chr16 from NCBI (~30 MB) … ", end="", flush=True)
    t0 = time.perf_counter()
    try:
        req = urllib.request.Request(
            CHR16_FASTA_URL,
            headers={"User-Agent": "hseeker-benchmark/1.0"},
        )
        with urllib.request.urlopen(req, timeout=300) as resp, \
                open(gz_path, "wb") as out:
            shutil.copyfileobj(resp, out)
    except Exception as exc:
        print(f"FAILED ({exc})")
        gz_path.unlink(missing_ok=True)
        return None
    print(f"{time.perf_counter() - t0:.1f}s")

    print("  decompressing … ", end="", flush=True)
    t1 = time.perf_counter()
    with gzip.open(gz_path, "rb") as f_in, open(fa_path, "wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    gz_path.unlink()
    print(f"{fa_path.stat().st_size / 1e6:.1f} MB ({time.perf_counter() - t1:.1f}s)")
    return fa_path


# ---------------------------------------------------------------------------
# Dataset caching with optional Zenodo fallback
# ---------------------------------------------------------------------------

def _ensure_dataset(
    name: str,
    n_records: int,
    record_length: int,
    profile: str,
    seed: int,
    force: bool = False,
    from_zenodo: bool = False,
) -> Path:
    """Return path to cached FASTA file.

    Resolution order:
    1. Local cache (benchmarks/data/<name>.fa) — reused if seed matches.
    2. If --from-zenodo: download from the Zenodo record listed in
       benchmarks/zenodo_record.json.
    3. Generate locally (always available as last resort).
    """
    path = DATA_DIR / f"{name}.fa"
    manifest = _load_seed_manifest()
    seed_changed = manifest.get(name) != seed

    if path.exists() and not force and not seed_changed:
        return path

    if seed_changed and path.exists() and not force:
        print(
            f"  seed changed for {name} "
            f"(was {manifest.get(name)!r} → {seed}) — regenerating",
            flush=True,
        )

    # ── Try Zenodo download ────────────────────────────────────────────────
    if from_zenodo and not force:
        record = _read_zenodo_record()
        if record is None:
            print(
                f"  WARNING: benchmarks/zenodo_record.json has no record_id; "
                f"falling back to local generation for {name}."
            )
        else:
            dl = _download_from_zenodo(name, record)
            if dl is not None:
                manifest[name] = seed
                _save_seed_manifest(manifest)
                return dl
            print(f"  Zenodo download failed for {name}; falling back to local generation.")

    # ── Generate locally ───────────────────────────────────────────────────
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
    manifest[name] = seed
    _save_seed_manifest(manifest)
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
# Seeds are stable — changing a seed triggers automatic regeneration via
# SEED_MANIFEST.json.   Small: 1xxx  |  Medium: 2xxx  |  Large: 3xxx
SMALL_DATASETS: list[DatasetSpec] = [
    # 5 records × 6 MB ≈ 30 MB — baseline throughput per profile
    DatasetSpec("small_uniform_5rec",      5,  6_000_000, "uniform",   seed=1001),
    # Dense GA purine content — stresses hit buffer and overlap removal
    DatasetSpec("small_ga_biased_5rec",    5,  6_000_000, "ga_biased", seed=1002),
    # Dense CT pyrimidine content — mirrors ga_biased for CT-mirror hits
    DatasetSpec("small_ct_biased_5rec",    5,  6_000_000, "ct_biased", seed=1003),
    # 24 realistic records × 1.25 MB — tests multi-record parallelism benefit
    DatasetSpec("small_realistic_24rec",  24,  1_250_000, "realistic", seed=1004),
    # Single long record (30 MB) — tests per-record latency (no parallel benefit)
    DatasetSpec("small_realistic_1rec",    1, 30_000_000, "realistic", seed=1005),
]

MEDIUM_DATASETS: list[DatasetSpec] = [
    DatasetSpec("medium_uniform_5rec",     5,  60_000_000, "uniform",   seed=2001),
    DatasetSpec("medium_ga_biased_5rec",   5,  60_000_000, "ga_biased", seed=2002),
    DatasetSpec("medium_ct_biased_5rec",   5,  60_000_000, "ct_biased", seed=2003),
    DatasetSpec("medium_realistic_24rec", 24,  12_500_000, "realistic", seed=2004),
    DatasetSpec("medium_realistic_1rec",   1, 300_000_000, "realistic", seed=2005),
]

LARGE_DATASETS: list[DatasetSpec] = [
    DatasetSpec("large_uniform_5rec",      5,  600_000_000, "uniform",   seed=3001),
    DatasetSpec("large_ga_biased_5rec",    5,  600_000_000, "ga_biased", seed=3002),
    DatasetSpec("large_ct_biased_5rec",    5,  600_000_000, "ct_biased", seed=3003),
    DatasetSpec("large_realistic_24rec",  24,  125_000_000, "realistic", seed=3004),
    DatasetSpec("large_realistic_1rec",    1, 3_000_000_000, "realistic", seed=3005),
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
# System info
# ---------------------------------------------------------------------------

def _get_system_info() -> dict[str, Any]:
    """Collect hardware/software facts for the benchmark report header."""
    info: dict[str, Any] = {
        "os":               platform.platform(),
        "python":           sys.version.split()[0],
        "hseeker":          hseeker.__version__,
        "numpy":            np.__version__,
        "cpu_cores_logical": os.cpu_count(),
    }

    # CPU brand string (macOS → sysctl; Linux → /proc/cpuinfo)
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
            info["disk_free_gb"] = round(psutil.disk_usage(str(DATA_DIR.parent)).free / 1e9, 1)
        except Exception:
            pass

    return info


# ---------------------------------------------------------------------------
# minrep sensitivity sweep
# ---------------------------------------------------------------------------

def _run_minrep_sweep(
    spec: DatasetSpec,
    path: Path,
    minrep_values: list[int],
) -> list[dict[str, Any]]:
    """Run scan_fasta with multiple minrep values; return timing + hit count rows."""
    file_mb = path.stat().st_size / 1e6
    rows: list[dict[str, Any]] = []
    for mr in minrep_values:
        print(f"    minrep={mr:>2} … ", end="", flush=True)
        with _measure() as r:
            hits = hseeker.scan_fasta(str(path), minrep=mr)
        r["minrep"] = mr
        r["hits"]   = len(hits)
        r["throughput_mbps"] = file_mb / r["wall_s"] if r["wall_s"] > 0 else 0.0
        rows.append(r)
        print(f"{r['wall_s']:.2f}s  {r['hits']:,} hits")
    return rows


# ---------------------------------------------------------------------------
# Markdown report generation
# ---------------------------------------------------------------------------

def _generate_markdown_report(
    sys_info:      dict[str, Any],
    all_output:    list[dict[str, Any]],
    completed:     list[tuple[DatasetSpec, Path, dict]],
    workers:       int,
    minrep:        int,
    minrep_sweep:  dict[str, list[dict[str, Any]]] | None = None,
) -> str:
    lines: list[str] = []
    now = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # ── Title ────────────────────────────────────────────────────────────
    lines += [
        "# HSeeker Performance Benchmark Report",
        "",
        f"**Date**: {now}  ",
        f"**hseeker version**: `{sys_info.get('hseeker', 'unknown')}`  ",
        f"**Workers**: {workers}  ",
        f"**minrep**: {minrep}  ",
        "",
    ]

    # ── System specs ──────────────────────────────────────────────────────
    lines += ["## System Specifications", "", "| Property | Value |", "|---|---|"]
    if "cpu_brand" in sys_info:
        lines.append(f"| CPU | {sys_info['cpu_brand']} |")
    phys = sys_info.get("cpu_cores_physical", "—")
    logi = sys_info.get("cpu_cores_logical",  "—")
    lines.append(f"| CPU cores | {phys} physical / {logi} logical |")
    if sys_info.get("cpu_freq_max_ghz"):
        lines.append(f"| CPU max clock | {sys_info['cpu_freq_max_ghz']} GHz |")
    if "ram_total_gb" in sys_info:
        lines.append(f"| RAM (total) | {sys_info['ram_total_gb']} GB |")
    if "ram_available_gb" in sys_info:
        lines.append(f"| RAM (available) | {sys_info['ram_available_gb']} GB |")
    lines.append(f"| OS | {sys_info.get('os', '—')} |")
    lines.append(f"| Python | {sys_info.get('python', '—')} |")
    lines.append(f"| NumPy | {sys_info.get('numpy', '—')} |")
    if "disk_free_gb" in sys_info:
        lines.append(f"| Disk free (benchmark dir) | {sys_info['disk_free_gb']} GB |")
    lines.append("")

    # ── Methodology ───────────────────────────────────────────────────────
    lines += [
        "## Methodology",
        "",
        "Synthetic FASTA files are generated deterministically with a fixed NumPy/Python "
        "random seed stored in `benchmarks/data/SEED_MANIFEST.json`. "
        "Re-running with `--no-cache` regenerates the same files bit-for-bit.",
        "",
        "### Sequence profiles",
        "",
        "| Profile | Composition | Purpose |",
        "|---|---|---|",
        "| `uniform` | 25% each ACGT | Sparse H-DNA; baseline throughput |",
        "| `ga_biased` | 45% G + 45% A (90% purine) | Dense GA-mirror hits; stresses hit buffer |",
        "| `ct_biased` | 45% C + 45% T (90% pyrimidine) | Dense CT-mirror hits; mirror of ga_biased |",
        "| `realistic` | ~41% GC, AT-biased | Human-chromosome-like density and composition |",
        "",
        "### Record layout scenarios",
        "",
        "| Scenario | Records | Record length | Purpose |",
        "|---|---|---|---|",
        "| `5rec` | 5 | 6 MB each | Multi-record parallelism baseline |",
        "| `24rec` | 24 | 1.25 MB each | Maximum parallelism benefit |",
        "| `1rec` | 1 | 30 MB | Per-record latency; validates chunk-level parallelism |",
        "",
        "### API paths benchmarked",
        "",
        "| Method | Description |",
        "|---|---|",
        "| `scan_fasta` | Sequential scan returning a full list |",
        "| `scan_fasta_iter` | Sequential streaming generator (lower peak RAM) |",
        "| `parallel_1w` | `scan_fasta_parallel(workers=1)` — overhead baseline |",
        f"| `parallel_{workers}w` | `scan_fasta_parallel(workers={workers})` — full parallelism |",
        "| `cli` | End-to-end `python -m hseeker` including TSV disk write |",
        "",
        "Timing: `time.perf_counter()` (wall), `time.process_time()` (CPU).  ",
        "Peak RAM: `tracemalloc` new allocations; RSS delta via `psutil`.  ",
        "",
    ]

    # ── Dataset summary ───────────────────────────────────────────────────
    lines += [
        "## Dataset Summary",
        "",
        "| Dataset | Profile | Records | Record length | Disk size | Seed |",
        "|---|---|---|---|---|---|",
    ]
    for spec, path, _ in completed:
        size_mb = path.stat().st_size / 1e6
        lines.append(
            f"| `{spec.name}` | {spec.profile} | {spec.n_records} | "
            f"{spec.record_length / 1e6:.1f} MB | {size_mb:.1f} MB | {spec.seed} |"
        )
    lines += [
        "",
        "> Seeds are fixed per dataset (see `SEED_MANIFEST.json`). "
        "Re-run with `--no-cache` to reproduce the exact same files.",
        "",
    ]

    # ── Per-dataset results ───────────────────────────────────────────────
    lines += ["## Benchmark Results by Dataset", ""]
    for spec, path, results in completed:
        file_mb = path.stat().st_size / 1e6
        lines += [
            f"### {spec.name}",
            "",
            f"**Profile**: `{spec.profile}` &nbsp;|&nbsp; "
            f"**Records**: {spec.n_records} × {spec.record_length / 1e6:.1f} MB &nbsp;|&nbsp; "
            f"**On disk**: {file_mb:.1f} MB",
            "",
            "| Method | Wall (s) | CPU (s) | Peak RAM (MB) | Hits | Throughput (MB/s) |",
            "|---|---|---|---|---|---|",
        ]
        for method, r in results.items():
            if "wall_s" not in r:
                continue
            tp = file_mb / r["wall_s"] if r["wall_s"] > 0 else 0.0
            hits_s = f"{r['hits']:,}" if isinstance(r.get("hits"), int) else "—"
            lines.append(
                f"| `{method}` | {r['wall_s']:.3f} | {r['cpu_s']:.3f} | "
                f"{r['peak_ram_mb']:.1f} | {hits_s} | {tp:.1f} |"
            )
        if "cli" in results:
            c = results["cli"]
            if c.get("returncode") == 0:
                lines.append(
                    f"| `cli` | {c['wall_s']:.3f} | — | — | — | "
                    f"{c['throughput_mbps']:.1f} |"
                )
        lines.append("")

    # ── Parallelism analysis ──────────────────────────────────────────────
    lines += [
        "## Parallelism Analysis",
        "",
        f"Comparing sequential `scan_fasta` vs `scan_fasta_parallel(workers={workers})`.",
        "",
        f"| Dataset | Sequential (s) | Parallel {workers}w (s) | Speedup | Efficiency |",
        "|---|---|---|---|---|",
    ]
    for spec, path, results in completed:
        seq = results.get("scan_fasta")
        par = results.get(f"parallel_{workers}w")
        if not seq or not par or par["wall_s"] <= 0:
            continue
        speedup    = seq["wall_s"] / par["wall_s"]
        efficiency = speedup / workers * 100
        lines.append(
            f"| `{spec.name}` | {seq['wall_s']:.3f} | {par['wall_s']:.3f} | "
            f"{speedup:.2f}× | {efficiency:.0f}% |"
        )
    lines.append("")

    # ── Memory analysis ───────────────────────────────────────────────────
    lines += [
        "## Memory Efficiency: `scan_fasta` vs `scan_fasta_iter`",
        "",
        "| Dataset | scan_fasta peak (MB) | scan_fasta_iter peak (MB) | Saving (MB) |",
        "|---|---|---|---|",
    ]
    for spec, path, results in completed:
        seq = results.get("scan_fasta")
        it  = results.get("scan_fasta_iter")
        if not seq or not it:
            continue
        saving = seq["peak_ram_mb"] - it["peak_ram_mb"]
        lines.append(
            f"| `{spec.name}` | {seq['peak_ram_mb']:.1f} | "
            f"{it['peak_ram_mb']:.1f} | {saving:+.1f} |"
        )
    lines.append("")

    # ── CLI disk I/O ──────────────────────────────────────────────────────
    cli_rows = [
        (s, p, r) for s, p, r in completed
        if "cli" in r and r["cli"].get("returncode") == 0
    ]
    if cli_rows:
        lines += [
            "## CLI Throughput (FASTA → TSV)",
            "",
            "| Dataset | FASTA (MB) | TSV output (MB) | FASTA/TSV ratio | Throughput (MB/s) |",
            "|---|---|---|---|---|",
        ]
        for spec, path, results in cli_rows:
            c = results["cli"]
            lines.append(
                f"| `{spec.name}` | {path.stat().st_size / 1e6:.1f} | "
                f"{c['output_mb']:.1f} | {c['compression_ratio']:.1f}× | "
                f"{c['throughput_mbps']:.1f} |"
            )
        lines.append("")

    # ── minrep sweep ──────────────────────────────────────────────────────
    if minrep_sweep:
        lines += [
            "## Parameter Sensitivity: `minrep` Sweep",
            "",
            "Effect of minimum arm length (`minrep`) on hit count and scan throughput.",
            "",
        ]
        for dataset_name, rows in minrep_sweep.items():
            # find the file size for throughput
            file_mb_map = {s.name: p.stat().st_size / 1e6 for s, p, _ in completed}
            fm = file_mb_map.get(dataset_name, 0.0)
            lines += [
                f"### {dataset_name}",
                "",
                "| minrep | Hits | Wall (s) | Throughput (MB/s) |",
                "|---|---|---|---|",
            ]
            for row in rows:
                lines.append(
                    f"| {row['minrep']} | {row['hits']:,} | "
                    f"{row['wall_s']:.3f} | {row['throughput_mbps']:.1f} |"
                )
            lines.append("")

    # ── Profile comparison ────────────────────────────────────────────────
    lines += [
        "## Profile Comparison (Sequential Throughput)",
        "",
        "Throughput (MB/s) of `scan_fasta` across sequence profiles for 5-record small datasets.",
        "",
        "| Profile | Wall (s) | Hits | Throughput (MB/s) |",
        "|---|---|---|---|",
    ]
    profile_order = ["small_uniform_5rec", "small_ga_biased_5rec",
                     "small_ct_biased_5rec", "small_realistic_24rec"]
    for spec, path, results in completed:
        if spec.name not in profile_order:
            continue
        r = results.get("scan_fasta")
        if not r:
            continue
        tp = path.stat().st_size / 1e6 / r["wall_s"] if r["wall_s"] > 0 else 0.0
        lines.append(
            f"| `{spec.profile}` ({spec.name}) | {r['wall_s']:.3f} | "
            f"{r.get('hits', '—'):,} | {tp:.1f} |"
        )
    lines.append("")

    # ── Key observations ──────────────────────────────────────────────────
    lines += ["## Key Observations", ""]

    # Throughput range
    throughputs = {
        s.name: p.stat().st_size / 1e6 / r["scan_fasta"]["wall_s"]
        for s, p, r in completed
        if "scan_fasta" in r and r["scan_fasta"]["wall_s"] > 0
    }
    if throughputs:
        hi = max(throughputs, key=throughputs.__getitem__)
        lo = min(throughputs, key=throughputs.__getitem__)
        lines += [
            f"- **Highest sequential throughput**: `{hi}` at "
            f"**{throughputs[hi]:.1f} MB/s** "
            f"(sparse hits → minimal post-processing).",
            f"- **Lowest sequential throughput**: `{lo}` at "
            f"**{throughputs[lo]:.1f} MB/s** "
            f"(dense hits require more overlap-removal work).",
        ]

    # Parallelism best case
    speedups = {}
    for spec, path, results in completed:
        seq = results.get("scan_fasta")
        par = results.get(f"parallel_{workers}w")
        if seq and par and par["wall_s"] > 0:
            speedups[spec.name] = seq["wall_s"] / par["wall_s"]
    if speedups:
        best = max(speedups, key=speedups.__getitem__)
        lines.append(
            f"- **Best parallelism speedup**: `{best}` at "
            f"**{speedups[best]:.2f}×** with {workers} workers."
        )
        # Single-record case
        single = [n for n in speedups if "1rec" in n]
        if single:
            lines.append(
                f"- **Single-record parallelism**: `{single[0]}` shows "
                f"**{speedups[single[0]]:.2f}×** speedup via chunk-level parallel scan "
                f"within one record."
            )

    # Memory saving
    for spec, path, results in completed:
        seq = results.get("scan_fasta")
        it  = results.get("scan_fasta_iter")
        if seq and it:
            saving = seq["peak_ram_mb"] - it["peak_ram_mb"]
            if saving > 0.5:
                lines.append(
                    f"- **`scan_fasta_iter` saves {saving:.1f} MB** peak RAM vs "
                    f"`scan_fasta` on `{spec.name}`."
                )
                break

    # Hit density comment
    ga_hits = next(
        (r["scan_fasta"].get("hits", 0) for s, p, r in completed
         if "ga_biased" in s.name and "5rec" in s.name and "scan_fasta" in r), 0
    )
    uni_hits = next(
        (r["scan_fasta"].get("hits", 0) for s, p, r in completed
         if "uniform" in s.name and "5rec" in s.name and "scan_fasta" in r), 0
    )
    if ga_hits and uni_hits:
        ratio = ga_hits / uni_hits if uni_hits else 0
        lines.append(
            f"- **GA-biased sequences produce ~{ratio:.0f}× more hits** than uniform "
            f"sequences ({ga_hits:,} vs {uni_hits:,} hits at minrep={minrep}), "
            f"directly impacting scan time."
        )

    lines += [
        "",
        "---",
        f"*Generated by HSeeker benchmark suite — {now}*",
        "",
    ]
    return "\n".join(lines)


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
        "--chr16", action="store_true",
        help=(
            "Download hg38 chromosome 16 from UCSC (~90 MB) and run the full "
            "benchmark suite on it as a real-genome validation dataset. "
            "The file is cached in benchmarks/data/hg38_chr16.fa."
        ),
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
    p.add_argument(
        "--minrep-sweep", action="store_true", dest="minrep_sweep",
        help="Run minrep sensitivity sweep on select datasets",
    )
    p.add_argument(
        "--minrep-values", type=int, nargs="+", metavar="N",
        dest="minrep_values", default=None,
        help="minrep values for the sweep (default: 6 8 10 15 20)",
    )
    p.add_argument(
        "--report", metavar="FILE",
        help="Write a detailed Markdown performance report to FILE",
    )
    p.add_argument(
        "--from-zenodo", action="store_true", dest="from_zenodo",
        help=(
            "Download benchmark FASTA files from Zenodo instead of generating "
            "them locally. Requires benchmarks/zenodo_record.json to contain a "
            "valid record_id (populated by benchmarks/zenodo_upload.py)."
        ),
    )
    return p.parse_args()


def _run_dataset(
    spec: DatasetSpec,
    workers: int,
    minrep: int,
    force: bool,
    cli: bool,
    from_zenodo: bool = False,
) -> tuple[Path, dict[str, dict[str, Any]]]:
    path    = _ensure_dataset(
        spec.name, spec.n_records, spec.record_length,
        spec.profile, spec.seed, force=force, from_zenodo=from_zenodo,
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
    sys_info = _get_system_info()
    _print_separator("=")
    print("  HSeeker Performance Benchmark".center(_W))
    print(datetime.datetime.now().strftime("  %Y-%m-%d %H:%M:%S").center(_W))
    _print_separator("=")
    print(f"  hseeker    : {sys_info['hseeker']}")
    print(f"  Python     : {sys_info['python']}")
    if "cpu_brand" in sys_info:
        print(f"  CPU        : {sys_info['cpu_brand']}")
    phys = sys_info.get("cpu_cores_physical", "?")
    logi = sys_info.get("cpu_cores_logical",  "?")
    print(f"  CPU cores  : {logi} logical / {phys} physical")
    if "cpu_freq_max_ghz" in sys_info:
        print(f"  CPU freq   : {sys_info['cpu_freq_max_ghz']} GHz (max)")
    print(f"  Workers    : {workers}")
    print(f"  minrep     : {args.minrep}")
    if "ram_total_gb" in sys_info:
        print(
            f"  RAM        : {sys_info['ram_total_gb']} GB total, "
            f"{sys_info.get('ram_available_gb', '?')} GB available"
        )
    else:
        print("  RAM        : install psutil for memory measurements")
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
            spec, workers, args.minrep, args.no_cache, not args.no_cli,
            from_zenodo=args.from_zenodo,
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

    # ── optional: hg38 chr16 ─────────────────────────────────────────────
    if args.chr16:
        print("\n── Downloading hg38 chr16 ──────────────────────────────────────")
        chr16_path = _download_chr16(force=args.no_cache)
        if chr16_path:
            chr16_mb = chr16_path.stat().st_size / 1e6
            print(f"\n[hg38_chr16]  {chr16_mb:.1f} MB")
            for label, fn in [
                ("scan_fasta",
                 lambda: bench_scan_fasta(chr16_path, args.minrep)),
                ("scan_fasta_iter",
                 lambda: bench_scan_fasta_iter(chr16_path, args.minrep)),
                ("parallel_1w",
                 lambda: bench_scan_fasta_parallel(chr16_path, 1, args.minrep)),
                (f"parallel_{workers}w",
                 lambda: bench_scan_fasta_parallel(chr16_path, workers, args.minrep)),
            ]:
                print(f"    {label} … ", end="", flush=True)
                r = fn()
                print(f"{r['wall_s']:.2f}s  ({r.get('hits', '?'):,} hits)")

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

    # ── minrep sweep ──────────────────────────────────────────────────────
    minrep_sweep: dict[str, list[dict[str, Any]]] = {}
    if args.minrep_sweep:
        print("\n\n── minrep Sensitivity Sweep ─────────────────────────────────────")
        sweep_target_names = {
            "small_realistic_24rec", "small_ga_biased_5rec", "small_ct_biased_5rec",
        }
        mv = args.minrep_values if args.minrep_values else MINREP_SWEEP_VALUES
        for spec, path, _results in completed:
            if spec.name in sweep_target_names:
                print(f"\n[{spec.name}]")
                minrep_sweep[spec.name] = _run_minrep_sweep(spec, path, mv)

    # ── JSON output ───────────────────────────────────────────────────────
    if args.json:
        out = Path(args.json)
        payload = {
            "system_info": sys_info,
            "workers":     workers,
            "minrep":      args.minrep,
            "datasets":    all_output,
            "minrep_sweep": {
                k: [
                    {kk: vv for kk, vv in row.items() if kk != "hits" or True}
                    for row in rows
                ]
                for k, rows in minrep_sweep.items()
            },
        }
        out.write_text(json.dumps(payload, indent=2, default=str))
        print(f"\n  Results saved → {out}")

    # ── Markdown report ───────────────────────────────────────────────────
    if args.report:
        report_path = Path(args.report)
        report_text = _generate_markdown_report(
            sys_info, all_output, completed, workers, args.minrep,
            minrep_sweep if minrep_sweep else None,
        )
        report_path.write_text(report_text, encoding="utf-8")
        print(f"\n  Markdown report → {report_path}")

    _print_separator("=")
    print("  Done.")
    _print_separator("=")


if __name__ == "__main__":
    main()
