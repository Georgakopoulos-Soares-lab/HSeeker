#!/usr/bin/env python3
"""Run HSeeker on hg38 chr1 with scoring."""

from __future__ import annotations

import csv
import json
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "benchmarks"))

import hseeker  # noqa: E402


OUT = ROOT / "analysis" / "results" / "chr1"
CHR1 = ROOT / "benchmarks" / "data" / "chr1.fa"
PARAMS = dict(minrep=10, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10)


def ensure_chr1() -> Path:
    if CHR1.exists():
        return CHR1
    import benchmark  # noqa: E402
    path = benchmark._download_chromosome("chr1")
    if path is None or not path.exists():
        raise RuntimeError("Could not obtain chr1 FASTA")
    return path


def write_tsv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    keys = [
        "seq_id", "source", "start", "end", "arm_length", "spacer_length",
        "total_length", "ga_pct", "ct_pct", "mirror_identity", "is_perfect",
        "left_arm", "spacer", "right_arm", "full_sequence",
        "stacking_score", "pairing_score", "total_score", "putative_triplex",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def run_hseeker_chr1(fasta: Path) -> tuple[list[dict[str, Any]], float]:
    t0 = time.perf_counter()
    hits = hseeker.scan_fasta_parallel(
        fasta,
        workers=16,
        chunk_size=1_000_000,
        score=True,
        purity_rmq=True,
        **PARAMS,
    )
    elapsed = time.perf_counter() - t0
    for h in hits:
        h["source"] = "hseeker_purity_rmq"
    hits.sort(key=lambda h: (h["seq_id"], h["start"]))
    return hits, elapsed


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fasta = ensure_chr1()
    h_hits, h_sec = run_hseeker_chr1(fasta)
    h_out = OUT / "chr1_hseeker_16threads_scored.tsv"
    write_tsv(h_out, h_hits)

    summary = {
        "fasta": str(fasta),
        "hseeker_output": str(h_out),
        "hseeker_hits": len(h_hits),
        "hseeker_runtime_sec": h_sec,
        "workers": 16,
        "scoring": True,
        "purity_rmq": True,
        "params": PARAMS,
        "triplex_status": "not_run_for_full_chr1",
        "triplex_note": "Full-chromosome Triplex was not included in this chr1 run.",
    }
    (OUT / "chr1_run_summary.json").write_text(json.dumps(summary, indent=2))
    (OUT / "summary.md").write_text(
        "# chr1 hg38 HSeeker Summary\n\n"
        f"FASTA: `{fasta}`\n\n"
        f"HSeeker output: `{h_out}`\n\n"
        f"HSeeker hits: {len(h_hits)}\n\n"
        "HSeeker settings: 16 worker threads, scoring enabled, "
        "`purity_rmq=True`, minrep=10, maxrep=1000, maxspacer=10, "
        "purity=0.90, mismatch=0.10.\n\n"
        f"HSeeker runtime: {h_sec:.2f} seconds.\n\n"
        "Triplex status: not run for full chr1 in this analysis.\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
