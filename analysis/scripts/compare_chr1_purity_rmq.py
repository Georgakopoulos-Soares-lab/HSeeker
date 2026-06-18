#!/usr/bin/env python3
"""Compare full chr1 HSeeker runtime with purity_rmq off vs on."""

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


def hit_signature(hit: dict[str, Any]) -> tuple[Any, ...]:
    return (
        hit.get("seq_id"),
        hit["start"],
        hit["end"],
        hit["arm_length"],
        hit["spacer_length"],
        round(float(hit["ga_pct"]), 4),
        round(float(hit["ct_pct"]), 4),
        round(float(hit["mirror_identity"]), 4),
    )


def score_summary(hits: list[dict[str, Any]]) -> dict[str, float | None]:
    scores = sorted(
        float(h["total_score"])
        for h in hits
        if h.get("total_score") not in (None, "", "None")
    )
    if not scores:
        return {"score_min": None, "score_median": None, "score_max": None}
    return {
        "score_min": scores[0],
        "score_median": scores[len(scores) // 2],
        "score_max": scores[-1],
    }


def run_mode(fasta: Path, purity_rmq: bool) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    start = time.perf_counter()
    hits = hseeker.scan_fasta_parallel(
        fasta,
        workers=16,
        chunk_size=1_000_000,
        score=True,
        purity_rmq=purity_rmq,
        **PARAMS,
    )
    elapsed = time.perf_counter() - start
    hits.sort(key=lambda h: (h["seq_id"], h["start"], h["end"], h["arm_length"], h["spacer_length"]))
    summary = {
        "purity_rmq": purity_rmq,
        "runtime_sec": elapsed,
        "hits": len(hits),
        **score_summary(hits),
    }
    return hits, summary


def write_tsv(path: Path, hits: list[dict[str, Any]]) -> None:
    fieldnames = [
        "seq_id", "source", "start", "end",
        "arm_length", "spacer_length", "total_length",
        "ga_pct", "ct_pct", "mirror_identity", "is_perfect",
        "left_arm", "spacer", "right_arm", "full_sequence",
        "stacking_score", "pairing_score", "total_score",
        "putative_triplex",
    ]
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(hits)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fasta = ensure_chr1()

    hits_false, summary_false = run_mode(fasta, False)
    for h in hits_false:
        h["source"] = "hseeker"
    write_tsv(OUT / "chr1_hseeker_16threads_scored_purityrmq_false.tsv", hits_false)

    hits_true, summary_true = run_mode(fasta, True)
    for h in hits_true:
        h["source"] = "hseeker_purity_rmq"
    write_tsv(OUT / "chr1_hseeker_16threads_scored_purityrmq_true.tsv", hits_true)

    sig_false = {hit_signature(h) for h in hits_false}
    sig_true = {hit_signature(h) for h in hits_true}
    comparison = {
        "fasta": str(fasta),
        "workers": 16,
        "scoring": True,
        "params": PARAMS,
        "purity_rmq_false": summary_false,
        "purity_rmq_true": summary_true,
        "speedup_true_vs_false": (
            summary_false["runtime_sec"] / summary_true["runtime_sec"]
            if summary_true["runtime_sec"] > 0
            else None
        ),
        "lost_hits_with_rmq": len(sig_false - sig_true),
        "extra_hits_with_rmq": len(sig_true - sig_false),
        "outputs": {
            "purity_rmq_false_tsv": str(OUT / "chr1_hseeker_16threads_scored_purityrmq_false.tsv"),
            "purity_rmq_true_tsv": str(OUT / "chr1_hseeker_16threads_scored_purityrmq_true.tsv"),
        },
    }

    (OUT / "chr1_purity_rmq_runtime_comparison.json").write_text(
        json.dumps(comparison, indent=2) + "\n"
    )
    (OUT / "chr1_purity_rmq_runtime_comparison.md").write_text(
        "# chr1 HSeeker purity_rmq Runtime Comparison\n\n"
        "Settings: 16 worker threads, scoring enabled, minrep=10, maxrep=1000, "
        "maxspacer=10, purity=0.90, mismatch=0.10.\n\n"
        f"- `purity_rmq=False`: {summary_false['runtime_sec']:.2f}s, "
        f"{summary_false['hits']:,} hits\n"
        f"- `purity_rmq=True`: {summary_true['runtime_sec']:.2f}s, "
        f"{summary_true['hits']:,} hits\n"
        f"- Speedup (`True` vs `False`): {comparison['speedup_true_vs_false']:.2f}x\n"
        f"- Lost hits with RMQ: {comparison['lost_hits_with_rmq']}\n"
        f"- Extra hits with RMQ: {comparison['extra_hits_with_rmq']}\n"
    )
    print(json.dumps(comparison, indent=2))


if __name__ == "__main__":
    main()
