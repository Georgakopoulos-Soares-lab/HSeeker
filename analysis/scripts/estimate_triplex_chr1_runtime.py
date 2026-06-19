#!/usr/bin/env python3
"""Estimate full hg38 chr1 Triplex runtime from smaller chr1 prefixes."""

from __future__ import annotations

import csv
import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "benchmarks"))

OUT = ROOT / "analysis" / "results" / "chr1"
TRIPLEX_R = ROOT / "analysis" / "scripts" / "triplex_validate.R"
CHR1 = ROOT / "benchmarks" / "data" / "chr1.fa"
SAMPLE_LENGTHS = (1_000, 5_000, 10_000, 25_000, 50_000, 100_000)


def ensure_chr1() -> tuple[Path, bool]:
    """Return chr1 FASTA path and whether this script created the cache."""
    if CHR1.exists():
        return CHR1, False

    import benchmark  # noqa: E402

    path = benchmark._download_chromosome("chr1")
    if path is None or not path.exists():
        raise RuntimeError("Could not obtain chr1 FASTA for Triplex runtime sampling")
    return path, True


def read_fasta_sequence(path: Path) -> str:
    chunks: list[str] = []
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                continue
            chunks.append(line.strip().upper())
    return "".join(chunks)


def write_triplex_input(seq: str, path: Path) -> list[dict[str, Any]]:
    rows = []
    for length in SAMPLE_LENGTHS:
        if length <= len(seq):
            rows.append(
                {
                    "sequence_id": f"chr1_prefix_{length}",
                    "label": "unknown",
                    "sequence": seq[:length],
                    "sample_length": length,
                }
            )

    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=["sequence_id", "label", "sequence"])
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "sequence_id": row["sequence_id"],
                    "label": row["label"],
                    "sequence": row["sequence"],
                }
            )
    return rows


def read_triplex_predictions(path: Path) -> dict[str, dict[str, str]]:
    with path.open(newline="") as fh:
        return {row["sequence_id"]: row for row in csv.DictReader(fh)}


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def estimate_runtime(rows: list[dict[str, Any]], chr1_length: int) -> dict[str, Any]:
    fit_rows = rows[-3:] if len(rows) >= 3 else rows
    denom = sum(float(r["sample_length"]) ** 2 for r in fit_rows)
    numer = sum(float(r["sample_length"]) * float(r["triplex_runtime_sec"]) for r in fit_rows)
    sec_per_base = numer / denom if denom else 0.0
    predicted_sec = sec_per_base * chr1_length
    return {
        "fit_method": "least_squares_through_origin_using_largest_3_samples",
        "fit_sample_lengths": [int(r["sample_length"]) for r in fit_rows],
        "sec_per_base": sec_per_base,
        "predicted_chr1_runtime_sec": predicted_sec,
        "predicted_chr1_runtime_min": predicted_sec / 60.0,
        "predicted_chr1_runtime_hr": predicted_sec / 3600.0,
    }


def update_chr1_summary(estimate: dict[str, Any], samples_csv: Path, estimate_json: Path) -> None:
    summary_path = OUT / "summary.md"
    existing = summary_path.read_text() if summary_path.exists() else "# chr1 hg38 HSeeker Summary\n"
    marker = "## Triplex Runtime Estimate"
    if marker in existing:
        existing = existing.split(marker, 1)[0].rstrip() + "\n\n"
    existing = existing.replace(
        "\nNo Triplex runtime extrapolation or subset-sample results are included in this cleaned output.\n",
        "\n",
    )
    hseeker_sec = None
    run_json = OUT / "chr1_run_summary.json"
    if run_json.exists():
        data = json.loads(run_json.read_text())
        hseeker_sec = float(data.get("hseeker_runtime_sec") or 0.0)
    ratio_line = ""
    if hseeker_sec and hseeker_sec > 0:
        ratio = estimate["predicted_chr1_runtime_sec"] / hseeker_sec
        ratio_line = f"Relative to the measured HSeeker chr1 runtime ({hseeker_sec:.2f} sec), this is about {ratio:.1f}x slower.\n\n"

    block = (
        f"{marker}\n\n"
        "Full chr1 Triplex was not run. Runtime was extrapolated from Triplex runs "
        "on hg38 chr1 prefix samples using the same default Triplex settings as the "
        "validation comparison (`min_score=15`, `p_value=0.05`, `min_len=6`, "
        "`max_len=25`, `min_loop=3`, `max_loop=10`).\n\n"
        f"Sample timing table: `{samples_csv}`\n\n"
        f"Estimate JSON: `{estimate_json}`\n\n"
        f"Fit method: {estimate['fit_method']}.\n\n"
        f"Predicted full chr1 Triplex runtime: "
        f"{estimate['predicted_chr1_runtime_sec']:.1f} seconds "
        f"({estimate['predicted_chr1_runtime_min']:.1f} minutes; "
        f"{estimate['predicted_chr1_runtime_hr']:.2f} hours).\n\n"
        f"{ratio_line}"
        "This is an extrapolation, not a completed full-chromosome Triplex run.\n"
    )
    summary_path.write_text(existing + block)

    if run_json.exists():
        data["triplex_chr1_runtime_estimate"] = estimate
        data["triplex_chr1_runtime_samples_csv"] = str(samples_csv)
        run_json.write_text(json.dumps(data, indent=2) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    fasta, created_cache = ensure_chr1()
    try:
        seq = read_fasta_sequence(fasta)
        triplex_out = OUT / "triplex_chr1_runtime_samples.csv"
        with tempfile.TemporaryDirectory() as tmpdir:
            triplex_in = Path(tmpdir) / "triplex_chr1_runtime_samples_input.csv"
            metadata_rows = write_triplex_input(seq, triplex_in)

            t0 = time.perf_counter()
            subprocess.run(["Rscript", str(TRIPLEX_R), str(triplex_in), str(triplex_out)], cwd=ROOT, check=True)
            wall_sec = time.perf_counter() - t0

        pred_by_id = read_triplex_predictions(triplex_out)
        rows = []
        for row in metadata_rows:
            pred = pred_by_id[row["sequence_id"]]
            rows.append(
                {
                    "sequence_id": row["sequence_id"],
                    "sample_length": row["sample_length"],
                    "triplex_runtime_sec": float(pred["triplex_runtime_sec"]),
                    "triplex_hits": int(pred["triplex_count"]),
                    "triplex_best_score": pred["triplex_best_score"],
                    "triplex_best_pvalue": pred["triplex_best_pvalue"],
                }
            )
        write_csv(triplex_out, rows)

        estimate = estimate_runtime(rows, len(seq))
        estimate.update(
            {
                "chr1_length_bp": len(seq),
                "samples_csv": str(triplex_out),
                "triplex_subprocess_wall_sec": wall_sec,
                "triplex_runtime_note": "Per-sample Triplex timings exclude R package startup; subprocess wall time includes it.",
            }
        )
        estimate_json = OUT / "triplex_chr1_runtime_estimate.json"
        estimate_md = OUT / "triplex_chr1_runtime_estimate.md"
        estimate_json.write_text(json.dumps(estimate, indent=2) + "\n")
        estimate_md.write_text(
            "# Triplex chr1 Runtime Estimate\n\n"
            f"chr1 length: {len(seq):,} bp\n\n"
            f"Predicted runtime: {estimate['predicted_chr1_runtime_sec']:.1f} seconds "
            f"({estimate['predicted_chr1_runtime_min']:.1f} minutes; "
            f"{estimate['predicted_chr1_runtime_hr']:.2f} hours)\n\n"
            f"Sample timings: `{triplex_out}`\n\n"
            f"Fit method: {estimate['fit_method']}\n"
        )
        update_chr1_summary(estimate, triplex_out, estimate_json)
        print(json.dumps(estimate, indent=2))
    finally:
        if created_cache:
            fasta.unlink(missing_ok=True)
            fasta.with_suffix(".fna.gz").unlink(missing_ok=True)


if __name__ == "__main__":
    main()
