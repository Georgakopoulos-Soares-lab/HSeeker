#!/usr/bin/env python3
"""Overlap-threshold sweep and background rate for the E. coli insertion analysis.

Reuses the matching logic of run_injected_sensitivity_analysis.py so that the
sweep is computed exactly as the primary result is, and additionally scans the
unmodified E. coli genome to establish the background locus rate.

Outputs (into the analysis directory given by --dir):
  overlap_sweep.csv   metrics for HSeeker and triplex at each overlap threshold
  background.json     loci found in the unmodified genome at default settings
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
SCRIPT = ROOT / "analysis" / "scripts" / "run_injected_sensitivity_analysis.py"

spec = importlib.util.spec_from_file_location("injected", SCRIPT)
injected = importlib.util.module_from_spec(spec)
sys.modules["injected"] = injected  # dataclasses resolve types via sys.modules
spec.loader.exec_module(injected)

THRESHOLDS = [0.50, 0.60, 0.70, 0.80, 0.90, 0.95]


def read_manifest(path: Path) -> list:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    return [
        injected.InsertRecord(
            sequence_id=r["sequence_id"],
            sequence_name=r["sequence_name"],
            label=r["label"],
            sequence=r["sequence"],
            insertion_index=int(r["insertion_index"]),
            original_insert_pos_0based=int(r["original_insert_pos_0based"]),
            start=int(r["start"]),
            end=int(r["end"]),
        )
        for r in rows
    ]


def read_triplex(path: Path) -> list[dict]:
    with path.open(newline="") as fh:
        rows = list(csv.DictReader(fh))
    for r in rows:
        r["start"] = int(r["start"])
        r["end"] = int(r["end"])
        r["triplex_score"] = float(r.get("score") or r.get("triplex_score") or 0)
    return rows


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", type=Path,
                    default=ROOT / "analysis" / "results" / "sensitivity_injected_2026-10-15")
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()
    d = args.dir

    manifest = read_manifest(d / "insertion_manifest.csv")
    hs_hits = injected.read_hseeker_hits(d / "hdna_benchmark_HDNA.tsv")
    tx_hits = read_triplex(d / "triplex_genome_hits.csv")

    rows = []
    for frac in THRESHOLDS:
        hs = injected.best_matches(manifest, hs_hits, "total_score", frac)
        for r in hs:
            r["pred"] = bool(r["matched"])
        m = injected.metrics_from_binary(hs, "pred")
        rows.append({"tool": "HSeeker", "overlap_fraction": frac, **m})

        tx = injected.best_matches(manifest, tx_hits, "triplex_score", frac)
        for r in tx:
            r["pred"] = bool(r["matched"])
        m = injected.metrics_from_binary(tx, "pred")
        rows.append({"tool": "triplex", "overlap_fraction": frac, **m})

    out = d / "overlap_sweep.csv"
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")
    for r in rows:
        print(f"  {r['tool']:8s} overlap>={r['overlap_fraction']:.2f}  "
              f"TP {r['TP']:3d} FN {r['FN']:3d} TN {r['TN']:3d} FP {r['FP']:3d}  "
              f"sens {r['sensitivity']:.3f}  MCC {r['MCC']:.3f}")

    # ---- background rate in the unmodified genome --------------------------
    ref = d / "data" / "NC_000913.3.fna"
    prefix = d / "background_unmodified"
    cmd = [sys.executable, "-m", "hseeker", "-seq", str(ref), "-minrep", "8",
           "-maxspacer", "10", "-mismatch", "0.1", "-purity", "0.9",
           "-maxrep", "1000", "-out", str(prefix), "-workers", str(args.workers),
           "-purity-rmq"]
    print("\nScanning unmodified genome:", " ".join(cmd))
    subprocess.run(cmd, cwd=ROOT, check=True, capture_output=True, text=True)
    bg = injected.read_hseeker_hits(Path(str(prefix) + "_HDNA.tsv"))

    ref_len = sum(len(l.strip()) for l in ref.read_text().splitlines()[1:])
    info = {
        "reference_accession": "NC_000913.3",
        "reference_length_bp": ref_len,
        "background_loci_unmodified_genome": len(bg),
        "loci_per_mb": len(bg) / (ref_len / 1e6),
        "injected_genome_total_hits": len(hs_hits),
        "inserted_sequences": len(manifest),
        "settings": "minrep 8, maxspacer 10, mismatch 0.10, purity 0.90, "
                    "AT 0.80 + homopolymer filter, greedy overlap removal",
    }
    (d / "background.json").write_text(json.dumps(info, indent=2))
    print(json.dumps(info, indent=2))


if __name__ == "__main__":
    main()
