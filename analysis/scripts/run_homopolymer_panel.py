#!/usr/bin/env python3
"""Homopolymer panel, reported separately from the benchmark (reviewer comment R2.1).

The six single-base homopolymer records were removed from the benchmark on 2026-10-02. They are
reported here instead, because (1) the composition filters decide their calls by construction,
and those filters were designed with these sequences in view, so including them would score the
filters on the data used to set them; (2) a run of one base is mirror-symmetric at every
position, so it does not test mirror-repeat detection; (3) their labels depend on the assay
conditions (poly-dC forms an Mg2+-dependent triple helix; A20 stays B-form in the same plasmid
series where its neighbours form H-DNA).

Each record is scanned at the shipped configuration under the four filter states.

Usage:  python3 analysis/scripts/run_homopolymer_panel.py
Output: analysis/results/homopolymer_panel/homopolymer_panel.csv, metadata.json
"""

from __future__ import annotations

import csv
import hashlib
import json
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import hseeker  # noqa: E402

OUT = ROOT / "analysis" / "results" / "homopolymer_panel"
BASE = dict(minrep=8, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10, score=False)
FILTERS = {"none": dict(at_threshold=None, filter_homopolymers=False),
           "AT_only": dict(at_threshold=0.80, filter_homopolymers=False),
           "homopolymer_only": dict(at_threshold=None, filter_homopolymers=True),
           "both (shipped CLI)": dict(at_threshold=0.80, filter_homopolymers=True)}
SYNTHETIC = [  # generator class E (RTR/generate_synthetic_negatives.py), removed 2026-10-02
    ("SYN0053", "poly-A15", "A" * 15), ("SYN0054", "poly-A30", "A" * 30),
    ("SYN0055", "poly-T25", "T" * 25), ("SYN0056", "poly-T35", "T" * 35)]
EVIDENCE = {
    "HDNA0053": "Hanvey 1988 (PMID 2835375): A20.T20 remains B-form under all conditions tested",
    "HDNA0054": "Kohwi & Kohwi-Shigematsu 1988 (PMID 3375241): Mg2+-dependent triple helix in "
                "homopurine-homopyrimidine (dG)n.(dC)n",
    "SYN0053": "synthetic design; poly-A extends the Hanvey A20 result (not experimentally tested)",
    "SYN0054": "synthetic design; poly-A extends the Hanvey A20 result (not experimentally tested)",
    "SYN0055": "synthetic design; poly-T extends the Hanvey A20 result (not experimentally tested)",
    "SYN0056": "synthetic design; poly-T extends the Hanvey A20 result (not experimentally tested)",
}


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    exp = {r["record_id"]: r for r in csv.DictReader(open(ROOT / "hdna_benchmark_experimental_v3.csv"))}
    panel = [(k, exp[k]["sequence_name"], exp[k]["sequence_5to3"].upper(), exp[k]["label"], "experimental")
             for k in ("HDNA0053", "HDNA0054")]
    panel += [(k, n, s, "non_forming", "synthetic") for k, n, s in SYNTHETIC]
    rows = []
    for rid, name, seq, label, kind in panel:
        row = dict(record_id=rid, sequence_name=name, length=len(seq), base=seq[0], label=label,
                   record_type=kind, evidence=EVIDENCE[rid])
        for fname, f in FILTERS.items():
            hits = hseeker.scan_sequence(seq, **BASE, **f)
            row[f"call[{fname}]"] = "hit" if hits else "no hit"
            row[f"correct[{fname}]"] = bool(hits) == (label == "forming")
        rows.append(row)
    with open(OUT / "homopolymer_panel.csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    summary = {f: f"{sum(r[f'correct[{f}]'] for r in rows)}/{len(rows)} correct" for f in FILTERS}
    meta = dict(date=time.strftime("%Y-%m-%d %H:%M:%S"), base=BASE, filters=FILTERS, summary=summary,
                experimental_csv_sha256=hashlib.sha256(
                    (ROOT / "hdna_benchmark_experimental_v3.csv").read_bytes()).hexdigest())
    (OUT / "metadata.json").write_text(json.dumps(meta, indent=1))
    for r in rows:
        print(r["record_id"], r["sequence_name"], r["label"], [r[f"call[{f}]"] for f in FILTERS])
    print(summary)


if __name__ == "__main__":
    main()
