"""Small, explicit helpers for evaluating HSeeker ablation results."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import hseeker


def load_matched_cohort(table_path: str | Path, fasta_path: str | Path) -> list[dict]:
    """Load labeled records only when table IDs and sequences match the FASTA."""
    with open(table_path, newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.DictReader(handle))
    fasta_records = [
        (seq_id.rstrip(","), sequence)
        for seq_id, sequence, _ in hseeker.parse_fasta(fasta_path)
    ]
    fasta = dict(fasta_records)
    if len(fasta) != len(fasta_records):
        raise ValueError("FASTA contains repeated record IDs")
    ids = [row["record_id"] for row in rows]
    if len(ids) != len(set(ids)) or set(ids) != set(fasta):
        raise ValueError("Supplementary table and FASTA record IDs do not match")

    records = []
    for row in rows:
        seq_id = row["record_id"]
        sequence = row["sequence_5to3"].upper()
        if sequence != fasta[seq_id]:
            raise ValueError(f"FASTA sequence differs from table for {seq_id}")
        label = row["label"].strip().replace("-", "_")
        if label not in ("forming", "non_forming"):
            raise ValueError(f"Unknown label {label!r} for {seq_id}")
        records.append({"record_id": seq_id, "sequence": sequence, "label": label})
    return records


def confusion_summary(rows: list[dict]) -> dict:
    """Summarize record-level any-hit predictions without fitting a threshold."""
    tp = sum(row["label"] == "forming" and row["has_hit"] for row in rows)
    fn = sum(row["label"] == "forming" and not row["has_hit"] for row in rows)
    tn = sum(row["label"] == "non_forming" and not row["has_hit"] for row in rows)
    fp = sum(row["label"] == "non_forming" and row["has_hit"] for row in rows)
    return {
        "TP": tp, "FN": fn, "TN": tn, "FP": fp,
        "sensitivity": tp / (tp + fn) if tp + fn else float("nan"),
        "specificity": tn / (tn + fp) if tn + fp else float("nan"),
    }


def rank_auc(rows: list[dict], score_key: str = "max_score") -> float:
    """Pairwise AUROC; tied positive/negative scores receive half credit."""
    positives = [row[score_key] for row in rows if row["label"] == "forming"]
    negatives = [row[score_key] for row in rows if row["label"] == "non_forming"]
    if not positives or not negatives:
        return float("nan")
    wins = sum(p > n for p in positives for n in negatives)
    ties = sum(p == n for p in positives for n in negatives)
    return (wins + 0.5 * ties) / (len(positives) * len(negatives))


def max_insertion_overlap_fraction(hits: list[dict], start: int, end: int) -> float:
    """Maximum fraction of an inclusive insertion interval covered by one hit."""
    if end < start:
        raise ValueError("Insertion end precedes start")
    length = end - start + 1
    return max(
        (max(0, min(end, hit["end"]) - max(start, hit["start"]) + 1) / length
         for hit in hits),
        default=0.0,
    )


def load_verified_insertions(
    manifest_path: str | Path, genome: str, records: list[dict]
) -> list[dict]:
    """Verify the manifest's zero-based positions against the injected genome."""
    insertions = json.loads(Path(manifest_path).read_text())
    by_id = {row["record_id"]: row for row in records}
    if len(insertions) != len(by_id):
        raise ValueError("Insertion manifest and direct cohort have different sizes")
    seen = set()
    result = []
    for item in insertions:
        seq_id = item["record_id"]
        if seq_id not in by_id or seq_id in seen:
            raise ValueError(f"Missing or repeated insertion ID {seq_id}")
        seen.add(seq_id)
        sequence = item["sequence"].upper()
        position = item["position"]  # zero-based in this manifest
        if sequence != by_id[seq_id]["sequence"]:
            raise ValueError(f"Insertion sequence differs from cohort for {seq_id}")
        if item["label"].strip().replace("-", "_") != by_id[seq_id]["label"]:
            raise ValueError(f"Insertion label differs from cohort for {seq_id}")
        if position < 0 or position + len(sequence) > len(genome):
            raise ValueError(f"Insertion position outside genome for {seq_id}")
        if genome[position:position + len(sequence)] != sequence:
            raise ValueError(f"Inserted sequence not found at position for {seq_id}")
        result.append({
            "record_id": seq_id,
            "label": by_id[seq_id]["label"],
            "start": position + 1,
            "end": position + len(sequence),
        })
    return result
