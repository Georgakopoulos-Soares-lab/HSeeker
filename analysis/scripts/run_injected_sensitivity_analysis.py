#!/usr/bin/env python3
"""Secondary HSeeker/Triplex sensitivity benchmark on injected E. coli K-12.

This pipeline:
1. Downloads NC_000913.3 if needed.
2. Inserts the 54 experimental H-DNA sequences exactly as provided, with
   seed=42 and >=500 bp between insertion sites.
3. Runs HSeeker on the injected genome with the paper/default parameters.
4. Runs Triplex on the same injected genome with package defaults.
5. Matches detected hits to inserted intervals using the 80% overlap rule.
6. Computes score ROC/AUC, Youden threshold, MCC, and summary tables/plots.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import os
import random
import re
import shutil
import subprocess
import sys
import time
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import hseeker  # noqa: E402

ACCESSION = "NC_000913.3"
NCBI_FASTA_URL = (
    "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
    f"?db=nuccore&id={ACCESSION}&rettype=fasta&retmode=text"
)
DEFAULT_CSV = ROOT / "hdna_experimental_sequences_final.csv"
DEFAULT_OUT = ROOT / "analysis" / "results" / "sensitivity_injected"
TRIPLEX_R = ROOT / "analysis" / "scripts" / "triplex_fasta_search_defaults.R"


def normalize_label(label: str) -> str:
    value = label.strip().lower().replace("_", "-")
    if value in {"forming", "non-forming"}:
        return value
    return value


@dataclass(frozen=True)
class InsertRecord:
    sequence_id: str
    sequence_name: str
    label: str
    sequence: str
    insertion_index: int
    original_insert_pos_0based: int
    start: int
    end: int


def read_fasta(path: Path) -> tuple[str, str]:
    header = ""
    parts: list[str] = []
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if not header:
                    header = line[1:]
                continue
            parts.append(line.upper())
    return header, "".join(parts)


def write_fasta(path: Path, header: str, seq: str, width: int = 80) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as handle:
        handle.write(f">{header}\n")
        for i in range(0, len(seq), width):
            handle.write(seq[i:i + width] + "\n")


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def download_reference(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    print(f"Downloading {ACCESSION} from NCBI -> {path}", flush=True)
    with urllib.request.urlopen(NCBI_FASTA_URL, timeout=120) as response:
        tmp.write_bytes(response.read())
    header, seq = read_fasta(tmp)
    if ACCESSION not in header or len(seq) < 4_000_000:
        raise RuntimeError(f"Downloaded FASTA does not look like {ACCESSION}: header={header!r}, length={len(seq)}")
    tmp.replace(path)


def load_experimental_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    required = {"record_id", "sequence_name", "sequence_5to3", "label"}
    missing = required.difference(reader.fieldnames or [])
    if missing:
        raise RuntimeError(f"Missing required CSV columns: {sorted(missing)}")
    bad: list[dict[str, str]] = []
    for row in rows:
        seq = (row["sequence_5to3"] or "").strip().upper()
        if not seq or re.search(r"[^ACGT]", seq):
            bad.append({
                "record_id": row.get("record_id", ""),
                "sequence_name": row.get("sequence_name", ""),
                "sequence_5to3": row.get("sequence_5to3", ""),
                "invalid_characters": "".join(sorted(set(re.sub(r"[ACGT]", "", seq)))),
            })
    if bad:
        raise RuntimeError(
            "Input CSV contains non-ACGT or empty sequences. "
            "See invalid_sequences.csv in the output directory after rerun with --write-invalid-report."
        )
    labels = {normalize_label(row["label"]) for row in rows}
    if labels != {"forming", "non-forming"}:
        raise RuntimeError(f"Expected labels forming/non-forming only; observed {sorted(labels)}")
    return rows


def select_insertion_positions(genome_len: int, n: int, min_distance: int, seed: int) -> list[int]:
    rng = random.Random(seed)
    positions: list[int] = []
    attempts = 0
    max_attempts = 1_000_000
    while len(positions) < n and attempts < max_attempts:
        attempts += 1
        pos = rng.randrange(1_000, genome_len - 1_000)
        if all(abs(pos - other) > min_distance for other in positions):
            positions.append(pos)
    if len(positions) != n:
        raise RuntimeError(f"Could not place {n} insertions with min distance {min_distance}")
    return sorted(positions)


def inject_sequences(genome: str, rows: list[dict[str, str]], seed: int, min_distance: int) -> tuple[str, list[InsertRecord]]:
    positions = select_insertion_positions(len(genome), len(rows), min_distance, seed)
    output_parts: list[str] = []
    manifest: list[InsertRecord] = []
    cursor = 0
    cumulative_inserted = 0
    for idx, (pos, row) in enumerate(zip(positions, rows), start=1):
        seq = row["sequence_5to3"].strip().upper()
        output_parts.append(genome[cursor:pos])
        start = pos + cumulative_inserted + 1
        end = start + len(seq) - 1
        output_parts.append(seq)
        manifest.append(
            InsertRecord(
                sequence_id=row["record_id"],
                sequence_name=row["sequence_name"],
                label=normalize_label(row["label"]),
                sequence=seq,
                insertion_index=idx,
                original_insert_pos_0based=pos,
                start=start,
                end=end,
            )
        )
        cursor = pos
        cumulative_inserted += len(seq)
    output_parts.append(genome[cursor:])
    return "".join(output_parts), manifest


def overlap_len(a_start: int, a_end: int, b_start: int, b_end: int) -> int:
    return max(0, min(a_end, b_end) - max(a_start, b_start) + 1)


def is_overlap_match(hit_start: int, hit_end: int, ins_start: int, ins_end: int, min_fraction: float) -> tuple[bool, int, float, float]:
    ov = overlap_len(hit_start, hit_end, ins_start, ins_end)
    hit_len = hit_end - hit_start + 1
    ins_len = ins_end - ins_start + 1
    frac_hit = ov / hit_len if hit_len > 0 else 0.0
    frac_insert = ov / ins_len if ins_len > 0 else 0.0
    return frac_hit >= min_fraction or frac_insert >= min_fraction, ov, frac_hit, frac_insert


def run_hseeker(fasta: Path, out_prefix: Path, workers: int) -> tuple[Path, float]:
    cmd = [
        sys.executable, "-m", "hseeker",
        "-seq", str(fasta),
        "-minrep", "8",
        "-maxspacer", "10",
        "-mismatch", "0.1",
        "-purity", "0.9",
        "-maxrep", "1000",
        "-out", str(out_prefix),
        "-workers", str(workers),
        "-purity-rmq",
    ]
    env = os.environ.copy()
    env["PYTHONPATH"] = str(SRC) + os.pathsep + env.get("PYTHONPATH", "")
    print("Running HSeeker:", " ".join(cmd), flush=True)
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=ROOT, env=env, capture_output=True, text=True)
    runtime = time.perf_counter() - t0
    (out_prefix.parent / "hseeker.stdout.txt").write_text(proc.stdout)
    (out_prefix.parent / "hseeker.stderr.txt").write_text(proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"HSeeker failed with return code {proc.returncode}; see hseeker.stderr.txt")
    return Path(str(out_prefix) + "_HDNA.tsv"), runtime


def read_hseeker_hits(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        rows = []
        for row in reader:
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            row["total_score"] = float(row["total_score"]) if row.get("total_score") not in ("", None, "None") else math.nan
            rows.append(row)
    return rows


def run_triplex(fasta: Path, out_csv: Path) -> tuple[list[dict[str, Any]], float]:
    cmd = ["Rscript", str(TRIPLEX_R), str(fasta), str(out_csv)]
    print("Running Triplex package defaults:", " ".join(cmd), flush=True)
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    runtime_wall = time.perf_counter() - t0
    (out_csv.parent / "triplex.stdout.txt").write_text(proc.stdout)
    (out_csv.parent / "triplex.stderr.txt").write_text(proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"Triplex failed with return code {proc.returncode}; see triplex.stderr.txt")
    runtime_file = Path(str(out_csv) + ".runtime.txt")
    runtime = runtime_wall
    if runtime_file.exists():
        text = runtime_file.read_text().strip().split(",")
        if len(text) == 2:
            runtime = float(text[1])
    rows: list[dict[str, Any]] = []
    with out_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            row["triplex_score"] = float(row["triplex_score"])
            rows.append(row)
    return rows, runtime


def best_matches(
    manifest: list[InsertRecord],
    hits: list[dict[str, Any]],
    score_key: str,
    min_fraction: float,
) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for ins in manifest:
        candidates: list[dict[str, Any]] = []
        for hit in hits:
            ok, ov, frac_hit, frac_insert = is_overlap_match(
                int(hit["start"]), int(hit["end"]), ins.start, ins.end, min_fraction
            )
            if ok:
                score = float(hit.get(score_key, math.nan))
                candidates.append({
                    "hit_start": int(hit["start"]),
                    "hit_end": int(hit["end"]),
                    "hit_length": int(hit["end"]) - int(hit["start"]) + 1,
                    "overlap_bp": ov,
                    "overlap_fraction_hit": frac_hit,
                    "overlap_fraction_insert": frac_insert,
                    "score": score,
                    "raw_hit": hit,
                })
        candidates.sort(key=lambda x: (
            -1 if math.isnan(x["score"]) else -x["score"],
            -x["overlap_bp"],
            x["hit_start"],
        ))
        best = candidates[0] if candidates else None
        rows.append({
            "sequence_id": ins.sequence_id,
            "sequence_name": ins.sequence_name,
            "label": ins.label,
            "insert_start": ins.start,
            "insert_end": ins.end,
            "insert_length": len(ins.sequence),
            "matched": best is not None,
            "matched_hit_count": len(candidates),
            "best_hit_start": best["hit_start"] if best else "",
            "best_hit_end": best["hit_end"] if best else "",
            "best_hit_length": best["hit_length"] if best else "",
            "best_overlap_bp": best["overlap_bp"] if best else 0,
            "best_overlap_fraction_hit": best["overlap_fraction_hit"] if best else 0.0,
            "best_overlap_fraction_insert": best["overlap_fraction_insert"] if best else 0.0,
            "best_score": best["score"] if best else math.nan,
        })
    return rows


def roc_points(labels: list[int], scores: list[float]) -> list[dict[str, float]]:
    finite = sorted({s for s in scores if not math.isnan(s)}, reverse=True)
    thresholds = [math.inf, *finite, -math.inf]
    points = []
    for thr in thresholds:
        preds = [0 if math.isnan(s) else int(s >= thr) for s in scores]
        tp = sum(1 for y, p in zip(labels, preds) if y == 1 and p == 1)
        fn = sum(1 for y, p in zip(labels, preds) if y == 1 and p == 0)
        tn = sum(1 for y, p in zip(labels, preds) if y == 0 and p == 0)
        fp = sum(1 for y, p in zip(labels, preds) if y == 0 and p == 1)
        sensitivity = tp / (tp + fn) if tp + fn else 0.0
        specificity = tn / (tn + fp) if tn + fp else 0.0
        fpr = 1.0 - specificity
        precision = tp / (tp + fp) if tp + fp else 0.0
        accuracy = (tp + tn) / len(labels) if labels else 0.0
        denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        mcc = ((tp * tn) - (fp * fn)) / denom if denom else 0.0
        points.append({
            "threshold": thr,
            "TP": tp,
            "FN": fn,
            "TN": tn,
            "FP": fp,
            "sensitivity": sensitivity,
            "specificity": specificity,
            "precision": precision,
            "accuracy": accuracy,
            "FPR": fpr,
            "TPR": sensitivity,
            "youden_j": sensitivity + specificity - 1.0,
            "MCC": mcc,
        })
    return points


def auc_roc_from_scores(labels: list[int], scores: list[float]) -> float:
    """Return ROC AUC for higher-is-more-positive scores.

    Missing scores are treated as -Inf because a sequence with no overlapping
    HSeeker hit is negative for every finite score threshold.
    """
    pos = [(-math.inf if math.isnan(s) else s) for y, s in zip(labels, scores) if y == 1]
    neg = [(-math.inf if math.isnan(s) else s) for y, s in zip(labels, scores) if y == 0]
    if not pos or not neg:
        return math.nan
    wins = 0.0
    for p in pos:
        for n in neg:
            if p > n:
                wins += 1.0
            elif p == n:
                wins += 0.5
    return wins / (len(pos) * len(neg))


def best_youden(points: list[dict[str, float]]) -> dict[str, float]:
    return max(points, key=lambda p: (p["youden_j"], p["sensitivity"], p["specificity"]))


def _pdf_escape(text: str) -> str:
    return text.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")


def _pdf_stream(commands: list[str]) -> bytes:
    body = "\n".join(commands).encode("ascii")
    return b"<< /Length " + str(len(body)).encode("ascii") + b" >>\nstream\n" + body + b"\nendstream"


def simple_pdf_roc(path: Path, points: list[dict[str, float]], auc: float) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width, height = 720, 520
    margin = 70
    plot_w = width - 2 * margin
    plot_h = height - 2 * margin
    xy = sorted({(p["FPR"], p["TPR"]) for p in points})
    coords = [
        (margin + x * plot_w, margin + y * plot_h)
        for x, y in xy
    ]
    commands = [
        "1 1 1 rg 0 0 720 520 re f",
        "0.13 0.13 0.13 RG 1 w",
        f"{margin} {margin} m {width - margin} {margin} l S",
        f"{margin} {margin} m {margin} {height - margin} l S",
        "0.6 0.6 0.6 RG 1 w [6 6] 0 d",
        f"{margin} {margin} m {width - margin} {height - margin} l S",
        "[] 0 d",
        "0.13 0.40 0.67 RG 3 w",
    ]
    if coords:
        x0, y0 = coords[0]
        commands.append(f"{x0:.1f} {y0:.1f} m")
        for x, y in coords[1:]:
            commands.append(f"{x:.1f} {y:.1f} l")
        commands.append("S")
    commands.extend([
        "0 0 0 rg",
        "BT /F1 20 Tf 70 485 Td (HSeeker score ROC on injected E. coli benchmark) Tj ET",
        f"BT /F1 16 Tf 70 460 Td (AUC = {auc:.3f}) Tj ET",
        "BT /F1 18 Tf 278 28 Td (False positive rate) Tj ET",
        "q 0 1 -1 0 24 205 cm BT /F1 18 Tf 0 0 Td (True positive rate) Tj ET Q",
        "BT /F1 12 Tf 63 50 Td (0) Tj ET",
        "BT /F1 12 Tf 640 50 Td (1) Tj ET",
        "BT /F1 12 Tf 45 67 Td (0) Tj ET",
        "BT /F1 12 Tf 45 445 Td (1) Tj ET",
    ])

    objects: list[bytes] = [
        b"<< /Type /Catalog /Pages 2 0 R >>",
        b"<< /Type /Pages /Kids [3 0 R] /Count 1 >>",
        b"<< /Type /Page /Parent 2 0 R /MediaBox [0 0 720 520] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>",
        b"<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>",
        _pdf_stream(commands),
    ]
    chunks = [b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n"]
    offsets = [0]
    for idx, obj in enumerate(objects, start=1):
        offsets.append(sum(len(c) for c in chunks))
        chunks.append(f"{idx} 0 obj\n".encode("ascii") + obj + b"\nendobj\n")
    xref_offset = sum(len(c) for c in chunks)
    chunks.append(f"xref\n0 {len(objects) + 1}\n".encode("ascii"))
    chunks.append(b"0000000000 65535 f \n")
    for offset in offsets[1:]:
        chunks.append(f"{offset:010d} 00000 n \n".encode("ascii"))
    chunks.append(
        f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref_offset}\n%%EOF\n".encode("ascii")
    )
    path.write_bytes(b"".join(chunks))


def metrics_from_binary(rows: list[dict[str, Any]], pred_key: str) -> dict[str, float]:
    tp = sum(1 for r in rows if r["label"] == "forming" and bool(r[pred_key]))
    fn = sum(1 for r in rows if r["label"] == "forming" and not bool(r[pred_key]))
    tn = sum(1 for r in rows if r["label"] == "non-forming" and not bool(r[pred_key]))
    fp = sum(1 for r in rows if r["label"] == "non-forming" and bool(r[pred_key]))
    sensitivity = tp / (tp + fn) if tp + fn else 0.0
    specificity = tn / (tn + fp) if tn + fp else 0.0
    precision = tp / (tp + fp) if tp + fp else 0.0
    accuracy = (tp + tn) / len(rows) if rows else 0.0
    f1 = 2 * precision * sensitivity / (precision + sensitivity) if precision + sensitivity else 0.0
    denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    mcc = ((tp * tn) - (fp * fn)) / denom if denom else 0.0
    return {
        "TP": tp, "FN": fn, "TN": tn, "FP": fp,
        "sensitivity": sensitivity, "specificity": specificity,
        "precision": precision, "F1": f1, "accuracy": accuracy, "MCC": mcc,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--min-distance", type=int, default=500)
    parser.add_argument("--overlap-fraction", type=float, default=0.80)
    parser.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    parser.add_argument("--skip-triplex", action="store_true")
    parser.add_argument("--write-invalid-report", action="store_true")
    args = parser.parse_args()

    out = args.outdir
    data_dir = out / "data"
    plots_dir = out / "plots"
    out.mkdir(parents=True, exist_ok=True)
    data_dir.mkdir(parents=True, exist_ok=True)
    plots_dir.mkdir(parents=True, exist_ok=True)

    rows = load_experimental_rows(args.csv)
    ref_fasta = data_dir / f"{ACCESSION}.fna"
    if not ref_fasta.exists():
        download_reference(ref_fasta)

    ref_header, ref_seq = read_fasta(ref_fasta)
    injected_seq, manifest = inject_sequences(ref_seq, rows, args.seed, args.min_distance)
    injected_fasta = data_dir / "hdna_benchmark_injected.fna"
    write_fasta(injected_fasta, f"{ACCESSION}_hdna_injected seed={args.seed}", injected_seq)

    manifest_rows = [record.__dict__ for record in manifest]
    write_csv(out / "insertion_manifest.csv", manifest_rows)

    hseeker_prefix = out / "hdna_benchmark"
    hseeker_tsv, hseeker_runtime = run_hseeker(injected_fasta, hseeker_prefix, args.workers)
    hseeker_hits = read_hseeker_hits(hseeker_tsv)
    hseeker_match_rows = best_matches(manifest, hseeker_hits, "total_score", args.overlap_fraction)
    for row in hseeker_match_rows:
        row["hseeker_pred_by_overlap"] = bool(row["matched"])
        row["hseeker_score"] = row["best_score"]
    write_csv(out / "hseeker_insert_matches.csv", hseeker_match_rows)

    triplex_match_rows: list[dict[str, Any]] = []
    triplex_runtime = math.nan
    triplex_hits: list[dict[str, Any]] = []
    if not args.skip_triplex:
        triplex_csv = out / "triplex_genome_hits.csv"
        triplex_hits, triplex_runtime = run_triplex(injected_fasta, triplex_csv)
        triplex_match_rows = best_matches(manifest, triplex_hits, "triplex_score", args.overlap_fraction)
        for row in triplex_match_rows:
            row["triplex_pred_by_overlap"] = bool(row["matched"])
            row["triplex_score"] = row["best_score"]
        write_csv(out / "triplex_insert_matches.csv", triplex_match_rows)

    labels = [1 if r.label == "forming" else 0 for r in manifest]
    scores = [float(r["hseeker_score"]) if not math.isnan(float(r["hseeker_score"])) else math.nan for r in hseeker_match_rows]
    roc = roc_points(labels, scores)
    auc = auc_roc_from_scores(labels, scores)
    youden = best_youden(roc)
    roc_rows = []
    for p in roc:
        row = dict(p)
        if row["threshold"] == math.inf:
            row["threshold"] = "Inf"
        elif row["threshold"] == -math.inf:
            row["threshold"] = "-Inf"
        roc_rows.append(row)
    write_csv(out / "hseeker_score_roc.csv", roc_rows)
    simple_pdf_roc(plots_dir / "hseeker_score_roc.pdf", roc, auc)

    comparison: list[dict[str, Any]] = []
    triplex_by_id = {r["sequence_id"]: r for r in triplex_match_rows}
    for hr in hseeker_match_rows:
        tr = triplex_by_id.get(hr["sequence_id"], {})
        comparison.append({
            "sequence_id": hr["sequence_id"],
            "sequence_name": hr["sequence_name"],
            "label": hr["label"],
            "insert_start": hr["insert_start"],
            "insert_end": hr["insert_end"],
            "insert_length": hr["insert_length"],
            "hseeker_overlap_pred": bool(hr["matched"]),
            "hseeker_score": hr["hseeker_score"],
            "hseeker_score_pred_youden": (
                False if math.isnan(float(hr["hseeker_score"])) else float(hr["hseeker_score"]) >= float(youden["threshold"])
            ),
            "triplex_overlap_pred": bool(tr.get("matched", False)),
            "triplex_score": tr.get("best_score", math.nan),
            "hseeker_best_hit_start": hr["best_hit_start"],
            "hseeker_best_hit_end": hr["best_hit_end"],
            "triplex_best_hit_start": tr.get("best_hit_start", ""),
            "triplex_best_hit_end": tr.get("best_hit_end", ""),
        })
    write_csv(out / "method_comparison_by_insertion.csv", comparison)

    hseeker_overlap_metrics = metrics_from_binary(comparison, "hseeker_overlap_pred")
    hseeker_score_metrics = metrics_from_binary(comparison, "hseeker_score_pred_youden")
    triplex_metrics = metrics_from_binary(comparison, "triplex_overlap_pred") if triplex_match_rows else {}

    metrics_rows = [
        {"method": "HSeeker overlap", **hseeker_overlap_metrics},
        {"method": "HSeeker score Youden", **hseeker_score_metrics},
    ]
    if triplex_metrics:
        metrics_rows.append({"method": "Triplex overlap", **triplex_metrics})
    write_csv(out / "summary_metrics.csv", metrics_rows)

    metadata = {
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "reference_accession": ACCESSION,
        "reference_url": NCBI_FASTA_URL,
        "reference_header": ref_header,
        "reference_length": len(ref_seq),
        "injected_length": len(injected_seq),
        "input_csv": str(args.csv),
        "dataset_size": len(manifest),
        "forming_count": sum(1 for r in manifest if r.label == "forming"),
        "nonforming_count": sum(1 for r in manifest if r.label == "non-forming"),
        "seed": args.seed,
        "min_inter_insertion_distance_bp": args.min_distance,
        "overlap_match_rule": f"overlap >= {args.overlap_fraction:.2f} of hit length OR insertion length",
        "hseeker_command": (
            "hseeker -seq hdna_benchmark_injected.fna -minrep 8 -maxspacer 10 "
            "-mismatch 0.1 -purity 0.9 -maxrep 1000 -out hdna_benchmark"
        ),
        "hseeker_purity_rmq": True,
        "hseeker_runtime_sec": hseeker_runtime,
        "hseeker_total_hits": len(hseeker_hits),
        "triplex_call": "triplex.search(DNAString(injected_genome))",
        "triplex_package_defaults": str(subprocess.run(
            ["Rscript", "-e", ".libPaths(c('.r-lib', .libPaths())); library(triplex); cat(capture.output(args(triplex.search)), sep='\\n')"],
            cwd=ROOT, capture_output=True, text=True
        ).stdout.strip()),
        "triplex_runtime_sec": triplex_runtime,
        "triplex_total_hits": len(triplex_hits),
        "hseeker_score_auc_roc": auc,
        "hseeker_youden_threshold": youden["threshold"],
        "hseeker_youden_j": youden["youden_j"],
        "hseeker_youden_mcc": youden["MCC"],
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    lines = [
        "# Secondary Injected E. coli H-DNA Sensitivity Analysis",
        "",
        "## Dataset",
        "",
        "- This injected-genome benchmark is a secondary/context validation. The primary validation is the direct-on-sequences benchmark.",
        f"- Experimental sequences: {len(manifest)} ({metadata['forming_count']} forming, {metadata['nonforming_count']} non-forming).",
        f"- Reference genome: {ACCESSION}; length {len(ref_seq):,} bp.",
        f"- Injected genome length: {len(injected_seq):,} bp.",
        f"- Insertion seed: {args.seed}; minimum inter-insertion distance: {args.min_distance} bp.",
        f"- Insertions used sequences exactly as provided in `{args.csv.name}`.",
        "",
        "## Methods",
        "",
        "- HSeeker command: `hseeker -seq hdna_benchmark_injected.fna -minrep 8 -maxspacer 10 -mismatch 0.1 -purity 0.9 -maxrep 1000 -out hdna_benchmark`.",
        "- This run used the exact purity RMQ prefilter (`-purity-rmq`) to preserve HSeeker output while reducing runtime.",
        "- Triplex was run on the same injected FASTA using `triplex.search(DNAString(seq))` with the installed package defaults.",
        f"- A detected hit matched an insertion if overlap was at least {args.overlap_fraction:.0%} of either hit length or insertion length.",
        "",
        "## Runtime",
        "",
        f"- HSeeker runtime: {hseeker_runtime:.3f} sec; hits: {len(hseeker_hits):,}.",
        f"- Triplex runtime: {triplex_runtime:.3f} sec; hits: {len(triplex_hits):,}." if not math.isnan(triplex_runtime) else "- Triplex was skipped.",
        "",
        "## Classification Metrics",
        "",
        "| method | TP | FN | TN | FP | sensitivity | specificity | precision | F1 | accuracy | MCC |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in metrics_rows:
        lines.append(
            f"| {row['method']} | {row['TP']} | {row['FN']} | {row['TN']} | {row['FP']} | "
            f"{row['sensitivity']:.3f} | {row['specificity']:.3f} | {row['precision']:.3f} | "
            f"{row['F1']:.3f} | {row['accuracy']:.3f} | {row['MCC']:.3f} |"
        )
    lines.extend([
        "",
        "## HSeeker Score ROC",
        "",
        f"- AUC-ROC: {auc:.3f}.",
        f"- Youden-optimal threshold: {float(youden['threshold']):.6g}.",
        f"- Youden J: {youden['youden_j']:.3f}; MCC at this threshold: {youden['MCC']:.3f}.",
        "",
        "## Output Files",
        "",
        "- `data/NC_000913.3.fna`: downloaded E. coli K-12 MG1655 reference.",
        "- `data/hdna_benchmark_injected.fna`: injected benchmark genome.",
        "- `insertion_manifest.csv`: inserted sequence coordinates.",
        "- `hseeker_insert_matches.csv`: best HSeeker match per insertion.",
        "- `triplex_insert_matches.csv`: best Triplex match per insertion.",
        "- `method_comparison_by_insertion.csv`: per-insertion HSeeker/Triplex comparison.",
        "- `hseeker_score_roc.csv`: ROC thresholds and metrics.",
        "- `plots/hseeker_score_roc.pdf`: ROC plot.",
    ])
    (out / "summary.md").write_text("\n".join(lines) + "\n")

    print("\nDone.")
    print(f"Summary: {out / 'summary.md'}")
    print(f"Metrics: {out / 'summary_metrics.csv'}")
    print(f"Comparison: {out / 'method_comparison_by_insertion.csv'}")


if __name__ == "__main__":
    main()
