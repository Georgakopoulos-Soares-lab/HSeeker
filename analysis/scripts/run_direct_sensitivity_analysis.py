#!/usr/bin/env python3
"""Primary HSeeker/Triplex sensitivity benchmark on experimental sequences.

This is the main validation analysis: HSeeker and Triplex are run directly on
the 54 experimentally curated sequences, without genomic insertion. The E. coli
injection benchmark is kept separately as a secondary/context validation.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import hseeker  # noqa: E402

DEFAULT_CSV = ROOT / "hdna_experimental_sequences_final.csv"
DEFAULT_OUT = ROOT / "analysis" / "results" / "sensitivity_direct"
TRIPLEX_R = ROOT / "analysis" / "scripts" / "triplex_fasta_search_defaults.R"


def normalize_label(label: str) -> str:
    value = label.strip().lower().replace("_", "-")
    if value in {"forming", "non-forming"}:
        return value
    return value


def write_csv(path: Path, rows: list[dict[str, Any]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if fieldnames is None:
        fieldnames = list(rows[0].keys()) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def load_rows(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    required = {"record_id", "sequence_name", "sequence_5to3", "label"}
    missing = required.difference(reader.fieldnames or [])
    if missing:
        raise RuntimeError(f"Missing required CSV columns: {sorted(missing)}")

    out: list[dict[str, Any]] = []
    invalid: list[dict[str, str]] = []
    for idx, row in enumerate(rows, start=1):
        sequence = (row["sequence_5to3"] or "").strip().upper()
        label = normalize_label(row["label"])
        bad = "".join(sorted(set(re.sub(r"[ACGT]", "", sequence))))
        if not sequence or bad:
            invalid.append({
                "sequence_id": row.get("record_id") or f"seq_{idx}",
                "sequence_name": row.get("sequence_name", ""),
                "invalid_characters": bad,
                "sequence": sequence,
            })
            continue
        if label not in {"forming", "non-forming"}:
            raise RuntimeError(f"Unexpected label for {row.get('record_id')}: {row['label']}")
        out.append({
            "sequence_id": row.get("record_id") or f"seq_{idx}",
            "sequence_name": row.get("sequence_name", ""),
            "label": label,
            "sequence": sequence,
            "sequence_length": len(sequence),
        })
    if invalid:
        raise RuntimeError(f"Found {len(invalid)} empty/non-ACGT sequences; fix the CSV before benchmarking.")
    return out


def write_fasta(path: Path, rows: list[dict[str, Any]]) -> None:
    with path.open("w") as handle:
        for row in rows:
            handle.write(f">{row['sequence_id']}\n")
            seq = row["sequence"]
            for i in range(0, len(seq), 80):
                handle.write(seq[i:i + 80] + "\n")


def run_hseeker_direct(rows: list[dict[str, Any]]) -> tuple[list[dict[str, Any]], float]:
    results: list[dict[str, Any]] = []
    t0 = time.perf_counter()
    for row in rows:
        hits = hseeker.scan_sequence(
            row["sequence"],
            minrep=8,
            maxrep=1000,
            maxspacer=10,
            purity=0.90,
            mismatch=0.10,
            remove_overlaps=True,
            seq_offset=1,
            score=True,
            purity_rmq=True,
        )
        best = None
        if hits:
            best = max(
                hits,
                key=lambda h: (
                    -math.inf if h.get("total_score") is None else float(h.get("total_score")),
                    h.get("arm_length", 0),
                    -h.get("spacer_length", 0),
                ),
            )
        results.append({
            "sequence_id": row["sequence_id"],
            "sequence_name": row["sequence_name"],
            "label": row["label"],
            "sequence_length": row["sequence_length"],
            "hseeker_pred_by_hit": best is not None,
            "hseeker_hit_count": len(hits),
            "hseeker_score": float(best["total_score"]) if best and best.get("total_score") is not None else math.nan,
            "hseeker_best_start": best["start"] if best else "",
            "hseeker_best_end": best["end"] if best else "",
            "hseeker_best_arm_length": best["arm_length"] if best else "",
            "hseeker_best_spacer_length": best["spacer_length"] if best else "",
            "hseeker_best_ga_pct": best["ga_pct"] if best else "",
            "hseeker_best_ct_pct": best["ct_pct"] if best else "",
            "hseeker_best_mirror_identity": best["mirror_identity"] if best else "",
        })
    return results, time.perf_counter() - t0


def run_triplex_direct(fasta: Path, out_csv: Path) -> tuple[list[dict[str, Any]], float]:
    cmd = ["Rscript", str(TRIPLEX_R), str(fasta), str(out_csv)]
    t0 = time.perf_counter()
    proc = subprocess.run(cmd, cwd=ROOT, capture_output=True, text=True)
    wall = time.perf_counter() - t0
    (out_csv.parent / "triplex.stdout.txt").write_text(proc.stdout)
    (out_csv.parent / "triplex.stderr.txt").write_text(proc.stderr)
    if proc.returncode != 0:
        raise RuntimeError(f"Triplex failed with return code {proc.returncode}; see triplex.stderr.txt")
    runtime = wall
    runtime_file = Path(str(out_csv) + ".runtime.txt")
    if runtime_file.exists():
        parts = runtime_file.read_text().strip().split(",")
        if len(parts) == 2:
            runtime = float(parts[1])
    hits: list[dict[str, Any]] = []
    with out_csv.open(newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            row["start"] = int(row["start"])
            row["end"] = int(row["end"])
            row["triplex_score"] = float(row["triplex_score"])
            row["triplex_pvalue"] = float(row["triplex_pvalue"])
            hits.append(row)
    return hits, runtime


def merge_triplex(rows: list[dict[str, Any]], hits: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_id: dict[str, list[dict[str, Any]]] = {}
    for hit in hits:
        by_id.setdefault(hit["seq_id"], []).append(hit)
    merged = []
    for row in rows:
        candidates = by_id.get(row["sequence_id"], [])
        best = max(candidates, key=lambda h: (h["triplex_score"], -h["triplex_pvalue"])) if candidates else None
        merged.append({
            "sequence_id": row["sequence_id"],
            "triplex_pred_by_hit": best is not None,
            "triplex_hit_count": len(candidates),
            "triplex_score": best["triplex_score"] if best else math.nan,
            "triplex_pvalue": best["triplex_pvalue"] if best else math.nan,
            "triplex_type": best["triplex_type"] if best else "",
            "triplex_best_start": best["start"] if best else "",
            "triplex_best_end": best["end"] if best else "",
        })
    return merged


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
        precision = tp / (tp + fp) if tp + fp else 0.0
        denom = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        points.append({
            "threshold": thr,
            "TP": tp, "FN": fn, "TN": tn, "FP": fp,
            "sensitivity": sensitivity,
            "specificity": specificity,
            "precision": precision,
            "accuracy": (tp + tn) / len(labels) if labels else 0.0,
            "FPR": 1.0 - specificity,
            "TPR": sensitivity,
            "youden_j": sensitivity + specificity - 1.0,
            "MCC": ((tp * tn) - (fp * fn)) / denom if denom else 0.0,
        })
    return points


def auc_roc_from_scores(labels: list[int], scores: list[float]) -> float:
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


def _pdf_stream(commands: list[str]) -> bytes:
    body = "\n".join(commands).encode("ascii")
    return b"<< /Length " + str(len(body)).encode("ascii") + b" >>\nstream\n" + body + b"\nendstream"


def simple_pdf_roc(path: Path, points: list[dict[str, float]], auc: float, title: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    width, height = 720, 520
    margin = 70
    plot_w = width - 2 * margin
    plot_h = height - 2 * margin
    coords = [(margin + p["FPR"] * plot_w, margin + p["TPR"] * plot_h) for p in sorted(points, key=lambda p: (p["FPR"], p["TPR"]))]
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
    safe_title = title.replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")
    commands.extend([
        "0 0 0 rg",
        f"BT /F1 20 Tf 70 485 Td ({safe_title}) Tj ET",
        f"BT /F1 16 Tf 70 460 Td (AUC = {auc:.3f}) Tj ET",
        "BT /F1 18 Tf 278 28 Td (False positive rate) Tj ET",
        "q 0 1 -1 0 24 205 cm BT /F1 18 Tf 0 0 Td (True positive rate) Tj ET Q",
        "BT /F1 12 Tf 63 50 Td (0) Tj ET",
        "BT /F1 12 Tf 640 50 Td (1) Tj ET",
        "BT /F1 12 Tf 45 67 Td (0) Tj ET",
        "BT /F1 12 Tf 45 445 Td (1) Tj ET",
    ])
    objects = [
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
    chunks.append(f"trailer\n<< /Size {len(objects) + 1} /Root 1 0 R >>\nstartxref\n{xref_offset}\n%%EOF\n".encode("ascii"))
    path.write_bytes(b"".join(chunks))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    args = parser.parse_args()

    out = args.outdir
    plots = out / "plots"
    out.mkdir(parents=True, exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)

    rows = load_rows(args.csv)
    fasta = out / "direct_sequences.fasta"
    write_fasta(fasta, rows)

    h_rows, h_runtime = run_hseeker_direct(rows)
    write_csv(out / "hseeker_direct_predictions.csv", h_rows)

    triplex_hits, t_runtime = run_triplex_direct(fasta, out / "triplex_direct_hits.csv")
    t_rows = merge_triplex(rows, triplex_hits)
    write_csv(out / "triplex_direct_predictions.csv", t_rows)

    t_by_id = {r["sequence_id"]: r for r in t_rows}
    comparison: list[dict[str, Any]] = []
    for h in h_rows:
        t = t_by_id[h["sequence_id"]]
        comparison.append({**h, **{k: v for k, v in t.items() if k != "sequence_id"}})

    labels = [1 if r["label"] == "forming" else 0 for r in comparison]
    scores = [float(r["hseeker_score"]) if not math.isnan(float(r["hseeker_score"])) else math.nan for r in comparison]
    roc = roc_points(labels, scores)
    auc = auc_roc_from_scores(labels, scores)
    youden = best_youden(roc)
    for row in comparison:
        score = float(row["hseeker_score"])
        row["hseeker_score_pred_youden"] = False if math.isnan(score) else score >= float(youden["threshold"])

    write_csv(out / "method_comparison_by_sequence.csv", comparison)

    roc_rows = []
    for p in roc:
        row = dict(p)
        row["threshold"] = "Inf" if row["threshold"] == math.inf else "-Inf" if row["threshold"] == -math.inf else row["threshold"]
        roc_rows.append(row)
    write_csv(out / "hseeker_score_roc.csv", roc_rows)
    simple_pdf_roc(plots / "hseeker_score_roc.pdf", roc, auc, "HSeeker score ROC on direct sequence validation")

    metrics_rows = [
        {"method": "HSeeker direct hit", **metrics_from_binary(comparison, "hseeker_pred_by_hit")},
        {"method": "HSeeker score Youden", **metrics_from_binary(comparison, "hseeker_score_pred_youden")},
        {"method": "Triplex direct hit", **metrics_from_binary(comparison, "triplex_pred_by_hit")},
    ]
    write_csv(out / "summary_metrics.csv", metrics_rows)

    metadata = {
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "input_csv": str(args.csv),
        "dataset_size": len(rows),
        "forming_count": sum(1 for r in rows if r["label"] == "forming"),
        "nonforming_count": sum(1 for r in rows if r["label"] == "non-forming"),
        "hseeker_parameters": {
            "minrep": 8, "maxrep": 1000, "maxspacer": 10, "purity": 0.90, "mismatch": 0.10,
            "score": True, "purity_rmq": True,
        },
        "triplex_call": "triplex.search(DNAString(seq))",
        "triplex_package_defaults": "type=0:7, min_score=15, p_value=0.05, min_len=6, max_len=25, min_loop=3, max_loop=10, seq_type='eukaryotic', scoring/penalty/statistical constants='default'",
        "hseeker_runtime_sec": h_runtime,
        "triplex_runtime_sec": t_runtime,
        "hseeker_score_auc_roc": auc,
        "hseeker_youden_threshold": youden["threshold"],
        "hseeker_youden_j": youden["youden_j"],
        "hseeker_youden_mcc": youden["MCC"],
    }
    (out / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

    lines = [
        "# Primary Direct H-DNA Sensitivity Analysis",
        "",
        "## Dataset",
        "",
        f"- Experimental sequences: {len(rows)} ({metadata['forming_count']} forming, {metadata['nonforming_count']} non-forming).",
        "- HSeeker and Triplex were run directly on each curated sequence. No genomic insertion or flanking sequence was used in this primary analysis.",
        "",
        "## Methods",
        "",
        "- HSeeker parameters: `minrep=8`, `maxrep=1000`, `maxspacer=10`, `purity=0.90`, `mismatch=0.10`, scoring enabled, exact `purity_rmq=True`.",
        "- Triplex was run with installed package defaults via `triplex.search(DNAString(seq))`.",
        "- HSeeker score discrimination was evaluated by ROC/AUC, and the classification threshold was selected by Youden's J.",
        "",
        "## Runtime",
        "",
        f"- HSeeker direct runtime: {h_runtime:.6f} sec.",
        f"- Triplex direct runtime: {t_runtime:.6f} sec.",
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
        "## Secondary Analysis",
        "",
        "- The E. coli injected-genome benchmark is retained separately in `analysis/results/sensitivity_injected/` as a secondary/context validation.",
        "",
        "## Output Files",
        "",
        "- `hseeker_direct_predictions.csv`: best HSeeker hit/score per sequence.",
        "- `triplex_direct_predictions.csv`: best Triplex hit/score per sequence.",
        "- `method_comparison_by_sequence.csv`: per-sequence HSeeker/Triplex comparison.",
        "- `summary_metrics.csv`: classification metrics.",
        "- `hseeker_score_roc.csv`: ROC thresholds and metrics.",
        "- `plots/hseeker_score_roc.pdf`: ROC plot.",
    ])
    (out / "summary.md").write_text("\n".join(lines) + "\n")

    print("Done.")
    print(f"Summary: {out / 'summary.md'}")
    print(f"Metrics: {out / 'summary_metrics.csv'}")


if __name__ == "__main__":
    main()
