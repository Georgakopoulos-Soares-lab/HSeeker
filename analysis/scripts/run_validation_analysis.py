#!/usr/bin/env python3
"""Validation/sensitivity analysis for experimental H-DNA sequences."""

from __future__ import annotations

import csv
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

import hseeker  # noqa: E402


INPUT = ROOT / "hdna_experimental_sequences_final.csv"
OUT = ROOT / "analysis" / "results" / "validation"
PLOTS = ROOT / "analysis" / "results" / "plots"
TRIPLEX_R = ROOT / "analysis" / "scripts" / "triplex_validate.R"
THRESHOLD = 60.0
DEFAULT_PARAMS = dict(minrep=10, maxrep=100, maxspacer=10, purity=0.90, mismatch=0.10)
USE_PURITY_RMQ = True


def load_rows() -> list[dict[str, Any]]:
    rows = []
    with INPUT.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        for i, row in enumerate(reader, start=1):
            label = row["label"].strip().lower()
            if label not in {"forming", "non-forming", "nonforming"}:
                continue
            rows.append(
                {
                    "sequence_id": row.get("record_id") or f"seq_{i}",
                    "sequence_name": row.get("sequence_name", ""),
                    "sequence": row["sequence_5to3"].strip().upper(),
                    "label": "forming" if label == "forming" else "non-forming",
                    "source": row.get("first_author", ""),
                    "paper_title": row.get("paper_title", ""),
                    "notes": row.get("notes", ""),
                }
            )
    return rows


def best_hit(hits: list[dict[str, Any]]) -> dict[str, Any]:
    if not hits:
        return {}
    return max(hits, key=lambda h: ((h.get("total_score") or 0.0), h["arm_length"], -h["spacer_length"]))


def run_hseeker(rows: list[dict[str, Any]], params: dict[str, Any], purity_rmq: bool = False) -> list[dict[str, Any]]:
    out = []
    for row in rows:
        t0 = time.perf_counter()
        hits = hseeker.scan_sequence(row["sequence"], score=True, purity_rmq=purity_rmq, **params)
        runtime = time.perf_counter() - t0
        best = best_hit(hits)
        score = float(best.get("total_score") or 0.0)
        out.append(
            {
                "sequence_id": row["sequence_id"],
                "sequence_name": row["sequence_name"],
                "label": row["label"],
                "sequence_length": len(row["sequence"]),
                "hseeker_score": score,
                "hseeker_pred": score >= THRESHOLD,
                "num_hits": len(hits),
                "runtime_sec": runtime,
                "arm_length": best.get("arm_length", 0),
                "spacer_length": best.get("spacer_length", 0),
                "ga_pct": best.get("ga_pct", 0.0),
                "ct_pct": best.get("ct_pct", 0.0),
                "mirror_identity": best.get("mirror_identity", 0.0),
                "left_arm": best.get("left_arm", ""),
                "spacer": best.get("spacer", ""),
                "right_arm": best.get("right_arm", ""),
                "full_sequence": best.get("full_sequence", ""),
            }
        )
    return out


def metrics(rows: list[dict[str, Any]], pred_key: str = "hseeker_pred") -> dict[str, Any]:
    tp = sum(r["label"] == "forming" and bool(r[pred_key]) for r in rows)
    fn = sum(r["label"] == "forming" and not bool(r[pred_key]) for r in rows)
    tn = sum(r["label"] != "forming" and not bool(r[pred_key]) for r in rows)
    fp = sum(r["label"] != "forming" and bool(r[pred_key]) for r in rows)
    sensitivity = tp / (tp + fn) if tp + fn else 0.0
    specificity = tn / (tn + fp) if tn + fp else 0.0
    precision = tp / (tp + fp) if tp + fp else 0.0
    f1 = 2 * precision * sensitivity / (precision + sensitivity) if precision + sensitivity else 0.0
    accuracy = (tp + tn) / len(rows) if rows else 0.0
    return dict(TP=tp, FN=fn, TN=tn, FP=fp, sensitivity=sensitivity, specificity=specificity,
                precision=precision, F1=f1, accuracy=accuracy, n=len(rows))


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    if not rows:
        path.write_text("")
        return
    keys = []
    seen = set()
    for row in rows:
        for key in row:
            if key not in seen:
                keys.append(key)
                seen.add(key)
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def parameter_grid(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    grid = []
    for minrep in (6, 8, 10, 12):
        for maxspacer in (3, 5, 7, 10, 20, 50):
            for purity in (0.85, 0.90, 0.95):
                for mismatch in (0.00, 0.05, 0.10, 0.15):
                    params = dict(minrep=minrep, maxrep=100, maxspacer=maxspacer, purity=purity, mismatch=mismatch)
                    preds = run_hseeker(rows, params, purity_rmq=USE_PURITY_RMQ)
                    m = metrics(preds)
                    grid.append({"minrep": minrep, "maxrep": 100, "maxspacer": maxspacer,
                                 "purity": purity, "mismatch": mismatch, **m})
    grid.sort(key=lambda r: (r["F1"], r["sensitivity"], r["specificity"]), reverse=True)
    return grid


def write_triplex_input(rows: list[dict[str, Any]], path: Path) -> None:
    write_csv(path, [{"sequence_id": r["sequence_id"], "label": r["label"], "sequence": r["sequence"]} for r in rows])


def run_triplex(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    triplex_in = OUT / "triplex_input.csv"
    triplex_out = OUT / "triplex_predictions.csv"
    write_triplex_input(rows, triplex_in)
    subprocess.run(["Rscript", str(TRIPLEX_R), str(triplex_in), str(triplex_out)],
                   cwd=ROOT, check=True)
    with triplex_out.open(newline="") as fh:
        return list(csv.DictReader(fh))


def merge_triplex(hseeker_rows: list[dict[str, Any]], triplex_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    t = {r["sequence_id"]: r for r in triplex_rows}
    merged = []
    for h in hseeker_rows:
        tr = t.get(h["sequence_id"], {})
        triplex_pred = str(tr.get("triplex_pred", "")).lower() == "true"
        merged.append(
            {
                "sequence_id": h["sequence_id"],
                "label": h["label"],
                "hseeker_score": h["hseeker_score"],
                "hseeker_pred": h["hseeker_pred"],
                "triplex_pred": triplex_pred,
                "triplex_count": tr.get("triplex_count", ""),
                "triplex_best_score": tr.get("triplex_best_score", ""),
                "triplex_best_pvalue": tr.get("triplex_best_pvalue", ""),
                "hseeker_correct": (h["label"] == "forming") == bool(h["hseeker_pred"]),
                "triplex_correct": (h["label"] == "forming") == triplex_pred,
            }
        )
    return merged


def simple_svg_bar(path: Path, title: str, labels: list[str], values: list[float], ylabel: str) -> None:
    width, height = 760, 420
    left, top, plot_w, plot_h = 70, 50, 640, 285
    maxv = max(values) if values else 1.0
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:20px;font-weight:bold}.small{font-size:12px}</style>',
             f'<text class="title" x="{left}" y="28">{title}</text>',
             f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333"/>',
             f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333"/>']
    bar_w = plot_w / max(1, len(labels)) * 0.7
    for i, (label, value) in enumerate(zip(labels, values)):
        x = left + (i + 0.15) * (plot_w / len(labels))
        h = 0 if maxv == 0 else value / maxv * plot_h
        y = top + plot_h - h
        parts.append(f'<rect x="{x:.1f}" y="{y:.1f}" width="{bar_w:.1f}" height="{h:.1f}" fill="#3b73b9"/>')
        parts.append(f'<text class="small" x="{x + bar_w/2:.1f}" y="{y-5:.1f}" text-anchor="middle">{value:.2f}</text>')
        parts.append(f'<text class="small" x="{x + bar_w/2:.1f}" y="{top+plot_h+18}" text-anchor="middle">{label}</text>')
    parts.append(f'<text class="small" transform="translate(20 {top+plot_h/2}) rotate(-90)" text-anchor="middle">{ylabel}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts))


def score_distribution_svg(path: Path, rows: list[dict[str, Any]]) -> None:
    forming = [float(r["hseeker_score"]) for r in rows if r["label"] == "forming"]
    non = [float(r["hseeker_score"]) for r in rows if r["label"] != "forming"]
    bins = [0, 20, 40, 60, 80, 100, 140, 180, 240]
    labels = [f"{bins[i]}-{bins[i+1]}" for i in range(len(bins)-1)]
    def counts(vals):
        out = []
        for lo, hi in zip(bins[:-1], bins[1:]):
            out.append(sum(lo <= v < hi for v in vals))
        return out
    width, height = 860, 440
    left, top, plot_w, plot_h = 70, 50, 720, 300
    c1, c2 = counts(forming), counts(non)
    maxv = max(c1 + c2 + [1])
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
             '<rect width="100%" height="100%" fill="white"/>',
             '<style>text{font-family:Arial,sans-serif;fill:#222}.title{font-size:20px;font-weight:bold}.small{font-size:12px}</style>',
             f'<text class="title" x="{left}" y="28">HSeeker Score Distribution</text>',
             f'<line x1="{left}" y1="{top+plot_h}" x2="{left+plot_w}" y2="{top+plot_h}" stroke="#333"/>',
             f'<line x1="{left}" y1="{top}" x2="{left}" y2="{top+plot_h}" stroke="#333"/>',
             f'<line x1="{left}" y1="{top}" x2="{left+plot_w}" y2="{top}" stroke="#ddd"/>']
    group_w = plot_w / len(labels)
    for i, label in enumerate(labels):
        x = left + i * group_w + 8
        h1 = c1[i] / maxv * plot_h
        h2 = c2[i] / maxv * plot_h
        parts.append(f'<rect x="{x:.1f}" y="{top+plot_h-h1:.1f}" width="{group_w*0.36:.1f}" height="{h1:.1f}" fill="#2e7d32"/>')
        parts.append(f'<rect x="{x+group_w*0.39:.1f}" y="{top+plot_h-h2:.1f}" width="{group_w*0.36:.1f}" height="{h2:.1f}" fill="#c62828"/>')
        parts.append(f'<text class="small" x="{x+group_w*0.36:.1f}" y="{top+plot_h+18}" text-anchor="middle">{label}</text>')
    tx = left + plot_w - 170
    parts.append(f'<rect x="{tx}" y="55" width="14" height="14" fill="#2e7d32"/><text class="small" x="{tx+20}" y="67">forming</text>')
    parts.append(f'<rect x="{tx}" y="78" width="14" height="14" fill="#c62828"/><text class="small" x="{tx+20}" y="90">non-forming</text>')
    parts.append(f'<line x1="{left + plot_w*(60/bins[-1]):.1f}" y1="{top}" x2="{left + plot_w*(60/bins[-1]):.1f}" y2="{top+plot_h}" stroke="#111" stroke-dasharray="4 4"/>')
    parts.append("</svg>")
    path.write_text("\n".join(parts))


def write_summary(rows: list[dict[str, Any]], grid: list[dict[str, Any]], comparison: list[dict[str, Any]]) -> None:
    m_h = metrics(rows)
    m_t = metrics([{**r, "triplex_pred": str(r["triplex_pred"]).lower() == "true"} for r in comparison], "triplex_pred")
    failures = [r for r in rows if (r["label"] == "forming") != bool(r["hseeker_pred"])]
    best = grid[0] if grid else {}
    lines = [
        "# HSeeker Validation / Sensitivity Summary",
        "",
        f"Dataset: `{INPUT}`",
        f"Dataset size: {len(rows)}",
        f"Forming: {sum(r['label']=='forming' for r in rows)}",
        f"Non-forming: {sum(r['label']!='forming' for r in rows)}",
        "",
        "## Default HSeeker",
        "",
        f"Score threshold: `{THRESHOLD}`",
        f"Purity RMQ enabled: `{USE_PURITY_RMQ}`",
        f"TP={m_h['TP']} FN={m_h['FN']} TN={m_h['TN']} FP={m_h['FP']}",
        f"Sensitivity={m_h['sensitivity']:.3f}, specificity={m_h['specificity']:.3f}, precision={m_h['precision']:.3f}, F1={m_h['F1']:.3f}, accuracy={m_h['accuracy']:.3f}",
        "",
        "## Best Grid Row",
        "",
        str(best),
        "",
        "## Triplex Presence/Absence",
        "",
        f"TP={m_t['TP']} FN={m_t['FN']} TN={m_t['TN']} FP={m_t['FP']}",
        f"Sensitivity={m_t['sensitivity']:.3f}, specificity={m_t['specificity']:.3f}, precision={m_t['precision']:.3f}, F1={m_t['F1']:.3f}, accuracy={m_t['accuracy']:.3f}",
        "",
        "Triplex is treated as a binary presence/absence caller here; its score is not directly comparable to HSeeker's thermodynamic score.",
        "",
        "## HSeeker Failure Cases",
        "",
    ]
    if not failures:
        lines.append("None at the default threshold.")
    else:
        for r in failures:
            lines.append(f"- {r['sequence_id']} ({r['label']}): score={float(r['hseeker_score']):.2f}, arm={r['arm_length']}, spacer={r['spacer_length']}, GA={float(r['ga_pct']):.1f}, CT={float(r['ct_pct']):.1f}, mirror={float(r['mirror_identity']):.1f}")
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    hseeker_rows = run_hseeker(rows, DEFAULT_PARAMS, purity_rmq=USE_PURITY_RMQ)
    write_csv(OUT / "hseeker_default_predictions.csv", hseeker_rows)
    grid = parameter_grid(rows)
    write_csv(OUT / "hseeker_parameter_grid.csv", grid)
    triplex_rows = run_triplex(rows)
    comparison = merge_triplex(hseeker_rows, triplex_rows)
    write_csv(OUT / "triplex_comparison.csv", comparison)
    write_csv(OUT / "validation_failure_inspection.csv", [r for r in hseeker_rows if (r["label"] == "forming") != bool(r["hseeker_pred"])])
    score_distribution_svg(PLOTS / "validation_score_distribution.svg", hseeker_rows)
    m_h = metrics(hseeker_rows)
    m_t = metrics([{**r, "triplex_pred": str(r["triplex_pred"]).lower() == "true"} for r in comparison], "triplex_pred")
    simple_svg_bar(PLOTS / "validation_f1_comparison.svg", "F1 Comparison", ["HSeeker", "Triplex"], [m_h["F1"], m_t["F1"]], "F1")
    write_summary(hseeker_rows, grid, comparison)
    print("Validation complete")
    print(f"Rows: {len(rows)}")
    print(f"HSeeker F1: {m_h['F1']:.3f}; Triplex F1: {m_t['F1']:.3f}")
    print(f"Summary: {OUT / 'summary.md'}")


if __name__ == "__main__":
    main()
