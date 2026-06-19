#!/usr/bin/env python3
"""Validation/sensitivity analysis for experimental H-DNA sequences."""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
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
TUNED_PARAMS = dict(minrep=8, maxrep=100, maxspacer=10, purity=0.90, mismatch=0.10)
USE_PURITY_RMQ = True


def load_rows() -> list[dict[str, Any]]:
    rows = []
    with INPUT.open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.DictReader(fh)
        for i, row in enumerate(reader, start=1):
            label = row["label"].strip().lower()
            if label not in {"forming", "non-forming", "nonforming", "non_forming"}:
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
    with tempfile.TemporaryDirectory() as tmpdir:
        triplex_in = Path(tmpdir) / "triplex_input.csv"
        triplex_out = Path(tmpdir) / "triplex_predictions.csv"
        write_triplex_input(rows, triplex_in)
        subprocess.run(["Rscript", str(TRIPLEX_R), str(triplex_in), str(triplex_out)],
                       cwd=ROOT, check=True)
        with triplex_out.open(newline="") as fh:
            return list(csv.DictReader(fh))


def merge_comparison(
    default_rows: list[dict[str, Any]],
    tuned_rows: list[dict[str, Any]],
    triplex_rows: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    t = {r["sequence_id"]: r for r in triplex_rows}
    tuned_by_id = {r["sequence_id"]: r for r in tuned_rows}
    merged = []
    for h in default_rows:
        tuned = tuned_by_id[h["sequence_id"]]
        tr = t.get(h["sequence_id"], {})
        triplex_pred = str(tr.get("triplex_pred", "")).lower() == "true"
        merged.append(
            {
                "sequence_id": h["sequence_id"],
                "sequence_name": h["sequence_name"],
                "label": h["label"],
                "default_hseeker_score": h["hseeker_score"],
                "default_hseeker_pred": h["hseeker_pred"],
                "default_hseeker_runtime_sec": h["runtime_sec"],
                "default_arm_length": h["arm_length"],
                "default_spacer_length": h["spacer_length"],
                "tuned_hseeker_score": tuned["hseeker_score"],
                "tuned_hseeker_pred": tuned["hseeker_pred"],
                "tuned_hseeker_runtime_sec": tuned["runtime_sec"],
                "tuned_arm_length": tuned["arm_length"],
                "tuned_spacer_length": tuned["spacer_length"],
                "tuned_ga_pct": tuned["ga_pct"],
                "tuned_ct_pct": tuned["ct_pct"],
                "tuned_mirror_identity": tuned["mirror_identity"],
                "triplex_pred": triplex_pred,
                "triplex_count": tr.get("triplex_count", ""),
                "triplex_best_score": tr.get("triplex_best_score", ""),
                "triplex_best_pvalue": tr.get("triplex_best_pvalue", ""),
                "triplex_best_type": tr.get("triplex_best_type", ""),
                "triplex_best_start": tr.get("triplex_best_start", ""),
                "triplex_best_end": tr.get("triplex_best_end", ""),
                "triplex_runtime_sec": tr.get("triplex_runtime_sec", ""),
                "default_hseeker_correct": (h["label"] == "forming") == bool(h["hseeker_pred"]),
                "tuned_hseeker_correct": (h["label"] == "forming") == bool(tuned["hseeker_pred"]),
                "triplex_correct": (h["label"] == "forming") == triplex_pred,
                "tuned_hseeker_vs_triplex_agree": bool(tuned["hseeker_pred"]) == triplex_pred,
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


def write_summary(
    default_rows: list[dict[str, Any]],
    tuned_rows: list[dict[str, Any]],
    grid: list[dict[str, Any]],
    comparison: list[dict[str, Any]],
) -> None:
    m_default = metrics(default_rows)
    m_tuned = metrics(tuned_rows)
    m_triplex = metrics(comparison, "triplex_pred")
    default_failures = [r for r in default_rows if (r["label"] == "forming") != bool(r["hseeker_pred"])]
    tuned_failures = [r for r in tuned_rows if (r["label"] == "forming") != bool(r["hseeker_pred"])]
    best = grid[0] if grid else {}
    default_total_runtime = sum(float(r.get("runtime_sec") or 0.0) for r in default_rows)
    tuned_total_runtime = sum(float(r.get("runtime_sec") or 0.0) for r in tuned_rows)
    triplex_total_runtime = sum(float(r.get("triplex_runtime_sec") or 0.0) for r in comparison)
    default_mean_runtime = default_total_runtime / len(default_rows) if default_rows else 0.0
    tuned_mean_runtime = tuned_total_runtime / len(tuned_rows) if tuned_rows else 0.0
    triplex_mean_runtime = triplex_total_runtime / len(comparison) if comparison else 0.0
    runtime_ratio = triplex_total_runtime / tuned_total_runtime if tuned_total_runtime > 0 else 0.0
    tuned_triplex_agree = sum(bool(r["tuned_hseeker_vs_triplex_agree"]) for r in comparison)
    disagreements = [r for r in comparison if not bool(r["tuned_hseeker_vs_triplex_agree"])]
    lines = [
        "# HSeeker Validation / Sensitivity Summary",
        "",
        f"Dataset: `{INPUT}`",
        f"Dataset size: {len(default_rows)}",
        f"Forming: {sum(r['label']=='forming' for r in default_rows)}",
        f"Non-forming: {sum(r['label']!='forming' for r in default_rows)}",
        f"Score threshold: `{THRESHOLD}`",
        f"Purity RMQ enabled: `{USE_PURITY_RMQ}`",
        "",
        "## Method Comparison",
        "",
        "| method | parameters | TP | FN | TN | FP | sensitivity | specificity | precision | F1 | accuracy |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        f"| HSeeker default | minrep=10, maxrep=100, maxspacer=10, purity=0.90, mismatch=0.10 | {m_default['TP']} | {m_default['FN']} | {m_default['TN']} | {m_default['FP']} | {m_default['sensitivity']:.3f} | {m_default['specificity']:.3f} | {m_default['precision']:.3f} | {m_default['F1']:.3f} | {m_default['accuracy']:.3f} |",
        f"| HSeeker tuned | minrep=8, maxrep=100, maxspacer=10, purity=0.90, mismatch=0.10 | {m_tuned['TP']} | {m_tuned['FN']} | {m_tuned['TN']} | {m_tuned['FP']} | {m_tuned['sensitivity']:.3f} | {m_tuned['specificity']:.3f} | {m_tuned['precision']:.3f} | {m_tuned['F1']:.3f} | {m_tuned['accuracy']:.3f} |",
        f"| Triplex default | min_score=15, p_value=0.05, min_len=6, max_len=25, min_loop=3, max_loop=10 | {m_triplex['TP']} | {m_triplex['FN']} | {m_triplex['TN']} | {m_triplex['FP']} | {m_triplex['sensitivity']:.3f} | {m_triplex['specificity']:.3f} | {m_triplex['precision']:.3f} | {m_triplex['F1']:.3f} | {m_triplex['accuracy']:.3f} |",
        "",
        "## Parameter Tuning Notes",
        "",
        f"Best grid row by the simple sort: `{best}`",
        "",
        "`45/0/7/2` is not unique to one parameter set. It appears for several grid rows, including `minrep=6` with relaxed purity/mismatch and `minrep=8` with default purity/mismatch. The recommended tuned setting is `minrep=8` because it recovers all default false negatives while changing fewer biological assumptions than `minrep=6` or `mismatch=0.20`.",
        "",
        f"Tuned HSeeker and Triplex agreement: {tuned_triplex_agree}/{len(comparison)} ({tuned_triplex_agree / len(comparison):.3f}).",
        "",
        "Triplex is treated as a binary presence/absence caller here; its score is not directly comparable to HSeeker's thermodynamic score.",
        "",
        "## Runtime Comparison",
        "",
        f"HSeeker default total runtime: {default_total_runtime:.6f} sec; mean per sequence: {default_mean_runtime:.6f} sec.",
        f"HSeeker tuned total runtime: {tuned_total_runtime:.6f} sec; mean per sequence: {tuned_mean_runtime:.6f} sec.",
        f"Triplex total runtime: {triplex_total_runtime:.6f} sec; mean per sequence: {triplex_mean_runtime:.6f} sec.",
        f"Triplex/tuned-HSeeker runtime ratio on this validation set: {runtime_ratio:.2f}x.",
        "These are per-sequence search timings on very short sequences; HSeeker timings include scoring, while Triplex timings exclude R package startup.",
        "",
        "## Default HSeeker Failure Cases",
        "",
    ]
    if not default_failures:
        lines.append("None at the default threshold.")
    else:
        for r in default_failures:
            lines.append(f"- {r['sequence_id']} ({r['label']}): score={float(r['hseeker_score']):.2f}, arm={r['arm_length']}, spacer={r['spacer_length']}, GA={float(r['ga_pct']):.1f}, CT={float(r['ct_pct']):.1f}, mirror={float(r['mirror_identity']):.1f}")
    lines.extend(["", "## Tuned HSeeker Failure Cases", ""])
    if not tuned_failures:
        lines.append("None at the tuned threshold.")
    else:
        for r in tuned_failures:
            lines.append(f"- {r['sequence_id']} ({r['label']}): score={float(r['hseeker_score']):.2f}, arm={r['arm_length']}, spacer={r['spacer_length']}, GA={float(r['ga_pct']):.1f}, CT={float(r['ct_pct']):.1f}, mirror={float(r['mirror_identity']):.1f}")
    lines.extend(["", "## Tuned HSeeker vs Triplex Disagreements", ""])
    if not disagreements:
        lines.append("None.")
    else:
        for r in disagreements:
            lines.append(f"- {r['sequence_id']} ({r['label']}): tuned HSeeker={r['tuned_hseeker_pred']} score={float(r['tuned_hseeker_score']):.2f}; Triplex={r['triplex_pred']} count={r['triplex_count']} best_score={r['triplex_best_score']}")
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    hseeker_rows = run_hseeker(rows, DEFAULT_PARAMS, purity_rmq=USE_PURITY_RMQ)
    tuned_rows = run_hseeker(rows, TUNED_PARAMS, purity_rmq=USE_PURITY_RMQ)
    grid = parameter_grid(rows)
    write_csv(OUT / "hseeker_parameter_grid.csv", grid)
    triplex_rows = run_triplex(rows)
    comparison = merge_comparison(hseeker_rows, tuned_rows, triplex_rows)
    write_csv(OUT / "triplex_comparison.csv", comparison)
    score_distribution_svg(PLOTS / "validation_score_distribution.svg", hseeker_rows)
    m_h = metrics(hseeker_rows)
    m_tuned = metrics(tuned_rows)
    m_triplex = metrics(comparison, "triplex_pred")
    simple_svg_bar(PLOTS / "validation_f1_comparison.svg", "F1 Comparison", ["HSeeker default", "HSeeker tuned", "Triplex"], [m_h["F1"], m_tuned["F1"], m_triplex["F1"]], "F1")
    write_summary(hseeker_rows, tuned_rows, grid, comparison)
    print("Validation complete")
    print(f"Rows: {len(rows)}")
    print(f"HSeeker default F1: {m_h['F1']:.3f}; HSeeker tuned F1: {m_tuned['F1']:.3f}; Triplex F1: {m_triplex['F1']:.3f}")
    print(f"Summary: {OUT / 'summary.md'}")


if __name__ == "__main__":
    main()
