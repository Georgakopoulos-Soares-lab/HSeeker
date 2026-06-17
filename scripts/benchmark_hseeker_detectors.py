#!/usr/bin/env python3
"""Primary FastGatedHSeeker benchmark suite.

Runs the experimental validation CSV through:

- original HSeeker
- FastGatedHSeeker exact-gate
- FastGatedHSeeker fast 4x
- FastGatedHSeeker fast 8x

It also runs the biological-parameter grid and writes the canonical
FastGated baseline outputs under ``results/fast_gated/``.
"""

from __future__ import annotations

import csv
import math
import os
import random
import statistics
import struct
import subprocess
import sys
import time
import zlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
SRC = ROOT / "src"
OUT = ROOT / "results" / "fast_gated"
PLOTS = OUT / "plots"
INPUT = ROOT / "hdna_experimental_sequences_final.csv"
SCORE_THRESHOLD = 60.0
BASE_PARAMS = {
    "minrep": 10,
    "maxrep": 1000,
    "maxspacer": 10,
    "purity": 0.90,
    "mismatch": 0.10,
}
GRID = {
    "minrep": [6, 8, 10, 12],
    "maxspacer": [3, 5, 7, 10, 20, 50],
    "purity": [0.85, 0.90, 0.95],
    "mismatch": [0.00, 0.05, 0.10, 0.15],
}


sys.path.insert(0, str(SRC))
import hseeker  # noqa: E402


def clean_sequence(seq: str) -> str:
    return "".join(ch for ch in seq.upper() if ch.isalpha())


def norm_label(label: str) -> str:
    label = label.strip().lower().replace("_", "-")
    if label in {"forming", "positive", "yes"}:
        return "forming"
    if label in {"non-forming", "nonforming", "negative", "no"}:
        return "non-forming"
    raise ValueError(f"Unknown label: {label!r}")


def load_records() -> list[dict]:
    with INPUT.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        rows = []
        for row in reader:
            keys = {k.lstrip("\ufeff"): v for k, v in row.items()}
            rows.append({
                "sequence_id": keys.get("record_id") or keys.get("sequence_id"),
                "sequence_name": keys.get("sequence_name", ""),
                "sequence": clean_sequence(keys.get("sequence_5to3") or keys.get("sequence", "")),
                "label": norm_label(keys.get("label", "")),
            })
        return rows


def best_hit(hits: list[dict]) -> dict | None:
    scored = [h for h in hits if h.get("total_score") is not None]
    if not scored:
        return None
    return max(
        scored,
        key=lambda h: (
            float(h.get("total_score") or 0.0),
            int(h.get("arm_length") or 0),
            float(h.get("mirror_identity") or 0.0),
        ),
    )


def method_config(method: str, params: dict) -> dict:
    if method == "original":
        return {"detector": "original"}
    if method == "fast_gated_exact_gate":
        return {
            "detector": "fast_gated",
            "fast_mode": False,
            "gate_window": params["minrep"],
            "gate_search_limit": params["maxrep"],
            "use_purity_prefilter": False,
        }
    if method == "fast_gated_fast_4x":
        return {
            "detector": "fast_gated",
            "fast_mode": True,
            "gate_window": params["minrep"],
            "gate_search_limit": 4 * params["minrep"],
            "use_purity_prefilter": True,
        }
    if method == "fast_gated_fast_8x":
        return {
            "detector": "fast_gated",
            "fast_mode": True,
            "gate_window": params["minrep"],
            "gate_search_limit": 8 * params["minrep"],
            "use_purity_prefilter": True,
        }
    raise ValueError(method)


def run_one(record: dict, method: str, params: dict, gate_limit: int | None = None) -> dict:
    kwargs = dict(params)
    cfg = method_config(method, params)
    if gate_limit is not None and method.startswith("fast_gated"):
        cfg["gate_search_limit"] = gate_limit
        cfg["fast_mode"] = gate_limit != params["maxrep"]
        cfg["use_purity_prefilter"] = gate_limit != params["maxrep"]
    start = time.perf_counter()
    hits = hseeker.search(
        record["sequence"],
        **kwargs,
        **cfg,
        remove_overlaps=True,
        score=True,
    )
    runtime = time.perf_counter() - start
    prof = hseeker.profiling_info()
    best = best_hit(hits)
    score = float(best.get("total_score", 0.0)) if best else 0.0
    pred = score >= SCORE_THRESHOLD
    return {
        "sequence_id": record["sequence_id"],
        "sequence_name": record["sequence_name"],
        "label": record["label"],
        "method": method,
        "minrep": params["minrep"],
        "maxrep": params["maxrep"],
        "maxspacer": params["maxspacer"],
        "purity": params["purity"],
        "mismatch": params["mismatch"],
        "gate_window": cfg.get("gate_window", ""),
        "gate_search_limit": cfg.get("gate_search_limit", ""),
        "best_score": round(score, 6),
        "predicted_forming_at_60": pred,
        "best_arm_length": best.get("arm_length", "") if best else "",
        "best_spacer_length": best.get("spacer_length", "") if best else "",
        "best_ga_pct": round(float(best.get("ga_pct", 0.0)), 6) if best else "",
        "best_ct_pct": round(float(best.get("ct_pct", 0.0)), 6) if best else "",
        "best_mirror_identity": round(float(best.get("mirror_identity", 0.0)), 6) if best else "",
        "num_raw_hits": prof.get("raw_hits", 0),
        "num_final_hits": prof.get("final_hits", len(hits)),
        "inner_iterations": prof.get("inner_iters", 0),
        "ctr_sp_pairs": prof.get("ctr_sp_pairs", 0),
        "skipped_invalid_kmax": prof.get("skipped_invalid_kmax", 0),
        "skipped_purity_prefilter": prof.get("skipped_purity_prefilter", 0),
        "gate_passed": prof.get("gate_passed", 0),
        "runtime_sec": runtime,
        "correct": (pred and record["label"] == "forming") or ((not pred) and record["label"] == "non-forming"),
    }


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def metrics(rows: list[dict], method: str | None = None) -> dict:
    if method:
        rows = [r for r in rows if r["method"] == method]
    tp = sum(r["label"] == "forming" and boolish(r["predicted_forming_at_60"]) for r in rows)
    fn = sum(r["label"] == "forming" and not boolish(r["predicted_forming_at_60"]) for r in rows)
    tn = sum(r["label"] == "non-forming" and not boolish(r["predicted_forming_at_60"]) for r in rows)
    fp = sum(r["label"] == "non-forming" and boolish(r["predicted_forming_at_60"]) for r in rows)
    sens = tp / (tp + fn) if tp + fn else 0.0
    spec = tn / (tn + fp) if tn + fp else 0.0
    prec = tp / (tp + fp) if tp + fp else 0.0
    f1 = 2 * sens * prec / (sens + prec) if sens + prec else 0.0
    acc = (tp + tn) / (tp + tn + fp + fn) if tp + tn + fp + fn else 0.0
    return {
        "TP": tp, "FN": fn, "TN": tn, "FP": fp,
        "sensitivity": sens, "specificity": spec,
        "precision": prec, "F1": f1, "accuracy": acc,
    }


def boolish(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).lower() in {"true", "1", "yes"}


def runtime_summary(rows: list[dict]) -> list[dict]:
    out = []
    for method in sorted({r["method"] for r in rows}):
        subset = [r for r in rows if r["method"] == method]
        runtimes = [float(r["runtime_sec"]) for r in subset]
        out.append({
            "method": method,
            "total_runtime_sec": sum(runtimes),
            "mean_runtime_per_sequence": statistics.mean(runtimes),
            "median_runtime_per_sequence": statistics.median(runtimes),
            "total_inner_iterations": sum(int(r["inner_iterations"]) for r in subset),
            "total_ctr_sp_pairs": sum(int(r["ctr_sp_pairs"]) for r in subset),
            "total_raw_hits": sum(int(r["num_raw_hits"]) for r in subset),
            "total_final_hits": sum(int(r["num_final_hits"]) for r in subset),
        })
    return out


def method_summary(rows: list[dict]) -> list[dict]:
    runtime = {r["method"]: r for r in runtime_summary(rows)}
    out = []
    for method in sorted({r["method"] for r in rows}):
        row = {"method": method}
        row.update(metrics(rows, method))
        row.update(runtime[method])
        out.append(row)
    return sorted(out, key=lambda r: (-r["F1"], r["total_runtime_sec"]))


def run_triplex(records: list[dict]) -> dict[str, dict]:
    helper = ROOT / "run_triplex_benchmark.R"
    if not (ROOT / ".r-lib" / "triplex").exists() or not helper.exists():
        return {}
    subprocess.run(["Rscript", str(helper)], cwd=ROOT, check=True)
    raw = ROOT / "results" / "triplex_predictions_raw.csv"
    with raw.open(newline="", encoding="utf-8") as fh:
        return {r["sequence_id"]: r for r in csv.DictReader(fh)}


def triplex_comparison(primary_rows: list[dict], records: list[dict]) -> list[dict]:
    triplex = run_triplex(records)
    by_id_method = {(r["sequence_id"], r["method"]): r for r in primary_rows}
    rows = []
    for rec in records:
        sid = rec["sequence_id"]
        orig = by_id_method[(sid, "original")]
        fast = by_id_method[(sid, "fast_gated_fast_4x")]
        t = triplex.get(sid, {})
        t_pred = t.get("triplex_pred", "not_run")
        rows.append({
            "sequence_id": sid,
            "label": rec["label"],
            "hseeker_original_score": orig["best_score"],
            "fast_gated_score": fast["best_score"],
            "triplex_call": t.get("triplex_n_hits", ""),
            "hseeker_original_pred": boolish(orig["predicted_forming_at_60"]),
            "fast_gated_pred": boolish(fast["predicted_forming_at_60"]),
            "triplex_pred": t_pred,
            "hseeker_original_correct": orig["correct"],
            "fast_gated_correct": fast["correct"],
            "triplex_correct": (t_pred == rec["label"]) if t_pred in {"forming", "non-forming"} else "",
        })
    return rows


def png_write(path: Path, width: int, height: int, pixels: list[list[tuple[int, int, int]]]) -> None:
    def chunk(tag: bytes, data: bytes) -> bytes:
        return struct.pack("!I", len(data)) + tag + data + struct.pack("!I", zlib.crc32(tag + data) & 0xffffffff)
    raw = b"".join(b"\x00" + bytes(c for px in row for c in px) for row in pixels)
    data = b"\x89PNG\r\n\x1a\n"
    data += chunk(b"IHDR", struct.pack("!IIBBBBB", width, height, 8, 2, 0, 0, 0))
    data += chunk(b"IDAT", zlib.compress(raw, 9))
    data += chunk(b"IEND", b"")
    path.write_bytes(data)


def blank(w: int, h: int, color=(255, 255, 255)):
    return [[color for _ in range(w)] for _ in range(h)]


def rect(img, x0, y0, x1, y1, color):
    h, w = len(img), len(img[0])
    for y in range(max(0, y0), min(h, y1)):
        row = img[y]
        for x in range(max(0, x0), min(w, x1)):
            row[x] = color


def plot_runtime(summary: list[dict]) -> None:
    img = blank(900, 420)
    vals = [float(r["total_runtime_sec"]) for r in summary]
    maxv = max(vals) if vals else 1.0
    colors = [(52, 120, 180), (90, 160, 90), (220, 150, 60), (180, 80, 80)]
    for i, row in enumerate(summary):
        x0 = 80 + i * 190
        bar_h = int(320 * float(row["total_runtime_sec"]) / maxv)
        rect(img, x0, 360 - bar_h, x0 + 110, 360, colors[i % len(colors)])
    png_write(PLOTS / "runtime_original_vs_fast.png", 900, 420, img)


def plot_sens_spec(summary: list[dict]) -> None:
    img = blank(900, 420)
    for i, row in enumerate(summary):
        x0 = 80 + i * 190
        sens_h = int(320 * float(row["sensitivity"]))
        spec_h = int(320 * float(row["specificity"]))
        rect(img, x0, 360 - sens_h, x0 + 50, 360, (40, 140, 90))
        rect(img, x0 + 60, 360 - spec_h, x0 + 110, 360, (160, 70, 70))
    png_write(PLOTS / "sensitivity_specificity_by_method.png", 900, 420, img)


def plot_score_distribution(rows: list[dict]) -> None:
    img = blank(900, 420)
    methods = ["original", "fast_gated_fast_4x"]
    colors = {"forming": (40, 140, 90), "non-forming": (180, 70, 70)}
    for mi, method in enumerate(methods):
        subset = [r for r in rows if r["method"] == method]
        max_score = max([float(r["best_score"]) for r in subset] + [60.0])
        for idx, row in enumerate(subset):
            x = 80 + mi * 390 + int(280 * min(float(row["best_score"]), max_score) / max_score)
            y = 70 + (idx * 6) % 280
            rect(img, x - 3, y - 3, x + 4, y + 4, colors[row["label"]])
        threshold_x = 80 + mi * 390 + int(280 * 60.0 / max_score)
        rect(img, threshold_x, 50, threshold_x + 2, 370, (80, 80, 80))
    png_write(PLOTS / "score_distribution_original_vs_fast.png", 900, 420, img)


def plot_f1_heatmap(grid_rows: list[dict]) -> None:
    img = blank(900, 420)
    mm = GRID["mismatch"]
    sp = GRID["maxspacer"]
    best = {}
    for row in grid_rows:
        key = (float(row["mismatch"]), int(row["maxspacer"]))
        best[key] = max(best.get(key, 0.0), float(row["F1"]))
    for yi, m in enumerate(mm):
        for xi, s in enumerate(sp):
            v = best.get((m, s), 0.0)
            color = (int(245 - 120 * v), int(245 - 30 * v), int(245 - 155 * v))
            rect(img, 90 + xi * 120, 70 + yi * 70, 195 + xi * 120, 125 + yi * 70, color)
    png_write(PLOTS / "f1_heatmap_mismatch_maxspacer.png", 900, 420, img)


def run_stress_runtime() -> list[dict]:
    rng = random.Random(17)
    seq = "".join(rng.choice("ACGT") for _ in range(20000))
    params = dict(minrep=10, maxrep=1000, maxspacer=50, purity=0.90, mismatch=0.10)
    configs = [
        ("original", dict(detector="original")),
        ("fast_gated_fast_4x", dict(
            detector="fast_gated",
            gate_window=params["minrep"],
            gate_search_limit=4 * params["minrep"],
            fast_mode=True,
            use_purity_prefilter=True,
        )),
        ("fast_gated_exact_gate", dict(
            detector="fast_gated",
            gate_window=params["minrep"],
            gate_search_limit=params["maxrep"],
            fast_mode=False,
            use_purity_prefilter=False,
        )),
    ]
    rows = []
    for method, cfg in configs:
        start = time.perf_counter()
        hits = hseeker.search(seq, **params, **cfg, remove_overlaps=True, score=False)
        elapsed = time.perf_counter() - start
        prof = hseeker.profiling_info()
        rows.append({
            "method": method,
            "sequence_length": len(seq),
            "minrep": params["minrep"],
            "maxrep": params["maxrep"],
            "maxspacer": params["maxspacer"],
            "purity": params["purity"],
            "mismatch": params["mismatch"],
            "runtime_sec": elapsed,
            "inner_iterations": prof.get("inner_iters", 0),
            "ctr_sp_pairs": prof.get("ctr_sp_pairs", 0),
            "skipped_invalid_kmax": prof.get("skipped_invalid_kmax", 0),
            "skipped_purity_prefilter": prof.get("skipped_purity_prefilter", 0),
            "gate_passed": prof.get("gate_passed", 0),
            "raw_hits": prof.get("raw_hits", 0),
            "final_hits": prof.get("final_hits", len(hits)),
        })
    return rows


def write_summary(records, primary_rows, grid_rows, summary_rows, runtime_rows, triplex_rows, stress_rows) -> None:
    forming = sum(r["label"] == "forming" for r in records)
    nonforming = len(records) - forming
    orig = next(r for r in summary_rows if r["method"] == "original")
    fast4 = next(r for r in summary_rows if r["method"] == "fast_gated_fast_4x")
    exact = next(r for r in summary_rows if r["method"] == "fast_gated_exact_gate")
    speedup = float(orig["total_runtime_sec"]) / float(fast4["total_runtime_sec"]) if fast4["total_runtime_sec"] else 0.0
    stress_by_method = {r["method"]: r for r in stress_rows}
    stress_speedup = 0.0
    if stress_by_method.get("fast_gated_fast_4x", {}).get("runtime_sec"):
        stress_speedup = (
            float(stress_by_method["original"]["runtime_sec"]) /
            float(stress_by_method["fast_gated_fast_4x"]["runtime_sec"])
        )
    lost = []
    by_id = {(r["sequence_id"], r["method"]): r for r in primary_rows}
    for rec in records:
        o = by_id[(rec["sequence_id"], "original")]
        f = by_id[(rec["sequence_id"], "fast_gated_fast_4x")]
        if boolish(o["predicted_forming_at_60"]) and not boolish(f["predicted_forming_at_60"]):
            lost.append(f"- {rec['sequence_id']} {rec['sequence_name']}: original score {o['best_score']}, fast score {f['best_score']}")
    if not lost:
        lost.append("- None at score threshold 60.")
    failures = [
        f"- {r['sequence_id']} {r['method']} label={r['label']} score={r['best_score']}"
        for r in primary_rows if not boolish(r["correct"])
    ]
    triplex_done = any(r["triplex_pred"] in {"forming", "non-forming"} for r in triplex_rows)
    text = f"""# FastGatedHSeeker Benchmark

## Dataset

- Total sequences: {len(records)}
- Forming: {forming}
- Non-forming: {nonforming}
- Score threshold: {SCORE_THRESHOLD:g}

## Primary Metrics

| Method | TP | FN | TN | FP | Sensitivity | Specificity | Precision | F1 | Accuracy | Runtime sec | Inner iterations |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
"""
    for r in summary_rows:
        text += (
            f"| {r['method']} | {r['TP']} | {r['FN']} | {r['TN']} | {r['FP']} | "
            f"{float(r['sensitivity']):.3f} | {float(r['specificity']):.3f} | "
            f"{float(r['precision']):.3f} | {float(r['F1']):.3f} | {float(r['accuracy']):.3f} | "
            f"{float(r['total_runtime_sec']):.6f} | {r['total_inner_iterations']} |\n"
        )
    text += f"""
## Speedup

- Fast 4x gate vs original wall-clock speedup on this dataset: {speedup:.2f}x
- Fast 4x gate vs original on synthetic 20 kb balanced-DNA stress input: {stress_speedup:.2f}x
- Exact-gate F1: {float(exact['F1']):.3f}
- Fast 4x F1: {float(fast4['F1']):.3f}

## Stress Runtime

| Method | Runtime sec | Inner iterations | Center/spacer pairs | Gate passed | Raw hits | Final hits |
|---|---:|---:|---:|---:|---:|---:|
"""
    for r in stress_rows:
        text += (
            f"| {r['method']} | {float(r['runtime_sec']):.6f} | {r['inner_iterations']} | "
            f"{r['ctr_sp_pairs']} | {r['gate_passed']} | {r['raw_hits']} | {r['final_hits']} |\n"
        )
    text += f"""

## Lost Calls In Fast Mode

{chr(10).join(lost)}

## Failure Cases

{chr(10).join(failures) if failures else '- None.'}

## Triplex

- Triplex comparison run: {'yes' if triplex_done else 'no'}
- Output: `triplex_comparison.csv`

## Recommended Gate Defaults

- `gate_window = minrep`
- `gate_mirror_frac = 1 - mismatch`
- `gate_purity_frac = purity`
- Use `exact_gate` (`gate_search_limit=maxrep`, purity prefilter disabled) when recall/equivalence matters.
- Use `fast` with `gate_search_limit=4*minrep` for speed-sensitive scans, while reviewing `lost_calls` in this report.
"""
    (OUT / "summary.md").write_text(text, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)
    records = load_records()

    methods = ["original", "fast_gated_exact_gate", "fast_gated_fast_4x", "fast_gated_fast_8x"]
    primary_rows = []
    for method in methods:
        for rec in records:
            primary_rows.append(run_one(rec, method, BASE_PARAMS))
    write_csv(OUT / "per_sequence_predictions.csv", primary_rows)

    summary_rows = method_summary(primary_rows)
    write_csv(OUT / "method_summary_metrics.csv", summary_rows)
    write_csv(OUT / "runtime_summary.csv", runtime_summary(primary_rows))

    grid_rows = []
    for minrep in GRID["minrep"]:
        for maxspacer in GRID["maxspacer"]:
            for purity in GRID["purity"]:
                for mismatch in GRID["mismatch"]:
                    params = {
                        "minrep": minrep,
                        "maxrep": BASE_PARAMS["maxrep"],
                        "maxspacer": maxspacer,
                        "purity": purity,
                        "mismatch": mismatch,
                    }
                    for gate_limit_name, gate_limit in [
                        ("4x_minrep", 4 * minrep),
                        ("8x_minrep", 8 * minrep),
                        ("maxrep", params["maxrep"]),
                    ]:
                        rows = [
                            run_one(rec, "fast_gated_fast_4x", params, gate_limit=gate_limit)
                            for rec in records
                        ]
                        m = metrics(rows)
                        m.update({
                            "minrep": minrep,
                            "maxrep": params["maxrep"],
                            "maxspacer": maxspacer,
                            "purity": purity,
                            "mismatch": mismatch,
                            "gate_search_limit": gate_limit_name,
                            "gate_search_limit_value": gate_limit,
                            "total_runtime_sec": sum(float(r["runtime_sec"]) for r in rows),
                            "total_inner_iterations": sum(int(r["inner_iterations"]) for r in rows),
                            "total_ctr_sp_pairs": sum(int(r["ctr_sp_pairs"]) for r in rows),
                            "total_raw_hits": sum(int(r["num_raw_hits"]) for r in rows),
                            "total_final_hits": sum(int(r["num_final_hits"]) for r in rows),
                        })
                        grid_rows.append(m)
    grid_rows.sort(key=lambda r: (-r["F1"], -r["sensitivity"], -r["specificity"], r["total_runtime_sec"]))
    write_csv(OUT / "parameter_grid_results.csv", grid_rows)

    triplex_rows = triplex_comparison(primary_rows, records)
    write_csv(OUT / "triplex_comparison.csv", triplex_rows)

    stress_rows = run_stress_runtime()
    write_csv(OUT / "stress_runtime.csv", stress_rows)

    plot_score_distribution(primary_rows)
    plot_runtime(summary_rows)
    plot_sens_spec(summary_rows)
    plot_f1_heatmap(grid_rows)

    write_summary(records, primary_rows, grid_rows, summary_rows, runtime_summary(primary_rows), triplex_rows, stress_rows)
    print(f"Wrote FastGatedHSeeker benchmark outputs to {OUT}")


if __name__ == "__main__":
    main()
