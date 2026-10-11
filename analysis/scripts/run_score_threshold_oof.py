#!/usr/bin/env python3
"""Out-of-fold evaluation of a stability-score threshold (reviewer comment R3.3).

R3.3 objects that the published threshold was chosen with Youden's J on the same benchmark it was
evaluated on. Here the threshold is chosen inside each training fold only and applied to the
held-out records:
  - repeated stratified 5-fold CV (20 repeats),
  - leave-one-family-out,
  - leave-one-study-out.
Each record's score is the highest total_score of its hits at the shipped configuration (minimum arm
8, maxspacer 10, purity 0.90, mismatch 0.10, AT and homopolymer filters, greedy overlap removal);
a record without hits scores 0 and is called non-forming. Inside a training fold the threshold
maximises Youden's J over the observed scores (ties: the lowest threshold, i.e. the most inclusive
call). The comparison is the shipped rule (any hit, no threshold) and the published in-sample
threshold of 76 applied as a fixed rule.

Usage:  python3 analysis/scripts/run_score_threshold_oof.py
Output: analysis/results/score_threshold_oof/
"""

from __future__ import annotations

import csv
import hashlib
import json
import random
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(Path(__file__).resolve().parent))
import hseeker  # noqa: E402
import run_grouped_validation as gv  # noqa: E402
import run_parameter_robustness as rb  # noqa: E402

OUT = ROOT / "analysis" / "results" / "score_threshold_oof"
PARAMS = dict(minrep=8, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10,
              at_threshold=0.80, filter_homopolymers=True, score=True, overlap_strategy="greedy")
PUBLISHED_THRESHOLD = 76.0
SEED = 20260923


def record_scores(records: list[dict]) -> list[float]:
    out = []
    for r in records:
        s = [h["total_score"] for h in hseeker.scan_sequence(r["sequence"], **PARAMS)
             if h.get("total_score") is not None]
        out.append(max(s) if s else 0.0)
    return out


def youden_threshold(y: list[int], s: list[float]) -> float:
    """Most inclusive threshold maximising Youden's J; a score of 0 (no hit) is never 'forming'."""
    best, best_t = -2.0, None
    for t in sorted({x for x in s if x > 0}):
        tp = sum(1 for a, b in zip(y, s) if a == 1 and b >= t)
        tn = sum(1 for a, b in zip(y, s) if a == 0 and b < t)
        npos, nneg = sum(y), len(y) - sum(y)
        j = tp / npos + tn / nneg - 1
        if j > best + 1e-12:
            best, best_t = j, t
    return best_t if best_t is not None else float("inf")


def oof(records, scores, folds) -> tuple[list[int], list[float]]:
    pred, chosen = [None] * len(records), []
    for test in folds:
        test_set = set(test)
        train = [i for i in range(len(records)) if i not in test_set]
        t = youden_threshold([records[i]["y"] for i in train], [scores[i] for i in train])
        chosen.append(t)
        for i in test:
            pred[i] = int(scores[i] > 0 and scores[i] >= t)
    return pred, chosen


def auroc(y, s) -> float:
    pos = [b for a, b in zip(y, s) if a == 1]
    neg = [b for a, b in zip(y, s) if a == 0]
    w = sum((p > q) + 0.5 * (p == q) for p in pos for q in neg)
    return w / (len(pos) * len(neg))


def summary(name, y, p, rng, extra=None) -> dict:
    m = rb.metrics_of(y, p)
    lo, hi = rb.bootstrap_mcc_ci(y, p, rng)
    return dict(scheme=name, **{k: m[k] for k in ("TP", "FN", "TN", "FP")},
                sensitivity=round(m["sensitivity"], 4), specificity=round(m["specificity"], 4),
                MCC=round(m["MCC"], 4), MCC_ci95=[round(lo, 4), round(hi, 4)], **(extra or {}))


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = rb.DEFAULT_CSV
    records = rb.load_balanced(csv_path)
    y = [r["y"] for r in records]
    scores = record_scores(records)
    rng = random.Random(SEED)
    rows = [summary("shipped rule: any hit, no threshold", y, [int(s > 0) for s in scores], rng),
            summary(f"published threshold {PUBLISHED_THRESHOLD:g}, fixed", y,
                    [int(s >= PUBLISHED_THRESHOLD) for s in scores], rng),
            summary("in-sample Youden threshold (circular, for reference)", y,
                    [int(s > 0 and s >= youden_threshold(y, scores)) for s in scores], rng,
                    dict(threshold=youden_threshold(y, scores)))]
    # repeated stratified 5-fold
    pooled_y, pooled_p, ts = [], [], []
    for rep in range(20):
        folds = rb.stratified_kfold(records, 5, SEED + rep)
        p, chosen = oof(records, scores, folds)
        pooled_y += y
        pooled_p += p
        ts += chosen
    rows.append(summary("out-of-fold Youden threshold, 5-fold x 20", pooled_y, pooled_p, rng,
                        dict(thresholds_min_median_max=[min(ts), sorted(ts)[len(ts) // 2], max(ts)])))
    for level in ("family", "study"):
        groups = gv.assign_groups(records, level)
        folds = [[i for i, g in enumerate(groups) if g == u] for u in sorted(set(groups))]
        p, chosen = oof(records, scores, folds)
        rows.append(summary(f"out-of-fold Youden threshold, leave-one-{level}-out", y, p, rng,
                            dict(thresholds_min_median_max=[min(chosen), sorted(chosen)[len(chosen) // 2],
                                                            max(chosen)])))
    a = auroc(y, scores)
    with open(OUT / "score_threshold_oof.csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(dict.fromkeys(k for r in rows for k in r)))
        w.writeheader()
        w.writerows({k: (json.dumps(v) if isinstance(v, list) else v) for k, v in r.items()} for r in rows)
    with open(OUT / "record_scores.csv", "w", newline="") as fo:
        w = csv.writer(fo)
        w.writerow(["record_id", "label", "max_total_score"])
        for r, s in zip(records, scores):
            w.writerow([r["record_id"], r["label"], round(s, 3)])
    neg_with_hits = sum(1 for r, s in zip(records, scores) if r["y"] == 0 and s > 0)
    meta = dict(date=time.strftime("%Y-%m-%d %H:%M:%S"), params=PARAMS, input_csv=csv_path.name,
                input_csv_sha256=hashlib.sha256(csv_path.read_bytes()).hexdigest(),
                n_records=len(records), auroc_max_score=round(a, 4),
                negatives_with_any_hit=neg_with_hits,
                note="With no negative carrying any hit, the AUROC equals (1 + sensitivity) / 2 and "
                     "any threshold above the lowest positive score can only lose forming records.")
    (OUT / "metadata.json").write_text(json.dumps(meta, indent=1))
    for r in rows:
        print(r)
    print("AUROC", round(a, 4), "negatives with any hit", neg_with_hits)


if __name__ == "__main__":
    main()
