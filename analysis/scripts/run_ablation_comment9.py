#!/usr/bin/env python3
"""Component ablation and one-at-a-time parameter variation (reviewer comment R3.9).

Every stage runs on the same records with the same labels. A record is called "forming" when
at least one candidate survives the stage (no score threshold). Score-based stages also report
an exploratory AUROC of the record's highest total_score (records without a candidate score 0);
it is a within-cohort description, not an independent estimate of predictive accuracy.

Stages (cumulative, overlaps retained so that recall is compared fairly):
  detection         C mirror-repeat detector, purity 0, no filters
  + purity          detector purity threshold
  + AT filter       full-motif A/T fraction < 0.80
  pairing           + pairing score on the detected arms
  pairing_stacking  + stacking score
  complete          + arm-boundary optimisation and homopolymer filter (the shipped scorer)
Then the shipped pipeline with each overlap mode (retain all / greedy / score-based), and the
shipped pipeline with one parameter varied at a time. The joint grid of 576 configurations is
in analysis/results/robustness_v3/ (run_parameter_robustness.py).

Configuration: the shipped defaults (minrep 8, maxrep 1000, maxspacer 10, purity 0.90,
mismatch 0.10, AT filter 0.80, homopolymer filter, greedy overlap removal).

Usage:  python3 analysis/scripts/run_ablation_comment9.py [--csv hdna_benchmark_balanced_v3.csv]
Outputs: analysis/results/ablation_comment9/
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import hseeker  # noqa: E402
from hseeker import ablation  # noqa: E402

OUT = ROOT / "analysis" / "results" / "ablation_comment9"
BASE = dict(minrep=8, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10)
AT = 0.80
STAGES = ("detection", "+purity", "+AT_filter", "pairing", "pairing_stacking", "complete")
VARY = dict(purity=[0.80, 0.85, 0.90, 0.95, 1.00], mismatch=[0.0, 0.05, 0.10, 0.15, 0.20],
            minrep=[6, 7, 8, 9, 10, 11, 12], maxspacer=[3, 5, 8, 10, 12, 15])


def load(path: Path) -> list[dict]:
    rows = [r for r in csv.DictReader(open(path)) if r["curation_decision"] == "kept"]
    return [dict(record_id=r["record_id"], sequence_name=r["sequence_name"],
                 sequence=r["sequence_5to3"].upper(), label=r["label"],
                 record_type=r.get("record_type", "experimental")) for r in rows]


def stage_hits(seq: str, stage: str) -> list[dict]:
    if stage == "detection":
        return ablation.scan_sequence_ablation(seq, model="detection", **BASE)
    if stage == "+purity":
        return ablation.scan_sequence_ablation(seq, model="composition", at_threshold=None, **BASE)
    if stage == "+AT_filter":
        return ablation.scan_sequence_ablation(seq, model="composition", at_threshold=AT, **BASE)
    return ablation.scan_sequence_ablation(seq, model=stage, at_threshold=AT,
                                           filter_homopolymers_enabled=True, **BASE)



def wilson(k: int, n: int, z: float = 1.959963985) -> list[float]:
    if n == 0:
        return [math.nan, math.nan]
    p = k / n
    d = 1 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return [round(c - h, 4), round(c + h, 4)]


def summarise(rows: list[dict], scored: bool) -> dict:
    tp = sum(r["label"] == "forming" and r["call"] for r in rows)
    fn = sum(r["label"] == "forming" and not r["call"] for r in rows)
    tn = sum(r["label"] == "non_forming" and not r["call"] for r in rows)
    fp = sum(r["label"] == "non_forming" and r["call"] for r in rows)
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    out = dict(n=len(rows), TP=tp, FN=fn, TN=tn, FP=fp,
               sensitivity=round(tp / (tp + fn), 4) if tp + fn else None,
               sensitivity_ci95=wilson(tp, tp + fn),
               specificity=round(tn / (tn + fp), 4) if tn + fp else None,
               specificity_ci95=wilson(tn, tn + fp),
               MCC=round((tp * tn - fp * fn) / den, 4) if den else 0.0)
    if scored:
        pos = [r["max_score"] for r in rows if r["label"] == "forming"]
        neg = [r["max_score"] for r in rows if r["label"] == "non_forming"]
        wins = sum(p > q for p in pos for q in neg) + 0.5 * sum(p == q for p in pos for q in neg)
        out["AUROC_exploratory"] = round(wins / (len(pos) * len(neg)), 4) if pos and neg else None
    return out


def subsets(records):
    return {"balanced": records,
            "experimental": [r for r in records if r["record_type"] == "experimental"]}


def run_stages(records: list[dict]) -> tuple[list[dict], list[dict]]:
    calls, metrics = [], []
    for stage in STAGES:
        scored = stage in ("pairing", "pairing_stacking", "complete")
        per = []
        for r in records:
            hits = stage_hits(r["sequence"], stage)
            s = [h["total_score"] for h in hits if h.get("total_score") is not None]
            per.append(dict(r, stage=stage, candidates=len(hits), call=bool(hits),
                            max_score=max(s) if s else 0.0))
        calls.extend({k: v for k, v in p.items() if k != "sequence"} for p in per)
        for name, sub in subsets(per).items():
            metrics.append(dict(stage=stage, subset=name, **summarise(sub, scored)))
    return calls, metrics


def run_overlap_modes(records: list[dict]) -> list[dict]:
    out = []
    params = dict(BASE, at_threshold=AT, filter_homopolymers=True, score=True)
    for mode, kw in (("retain_all", dict(remove_overlaps=False)),
                     ("greedy", dict(overlap_strategy="greedy")),
                     ("score", dict(overlap_strategy="score"))):
        per = []
        for r in records:
            hits = hseeker.scan_sequence(r["sequence"], **kw, **params)
            s = [h["total_score"] for h in hits if h.get("total_score") is not None]
            per.append(dict(r, call=bool(hits), max_score=max(s) if s else 0.0))
        for name, sub in subsets(per).items():
            out.append(dict(overlap_mode=mode, subset=name, **summarise(sub, True)))
    return out


def run_vary(records: list[dict]) -> list[dict]:
    out = []
    for param, values in VARY.items():
        for v in values:
            p = dict(BASE, at_threshold=AT, filter_homopolymers=True, score=False)
            p[param] = v
            per = [dict(r, call=bool(hseeker.scan_sequence(r["sequence"], **p))) for r in records]
            for name, sub in subsets(per).items():
                out.append(dict(parameter=param, value=v, is_default=BASE[param] == v,
                                subset=name, **summarise(sub, False)))
    return out


def write_csv(path: Path, rows: list[dict]) -> None:
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(path, "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=keys)
        w.writeheader()
        for r in rows:
            w.writerow({k: (json.dumps(v) if isinstance(v, list) else v) for k, v in r.items()})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", type=Path, default=ROOT / "hdna_benchmark_balanced_v3.csv")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    records = load(args.csv)
    calls, stages = run_stages(records)
    overlap = run_overlap_modes(records)
    vary = run_vary(records)
    write_csv(OUT / "record_calls.csv", calls)
    write_csv(OUT / "stage_metrics.csv", stages)
    write_csv(OUT / "overlap_mode_metrics.csv", overlap)
    write_csv(OUT / "parameter_one_at_a_time.csv", vary)
    meta = dict(date=time.strftime("%Y-%m-%d %H:%M:%S"), base=BASE, at_threshold=AT,
                homopolymer_filter=True, input_csv=args.csv.name,
                input_csv_sha256=hashlib.sha256(args.csv.read_bytes()).hexdigest(),
                n_records=len(records), n_forming=sum(r["label"] == "forming" for r in records),
                git_commit=subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                                          capture_output=True, text=True).stdout.strip(),
                git_dirty=bool(subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain"],
                                              capture_output=True, text=True).stdout.strip()),
                notes=["AUROC values are within-cohort descriptions, not independent estimates.",
                       "The scorer exposes no separate mismatch-penalty switch; that penalty is "
                       "not ablated on its own.",
                       "Joint 576-configuration grid: analysis/results/robustness_v3/."])
    (OUT / "metadata.json").write_text(json.dumps(meta, indent=1))
    for r in stages:
        if r["subset"] == "balanced":
            print(r["stage"], {k: r[k] for k in ("TP", "FN", "TN", "FP", "MCC")}, r.get("AUROC_exploratory"))
    for r in overlap:
        if r["subset"] == "balanced":
            print(r["overlap_mode"], {k: r[k] for k in ("TP", "FN", "TN", "FP", "MCC")}, r.get("AUROC_exploratory"))


if __name__ == "__main__":
    main()
