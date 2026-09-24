#!/usr/bin/env python3
"""Parameter robustness and out-of-sample parameter selection on benchmark v3.

Answers reviewer R3.9 (parameter robustness + component ablation) without
reintroducing the circularity of R3.3: the parameter grid is searched on a TRAIN
split only, and the selected configuration is scored once on a HELD-OUT EVAL split
that took no part in selection.

Design (fixed by the authors, 2026-09-23):
  * Data      : hdna_benchmark_balanced_v3.csv, the 65:65 balanced set only.
  * Splitting : SEQUENCE level, stratified by (label, record_type). Explicitly NOT
                family- or study-aware -- see the leakage audit this script emits,
                and the caveat it prints, because reviewer R3.2 asks for the
                family-aware design that this deliberately does not use.
  * Selection : highest MCC on TRAIN, with a deterministic tie-break.
  * Reporting : the selected config on EVAL, against two reference configs, plus
                the full grid surface and a repeated-split stability analysis.

Everything is deterministic given --seed.

Outputs (written to --outdir):
  grid_train_metrics.csv      every grid cell scored on TRAIN (selection surface)
  grid_eval_metrics.csv       every grid cell scored on EVAL  (reported for shape only)
  selected_vs_reference.csv   selected config vs reference configs on EVAL, with CIs
  split_manifest.csv          which record went to which split
  leakage_audit.csv           eval records having a near-duplicate or source in train
  repeated_split_stability.csv  selection procedure repeated over many random splits
  metadata.json               inputs, checksums, parameters, provenance
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
import random
import sys
import time
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

import hseeker  # noqa: E402

DEFAULT_CSV = ROOT / "hdna_benchmark_balanced_v3.csv"
DEFAULT_OUT = ROOT / "analysis" / "results" / "robustness_v3"

# ---------------------------------------------------------------- grid
PURITY = [0.85, 0.90, 0.95, 1.0]
MISMATCH = [0.0, 0.05, 0.10, 0.15]
MINREP = [8, 10, 12]
MAXSPACER = [5, 10, 15]
# Composition-filter states double as the R3.9 ablation axis
# ("sequence detection alone; detection plus composition filtering").
FILTERS = {
    "none": {"at_threshold": None, "filter_homopolymers": False},
    "at": {"at_threshold": 0.8, "filter_homopolymers": False},
    "homopolymer": {"at_threshold": None, "filter_homopolymers": True},
    "both": {"at_threshold": 0.8, "filter_homopolymers": True},
}

# Reference configurations we must beat (or fail to beat) out of sample.
REFERENCES = {
    "script_default": dict(purity=0.90, mismatch=0.10, minrep=8, maxspacer=10, filters="none"),
    "production_cli": dict(purity=0.90, mismatch=0.10, minrep=10, maxspacer=10, filters="both"),
    "manuscript_stated": dict(purity=0.90, mismatch=0.15, minrep=8, maxspacer=10, filters="both"),
}


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_balanced(path: Path) -> list[dict[str, Any]]:
    with path.open(newline="", encoding="utf-8-sig") as fh:
        rows = [r for r in csv.DictReader(fh)]
    out = []
    for r in rows:
        if (r.get("curation_decision") or "").strip() != "kept":
            continue
        seq = (r["sequence_5to3"] or "").strip().upper()
        label = r["label"].strip().lower().replace("_", "-")
        if label not in {"forming", "non-forming"}:
            raise RuntimeError(f"unexpected label {r['label']!r} for {r['record_id']}")
        out.append({
            "record_id": r["record_id"],
            "sequence_name": r["sequence_name"],
            "sequence": seq,
            "label": label,
            "y": 1 if label == "forming" else 0,
            "record_type": r.get("record_type", ""),
            "family_id": r.get("family_id", ""),
            "study_id": r.get("study_id", ""),
            "source_record": r.get("source_record", ""),
        })
    return out


# ---------------------------------------------------------------- metrics
def confusion(y: list[int], p: list[int]) -> dict[str, int]:
    tp = sum(1 for a, b in zip(y, p) if a == 1 and b == 1)
    fn = sum(1 for a, b in zip(y, p) if a == 1 and b == 0)
    tn = sum(1 for a, b in zip(y, p) if a == 0 and b == 0)
    fp = sum(1 for a, b in zip(y, p) if a == 0 and b == 1)
    return {"TP": tp, "FN": fn, "TN": tn, "FP": fp}


def mcc_of(c: dict[str, int]) -> float:
    tp, fn, tn, fp = c["TP"], c["FN"], c["TN"], c["FP"]
    den = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return ((tp * tn) - (fp * fn)) / den if den else 0.0


def metrics_of(y: list[int], p: list[int]) -> dict[str, float]:
    c = confusion(y, p)
    tp, fn, tn, fp = c["TP"], c["FN"], c["TN"], c["FP"]
    sens = tp / (tp + fn) if tp + fn else 0.0
    spec = tn / (tn + fp) if tn + fp else 0.0
    prec = tp / (tp + fp) if tp + fp else 0.0
    npv = tn / (tn + fn) if tn + fn else 0.0
    f1 = 2 * prec * sens / (prec + sens) if prec + sens else 0.0
    return {**c, "n": len(y), "sensitivity": sens, "specificity": spec,
            "precision": prec, "NPV": npv, "F1": f1,
            "accuracy": (tp + tn) / len(y) if y else 0.0, "MCC": mcc_of(c)}


def wilson(k: int, n: int, z: float = 1.959963985) -> tuple[float, float]:
    """Wilson score interval for a binomial proportion (R2.7)."""
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return (max(0.0, centre - half), min(1.0, centre + half))


def bootstrap_mcc_ci(y: list[int], p: list[int], rng: random.Random,
                     n_boot: int = 2000) -> tuple[float, float]:
    """Percentile CI for MCC, resampling RECORDS (sequence-level, as specified)."""
    n = len(y)
    vals = []
    for _ in range(n_boot):
        idx = [rng.randrange(n) for _ in range(n)]
        yy = [y[i] for i in idx]
        pp = [p[i] for i in idx]
        if len(set(yy)) < 2:
            continue
        vals.append(mcc_of(confusion(yy, pp)))
    if not vals:
        return (float("nan"), float("nan"))
    vals.sort()
    lo = vals[int(0.025 * len(vals))]
    hi = vals[min(len(vals) - 1, int(0.975 * len(vals)))]
    return (lo, hi)


# ---------------------------------------------------------------- prediction matrix
def grid_configs() -> list[dict[str, Any]]:
    cfgs = []
    for f in FILTERS:
        for pu in PURITY:
            for mm in MISMATCH:
                for mr in MINREP:
                    for ms in MAXSPACER:
                        cfgs.append({"purity": pu, "mismatch": mm, "minrep": mr,
                                     "maxspacer": ms, "filters": f})
    return cfgs


def config_key(c: dict[str, Any]) -> tuple:
    return (c["filters"], c["purity"], c["mismatch"], c["minrep"], c["maxspacer"])


def predict_all(records: list[dict], configs: list[dict]) -> dict[tuple, list[int]]:
    """Predict every record under every config ONCE.

    Splitting afterwards is pure bookkeeping, so train/eval and the repeated-split
    analysis cost nothing extra and cannot accidentally use different predictions.
    """
    preds: dict[tuple, list[int]] = {}
    t0 = time.perf_counter()
    for i, c in enumerate(configs, 1):
        fk = FILTERS[c["filters"]]
        row = []
        for r in records:
            hits = hseeker.scan_sequence(
                r["sequence"], minrep=c["minrep"], maxrep=1000,
                maxspacer=c["maxspacer"], purity=c["purity"], mismatch=c["mismatch"],
                remove_overlaps=True, seq_offset=1, score=True, purity_rmq=True,
                at_threshold=fk["at_threshold"],
                filter_homopolymers=fk["filter_homopolymers"],
            )
            row.append(1 if hits else 0)
        preds[config_key(c)] = row
        if i % 100 == 0:
            print(f"  {i}/{len(configs)} configs  ({time.perf_counter()-t0:.1f}s)", flush=True)
    print(f"  prediction matrix complete: {len(configs)} configs x {len(records)} records "
          f"in {time.perf_counter()-t0:.1f}s")
    return preds


# ---------------------------------------------------------------- splitting
def stratified_split(records: list[dict], frac: float, seed: int) -> tuple[list[int], list[int]]:
    """Sequence-level stratified split on (label, record_type).

    Per-label train sizes are computed from the label total so the EVAL split stays
    exactly class-balanced, then divided across record_type in proportion.
    """
    rng = random.Random(seed)
    by_label: dict[int, dict[str, list[int]]] = {}
    for i, r in enumerate(records):
        by_label.setdefault(r["y"], {}).setdefault(r["record_type"], []).append(i)

    train: list[int] = []
    eval_: list[int] = []
    for y, by_type in sorted(by_label.items()):
        n_label = sum(len(v) for v in by_type.values())
        n_train_label = int(frac * n_label)
        allocated = 0
        types = sorted(by_type)
        for j, t in enumerate(types):
            idx = by_type[t][:]
            rng.shuffle(idx)
            if j == len(types) - 1:
                k = n_train_label - allocated
            else:
                k = round(n_train_label * len(idx) / n_label)
            k = max(0, min(len(idx), k))
            allocated += k
            train.extend(idx[:k])
            eval_.extend(idx[k:])
    return sorted(train), sorted(eval_)


def select_config(records, preds, configs, idx) -> dict[str, Any]:
    """Pick the best configuration on the given indices.

    Rule, fixed in advance and applied identically everywhere: maximise MCC; break
    ties by higher sensitivity, then by the *simpler* model -- fewer filters, then
    larger minrep, then smaller maxspacer, then lower mismatch. The tie-break
    prefers simplicity so a tie never silently selects a more elaborate model.
    """
    filt_rank = {"none": 0, "at": 1, "homopolymer": 1, "both": 2}
    y = [records[i]["y"] for i in idx]
    best = None
    for c in configs:
        p = [preds[config_key(c)][i] for i in idx]
        m = metrics_of(y, p)
        key = (m["MCC"], m["sensitivity"], -filt_rank[c["filters"]],
               c["minrep"], -c["maxspacer"], -c["mismatch"])
        if best is None or key > best[0]:
            best = (key, c, m)
    return {"config": best[1], "train_metrics": best[2]}


# ---------------------------------------------------------------- leakage audit
def hamming_ratio(a: str, b: str) -> float:
    if len(a) != len(b):
        return 0.0
    same = sum(1 for x, y in zip(a, b) if x == y)
    return same / len(a)


def leakage_audit(records, train_idx, eval_idx) -> list[dict[str, Any]]:
    """Quantify what a sequence-level split leaks that a family-aware one would not.

    Reviewer R3.2 asks for family/study-aware validation. This design is
    sequence-level by explicit instruction, so the cost is measured rather than
    assumed: for every EVAL record, report whether TRAIN contains a record of the
    same family or study, or the synthetic record's own source, or a near-identical
    sequence (>=90% identity at equal length).
    """
    tr = [records[i] for i in train_idx]
    tr_ids = {r["record_id"] for r in tr}
    tr_fam = {r["family_id"] for r in tr if r["family_id"]}
    tr_study = {r["study_id"] for r in tr if r["study_id"]}
    out = []
    for i in eval_idx:
        r = records[i]
        near = ""
        best = 0.0
        for t in tr:
            h = hamming_ratio(r["sequence"], t["sequence"])
            if h > best:
                best, near = h, t["record_id"]
        src_in_train = r["source_record"] in tr_ids if r["source_record"] else False
        out.append({
            "record_id": r["record_id"],
            "label": r["label"],
            "record_type": r["record_type"],
            "family_in_train": r["family_id"] in tr_fam,
            "study_in_train": r["study_id"] in tr_study,
            "synthetic_source_in_train": src_in_train,
            "nearest_train_record": near,
            "nearest_train_identity": round(best, 4),
            "near_duplicate_ge_0.90": best >= 0.90,
        })
    return out


def stratified_kfold(records: list[dict], k: int, seed: int) -> list[list[int]]:
    """Partition record indices into k folds, stratified on (label, record_type).

    Each stratum is shuffled and dealt round-robin across folds, so every fold
    carries the class and experimental/synthetic mix of the whole set.
    """
    rng = random.Random(seed)
    folds: list[list[int]] = [[] for _ in range(k)]
    strata: dict[tuple, list[int]] = {}
    for i, r in enumerate(records):
        strata.setdefault((r["y"], r["record_type"]), []).append(i)
    for key in sorted(strata):
        idx = strata[key][:]
        rng.shuffle(idx)
        for j, i in enumerate(idx):
            folds[j % k].append(i)
    return [sorted(f) for f in folds]


def nested_cv(records, preds, configs, k: int, repeats: int, seed: int):
    """Repeated stratified k-fold CV with per-fold parameter selection.

    For each outer fold: select the configuration on the OTHER folds only, then
    predict the held-out fold with it. Pooling the held-out predictions gives one
    out-of-fold confusion matrix per repeat covering every record exactly once.

    The design (k, repeats, stratification, selection rule) is fixed in advance and
    is NOT tuned against the resulting metrics -- choosing a partition because it
    flatters the result would reintroduce, at the split level, exactly the
    circularity reviewer R3.3 objects to.
    """
    per_repeat: list[dict[str, Any]] = []
    selections: list[dict[str, Any]] = []
    ref_per_repeat: dict[str, list[dict[str, Any]]] = {n: [] for n in REFERENCES}
    for rep in range(repeats):
        folds = stratified_kfold(records, k, seed + 7919 * rep)
        oof_pred = [None] * len(records)
        ref_oof = {n: [None] * len(records) for n in REFERENCES}
        for fi, test_idx in enumerate(folds):
            train_idx = [i for j, f in enumerate(folds) if j != fi for i in f]
            sel = select_config(records, preds, configs, train_idx)["config"]
            selections.append({"repeat": rep, "fold": fi, "config": sel})
            pk = preds[config_key(sel)]
            for i in test_idx:
                oof_pred[i] = pk[i]
            for n, c in REFERENCES.items():
                rk = preds[config_key(c)]
                for i in test_idx:
                    ref_oof[n][i] = rk[i]
        y = [r["y"] for r in records]
        m = metrics_of(y, oof_pred)
        per_repeat.append({"repeat": rep, **m})
        for n in REFERENCES:
            ref_per_repeat[n].append(metrics_of(y, ref_oof[n]))
        if rep == 0:
            pooled_first = {"NESTED-CV (per-fold tuned)": list(oof_pred),
                            **{f"fixed: {n}": list(ref_oof[n]) for n in REFERENCES}}
    return per_repeat, selections, ref_per_repeat, pooled_first


def write_csv(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=Path, default=DEFAULT_CSV)
    ap.add_argument("--outdir", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--train-frac", type=float, default=0.7)
    ap.add_argument("--seed", type=int, default=20260923)
    ap.add_argument("--repeats", type=int, default=200,
                    help="random splits for the holdout stability analysis")
    ap.add_argument("--folds", type=int, default=5,
                    help="outer folds for the nested-CV primary estimate")
    ap.add_argument("--cv-repeats", type=int, default=20,
                    help="times the whole k-fold CV is repeated")
    args = ap.parse_args()

    out = args.outdir
    out.mkdir(parents=True, exist_ok=True)

    records = load_balanced(args.csv)
    npos = sum(r["y"] for r in records)
    print(f"Loaded {len(records)} kept records: {npos} forming / {len(records)-npos} non-forming")
    if npos != len(records) - npos:
        print("  WARNING: input is not class-balanced; this analysis assumes the 65:65 set")

    configs = grid_configs()
    print(f"Grid: {len(configs)} configurations "
          f"({len(FILTERS)} filter states x {len(PURITY)} purity x {len(MISMATCH)} mismatch "
          f"x {len(MINREP)} minrep x {len(MAXSPACER)} maxspacer)")
    preds = predict_all(records, configs)

    train_idx, eval_idx = stratified_split(records, args.train_frac, args.seed)
    ytr = [records[i]["y"] for i in train_idx]
    yev = [records[i]["y"] for i in eval_idx]
    print(f"\nSplit (seed={args.seed}, train_frac={args.train_frac}): "
          f"train n={len(train_idx)} ({sum(ytr)} pos / {len(ytr)-sum(ytr)} neg), "
          f"eval n={len(eval_idx)} ({sum(yev)} pos / {len(yev)-sum(yev)} neg)")

    write_csv(out / "split_manifest.csv", [
        {"record_id": records[i]["record_id"], "sequence_name": records[i]["sequence_name"],
         "label": records[i]["label"], "record_type": records[i]["record_type"],
         "family_id": records[i]["family_id"],
         "split": ("train" if i in set(train_idx) else "eval")}
        for i in range(len(records))
    ])

    # ---- grid surfaces
    grid_train, grid_eval = [], []
    for c in configs:
        pk = preds[config_key(c)]
        mtr = metrics_of(ytr, [pk[i] for i in train_idx])
        mev = metrics_of(yev, [pk[i] for i in eval_idx])
        grid_train.append({**c, **mtr})
        grid_eval.append({**c, **mev})
    write_csv(out / "grid_train_metrics.csv", grid_train)
    write_csv(out / "grid_eval_metrics.csv", grid_eval)

    # ---- selection on TRAIN only
    sel = select_config(records, preds, configs, train_idx)
    cfg, mtr = sel["config"], sel["train_metrics"]
    print("\nSelected on TRAIN (eval untouched):")
    print(f"  {cfg}")
    print(f"  train MCC={mtr['MCC']:.4f} sens={mtr['sensitivity']:.4f} spec={mtr['specificity']:.4f}")

    rng = random.Random(args.seed)
    rows = []
    for name, c in [("SELECTED (train-tuned)", cfg)] + [(k, v) for k, v in REFERENCES.items()]:
        pk = preds[config_key(c)]
        pev = [pk[i] for i in eval_idx]
        m = metrics_of(yev, pev)
        sl, sh = wilson(m["TP"], m["TP"] + m["FN"])
        pl, ph = wilson(m["TN"], m["TN"] + m["FP"])
        ml, mh = bootstrap_mcc_ci(yev, pev, rng)
        rows.append({
            "config_name": name, "purity": c["purity"], "mismatch": c["mismatch"],
            "minrep": c["minrep"], "maxspacer": c["maxspacer"], "filters": c["filters"],
            **{k: m[k] for k in ("n", "TP", "FN", "TN", "FP")},
            "sensitivity": round(m["sensitivity"], 4),
            "sens_CI95_lo": round(sl, 4), "sens_CI95_hi": round(sh, 4),
            "specificity": round(m["specificity"], 4),
            "spec_CI95_lo": round(pl, 4), "spec_CI95_hi": round(ph, 4),
            "precision": round(m["precision"], 4), "NPV": round(m["NPV"], 4),
            "F1": round(m["F1"], 4), "accuracy": round(m["accuracy"], 4),
            "MCC": round(m["MCC"], 4),
            "MCC_CI95_lo": round(ml, 4), "MCC_CI95_hi": round(mh, 4),
        })
    write_csv(out / "selected_vs_reference.csv", rows)
    print("\nHELD-OUT EVAL:")
    for r in rows:
        print(f"  {r['config_name']:24s} sens={r['sensitivity']:.3f} spec={r['specificity']:.3f} "
              f"MCC={r['MCC']:.3f} [{r['MCC_CI95_lo']:.3f},{r['MCC_CI95_hi']:.3f}]")

    # ---- leakage audit
    la = leakage_audit(records, train_idx, eval_idx)
    write_csv(out / "leakage_audit.csv", la)
    n_fam = sum(1 for r in la if r["family_in_train"])
    n_src = sum(1 for r in la if r["synthetic_source_in_train"])
    n_dup = sum(1 for r in la if r["near_duplicate_ge_0.90"])
    print(f"\nLeakage audit (sequence-level split, as specified):")
    print(f"  eval records whose FAMILY also appears in train : {n_fam}/{len(la)}")
    print(f"  synthetic eval records whose SOURCE is in train : {n_src}/{len(la)}")
    print(f"  eval records with a >=90% identical train record: {n_dup}/{len(la)}")

    # ---- repeated-split stability of the SELECTION PROCEDURE
    stab = []
    for k in range(args.repeats):
        s = args.seed + 1000 + k
        tr_i, ev_i = stratified_split(records, args.train_frac, s)
        sc = select_config(records, preds, configs, tr_i)["config"]
        yv = [records[i]["y"] for i in ev_i]
        pv = [preds[config_key(sc)][i] for i in ev_i]
        m = metrics_of(yv, pv)
        ref = metrics_of(yv, [preds[config_key(REFERENCES["script_default"])][i] for i in ev_i])
        stab.append({"split_seed": s, **{f"sel_{k2}": v for k2, v in sc.items()},
                     "sel_eval_MCC": round(m["MCC"], 4),
                     "sel_eval_sens": round(m["sensitivity"], 4),
                     "sel_eval_spec": round(m["specificity"], 4),
                     "default_eval_MCC": round(ref["MCC"], 4),
                     "delta_MCC_sel_minus_default": round(m["MCC"] - ref["MCC"], 4)})
    write_csv(out / "repeated_split_stability.csv", stab)
    deltas = sorted(r["delta_MCC_sel_minus_default"] for r in stab)
    sel_mccs = sorted(r["sel_eval_MCC"] for r in stab)
    wins = sum(1 for d in deltas if d > 0)
    ties = sum(1 for d in deltas if d == 0)

    def pct(v, q):
        return v[min(len(v) - 1, int(q * len(v)))]

    print(f"\nRepeated-split stability ({args.repeats} random splits):")
    print(f"  selected-config eval MCC : median {pct(sel_mccs,0.5):.3f} "
          f"[{pct(sel_mccs,0.025):.3f}, {pct(sel_mccs,0.975):.3f}]")
    print(f"  tuned - default eval MCC : median {pct(deltas,0.5):+.4f} "
          f"[{pct(deltas,0.025):+.4f}, {pct(deltas,0.975):+.4f}]")
    print(f"  tuning beat default in {wins}/{len(deltas)} splits, tied in {ties}")

    # ---- PRIMARY ESTIMATE: repeated stratified k-fold nested CV
    print(f"\nNested CV: {args.folds}-fold stratified, {args.cv_repeats} repeats "
          f"(every record held out once per repeat)")
    per_rep, sels, ref_rep, pooled = nested_cv(records, preds, configs,
                                               args.folds, args.cv_repeats, args.seed)
    write_csv(out / "nested_cv_per_repeat.csv", per_rep)
    write_csv(out / "nested_cv_selections.csv",
              [{"repeat": s["repeat"], "fold": s["fold"], **s["config"]} for s in sels])

    def summarise(rows_, field):
        v = sorted(r[field] for r in rows_)
        return (v[len(v) // 2], v[max(0, int(0.025 * len(v)))],
                v[min(len(v) - 1, int(0.975 * len(v)))])

    cv_rows = []
    for name, rows_ in [("NESTED-CV (per-fold tuned)", per_rep)] + \
                       [(f"fixed: {n}", ref_rep[n]) for n in REFERENCES]:
        med, lo, hi = summarise(rows_, "MCC")
        s_med, *_ = summarise(rows_, "sensitivity")
        p_med, *_ = summarise(rows_, "specificity")
        # Record-level uncertainty on the pooled out-of-fold predictions (R2.7).
        # Distinct from the across-repeat spread above, which for a FIXED config is
        # zero by construction: every record gets the same call whatever fold it is in.
        yall = [r["y"] for r in records]
        pp = pooled[name]
        brng = random.Random(args.seed + 31)
        bl, bh = bootstrap_mcc_ci(yall, pp, brng)
        mm = metrics_of(yall, pp)
        sl, sh = wilson(mm["TP"], mm["TP"] + mm["FN"])
        kl, kh = wilson(mm["TN"], mm["TN"] + mm["FP"])
        cv_rows.append({"config_name": name, "n_records": len(records),
                        "repeats": len(rows_),
                        "MCC_median": round(med, 4), "MCC_p2.5": round(lo, 4),
                        "MCC_p97.5": round(hi, 4),
                        "sensitivity_median": round(s_med, 4),
                        "specificity_median": round(p_med, 4),
                        "TP": mm["TP"], "FN": mm["FN"], "TN": mm["TN"], "FP": mm["FP"],
                        "MCC_boot_lo": round(bl, 4), "MCC_boot_hi": round(bh, 4),
                        "sens_CI95_lo": round(sl, 4), "sens_CI95_hi": round(sh, 4),
                        "spec_CI95_lo": round(kl, 4), "spec_CI95_hi": round(kh, 4)})
    write_csv(out / "nested_cv_summary.csv", cv_rows)
    print("  out-of-fold metrics pooled over all 130 records, across repeats:")
    for r in cv_rows:
        print(f"    {r['config_name']:28s} MCC {r['MCC_median']:.3f} "
              f"boot95[{r['MCC_boot_lo']:.3f},{r['MCC_boot_hi']:.3f}]  "
              f"sens {r['sensitivity_median']:.3f} spec {r['specificity_median']:.3f}  "
              f"(TP{r['TP']} FN{r['FN']} TN{r['TN']} FP{r['FP']})")

    # ---- which parameter values does selection actually prefer?
    freq_rows = []
    total = len(sels)
    for axis in ("purity", "mismatch", "minrep", "maxspacer", "filters"):
        counts: dict[Any, int] = {}
        for s in sels:
            counts[s["config"][axis]] = counts.get(s["config"][axis], 0) + 1
        for val, n in sorted(counts.items(), key=lambda kv: -kv[1]):
            freq_rows.append({"axis": axis, "value": val, "times_selected": n,
                              "fraction_of_folds": round(n / total, 4)})
    write_csv(out / "nested_cv_selection_frequency.csv", freq_rows)
    print(f"  selection frequency across {total} outer folds (top value per axis):")
    for axis in ("purity", "mismatch", "minrep", "maxspacer", "filters"):
        top = max((r for r in freq_rows if r["axis"] == axis),
                  key=lambda r: r["times_selected"])
        print(f"    {axis:10s} -> {str(top['value']):12s} in {top['fraction_of_folds']*100:.0f}% of folds")

    meta = {
        "date": time.strftime("%Y-%m-%d %H:%M:%S"),
        "status": "PROVISIONAL - benchmark not frozen (WS0-A1) and DEC-5 unresolved",
        "input_csv": str(args.csv),
        "input_csv_sha256": sha256_of(args.csv),
        "n_records": len(records), "n_forming": npos, "n_non_forming": len(records) - npos,
        "split": {"level": "sequence (NOT family/study aware, by author instruction)",
                  "stratified_by": ["label", "record_type"],
                  "train_frac": args.train_frac, "seed": args.seed,
                  "n_train": len(train_idx), "n_eval": len(eval_idx)},
        "grid": {"purity": PURITY, "mismatch": MISMATCH, "minrep": MINREP,
                 "maxspacer": MAXSPACER, "filter_states": list(FILTERS),
                 "n_configs": len(configs)},
        "selection_rule": "max MCC on TRAIN; ties -> higher sensitivity, then fewer "
                          "filters, then larger minrep, smaller maxspacer, lower mismatch",
        "selected_config": cfg,
        "selected_train_MCC": round(mtr["MCC"], 4),
        "leakage": {"eval_n": len(la), "family_in_train": n_fam,
                    "synthetic_source_in_train": n_src, "near_duplicate_ge_0.90": n_dup},
        "repeats": args.repeats,
        "nested_cv": {"folds": args.folds, "cv_repeats": args.cv_repeats,
                      "note": "PRIMARY estimate. Design fixed a priori; not tuned "
                              "against the resulting metrics.",
                      "MCC_median": cv_rows[0]["MCC_median"],
                      "MCC_p2.5": cv_rows[0]["MCC_p2.5"],
                      "MCC_p97.5": cv_rows[0]["MCC_p97.5"]},
        "hseeker_version": getattr(hseeker, "__version__", "unknown"),
    }
    (out / "metadata.json").write_text(json.dumps(meta, indent=2) + "\n")
    print(f"\nWrote outputs to {out}")


if __name__ == "__main__":
    main()
