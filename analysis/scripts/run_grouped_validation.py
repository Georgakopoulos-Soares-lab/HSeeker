#!/usr/bin/env python3
"""Family- and study-aware validation of HSeeker on benchmark v3 (reviewer R3.2).

R3.2 asks, in the reviewer's own words, for "a family-aware validation strategy in
which closely related constructs from the same experimental series are grouped",
ideally "leave-one-family-out or leave-one-study-out validation", and "at minimum"
for results "recalculated after collapsing highly related sequence families".

This script implements all three:

  lofo       leave-one-FAMILY-out
  loso       leave-one-STUDY-out
  collapsed  one representative record per family, then plain evaluation

Grouping rule (the part that makes it honest):
  A synthetic negative generated FROM an experimental record inherits that record's
  family/study. A mirror-disrupted mutant of pGG32 carries pGG32's sequence, so
  leaving out the pGG32 family while training on its mutant would leak precisely
  what the family-aware design exists to prevent. Synthetic records with no source
  (G4, Z-DNA, B-DNA, homopolymer, perfect-mirror controls) form their own groups.

Parameters are selected on the retained groups only and the held-out group is
predicted with that selection; predictions are pooled across folds so one confusion
matrix covers every record. Fixed reference configurations are evaluated under the
identical scheme for comparison.

Usage:
    python3 analysis/scripts/run_grouped_validation.py
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import random
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]

_spec = importlib.util.spec_from_file_location(
    "rb", HERE / "run_parameter_robustness.py")
rb = importlib.util.module_from_spec(_spec)
sys.modules["rb"] = rb
_spec.loader.exec_module(rb)


def assign_groups(records: list[dict], level: str) -> list[str]:
    """Map every record to a validation group.

    level='family' -> family_id ; level='study' -> study_id.
    Synthetic records inherit the group of their source record when they have one.
    """
    key = "family_id" if level == "family" else "study_id"
    by_id = {r["record_id"]: r for r in records}
    groups = []
    for r in records:
        src = r.get("source_record") or ""
        if src and src in by_id:
            g = by_id[src][key] or f"ungrouped:{src}"
        else:
            g = r[key] or f"ungrouped:{r['record_id']}"
        groups.append(g)
    return groups


def leave_one_group_out(records, preds, configs, groups, seed: int):
    """Leave-one-group-out with per-fold parameter selection, pooled out-of-fold."""
    uniq = sorted(set(groups))
    oof = [None] * len(records)
    ref_oof = {n: [None] * len(records) for n in rb.REFERENCES}
    selections = []
    for g in uniq:
        test_idx = [i for i, gg in enumerate(groups) if gg == g]
        train_idx = [i for i, gg in enumerate(groups) if gg != g]
        ytr = {records[i]["y"] for i in train_idx}
        if len(ytr) < 2:
            raise RuntimeError(f"training portion for held-out group {g!r} is single-class")
        sel = rb.select_config(records, preds, configs, train_idx)["config"]
        selections.append({"group": g, "n_held_out": len(test_idx), **sel})
        pk = preds[rb.config_key(sel)]
        for i in test_idx:
            oof[i] = pk[i]
        for n, c in rb.REFERENCES.items():
            rk = preds[rb.config_key(c)]
            for i in test_idx:
                ref_oof[n][i] = rk[i]
    return oof, ref_oof, selections, uniq


def collapsed_families(records, groups, by_label: bool = False):
    """Reviewer's 'at minimum' option: one representative record per group.

    ``by_label=False`` takes a single representative per family. Because synthetic
    negatives inherit their source positive's family, that variant collapses a
    family's negatives away with it and leaves a mostly-positive set -- faithful to
    a literal reading, but it discards the balanced design.

    ``by_label=True`` takes one representative per (family, class), which removes
    within-family redundancy while keeping both classes represented. This is the
    more informative reading and is reported alongside.
    """
    seen: dict[Any, int] = {}
    for i, (r, g) in enumerate(zip(records, groups)):
        key = (g, r["y"]) if by_label else g
        if key not in seen:
            seen[key] = i
    return sorted(seen.values())


def cluster_bootstrap(y, p, groups, rng, n_boot: int = 4000):
    """Cluster (family/study) bootstrap CIs for sensitivity, specificity and MCC.

    Resamples whole GROUPS with replacement rather than individual records. This is
    the uncertainty R3.2 is really about: 16 pXY32 variants from one study are not
    16 independent observations, so record-level intervals understate the true
    spread. Reported alongside the record-level interval, never instead of it --
    the two answer different questions (R2.7 asks for the record-level interval on
    the observed sample; R3.2 asks what happens if the *families* had differed).
    """
    by_g: dict[str, list[int]] = {}
    for i, g in enumerate(groups):
        by_g.setdefault(g, []).append(i)
    keys = list(by_g)
    sens, spec, mccs = [], [], []
    for _ in range(n_boot):
        idx = [i for _ in keys for i in by_g[rng.choice(keys)]]
        yy = [y[i] for i in idx]
        pp = [p[i] for i in idx]
        if len(set(yy)) < 2:
            continue
        c = rb.confusion(yy, pp)
        if c["TP"] + c["FN"]:
            sens.append(c["TP"] / (c["TP"] + c["FN"]))
        if c["TN"] + c["FP"]:
            spec.append(c["TN"] / (c["TN"] + c["FP"]))
        mccs.append(rb.mcc_of(c))

    def ci(v):
        if not v:
            return (float("nan"), float("nan"))
        v = sorted(v)
        return (v[int(0.025 * len(v))], v[min(len(v) - 1, int(0.975 * len(v)))])

    return ci(sens), ci(spec), ci(mccs)


def report(name, y, p, groups, rng, extra=""):
    m = rb.metrics_of(y, p)
    # record-level (kept: this is the R2.7 interval on the observed sample)
    lo, hi = rb.bootstrap_mcc_ci(y, p, rng)
    sl, sh = rb.wilson(m["TP"], m["TP"] + m["FN"])
    kl, kh = rb.wilson(m["TN"], m["TN"] + m["FP"])
    # cluster-level (added: the R3.2 interval treating each family as the unit)
    crng = random.Random(20260923 + 77)
    (csl, csh), (ckl, ckh), (clo, chi) = cluster_bootstrap(y, p, groups, crng)
    print(f"  {name:34s} TP{m['TP']:3d} FN{m['FN']:3d} TN{m['TN']:3d} FP{m['FP']:3d} | "
          f"sens {m['sensitivity']:.3f} rec[{sl:.3f},{sh:.3f}] clu[{csl:.3f},{csh:.3f}] | "
          f"spec {m['specificity']:.3f} rec[{kl:.3f},{kh:.3f}] clu[{ckl:.3f},{ckh:.3f}] | "
          f"MCC {m['MCC']:.3f} rec[{lo:.3f},{hi:.3f}] clu[{clo:.3f},{chi:.3f}] {extra}")
    return {"scheme": name, **{k: m[k] for k in ("n", "TP", "FN", "TN", "FP")},
            "sensitivity": round(m["sensitivity"], 4),
            "sens_record_lo": round(sl, 4), "sens_record_hi": round(sh, 4),
            "sens_cluster_lo": round(csl, 4), "sens_cluster_hi": round(csh, 4),
            "specificity": round(m["specificity"], 4),
            "spec_record_lo": round(kl, 4), "spec_record_hi": round(kh, 4),
            "spec_cluster_lo": round(ckl, 4), "spec_cluster_hi": round(ckh, 4),
            "precision": round(m["precision"], 4), "NPV": round(m["NPV"], 4),
            "MCC": round(m["MCC"], 4),
            "MCC_record_lo": round(lo, 4), "MCC_record_hi": round(hi, 4),
            "MCC_cluster_lo": round(clo, 4), "MCC_cluster_hi": round(chi, 4),
            "MCC_record_width": round(hi - lo, 4),
            "MCC_cluster_width": round(chi - clo, 4),
            "cluster_over_record_width": round((chi - clo) / (hi - lo), 3) if hi > lo else None}


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--csv", type=Path, default=rb.DEFAULT_CSV)
    ap.add_argument("--outdir", type=Path,
                    default=ROOT / "analysis" / "results" / "grouped_validation_v3")
    ap.add_argument("--seed", type=int, default=20260923)
    args = ap.parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)

    records = rb.load_balanced(args.csv)
    y = [r["y"] for r in records]
    configs = rb.grid_configs()
    print(f"Records: {len(records)} ({sum(y)} forming / {len(y)-sum(y)} non-forming)")
    print(f"Grid: {len(configs)} configurations\n")
    preds = rb.predict_all(records, configs)

    rows_out: list[dict[str, Any]] = []
    sel_rows: list[dict[str, Any]] = []

    for level, label in (("family", "LEAVE-ONE-FAMILY-OUT"), ("study", "LEAVE-ONE-STUDY-OUT")):
        groups = assign_groups(records, level)
        uniq = sorted(set(groups))
        sizes = Counter(groups)
        inherited = sum(1 for r, g in zip(records, groups)
                        if r["record_type"] == "synthetic" and r.get("source_record"))
        print(f"\n===== {label} =====")
        print(f"  groups: {len(uniq)} | largest: {sizes.most_common(3)}")
        print(f"  singleton groups: {sum(1 for g,c in sizes.items() if c==1)}")
        print(f"  synthetic records inheriting their source's group: {inherited}")
        oof, ref_oof, sels, _ = leave_one_group_out(records, preds, configs, groups, args.seed)
        rng = random.Random(args.seed + 5)
        rows_out.append({"level": level, **report(f"{level}: per-fold tuned", y, oof, groups, rng)})
        for n in rb.REFERENCES:
            rng2 = random.Random(args.seed + 5)
            rows_out.append({"level": level,
                             **report(f"{level}: fixed {n}", y, ref_oof[n], groups, rng2)})
        freq = defaultdict(Counter)
        for s in sels:
            for axis in ("purity", "mismatch", "minrep", "maxspacer", "filters"):
                freq[axis][s[axis]] += 1
            sel_rows.append({"level": level, **s})
        print(f"  selection across {len(sels)} folds:")
        for axis in ("purity", "mismatch", "minrep", "maxspacer", "filters"):
            top, n = freq[axis].most_common(1)[0]
            print(f"    {axis:10s} -> {str(top):12s} {n}/{len(sels)} folds "
                  f"({', '.join(f'{k}:{v}' for k,v in freq[axis].most_common())})")

    # ---- reviewer's 'at minimum': collapse families, then evaluate fixed configs
    groups = assign_groups(records, "family")
    keep = collapsed_families(records, groups)
    yk = [records[i]["y"] for i in keep]
    print(f"\n===== COLLAPSED FAMILIES (reviewer's 'at minimum') =====")
    print(f"  {len(keep)} representative records ({sum(yk)} forming / {len(yk)-sum(yk)} non-forming)")
    for n, c in rb.REFERENCES.items():
        rng = random.Random(args.seed + 9)
        pk = preds[rb.config_key(c)]
        rows_out.append({"level": "collapsed",
                         **report(f"collapsed: fixed {n}", yk, [pk[i] for i in keep],
                                  [groups[i] for i in keep], rng)})

    keep2 = collapsed_families(records, groups, by_label=True)
    yk2 = [records[i]["y"] for i in keep2]
    print(f"\n===== COLLAPSED per (FAMILY, CLASS) — keeps both classes =====")
    print(f"  {len(keep2)} representatives ({sum(yk2)} forming / {len(yk2)-sum(yk2)} non-forming)")
    for n, c in rb.REFERENCES.items():
        rng = random.Random(args.seed + 13)
        pk = preds[rb.config_key(c)]
        rows_out.append({"level": "collapsed_by_class",
                         **report(f"collapsed2: fixed {n}", yk2, [pk[i] for i in keep2],
                                  [groups[i] for i in keep2], rng)})

    with (args.outdir / "grouped_validation_metrics.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows_out[0].keys()))
        w.writeheader()
        w.writerows(rows_out)
    with (args.outdir / "grouped_selections.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(sel_rows[0].keys()))
        w.writeheader()
        w.writerows(sel_rows)
    print(f"\nWrote {args.outdir}")


if __name__ == "__main__":
    main()
