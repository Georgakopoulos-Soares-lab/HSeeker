#!/usr/bin/env python3
"""Compare FastGatedHSeeker validation calls against legacy Triplex.

Requires the local R library cache created for this project, including the
Bioconductor ``triplex`` package. The helper R script writes raw Triplex calls;
this wrapper joins them to ``results/fast_gated/per_sequence_predictions.csv``.
"""

from __future__ import annotations

import csv
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "fast_gated"


def read_csv(path: Path) -> list[dict]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def write_csv(path: Path, rows: list[dict]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def boolish(value) -> bool:
    if isinstance(value, bool):
        return value
    return str(value).lower() in {"true", "1", "yes"}


def main() -> None:
    predictions = read_csv(OUT / "per_sequence_predictions.csv")
    by_id_method = {(r["sequence_id"], r["method"]): r for r in predictions}

    helper = ROOT / "scripts" / "run_triplex_benchmark.R"
    triplex_raw = OUT / "triplex_predictions_raw.csv"
    if (ROOT / ".r-lib" / "triplex").exists() and helper.exists():
        subprocess.run(["Rscript", str(helper)], cwd=ROOT, check=True)
    triplex = {}
    if triplex_raw.exists():
        triplex = {r["sequence_id"]: r for r in read_csv(triplex_raw)}

    sequence_ids = sorted({r["sequence_id"] for r in predictions})
    rows = []
    for sid in sequence_ids:
        orig = by_id_method[(sid, "original")]
        fast = by_id_method[(sid, "fast_gated_fast_4x")]
        t = triplex.get(sid, {})
        t_pred = t.get("triplex_pred", "not_run")
        rows.append({
            "sequence_id": sid,
            "label": orig["label"],
            "hseeker_original_score": orig["best_score"],
            "fast_gated_score": fast["best_score"],
            "triplex_call": t.get("triplex_n_hits", ""),
            "hseeker_original_pred": boolish(orig["predicted_forming_at_60"]),
            "fast_gated_pred": boolish(fast["predicted_forming_at_60"]),
            "triplex_pred": t_pred,
            "hseeker_original_correct": boolish(orig["correct"]),
            "fast_gated_correct": boolish(fast["correct"]),
            "triplex_correct": (t_pred == orig["label"]) if t_pred in {"forming", "non-forming"} else "",
        })
    write_csv(OUT / "triplex_comparison.csv", rows)
    print(f"Wrote {OUT / 'triplex_comparison.csv'}")


if __name__ == "__main__":
    main()
