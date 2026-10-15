#!/usr/bin/env python3
"""Supplementary Figure 1: confusion matrices for all five tools.

Each tool at its documented default settings on the balanced 130-record
benchmark (65 forming / 65 non-forming), with the Matthews correlation
coefficient on each panel.

Input: analysis/results/tool_comparison/per_record_calls.csv
"""
from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parents[2]
CALLS = ROOT / "analysis" / "results" / "tool_comparison" / "per_record_calls.csv"
OUT = ROOT / "analysis" / "results" / "plots"

# The documented-default setting for each tool, as pre-specified in the
# comparison design (NeSSie's -T triplex motif; nBMST's triplex-classified
# mirror repeats).
PRIMARY = {
    "HSeeker": "defaults",
    "triplex": "defaults",
    "NeSSie": "T_defaults",
    "nBMST": "triplex_defaults",
    "Triplexator": "defaults",
}
ORDER = ["HSeeker", "triplex", "NeSSie", "nBMST", "Triplexator"]


def mcc(tp: int, fn: int, tn: int, fp: int) -> float:
    d = math.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return ((tp * tn) - (fp * fn)) / d if d else 0.0


def main() -> None:
    calls: dict[str, dict[str, int]] = {}
    labels: dict[str, str] = {}
    for r in csv.DictReader(open(CALLS)):
        if PRIMARY.get(r["tool"]) != r["setting"]:
            continue
        calls.setdefault(r["tool"], {})[r["record_id"]] = int(r["call"])
        labels[r["record_id"]] = r["label"]

    tools = [t for t in ORDER if t in calls]
    fig, axes = plt.subplots(1, len(tools), figsize=(3.05 * len(tools), 3.5))
    fig.subplots_adjust(left=0.065, right=0.99, top=0.80, bottom=0.17, wspace=0.30)

    for ax, tool in zip(np.atleast_1d(axes), tools):
        tp = fn = tn = fp = 0
        for rid, lab in labels.items():
            c = calls[tool].get(rid, 0)
            if lab == "forming":
                tp += c == 1
                fn += c == 0
            else:
                fp += c == 1
                tn += c == 0
        m = np.array([[tp, fn], [fp, tn]], dtype=float)
        frac = m / m.sum(axis=1, keepdims=True)

        ax.imshow(frac, cmap="Blues", vmin=0, vmax=1)
        for i in range(2):
            for j in range(2):
                ax.text(j, i, f"{int(m[i, j])}\n{frac[i, j]*100:.1f}%",
                        ha="center", va="center", fontsize=11,
                        color="white" if frac[i, j] > 0.55 else "#222222")
        ax.set_xticks([0, 1]); ax.set_yticks([0, 1])
        ax.set_xticklabels(["Forming", "Non-forming"], fontsize=9)
        ax.set_yticklabels(["Forming", "Non-forming"], fontsize=9, rotation=90, va="center")
        ax.set_xlabel("Predicted", fontsize=9.5)
        if tool == tools[0]:
            ax.set_ylabel("Reported", fontsize=9.5)
        ax.set_title(f"{tool}\nMCC = {mcc(tp, fn, tn, fp):.3f}",
                     fontsize=11, fontweight="bold", pad=9)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
        print(f"  {tool:12s} TP {tp:2d} FN {fn:2d} TN {tn:2d} FP {fp:2d}  "
              f"MCC {mcc(tp, fn, tn, fp):.3f}")

    OUT.mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png", "svg"):
        fig.savefig(OUT / f"suppfig1_confusion_matrices.{ext}", dpi=300, bbox_inches="tight")
    print("wrote", OUT / "suppfig1_confusion_matrices.pdf")


if __name__ == "__main__":
    main()
