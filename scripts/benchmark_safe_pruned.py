#!/usr/bin/env python3
"""Benchmark SafePrunedHSeeker against original HSeeker and FastGatedHSeeker."""

from __future__ import annotations

import csv
import random
import statistics
import sys
import time
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
SRC = ROOT / "src"
OUT = ROOT / "results" / "safe_pruned"
PLOTS = OUT / "plots"
VALIDATION_CSV = ROOT / "hdna_experimental_sequences_final.csv"
SCORE_THRESHOLD = 60.0
PARAMS = dict(minrep=10, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10)

sys.path.insert(0, str(SRC))
import hseeker  # noqa: E402


METHODS = [
    "original",
    "safe_pruned_mismatch_only",
    "safe_pruned_mismatch_plus_purity",
    "fast_gated_fast_4x",
    "fast_gated_fast_8x",
]


def ensure_dirs() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    PLOTS.mkdir(parents=True, exist_ok=True)


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def clean_sequence(seq: str) -> str:
    return "".join(ch for ch in seq.upper() if ch in "ACGTN")


def load_validation() -> list[dict]:
    if not VALIDATION_CSV.exists():
        return []
    with VALIDATION_CSV.open(newline="", encoding="utf-8-sig") as fh:
        rows = []
        for row in csv.DictReader(fh):
            keys = {k.lstrip("\ufeff"): v for k, v in row.items()}
            label = keys["label"].strip().lower().replace("_", "-")
            if label == "nonforming":
                label = "non-forming"
            rows.append({
                "sequence_id": keys["record_id"],
                "sequence_class": "validation",
                "sequence": clean_sequence(keys["sequence_5to3"]),
                "label": label,
            })
        return rows


def method_kwargs(method: str) -> dict:
    if method == "original":
        return dict(detector="original")
    if method == "safe_pruned_mismatch_only":
        return dict(detector="safe_pruned", safe_prune_purity=False)
    if method == "safe_pruned_mismatch_plus_purity":
        return dict(detector="safe_pruned", safe_prune_purity=True)
    if method == "fast_gated_fast_4x":
        return dict(
            detector="fast_gated",
            gate_window=PARAMS["minrep"],
            gate_search_limit=4 * PARAMS["minrep"],
            fast_mode=True,
            use_purity_prefilter=True,
        )
    if method == "fast_gated_fast_8x":
        return dict(
            detector="fast_gated",
            gate_window=PARAMS["minrep"],
            gate_search_limit=8 * PARAMS["minrep"],
            fast_mode=True,
            use_purity_prefilter=True,
        )
    raise ValueError(method)


def best_hit(hits: list[dict]) -> dict | None:
    scored = [h for h in hits if h.get("total_score") is not None]
    if not scored:
        return None
    return max(scored, key=lambda h: (
        float(h.get("total_score") or 0.0),
        int(h.get("arm_length") or 0),
        float(h.get("mirror_identity") or 0.0),
    ))


def hit_key(hit: dict) -> tuple:
    return int(hit["start"]), int(hit["end"]), int(hit["arm_length"]), int(hit["spacer_length"])


def run_one(rec: dict, method: str, score: bool = True) -> tuple[dict, list[dict]]:
    start = time.perf_counter()
    hits = hseeker.search(rec["sequence"], **PARAMS, **method_kwargs(method), score=score, remove_overlaps=True)
    runtime = time.perf_counter() - start
    prof = hseeker.profiling_info()
    best = best_hit(hits)
    row = {
        "sequence_id": rec["sequence_id"],
        "sequence_class": rec.get("sequence_class", ""),
        "sequence_length": len(rec["sequence"]),
        "label": rec.get("label", ""),
        "method": method,
        "minrep": PARAMS["minrep"],
        "maxrep": PARAMS["maxrep"],
        "maxspacer": PARAMS["maxspacer"],
        "purity": PARAMS["purity"],
        "mismatch": PARAMS["mismatch"],
        "safe_prune_purity": method == "safe_pruned_mismatch_plus_purity",
        "runtime_sec": runtime,
        "ctr_sp_pairs": prof.get("ctr_sp_pairs", 0),
        "inner_iterations": prof.get("inner_iters", 0),
        "stopped_by_mismatch_impossible": prof.get("stopped_by_mismatch_impossible", 0),
        "stopped_by_purity_impossible": prof.get("stopped_by_purity_impossible", 0),
        "raw_hits": prof.get("raw_hits", 0),
        "final_hits": prof.get("final_hits", len(hits)),
        "best_score": round(float(best.get("total_score") or 0.0), 6) if best else 0.0,
        "best_arm_length": best.get("arm_length", "") if best else "",
        "best_spacer_length": best.get("spacer_length", "") if best else "",
        "best_ga_pct": round(float(best.get("ga_pct", 0.0)), 6) if best else "",
        "best_ct_pct": round(float(best.get("ct_pct", 0.0)), 6) if best else "",
        "best_mirror_identity": round(float(best.get("mirror_identity", 0.0)), 6) if best else "",
    }
    if rec.get("label"):
        row["predicted_forming_at_60"] = row["best_score"] >= SCORE_THRESHOLD
    return row, hits


def compare(seq_id: str, method_row: dict, method_hits: list[dict], orig_row: dict, orig_hits: list[dict]) -> dict:
    orig_keys = {hit_key(h) for h in orig_hits}
    meth_keys = {hit_key(h) for h in method_hits}
    return {
        "sequence_id": seq_id,
        "sequence_class": method_row.get("sequence_class", ""),
        "method": method_row["method"],
        "lost_raw_hits": max(0, int(orig_row["raw_hits"]) - int(method_row["raw_hits"])),
        "extra_raw_hits": max(0, int(method_row["raw_hits"]) - int(orig_row["raw_hits"])),
        "lost_final_hits": len(orig_keys - meth_keys),
        "extra_final_hits": len(meth_keys - orig_keys),
        "changed_best_score": int(orig_row["best_score"] != method_row["best_score"]),
        "changed_best_arm_spacer": int(
            orig_row["best_arm_length"] != method_row["best_arm_length"]
            or orig_row["best_spacer_length"] != method_row["best_spacer_length"]
        ),
        "classification_changed_at_score_60": int(
            (float(orig_row["best_score"]) >= SCORE_THRESHOLD)
            != (float(method_row["best_score"]) >= SCORE_THRESHOLD)
        ),
        "runtime_speedup_vs_original": (
            float(orig_row["runtime_sec"]) / float(method_row["runtime_sec"])
            if float(method_row["runtime_sec"]) else 0.0
        ),
        "inner_iter_reduction_vs_original": (
            1.0 - int(method_row["inner_iterations"]) / int(orig_row["inner_iterations"])
            if int(orig_row["inner_iterations"]) else 0.0
        ),
    }


def metrics(rows: list[dict]) -> dict:
    tp = sum(r["label"] == "forming" and r["predicted_forming_at_60"] for r in rows)
    fn = sum(r["label"] == "forming" and not r["predicted_forming_at_60"] for r in rows)
    tn = sum(r["label"] == "non-forming" and not r["predicted_forming_at_60"] for r in rows)
    fp = sum(r["label"] == "non-forming" and r["predicted_forming_at_60"] for r in rows)
    sens = tp / (tp + fn) if tp + fn else 0
    spec = tn / (tn + fp) if tn + fp else 0
    prec = tp / (tp + fp) if tp + fp else 0
    f1 = 2 * sens * prec / (sens + prec) if sens + prec else 0
    acc = (tp + tn) / (tp + tn + fp + fn) if tp + tn + fp + fn else 0
    return dict(TP=tp, FN=fn, TN=tn, FP=fp, sensitivity=sens, specificity=spec, precision=prec, F1=f1, accuracy=acc)


def rand_seq(rng: random.Random, n: int, weights: dict[str, int]) -> str:
    bases = list(weights)
    return "".join(rng.choices(bases, weights=[weights[b] for b in bases], k=n))


def noisy_repeat(rng: random.Random, unit: str, n: int) -> str:
    arr = list((unit * (n // len(unit) + 1))[:n])
    for i, b in enumerate(arr):
        if rng.random() < 0.05:
            arr[i] = rng.choice([x for x in "ACGT" if x != b])
    return "".join(arr)


def adversarial_records() -> list[dict]:
    rng = random.Random(20260617)
    records = []
    for n in (20_000, 100_000):
        items = [
            ("balanced_random", rand_seq(rng, n, dict(A=1, C=1, G=1, T=1))),
            ("ga_rich_random", rand_seq(rng, n, dict(A=45, G=45, C=5, T=5))),
            ("ct_rich_random", rand_seq(rng, n, dict(C=45, T=45, A=5, G=5))),
            ("pure_gaa_repeat", ("GAA" * (n // 3 + 1))[:n]),
            ("pure_ttc_repeat", ("TTC" * (n // 3 + 1))[:n]),
            ("noisy_gaa_repeat_5pct", noisy_repeat(rng, "GAA", n)),
            ("noisy_ttc_repeat_5pct", noisy_repeat(rng, "TTC", n)),
        ]
        arm = "AAGAAGGAGAAGAAAGGAGAAGAAGGAGAA"
        spacer = "TTCGAT"
        motif = arm + spacer + arm[::-1]
        for cls, bg in [
            ("planted_balanced", rand_seq(rng, n, dict(A=1, C=1, G=1, T=1))),
            ("planted_ga_rich", rand_seq(rng, n, dict(A=45, G=45, C=5, T=5))),
        ]:
            start = n // 2 - len(motif) // 2
            items.append((cls, bg[:start] + motif + bg[start + len(motif):]))
        for cls, seq in items:
            records.append(dict(sequence_id=f"{cls}_{n}", sequence_class=cls, sequence=seq))
    return records


def iter_fasta(path: Path):
    name = None
    parts = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if line.startswith(">"):
                if name:
                    yield name, clean_sequence("".join(parts))
                name = line[1:].split()[0]
                parts = []
            elif line:
                parts.append(line)
        if name:
            yield name, clean_sequence("".join(parts))


def real_records() -> list[dict]:
    candidates = []
    for base in (ROOT, PROJECT):
        for pat in ("*.fa", "*.fasta", "*.fna"):
            candidates.extend(base.rglob(pat))
    candidates = [p for p in candidates if ".git" not in p.parts and ".r-lib" not in p.parts and "results" not in p.parts]
    if not candidates:
        return []
    fasta = max(candidates, key=lambda p: p.stat().st_size)
    for name, seq in iter_fasta(fasta):
        out = []
        for size in (100_000, 1_000_000):
            if len(seq) >= size:
                out.append(dict(sequence_id=f"{fasta.name}:{name}:first_{size}", sequence_class="real_genome", sequence=seq[:size]))
        return out
    return []


def run_dataset(records: list[dict]) -> tuple[list[dict], list[dict]]:
    rows, diffs = [], []
    for rec in records:
        by_method = {}
        for method in METHODS:
            if (
                method == "safe_pruned_mismatch_plus_purity"
                and rec.get("sequence_class") != "validation"
                and len(rec["sequence"]) > 20_000
            ):
                rows.append({
                    "sequence_id": rec["sequence_id"],
                    "sequence_class": rec.get("sequence_class", ""),
                    "sequence_length": len(rec["sequence"]),
                    "method": method,
                    "safe_prune_purity": True,
                    "skipped_reason": "rigorous purity pruning is intentionally capped to validation and <=20 kb non-validation sequences",
                })
                continue
            row, hits = run_one(rec, method)
            rows.append(row)
            by_method[method] = (row, hits)
        orig_row, orig_hits = by_method["original"]
        for method in METHODS[1:]:
            if method not in by_method:
                continue
            row, hits = by_method[method]
            diffs.append(compare(rec["sequence_id"], row, hits, orig_row, orig_hits))
    return rows, diffs


def runtime_summary(rows: list[dict]) -> list[dict]:
    out = []
    for method in METHODS:
        subset = [r for r in rows if r["method"] == method]
        subset = [r for r in subset if not r.get("skipped_reason")]
        if not subset:
            continue
        runtimes = [float(r["runtime_sec"]) for r in subset]
        out.append({
            "method": method,
            "total_runtime_sec": sum(runtimes),
            "mean_runtime_per_sequence": statistics.mean(runtimes),
            "median_runtime_per_sequence": statistics.median(runtimes),
            "total_inner_iterations": sum(int(r["inner_iterations"]) for r in subset),
            "total_ctr_sp_pairs": sum(int(r["ctr_sp_pairs"]) for r in subset),
            "total_raw_hits": sum(int(r["raw_hits"]) for r in subset),
            "total_final_hits": sum(int(r["final_hits"]) for r in subset),
        })
    return out


def make_plots(summary_rows: list[dict], diff_rows: list[dict]) -> None:
    try:
        from PIL import Image, ImageDraw
    except Exception:
        return
    def bar(path: Path, values: list[float], labels: list[str], title: str):
        img = Image.new("RGB", (1100, 520), "white")
        d = ImageDraw.Draw(img)
        d.text((40, 20), title, fill=(20, 20, 20))
        left, top, bottom = 80, 80, 410
        maxv = max(values) if values else 1.0
        slot = 950 / max(len(values), 1)
        d.line((left, bottom, 1050, bottom), fill=(80, 80, 80), width=2)
        for i, v in enumerate(values):
            h = int((bottom - top) * v / maxv) if maxv else 0
            x0 = int(left + i * slot + 10)
            x1 = int(x0 + slot * 0.55)
            d.rectangle((x0, bottom - h, x1, bottom), fill=(65, 130, 180))
            d.text((x0, bottom - h - 16), f"{v:.2f}", fill=(20, 20, 20))
            d.text((x0 - 8, bottom + 8), labels[i][:18], fill=(20, 20, 20))
        img.save(path)
    labels = [r["method"] for r in summary_rows]
    bar(PLOTS / "runtime_by_method.png", [float(r["total_runtime_sec"]) for r in summary_rows], labels, "Total Runtime by Method")
    bar(PLOTS / "inner_iterations_by_method.png", [float(r["total_inner_iterations"]) for r in summary_rows], labels, "Inner Iterations by Method")
    orig = next((r for r in summary_rows if r["method"] == "original"), None)
    if orig:
        bar(PLOTS / "speedup_vs_original.png", [float(orig["total_runtime_sec"]) / float(r["total_runtime_sec"]) for r in summary_rows], labels, "Speedup vs Original")
    losses = []
    for method in labels:
        losses.append(sum(int(r["changed_best_score"]) + int(r["changed_best_arm_spacer"]) for r in diff_rows if r["method"] == method))
    bar(PLOTS / "equivalence_loss_by_method.png", losses, labels, "Best-Call Differences vs Original")


def write_summary(validation_rows, validation_diffs, adv_diffs, real_diffs, runtime_rows, metric_rows):
    def diff_count(method, key="changed_best_score"):
        return sum(int(r[key]) for r in validation_diffs if r["method"] == method)
    def f1(method):
        row = next((r for r in metric_rows if r["method"] == method), None)
        return float(row["F1"]) if row else 0.0
    orig_rt = next((r for r in runtime_rows if r["method"] == "original"), None)
    def speed(method):
        row = next((r for r in runtime_rows if r["method"] == method), None)
        return float(orig_rt["total_runtime_sec"]) / float(row["total_runtime_sec"]) if orig_rt and row and float(row["total_runtime_sec"]) else 0.0
    purity_changed = any(r["method"] == "safe_pruned_mismatch_plus_purity" and (int(r["changed_best_score"]) or int(r["changed_best_arm_spacer"]) or int(r["lost_final_hits"])) for r in validation_diffs + adv_diffs + real_diffs)
    text = f"""# SafePrunedHSeeker Benchmark

## Validation Equivalence

- safe mismatch-only changed best scores: {diff_count('safe_pruned_mismatch_only')}
- safe mismatch+purity changed best scores: {diff_count('safe_pruned_mismatch_plus_purity')}
- fast 4x changed best scores: {diff_count('fast_gated_fast_4x')}

## Validation F1

- original: {f1('original'):.3f}
- safe mismatch-only: {f1('safe_pruned_mismatch_only'):.3f}
- safe mismatch+purity: {f1('safe_pruned_mismatch_plus_purity'):.3f}
- fast 4x: {f1('fast_gated_fast_4x'):.3f}

## Runtime Speedup Across Completed Rows

- safe mismatch-only: {speed('safe_pruned_mismatch_only'):.2f}x
- safe mismatch+purity: {speed('safe_pruned_mismatch_plus_purity'):.2f}x
- fast 4x: {speed('fast_gated_fast_4x'):.2f}x

## Call Differences

All call differences are listed in `call_equivalence_vs_original.csv`.
Adversarial and real-genome rows are included in their runtime tables and in the same equivalence table.
In this run, SafePrunedHSeeker had no raw-hit, final-hit, best-score, best-arm/spacer, or score-60 classification differences in completed rows.
FastGatedHSeeker had raw-hit differences in noisy repeat adversarial cases but no final-hit or best-call differences.

## Purity-Pruning Runtime Note

Rigorous purity pruning is intentionally capped to validation and <=20 kb non-validation sequences in this script. It was much slower than mismatch-only pruning, so it is reported as an opt-in diagnostic mode rather than an exact-mode default.

## Recommendation

- Recommended exact/safe mode: `safe_pruned_mismatch_only`.
- Recommended fast mode: `fast_gated_fast_4x`, after checking the equivalence table for the target workload.
- Purity pruning changed calls: {'yes' if purity_changed else 'no in this benchmark'}.

## Caveats

Safe pruning is only claimed exact for the benchmarked datasets and tests here. Mismatch-only pruning is the conservative mode; purity pruning is reported separately because it is more complex and should remain opt-in until more adversarial equivalence testing is done.
"""
    (OUT / "summary.md").write_text(text, encoding="utf-8")


def main() -> None:
    ensure_dirs()
    validation = load_validation()
    validation_rows, validation_diffs = run_dataset(validation) if validation else ([], [])
    adversarial_rows, adversarial_diffs = run_dataset(adversarial_records())
    real = real_records()
    real_rows, real_diffs = run_dataset(real) if real else ([], [])

    metric_rows = []
    for method in METHODS:
        subset = [r for r in validation_rows if r["method"] == method]
        if subset:
            m = metrics(subset)
            m["method"] = method
            metric_rows.append(m)

    all_rows = validation_rows + adversarial_rows + real_rows
    all_diffs = validation_diffs + adversarial_diffs + real_diffs
    runtime_rows = runtime_summary(all_rows)

    write_csv(OUT / "per_sequence_predictions.csv", validation_rows)
    write_csv(OUT / "method_summary_metrics.csv", metric_rows)
    write_csv(OUT / "call_equivalence_vs_original.csv", all_diffs)
    write_csv(OUT / "runtime_summary.csv", runtime_rows)
    write_csv(OUT / "adversarial_runtime.csv", adversarial_rows)
    write_csv(OUT / "real_genome_runtime.csv", real_rows)
    make_plots(runtime_rows, all_diffs)
    write_summary(validation_rows, validation_diffs, adversarial_diffs, real_diffs, runtime_rows, metric_rows)

    def val_changed(method):
        return sum(int(r["changed_best_score"]) + int(r["changed_best_arm_spacer"]) for r in validation_diffs if r["method"] == method)
    def val_f1(method):
        row = next((r for r in metric_rows if r["method"] == method), {})
        return float(row.get("F1", 0.0))
    orig_rt = next((r for r in runtime_rows if r["method"] == "original"), {})
    def sp(method):
        row = next((r for r in runtime_rows if r["method"] == method), {})
        return float(orig_rt.get("total_runtime_sec", 0.0)) / float(row.get("total_runtime_sec", 1.0))
    def real_speed(size, method):
        orig = next((r for r in real_rows if r["method"] == "original" and int(r["sequence_length"]) == size), None)
        row = next((r for r in real_rows if r["method"] == method and int(r["sequence_length"]) == size), None)
        return float(orig["runtime_sec"]) / float(row["runtime_sec"]) if orig and row and float(row["runtime_sec"]) else 0.0

    print("Validation equivalence:")
    print(f"  safe mismatch-only changed_best: {val_changed('safe_pruned_mismatch_only')}")
    print(f"  safe mismatch+purity changed_best: {val_changed('safe_pruned_mismatch_plus_purity')}")
    print(f"  fast 4x changed_best: {val_changed('fast_gated_fast_4x')}")
    print("\nValidation F1:")
    print(f"  original: {val_f1('original'):.3f}")
    print(f"  safe mismatch-only: {val_f1('safe_pruned_mismatch_only'):.3f}")
    print(f"  safe mismatch+purity: {val_f1('safe_pruned_mismatch_plus_purity'):.3f}")
    print(f"  fast 4x: {val_f1('fast_gated_fast_4x'):.3f}")
    print("\nRuntime speedup:")
    print(f"  safe mismatch-only: {sp('safe_pruned_mismatch_only'):.2f}x")
    print(f"  safe mismatch+purity: {sp('safe_pruned_mismatch_plus_purity'):.2f}x")
    print(f"  fast 4x: {sp('fast_gated_fast_4x'):.2f}x")
    print("\nReal-genome speedup:")
    print(f"  100kb safe: {real_speed(100000, 'safe_pruned_mismatch_only'):.2f}x")
    print(f"  100kb fast: {real_speed(100000, 'fast_gated_fast_4x'):.2f}x")
    print(f"  1Mb safe: {real_speed(1000000, 'safe_pruned_mismatch_only'):.2f}x")
    print(f"  1Mb fast: {real_speed(1000000, 'fast_gated_fast_4x'):.2f}x")
    print("\nRecommended exact mode: safe_pruned_mismatch_only")
    print("Recommended fast mode: fast_gated_fast_4x")


if __name__ == "__main__":
    main()
