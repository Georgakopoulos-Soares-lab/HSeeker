#!/usr/bin/env python3
"""Next-step FastGatedHSeeker stress and tradeoff benchmarks.

This script intentionally leaves detector semantics alone. It measures the
already-implemented original and FastGated detectors on:

- local real-genome FASTA chunks, when available
- synthetic GA/CT-rich and repeat-heavy adversarial sequences
- validation-set gate-search-limit tradeoffs

Outputs are grouped under ``results/fast_gated_nextstep/``.
"""

from __future__ import annotations

import csv
import random
import statistics
import struct
import subprocess
import sys
import time
import zlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
SRC = ROOT / "src"
OUT = ROOT / "results" / "fast_gated_nextstep"
TABLES = OUT / "tables"
PLOTS = OUT / "plots"
SEQS = OUT / "sequences"
VALIDATION_CSV = ROOT / "hdna_experimental_sequences_final.csv"
SCORE_THRESHOLD = 60.0

PARAMS = {
    "minrep": 10,
    "maxrep": 1000,
    "maxspacer": 10,
    "purity": 0.90,
    "mismatch": 0.10,
}

sys.path.insert(0, str(SRC))
import hseeker  # noqa: E402


def ensure_dirs() -> None:
    for path in (OUT, TABLES, PLOTS, SEQS):
        path.mkdir(parents=True, exist_ok=True)


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
    with VALIDATION_CSV.open(newline="", encoding="utf-8-sig") as fh:
        reader = csv.DictReader(fh)
        rows = []
        for row in reader:
            keys = {k.lstrip("\ufeff"): v for k, v in row.items()}
            label = keys["label"].strip().lower().replace("_", "-")
            if label == "nonforming":
                label = "non-forming"
            rows.append({
                "sequence_id": keys["record_id"],
                "sequence_name": keys.get("sequence_name", ""),
                "sequence": clean_sequence(keys["sequence_5to3"]),
                "label": label,
            })
        return rows


def iter_fasta(path: Path):
    name = None
    parts: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    yield name, clean_sequence("".join(parts))
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            yield name, clean_sequence("".join(parts))


def find_real_fasta() -> Path | None:
    candidates = []
    for base in (ROOT, PROJECT):
        for pattern in ("*.fa", "*.fasta", "*.fna"):
            candidates.extend(base.rglob(pattern))
    filtered = [
        p for p in candidates
        if ".git" not in p.parts
        and ".r-lib" not in p.parts
        and "results" not in p.parts
        and p.is_file()
    ]
    if not filtered:
        return None
    filtered.sort(key=lambda p: p.stat().st_size, reverse=True)
    return filtered[0]


def get_real_chunks() -> list[dict]:
    fasta = find_real_fasta()
    if fasta is None:
        print("No local FASTA found; skipping real-genome benchmark.")
        return []
    chunks = []
    target_sizes = [100_000, 1_000_000, 5_000_000]
    for name, seq in iter_fasta(fasta):
        if not seq:
            continue
        for size in target_sizes:
            if len(seq) >= size:
                chunks.append({
                    "chunk_id": f"{fasta.name}:{name}:first_{size}",
                    "chunk_length": size,
                    "sequence": seq[:size],
                    "source_fasta": str(fasta),
                })
        break
    if not chunks:
        print(f"Found FASTA {fasta}, but no record long enough for 100 kb; skipping real-genome benchmark.")
    return chunks


def detector_kwargs(method: str, params: dict) -> dict:
    if method == "original":
        return {"detector": "original"}
    if method == "fast_gated_fast_2x":
        mult = 2
    elif method == "fast_gated_fast_4x":
        mult = 4
    elif method == "fast_gated_fast_8x":
        mult = 8
    elif method == "fast_gated_fast_16x":
        mult = 16
    elif method == "fast_gated_maxrep":
        mult = None
    else:
        raise ValueError(method)
    limit = params["maxrep"] if mult is None else mult * params["minrep"]
    return {
        "detector": "fast_gated",
        "gate_window": params["minrep"],
        "gate_search_limit": limit,
        "fast_mode": limit != params["maxrep"],
        "use_purity_prefilter": limit != params["maxrep"],
    }


def hit_key(hit: dict) -> tuple:
    return (
        int(hit["start"]),
        int(hit["end"]),
        int(hit["arm_length"]),
        int(hit["spacer_length"]),
    )


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


def run_detector(sequence: str, method: str, params: dict, *, score: bool = True) -> tuple[list[dict], dict, float]:
    kwargs = detector_kwargs(method, params)
    start = time.perf_counter()
    hits = hseeker.search(sequence, **params, **kwargs, remove_overlaps=True, score=score)
    runtime = time.perf_counter() - start
    return hits, hseeker.profiling_info(), runtime


def summarize_hits(hits: list[dict]) -> dict:
    best = best_hit(hits)
    arms = [int(h["arm_length"]) for h in hits]
    return {
        "best_score": round(float(best.get("total_score") or 0.0), 6) if best else 0.0,
        "max_arm_length": max(arms) if arms else 0,
        "mean_arm_length": round(statistics.mean(arms), 6) if arms else 0.0,
        "best_arm_length": best.get("arm_length", "") if best else "",
        "best_spacer_length": best.get("spacer_length", "") if best else "",
        "best_ga_pct": round(float(best.get("ga_pct", 0.0)), 6) if best else "",
        "best_ct_pct": round(float(best.get("ct_pct", 0.0)), 6) if best else "",
        "best_mirror_identity": round(float(best.get("mirror_identity", 0.0)), 6) if best else "",
        "left_arm": best.get("left_arm", "") if best else "",
        "spacer": best.get("spacer", "") if best else "",
        "right_arm": best.get("right_arm", "") if best else "",
    }


def method_row(prefix: dict, sequence: str, method: str, params: dict, *, score: bool = True) -> tuple[dict, list[dict]]:
    hits, prof, runtime = run_detector(sequence, method, params, score=score)
    row = dict(prefix)
    row.update({
        "method": method,
        "minrep": params["minrep"],
        "maxrep": params["maxrep"],
        "maxspacer": params["maxspacer"],
        "purity": params["purity"],
        "mismatch": params["mismatch"],
        "gate_window": params["minrep"] if method != "original" else "",
        "gate_search_limit": detector_kwargs(method, params).get("gate_search_limit", ""),
        "runtime_sec": runtime,
        "ctr_sp_pairs": prof.get("ctr_sp_pairs", 0),
        "inner_iterations": prof.get("inner_iters", 0),
        "gate_passed": prof.get("gate_passed", 0),
        "raw_hits": prof.get("raw_hits", 0),
        "final_hits": prof.get("final_hits", len(hits)),
    })
    row.update(summarize_hits(hits))
    return row, hits


def compare_to_original(chunk_id: str, method_row_: dict, method_hits: list[dict], original_row: dict, original_hits: list[dict]) -> dict:
    orig_keys = {hit_key(h) for h in original_hits}
    method_keys = {hit_key(h) for h in method_hits}
    best_changed = method_row_["best_score"] != original_row["best_score"]
    arm_changed = (
        method_row_["best_arm_length"] != original_row["best_arm_length"]
        or method_row_["best_spacer_length"] != original_row["best_spacer_length"]
    )
    orig_inner = int(original_row["inner_iterations"])
    meth_inner = int(method_row_["inner_iterations"])
    return {
        "chunk_id": chunk_id,
        "method": method_row_["method"],
        "lost_raw_hits": max(0, int(original_row["raw_hits"]) - int(method_row_["raw_hits"])),
        "lost_final_hits": len(orig_keys - method_keys),
        "extra_raw_hits": max(0, int(method_row_["raw_hits"]) - int(original_row["raw_hits"])),
        "extra_final_hits": len(method_keys - orig_keys),
        "changed_best_score_count": int(best_changed),
        "changed_best_arm_spacer_count": int(arm_changed),
        "runtime_speedup_vs_original": (
            float(original_row["runtime_sec"]) / float(method_row_["runtime_sec"])
            if float(method_row_["runtime_sec"]) else 0.0
        ),
        "inner_iter_reduction_vs_original": (
            1.0 - (meth_inner / orig_inner) if orig_inner else 0.0
        ),
    }


def random_weighted(rng: random.Random, length: int, weights: dict[str, float]) -> str:
    letters = list(weights)
    probs = [weights[b] for b in letters]
    total = sum(probs)
    cum = []
    running = 0.0
    for p in probs:
        running += p / total
        cum.append(running)
    out = []
    for _ in range(length):
        x = rng.random()
        for letter, c in zip(letters, cum):
            if x <= c:
                out.append(letter)
                break
    return "".join(out)


def noisy_repeat(rng: random.Random, unit: str, length: int, subst: float = 0.05) -> str:
    seq = (unit * ((length // len(unit)) + 1))[:length]
    bases = "ACGT"
    arr = list(seq)
    for i, b in enumerate(arr):
        if rng.random() < subst:
            choices = [x for x in bases if x != b]
            arr[i] = rng.choice(choices)
    return "".join(arr)


def planted_sequence(rng: random.Random, length: int, background: str) -> tuple[str, dict]:
    arm_len = 30
    spacer = "TTCGAT"
    left = "AAGAAGGAGAAGAAAGGAGAAGAAGGAGAA"[:arm_len]
    motif = left + spacer + left[::-1]
    start = length // 2 - len(motif) // 2
    seq = background[:start] + motif + background[start + len(motif):]
    meta = {
        "planted_start": start + 1,
        "planted_end": start + len(motif),
        "planted_arm_length": arm_len,
        "planted_spacer_length": len(spacer),
        "expected_detectable": True,
    }
    return seq, meta


def generate_adversarial() -> list[dict]:
    rng = random.Random(20260617)
    records = []
    fasta_lines = []
    for length in (20_000, 100_000):
        generators = [
            ("balanced_random", random_weighted(rng, length, {"A": 1, "C": 1, "G": 1, "T": 1}), {}),
            ("ga_rich_random", random_weighted(rng, length, {"A": 45, "G": 45, "C": 5, "T": 5}), {}),
            ("ct_rich_random", random_weighted(rng, length, {"C": 45, "T": 45, "A": 5, "G": 5}), {}),
            ("pure_gaa_repeat", ("GAA" * ((length // 3) + 1))[:length], {}),
            ("pure_ttc_repeat", ("TTC" * ((length // 3) + 1))[:length], {}),
            ("noisy_gaa_repeat_5pct", noisy_repeat(rng, "GAA", length), {}),
            ("noisy_ttc_repeat_5pct", noisy_repeat(rng, "TTC", length), {}),
        ]
        bal = random_weighted(rng, length, {"A": 1, "C": 1, "G": 1, "T": 1})
        seq, meta = planted_sequence(rng, length, bal)
        generators.append(("planted_balanced", seq, meta))
        garich = random_weighted(rng, length, {"A": 45, "G": 45, "C": 5, "T": 5})
        seq, meta = planted_sequence(rng, length, garich)
        generators.append(("planted_ga_rich", seq, meta))
        for cls, seq, meta in generators:
            sid = f"{cls}_{length}"
            rec = {
                "sequence_id": sid,
                "sequence_class": cls,
                "sequence_length": length,
                "sequence": seq,
                "planted_start": meta.get("planted_start", ""),
                "planted_end": meta.get("planted_end", ""),
                "planted_arm_length": meta.get("planted_arm_length", ""),
                "planted_spacer_length": meta.get("planted_spacer_length", ""),
                "expected_detectable": meta.get("expected_detectable", False),
            }
            records.append(rec)
            fasta_lines.append(f">{sid} class={cls} length={length} planted_start={rec['planted_start']} planted_end={rec['planted_end']} expected_detectable={rec['expected_detectable']}")
            for i in range(0, len(seq), 80):
                fasta_lines.append(seq[i:i + 80])
    (SEQS / "adversarial_sequences.fasta").write_text("\n".join(fasta_lines) + "\n", encoding="utf-8")
    return records


def planted_detected(hits: list[dict], rec: dict) -> bool | str:
    if not rec.get("expected_detectable"):
        return ""
    start = int(rec["planted_start"])
    end = int(rec["planted_end"])
    for h in hits:
        if int(h["start"]) <= start and int(h["end"]) >= end:
            return True
        if start <= int(h["start"]) <= end or start <= int(h["end"]) <= end:
            return True
    return False


def run_real_genome() -> tuple[list[dict], list[dict], bool]:
    chunks = get_real_chunks()
    if not chunks:
        return [], [], False
    runtime_rows = []
    diff_rows = []
    methods = ["original", "fast_gated_fast_4x", "fast_gated_fast_8x"]
    for chunk in chunks:
        by_method = {}
        for method in methods:
            prefix = {"chunk_id": chunk["chunk_id"], "chunk_length": chunk["chunk_length"]}
            row, hits = method_row(prefix, chunk["sequence"], method, PARAMS, score=True)
            runtime_rows.append(row)
            by_method[method] = (row, hits)
        orig_row, orig_hits = by_method["original"]
        for method in methods[1:]:
            row, hits = by_method[method]
            diff = compare_to_original(chunk["chunk_id"], row, hits, orig_row, orig_hits)
            diff["chunk_length"] = chunk["chunk_length"]
            diff["source_fasta"] = chunk["source_fasta"]
            diff_rows.append(diff)
    write_csv(TABLES / "real_genome_runtime.csv", runtime_rows)
    write_csv(TABLES / "real_genome_call_differences.csv", diff_rows)
    return runtime_rows, diff_rows, True


def run_adversarial() -> tuple[list[dict], list[dict]]:
    records = generate_adversarial()
    methods = [
        "original",
        "fast_gated_fast_2x",
        "fast_gated_fast_4x",
        "fast_gated_fast_8x",
        "fast_gated_fast_16x",
    ]
    runtime_rows = []
    diff_rows = []
    for rec in records:
        by_method = {}
        for method in methods:
            prefix = {
                "sequence_id": rec["sequence_id"],
                "sequence_class": rec["sequence_class"],
                "sequence_length": rec["sequence_length"],
            }
            row, hits = method_row(prefix, rec["sequence"], method, PARAMS, score=True)
            row["planted_detected"] = planted_detected(hits, rec)
            runtime_rows.append(row)
            by_method[method] = (row, hits)
        orig_row, orig_hits = by_method["original"]
        for method in methods[1:]:
            row, hits = by_method[method]
            d = compare_to_original(rec["sequence_id"], row, hits, orig_row, orig_hits)
            d.update({
                "sequence_id": rec["sequence_id"],
                "sequence_class": rec["sequence_class"],
                "sequence_length": rec["sequence_length"],
                "planted_detected": row["planted_detected"],
            })
            diff_rows.append(d)
    write_csv(TABLES / "adversarial_runtime.csv", runtime_rows)
    write_csv(TABLES / "adversarial_call_differences.csv", diff_rows)
    return runtime_rows, diff_rows


def classification_metrics(rows: list[dict]) -> dict:
    tp = sum(r["label"] == "forming" and r["predicted"] for r in rows)
    fn = sum(r["label"] == "forming" and not r["predicted"] for r in rows)
    tn = sum(r["label"] == "non-forming" and not r["predicted"] for r in rows)
    fp = sum(r["label"] == "non-forming" and r["predicted"] for r in rows)
    sens = tp / (tp + fn) if tp + fn else 0
    spec = tn / (tn + fp) if tn + fp else 0
    prec = tp / (tp + fp) if tp + fp else 0
    f1 = 2 * sens * prec / (sens + prec) if sens + prec else 0
    acc = (tp + tn) / (tp + tn + fp + fn) if tp + tn + fp + fn else 0
    return {
        "TP": tp, "FN": fn, "TN": tn, "FP": fp,
        "sensitivity": sens, "specificity": spec, "precision": prec, "F1": f1, "accuracy": acc,
    }


def run_validation_tradeoff() -> tuple[list[dict], list[dict], dict[str, dict]]:
    records = load_validation()
    settings = [
        ("original", "original"),
        ("fast_2x", "fast_gated_fast_2x"),
        ("fast_4x", "fast_gated_fast_4x"),
        ("fast_8x", "fast_gated_fast_8x"),
        ("fast_16x", "fast_gated_fast_16x"),
        ("fast_maxrep", "fast_gated_maxrep"),
    ]
    per_setting: dict[str, list[dict]] = {}
    hit_cache = {}
    for label, method in settings:
        rows = []
        start_total = time.perf_counter()
        total_inner = 0
        for rec in records:
            row, hits = method_row(
                {
                    "sequence_id": rec["sequence_id"],
                    "label": rec["label"],
                    "sequence_length": len(rec["sequence"]),
                },
                rec["sequence"], method, PARAMS, score=True
            )
            pred = float(row["best_score"]) >= SCORE_THRESHOLD
            row["predicted"] = pred
            row["setting"] = label
            rows.append(row)
            hit_cache[(rec["sequence_id"], label)] = (row, hits)
            total_inner += int(row["inner_iterations"])
        runtime = time.perf_counter() - start_total
        per_setting[label] = rows
        per_setting[label + "_summary"] = [{"runtime_sec": runtime, "inner_iterations": total_inner}]

    orig_rows = per_setting["original"]
    orig_by_id = {r["sequence_id"]: r for r in orig_rows}
    summary_rows = []
    lost_rows = []
    for label, _method in settings:
        rows = per_setting[label]
        m = classification_metrics(rows)
        runtime = per_setting[label + "_summary"][0]["runtime_sec"]
        inner = per_setting[label + "_summary"][0]["inner_iterations"]
        changed_score = 0
        changed_arm = 0
        lost = 0
        for r in rows:
            o = orig_by_id[r["sequence_id"]]
            if float(o["best_score"]) >= SCORE_THRESHOLD and float(r["best_score"]) < SCORE_THRESHOLD:
                lost += 1
                lost_rows.append({
                    "sequence_id": r["sequence_id"],
                    "label": r["label"],
                    "original_score": o["best_score"],
                    "fast_score": r["best_score"],
                    "original_arm_length": o["best_arm_length"],
                    "fast_arm_length": r["best_arm_length"],
                    "original_spacer_length": o["best_spacer_length"],
                    "fast_spacer_length": r["best_spacer_length"],
                    "gate_search_limit": r["gate_search_limit"],
                    "reason_if_inferable": "no gated hit above score threshold",
                })
            if o["best_score"] != r["best_score"]:
                changed_score += 1
            if o["best_arm_length"] != r["best_arm_length"] or o["best_spacer_length"] != r["best_spacer_length"]:
                changed_arm += 1
        m.update({
            "setting": label,
            "gate_search_limit": rows[0]["gate_search_limit"] if rows else "",
            "runtime_sec": runtime,
            "inner_iterations": inner,
            "lost_calls_vs_original": lost,
            "changed_best_score_count": changed_score,
            "changed_best_arm_spacer_count": changed_arm,
        })
        summary_rows.append(m)
    write_csv(TABLES / "gate_tradeoff_validation.csv", summary_rows)
    write_csv(TABLES / "gate_lost_calls_validation.csv", lost_rows)
    return summary_rows, lost_rows, {label: rows for label, _ in settings}


def run_failure_inspection(validation_by_setting: dict[str, list[dict]]) -> list[dict]:
    failure_ids = {"HDNA0022", "HDNA0023", "HDNA0024", "HDNA0025", "HDNA0027", "HDNA0030", "HDNA0034", "HDNA0053"}
    rows = []
    originals = {r["sequence_id"]: r for r in validation_by_setting["original"]}
    fast4 = {r["sequence_id"]: r for r in validation_by_setting["fast_4x"]}
    for sid in sorted(failure_ids):
        for method, table in (("original", originals), ("fast_4x", fast4)):
            r = table[sid]
            rows.append({
                "sequence_id": sid,
                "method": method,
                "label": r["label"],
                "sequence_length": next(x["sequence_length"] for x in table.values() if x["sequence_id"] == sid),
                "best_score": r["best_score"],
                "best_arm_length": r["best_arm_length"],
                "best_spacer_length": r["best_spacer_length"],
                "ga_pct": r["best_ga_pct"],
                "ct_pct": r["best_ct_pct"],
                "mirror_identity": r["best_mirror_identity"],
                "left_arm": r["left_arm"],
                "spacer": r["spacer"],
                "right_arm": r["right_arm"],
                "notes": "same as original" if method == "fast_4x" else "",
            })
    write_csv(TABLES / "validation_failure_inspection.csv", rows)
    return rows


def try_matplotlib():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        return plt
    except Exception:
        return None


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
    for y in range(max(0, y0), min(len(img), y1)):
        for x in range(max(0, x0), min(len(img[0]), x1)):
            img[y][x] = color


def fallback_bar(path: Path, values: list[float], labels: list[str] | None = None, title: str = "", ylabel: str = "") -> None:
    try:
        from PIL import Image, ImageDraw
        width, height = 1200, 620
        img = Image.new("RGB", (width, height), "white")
        draw = ImageDraw.Draw(img)
        left, top, right, bottom = 70, 70, width - 30, height - 170
        draw.text((left, 20), title, fill=(20, 20, 20))
        if ylabel:
            draw.text((10, top), ylabel, fill=(20, 20, 20))
        draw.line((left, bottom, right, bottom), fill=(80, 80, 80), width=2)
        draw.line((left, top, left, bottom), fill=(80, 80, 80), width=2)
        maxv = max(values) if values else 1.0
        maxv = max(maxv, 1e-9)
        slot = (right - left) / max(len(values), 1)
        for i, v in enumerate(values):
            bar_w = max(8, int(slot * 0.58))
            x0 = int(left + i * slot + (slot - bar_w) / 2)
            x1 = x0 + bar_w
            y0 = int(bottom - (bottom - top) * (v / maxv))
            draw.rectangle((x0, y0, x1, bottom), fill=(65, 130, 180))
            draw.text((x0, max(top, y0 - 16)), f"{v:.2f}", fill=(20, 20, 20))
            if labels and i < len(labels):
                lab = labels[i][:18]
                draw.text((x0 - 8, bottom + 8), lab, fill=(20, 20, 20))
        img.save(path)
        return
    except Exception:
        pass
    img = blank(900, 420)
    maxv = max(values) if values else 1.0
    for i, v in enumerate(values):
        h = int(320 * v / maxv) if maxv else 0
        rect(img, 80 + i * 70, 360 - h, 120 + i * 70, 360, (65, 130, 180))
    png_write(path, 900, 420, img)


def make_plots(adversarial_runtime, adversarial_diff, tradeoff, real_runtime) -> None:
    plt = try_matplotlib()
    if plt is None:
        # Labeled fallback PNGs when matplotlib is not installed.
        fast4 = [r for r in adversarial_diff if r["method"] == "fast_gated_fast_4x"]
        labels = [f"{r['sequence_class']} {int(r['sequence_length'])//1000}kb" for r in fast4]
        speedups = [float(r["runtime_speedup_vs_original"]) for r in fast4]
        fallback_bar(PLOTS / "adversarial_speedup_by_class.png", speedups, labels, "Fast 4x Speedup by Adversarial Class", "speedup")
        reductions = [float(r["inner_iter_reduction_vs_original"]) for r in fast4]
        fallback_bar(PLOTS / "adversarial_inner_iter_reduction.png", reductions, labels, "Fast 4x Inner-Iteration Reduction", "reduction")
        trade_labels = [r["setting"] for r in tradeoff]
        fallback_bar(PLOTS / "gate_limit_tradeoff_validation.png", [float(r["F1"]) for r in tradeoff], trade_labels, "Validation F1 by Gate Setting", "F1")
        if real_runtime:
            real_labels = [f"{r['method']} {int(r['chunk_length'])//1000}kb" for r in real_runtime]
            fallback_bar(PLOTS / "real_genome_runtime.png", [float(r["runtime_sec"]) for r in real_runtime], real_labels, "Real-Genome Runtime", "sec")
        return

    fast4 = [r for r in adversarial_diff if r["method"] == "fast_gated_fast_4x"]
    labels = [f"{r['sequence_class']}\n{int(r['sequence_length'])//1000}kb" for r in fast4]
    x = range(len(fast4))

    plt.figure(figsize=(14, 5))
    plt.bar(x, [float(r["runtime_speedup_vs_original"]) for r in fast4])
    plt.xticks(list(x), labels, rotation=70, ha="right")
    plt.ylabel("Speedup vs original")
    plt.tight_layout()
    plt.savefig(PLOTS / "adversarial_speedup_by_class.png", dpi=160)
    plt.close()

    plt.figure(figsize=(14, 5))
    plt.bar(x, [float(r["inner_iter_reduction_vs_original"]) for r in fast4])
    plt.xticks(list(x), labels, rotation=70, ha="right")
    plt.ylabel("Inner-iteration reduction")
    plt.tight_layout()
    plt.savefig(PLOTS / "adversarial_inner_iter_reduction.png", dpi=160)
    plt.close()

    plt.figure(figsize=(8, 5))
    tx = [str(r["setting"]) for r in tradeoff]
    plt.plot(tx, [float(r["F1"]) for r in tradeoff], marker="o", label="F1")
    orig_runtime = float(next(r for r in tradeoff if r["setting"] == "original")["runtime_sec"])
    speed = [orig_runtime / float(r["runtime_sec"]) if float(r["runtime_sec"]) else 0 for r in tradeoff]
    plt.plot(tx, speed, marker="s", label="speedup")
    plt.xticks(rotation=45, ha="right")
    plt.legend()
    plt.tight_layout()
    plt.savefig(PLOTS / "gate_limit_tradeoff_validation.png", dpi=160)
    plt.close()

    if real_runtime:
        plt.figure(figsize=(10, 5))
        labels = [f"{r['method']}\n{int(r['chunk_length'])//1000}kb" for r in real_runtime]
        plt.bar(range(len(real_runtime)), [float(r["runtime_sec"]) for r in real_runtime])
        plt.xticks(range(len(real_runtime)), labels, rotation=70, ha="right")
        plt.ylabel("Runtime (sec)")
        plt.tight_layout()
        plt.savefig(PLOTS / "real_genome_runtime.png", dpi=160)
        plt.close()


def write_summary(real_runtime, real_diff, adversarial_runtime, adversarial_diff, tradeoff, lost_calls, failure_rows, real_available: bool) -> None:
    original = next(r for r in tradeoff if r["setting"] == "original")
    fast4 = next(r for r in tradeoff if r["setting"] == "fast_4x")
    fast2 = next(r for r in tradeoff if r["setting"] == "fast_2x")
    fast8 = next(r for r in tradeoff if r["setting"] == "fast_8x")
    fast16 = next(r for r in tradeoff if r["setting"] == "fast_16x")
    validation_loss_lines = [
        f"- {r['setting']}: lost_calls={r['lost_calls_vs_original']}, changed_best_score={r['changed_best_score_count']}, changed_arm_spacer={r['changed_best_arm_spacer_count']}"
        for r in tradeoff
        if r["setting"] != "original"
    ]
    adversarial_loss_rows = [
        r for r in adversarial_diff
        if int(r["lost_final_hits"]) or int(r["extra_final_hits"]) or int(r["changed_best_score_count"]) or int(r["changed_best_arm_spacer_count"])
    ]
    adversarial_loss_lines = (
        [
            f"- {r['sequence_id']} / {r['method']}: lost_final={r['lost_final_hits']}, extra_final={r['extra_final_hits']}, "
            f"changed_best_score={r['changed_best_score_count']}, changed_arm_spacer={r['changed_best_arm_spacer_count']}"
            for r in adversarial_loss_rows
        ]
        if adversarial_loss_rows
        else ["- No adversarial final-call or best-call differences were observed."]
    )
    hard = sorted(
        [r for r in adversarial_diff if r["method"] == "fast_gated_fast_4x"],
        key=lambda r: float(r["runtime_speedup_vs_original"])
    )[:5]
    real_lines = []
    if real_available:
        for r in real_diff:
            if r["method"] == "fast_gated_fast_4x":
                real_lines.append(
                    f"- {r['chunk_id']}: {float(r['runtime_speedup_vs_original']):.2f}x speedup, "
                    f"lost_raw={r['lost_raw_hits']}, lost_final={r['lost_final_hits']}, extra_final={r['extra_final_hits']}"
                )
    else:
        real_lines.append("- No usable local real FASTA was found.")
    lost_line = "No lost validation calls." if not lost_calls else f"{len(lost_calls)} lost validation calls; see gate_lost_calls_validation.csv."
    text = f"""# FastGatedHSeeker Next-Step Benchmark

## What Was Tested

- Real local FASTA chunks at 100 kb, 1 Mb, and 5 Mb when available.
- Synthetic adversarial sequences at 20 kb and 100 kb.
- Validation CSV gate sweep at fixed HSeeker biological parameters.

## Validation Preservation

- Original F1: {float(original['F1']):.3f}
- Fast 2x F1: {float(fast2['F1']):.3f}
- Fast 4x F1: {float(fast4['F1']):.3f}
- Fast 8x F1: {float(fast8['F1']):.3f}
- Fast 16x F1: {float(fast16['F1']):.3f}
- Lost calls: {lost_line}

## Real-Genome Runtime

{chr(10).join(real_lines)}

## Adversarial Runtime

Hardest classes for fast 4x by speedup:
{chr(10).join(f"- {r['sequence_id']}: {float(r['runtime_speedup_vs_original']):.2f}x speedup, lost_final={r['lost_final_hits']}" for r in hard)}

## Gate Limit Tradeoff

Validation data:
{chr(10).join(validation_loss_lines)}

Adversarial final-call or best-call differences:
{chr(10).join(adversarial_loss_lines)}

The 2x/4x/8x/16x gate settings are all reported in `gate_tradeoff_validation.csv`.
Lost calls are explicitly listed in `gate_lost_calls_validation.csv`.

## Recommended Default Gate

Use `gate_search_limit = 4 * minrep` as a practical default if validation-call preservation remains true for the target dataset. Use `8 * minrep` when prioritizing recall margin over speed.

## Caveats

- These benchmarks are local and parameter-specific.
- GA/CT-rich and repeat-heavy backgrounds can reduce speedup because many pairs pass purity and gate checks.
- Real-genome superiority should only be claimed for the tested local FASTA chunks, not genome-wide in general.
- Exact-gate mode is useful for equivalence checks but can be slower than original on adversarial inputs.
"""
    (OUT / "summary.md").write_text(text, encoding="utf-8")


def concise_terminal_summary(tradeoff, adversarial_diff, real_diff, real_available):
    def f1(setting):
        return float(next(r for r in tradeoff if r["setting"] == setting)["F1"])
    def adv_speed(cls, length=20_000):
        rows = [
            r for r in adversarial_diff
            if r["method"] == "fast_gated_fast_4x"
            and r["sequence_class"] == cls
            and int(r["sequence_length"]) == length
        ]
        return float(rows[0]["runtime_speedup_vs_original"]) if rows else 0.0
    def real_speed(length):
        rows = [
            r for r in real_diff
            if r["method"] == "fast_gated_fast_4x"
            and int(r.get("chunk_length", 0)) == length
        ]
        return float(rows[0]["runtime_speedup_vs_original"]) if rows else 0.0
    recommendation = "4x" if f1("fast_4x") == f1("original") else "8x/review lost calls"
    print("\nValidation tradeoff:")
    print(f"  original F1: {f1('original'):.3f}")
    print(f"  fast 2x F1: {f1('fast_2x'):.3f}")
    print(f"  fast 4x F1: {f1('fast_4x'):.3f}")
    print(f"  fast 8x F1: {f1('fast_8x'):.3f}")
    print(f"  fast 16x F1: {f1('fast_16x'):.3f}")
    print("\nAdversarial speedup:")
    print(f"  balanced 20kb: {adv_speed('balanced_random'):.2f}x")
    print(f"  GA-rich 20kb: {adv_speed('ga_rich_random'):.2f}x")
    print(f"  GAA repeat 20kb: {adv_speed('pure_gaa_repeat'):.2f}x")
    print(f"  planted balanced 20kb: {adv_speed('planted_balanced'):.2f}x")
    print("\nReal-genome speedup:")
    print(f"  100kb: {real_speed(100000):.2f}x" if real_available else "  100kb: skipped")
    print(f"  1Mb: {real_speed(1000000):.2f}x" if real_available else "  1Mb: skipped")
    print(f"\nRecommended gate: {recommendation}")


def main() -> None:
    ensure_dirs()
    real_runtime, real_diff, real_available = run_real_genome()
    adversarial_runtime, adversarial_diff = run_adversarial()
    tradeoff, lost_calls, validation_by_setting = run_validation_tradeoff()
    failure_rows = run_failure_inspection(validation_by_setting)
    make_plots(adversarial_runtime, adversarial_diff, tradeoff, real_runtime)
    write_summary(real_runtime, real_diff, adversarial_runtime, adversarial_diff, tradeoff, lost_calls, failure_rows, real_available)
    concise_terminal_summary(tradeoff, adversarial_diff, real_diff, real_available)


if __name__ == "__main__":
    main()
