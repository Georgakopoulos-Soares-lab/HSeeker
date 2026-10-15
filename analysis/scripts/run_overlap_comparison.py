#!/usr/bin/env python3
"""Greedy (longest-arm) vs score-based overlap selection (reviewer comment R3.10).

Regenerates every number in the R3.10 response on the current code:
  1. the two illustrative inputs ("agreement", "difference") and the GAA x 20 overlap chain;
  2. the benchmark: records with several overlapping candidates, and whether greedy and
     score-based selection return the same intervals and the same calls;
  3. hg38 chr1 (optional, --chr1): loci per mode, identical coordinates, greedy loci covered by
     a score-based locus, score-only loci.

Configuration: the shipped CLI defaults (minrep 8, maxrep 1000, maxspacer 10, purity 0.90,
mismatch 0.10, AT filter 0.80, homopolymer filter, scoring on), unless stated otherwise.

Usage:
  python3 analysis/scripts/run_overlap_comparison.py [--chr1 path/to/chr1.fa] [--workers 16]
Outputs: analysis/results/overlap_selection/overlap_selection_comparison.csv, examples.json,
         chr1_comparison.json (with --chr1), metadata.json
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import subprocess
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import hseeker  # noqa: E402

OUT = ROOT / "analysis" / "results" / "overlap_selection"
DEFAULTS = dict(minrep=8, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10,
                at_threshold=0.80, filter_homopolymers=True, score=True)
EXAMPLE_PARAMS = dict(minrep=7, maxrep=25, maxspacer=5, purity=0.9, mismatch=0.2,
                      at_threshold=0.8, filter_homopolymers=True, score=True)
EXAMPLES = {"agreement": "GAGAAAAGAGGGAGAAAAAGGAA",
            "difference": "GGAAGAAGAAGGAGGAGAGAGGGAAGGAA"}
CHAIN_PARAMS = dict(minrep=10, maxrep=15, maxspacer=10, purity=0.9, mismatch=0.1,
                    at_threshold=0.8, filter_homopolymers=True, score=True)


def brief(h: dict) -> dict:
    return dict(start=h["start"], end=h["end"], arm=h["arm_length"], spacer=h["spacer_length"],
                full_sequence=h["full_sequence"],
                total_score=None if h.get("total_score") is None else round(h["total_score"], 3))


def modes(seq: str, params: dict) -> dict[str, list[dict]]:
    return {"all": hseeker.scan_sequence(seq, remove_overlaps=False, **params),
            "greedy": hseeker.scan_sequence(seq, overlap_strategy="greedy", **params),
            "score": hseeker.scan_sequence(seq, overlap_strategy="score", **params)}


def has_overlapping_pair(hits: list[dict]) -> bool:
    iv = sorted((h["start"], h["end"]) for h in hits)
    return any(b[0] <= a[1] for a, b in zip(iv, iv[1:]))


def benchmark_rows(path: Path, subset: str) -> list[dict]:
    rows = [r for r in csv.DictReader(open(path)) if r["curation_decision"] == "kept"]
    if subset == "experimental":
        rows = [r for r in rows if r.get("record_type", "experimental") == "experimental"]
    out = []
    for r in rows:
        m = modes(r["sequence_5to3"].upper(), DEFAULTS)
        g = [(h["start"], h["end"]) for h in m["greedy"]]
        s = [(h["start"], h["end"]) for h in m["score"]]
        out.append(dict(subset=subset, record_id=r["record_id"], sequence_name=r["sequence_name"],
                        label=r["label"], candidates=len(m["all"]),
                        overlapping_candidates=has_overlapping_pair(m["all"]),
                        greedy_hits=len(g), score_hits=len(s),
                        greedy_intervals=";".join(f"{a}-{b}" for a, b in sorted(g)),
                        score_intervals=";".join(f"{a}-{b}" for a, b in sorted(s)),
                        identical_intervals=sorted(g) == sorted(s),
                        call_greedy=bool(g), call_score=bool(s)))
    return out


def interval_index(loci):
    by = {}
    for c, s, e in loci:
        by.setdefault(c, []).append((s, e))
    for v in by.values():
        v.sort()
    return by


def any_overlap(by, c, s, e) -> bool:
    import bisect
    v = by.get(c, [])
    i = bisect.bisect_right(v, (e, float("inf")))
    j = i - 1
    while j >= 0:
        if v[j][1] >= s:
            return True
        if v[j][0] < s - 100_000:   # loci are far shorter than this
            break
        j -= 1
    return False


def chr1_comparison(fa: Path, workers: int) -> dict:
    res = {}
    for label, params in (("cli_defaults", DEFAULTS),
                          ("no_filters", {**DEFAULTS, "at_threshold": None,
                                          "filter_homopolymers": False})):
        runs = {}
        for strategy in ("greedy", "score"):
            t = time.perf_counter()
            hits = hseeker.scan_fasta_parallel(str(fa), workers=workers,
                                               overlap_strategy=strategy, **params)
            runs[strategy] = dict(seconds=round(time.perf_counter() - t, 1),
                                  loci=[(h["seq_id"], h["start"], h["end"]) for h in hits],
                                  motif={(h["seq_id"], h["start"], h["end"]): h["full_sequence"]
                                         for h in hits})
        g, s = set(runs["greedy"]["loci"]), set(runs["score"]["loci"])
        gi, si = interval_index(g), interval_index(s)
        covered = sum(any_overlap(si, *x) for x in g)
        score_only = sum(not any_overlap(gi, *x) for x in s)
        res[label] = dict(params={k: v for k, v in params.items()},
                          greedy_loci=len(g), score_loci=len(s),
                          identical_coordinates=len(g & s),
                          identical_fraction_of_greedy=round(len(g & s) / len(g), 4),
                          greedy_loci_overlapped_by_score=covered,
                          greedy_covered_fraction=round(covered / len(g), 4),
                          score_only_loci=score_only,
                          wall_seconds=dict(greedy=runs["greedy"]["seconds"],
                                            score=runs["score"]["seconds"]),
                          wall_note="single run on a shared node; not a runtime benchmark")
        print(label, json.dumps({k: v for k, v in res[label].items() if k != "params"}), flush=True)
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--chr1", type=Path)
    ap.add_argument("--workers", type=int, default=16)
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)

    ex = {name: {k: [brief(h) for h in v] for k, v in modes(seq, EXAMPLE_PARAMS).items()}
          for name, seq in EXAMPLES.items()}
    ex["chain_GAAx20"] = {k: [brief(h) for h in v]
                          for k, v in modes("GAA" * 20, CHAIN_PARAMS).items()}
    examples = dict(example_params=EXAMPLE_PARAMS, chain_params=CHAIN_PARAMS,
                    inputs={**EXAMPLES, "chain_GAAx20": "GAA" * 20}, results=ex)
    (OUT / "examples.json").write_text(json.dumps(examples, indent=1))

    csv_path = ROOT / "hdna_benchmark_balanced_v3.csv"
    rows = benchmark_rows(csv_path, "balanced") + benchmark_rows(csv_path, "experimental")
    with open(OUT / "overlap_selection_comparison.csv", "w", newline="") as fo:
        w = csv.DictWriter(fo, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    summary = {}
    for subset in ("balanced", "experimental"):
        rs = [r for r in rows if r["subset"] == subset]
        multi = [r for r in rs if r["overlapping_candidates"]]
        summary[subset] = dict(records=len(rs), with_overlapping_candidates=len(multi),
                               identical_intervals_among_those=sum(r["identical_intervals"] for r in multi),
                               differing_intervals=[r["record_id"] for r in multi if not r["identical_intervals"]],
                               calls_differ=[r["record_id"] for r in rs if r["call_greedy"] != r["call_score"]])
    print(json.dumps(summary, indent=1))
    meta = dict(date=time.strftime("%Y-%m-%d %H:%M:%S"), defaults=DEFAULTS, benchmark=str(csv_path.name),
                benchmark_sha256=hashlib.sha256(csv_path.read_bytes()).hexdigest(),
                git_commit=subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"],
                                          capture_output=True, text=True).stdout.strip(),
                git_dirty=bool(subprocess.run(["git", "-C", str(ROOT), "status", "--porcelain"],
                                              capture_output=True, text=True).stdout.strip()),
                benchmark_summary=summary)
    if args.chr1:
        meta["chr1"] = chr1_comparison(args.chr1, args.workers)
        meta["chr1_fasta_sha256"] = hashlib.sha256(args.chr1.read_bytes()).hexdigest()
        (OUT / "chr1_comparison.json").write_text(json.dumps(meta["chr1"], indent=1))
    (OUT / "metadata.json").write_text(json.dumps(meta, indent=1))


if __name__ == "__main__":
    main()
