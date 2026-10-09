#!/usr/bin/env python3
"""Which of the authors' ssDNA+ H-DNA structures (Kouzine 2017) does HSeeker miss, and why?

For every ssDNA+ H-DNA motif not overlapped by an HSeeker locus (CLI defaults), the motif
+-50 bp is re-scanned with one setting relaxed at a time; the first relaxation that yields a
hit overlapping the motif is recorded as the reason. Sequence features of caught and missed
motifs are compared, and each motif is checked against the authors' own SIDD calls
(stress-induced duplex destabilisation, i.e. melting-prone DNA) to ask whether the missed
motifs' ssDNA signal could come from duplex melting rather than a triplex.

Usage:  python3 analysis/scripts/kouzine_missed_hdna.py --genome mm9|hg19
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_kouzine_comparison as rk  # noqa: E402
import run_kouzine_published as rp  # noqa: E402

import hseeker  # noqa: E402

FLANK = 50
NEAR = 100  # a missed motif with an HSeeker locus this close sits in a tract HSeeker does call
CLI = dict(minrep=10, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10,
           at_threshold=0.8, filter_homopolymers=True, overlap_strategy="greedy")
LADDER = [  # (reason if this step rescues, cumulative relaxation)
    ("AT filter (motif >= 80% A/T)", dict(at_threshold=None)),
    ("homopolymer filter", dict(filter_homopolymers=False)),
    ("arm length 8-9 (< minrep 10)", dict(minrep=8)),
    ("arm length 6-7", dict(minrep=6)),
    ("mirror mismatch 10-20%", dict(mismatch=0.20)),
    ("purine/pyrimidine purity 80-90%", dict(purity=0.80)),
]


def features(s: str) -> dict:
    s = s.upper()
    n = len(s)
    at = (s.count("A") + s.count("T")) / n
    pur = (s.count("A") + s.count("G")) / n
    run, best = 1, 1
    for i in range(1, n):
        run = run + 1 if s[i] == s[i - 1] else 1
        best = max(best, run)
    unit = None
    for k in (1, 2, 3, 4):
        if n > k and sum(s[i] == s[i + k] for i in range(n - k)) / (n - k) >= 0.9:
            unit = k
            break
    return dict(length=n, at=at, purine_skew=max(pur, 1 - pur), longest_run_frac=best / n,
                period=unit or 0)


def composition_group(s: str) -> str:
    s = s.upper()
    if len(set(s)) == 1:
        return "homopolymer"
    return "A/T>=0.8" if (s.count("A") + s.count("T")) / len(s) >= 0.8 else "A/T<0.8"


def overlaps(hits, a, b) -> bool:
    return any(h["start"] - 1 < b and h["end"] > a for h in hits)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", choices=rp.SETS, default="mm9")
    g = ap.parse_args().genome
    cfg = rp.SETS[g]
    _, hidx = rp.load_hseeker(g)
    called = sorted(m for m in rp.read_bed(rp.AUTH / "nonB_DNA_ssDNA_enriched" / cfg["species"]
                                           / f"{cfg['called']}.bed.gz") if m[0] in set(rk.GENOMES[g]["primary"]))
    genome = dict(rk.read_fasta(rk.KOUZINE / "work" / g / f"{g}.primary.fa", {m[0] for m in called}))
    caught, missed = [], []
    for m in called:
        (caught if hidx.any(*m) else missed).append(m)

    auth = rp.AUTH
    sidd_pred = rp.Index(rp.read_bed(auth / "nonB_DNA_predicted" / cfg["species"] / "SIDD.bed.gz"))
    sidd_call = rp.Index(rp.read_bed(auth / "nonB_DNA_ssDNA_enriched" / cfg["species"]
                                     / f"{cfg['called'].replace('H-DNA', 'SIDD')}.bed.gz"))
    sidd = collections.defaultdict(lambda: [0, 0, 0])  # n, overlaps SIDD motif, overlaps ssDNA+ SIDD

    reasons = collections.Counter()
    near = collections.Counter()
    examples = collections.defaultdict(list)
    feats = {"caught": [], "missed": []}
    for label, group in (("caught", caught), ("missed", missed)):
        for c, s, e in group:
            seq = genome[c][s:e].decode()
            feats[label].append(features(seq))
            row = sidd[f"{label} | {composition_group(seq)}"]
            row[0] += 1
            row[1] += sidd_pred.any(c, s, e)
            row[2] += sidd_call.any(c, s, e)
            if label == "caught":
                continue
            a0 = max(0, s - FLANK)
            window = genome[c][a0:e + FLANK].decode()
            ms, me = s - a0, e - a0  # motif in window, 0-based half-open
            params = dict(CLI)
            reason = "not recovered by any relaxation"
            if overlaps(hseeker.scan_sequence(window, **params), ms, me):
                # with +-50 bp context the motif's candidate survives; in the genome scan greedy
                # overlap removal keeps a longer neighbouring candidate that does not cover the motif
                reason = "overlap removal (a longer neighbouring locus was kept)"
            else:
                for name, relax in LADDER:
                    params.update(relax)
                    if overlaps(hseeker.scan_sequence(window, **params), ms, me):
                        reason = name
                        break
            reasons[reason] += 1
            near[reason] += hidx.any(c, s - NEAR, e + NEAR)
            if len(examples[reason]) < 6:
                examples[reason].append(f"{c}:{s}-{e} {seq}")

    def summary(rows):
        arr = {k: np.array([r[k] for r in rows]) for k in rows[0]}
        out = {k: dict(median=float(np.median(v)), p10=float(np.percentile(v, 10)), p90=float(np.percentile(v, 90)))
               for k, v in arr.items() if k != "period"}
        per = collections.Counter(int(x) for x in arr["period"])
        out["repeat_period_share"] = {("none" if k == 0 else f"{k}-mer"): round(v / len(rows), 3)
                                      for k, v in sorted(per.items())}
        out["at_ge_0.8_share"] = round(float((arr["at"] >= 0.8).mean()), 3)
        return out

    res = dict(genome=g, hseeker_config="cli", called=len(called), caught=len(caught), missed=len(missed),
               missed_fraction=len(missed) / len(called),
               missed_by_reason={k: dict(n=v, share=round(v / len(missed), 3),
                                         hseeker_locus_within_100bp=round(near[k] / v, 3))
                                 for k, v in reasons.most_common()},
               sidd_overlap={k: dict(n=n, sidd_motif=round(p / n, 3), ssDNA_plus_sidd=round(q / n, 3))
                             for k, (n, p, q) in sorted(sidd.items())},
               examples=examples, features=dict(caught=summary(feats["caught"]), missed=summary(feats["missed"])))
    out = rk.ROOT / "analysis" / "results" / "kouzine_v1" / "published"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{g}_missed_analysis.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "examples"}, indent=1))
    for k, v in examples.items():
        print(f"\n[{k}]")
        for x in v[:4]:
            print("  ", x)


if __name__ == "__main__":
    main()
