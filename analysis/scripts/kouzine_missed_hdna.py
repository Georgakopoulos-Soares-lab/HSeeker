#!/usr/bin/env python3
"""Which of the authors' ssDNA+ H-DNA structures (Kouzine 2017) does HSeeker miss, and why?

For every ssDNA+ H-DNA motif not overlapped by an HSeeker locus (configuration --config), the
motif +-50 bp is re-scanned with one setting relaxed at a time; the first relaxation that yields
a hit overlapping the motif is recorded as the reason. Sequence features of caught and missed
motifs are compared, and each motif is checked against the authors' own SIDD calls
(stress-induced duplex destabilisation, i.e. melting-prone DNA) to ask whether the missed
motifs' ssDNA signal could come from duplex melting rather than a triplex.

AT stratification: on the authors' analysed H-DNA motifs (SINE filter re-applied, as in P2),
the ssDNA+ rate by composition (A/T < 0.8, A/T >= 0.8, single-base) and by overlap with a
*predicted* SIDD motif (sequence-based, so not circular with the ssDNA calls).

Usage:  python3 analysis/scripts/kouzine_missed_hdna.py --genome mm9|hg19 [--config cli]
"""

from __future__ import annotations

import argparse
import collections
import json
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_kouzine_comparison as rk  # noqa: E402
import run_kouzine_published as rp  # noqa: E402

import hseeker  # noqa: E402

FLANK = 50
NEAR = 100  # a missed motif with an HSeeker locus this close sits in a tract HSeeker does call
LADDER = [  # (reason if this step rescues, cumulative relaxation); no-op steps are skipped
    ("AT filter (motif >= 80% A/T)", dict(at_threshold=None)),
    ("homopolymer filter", dict(filter_homopolymers=False)),
    ("arm length 8-9 (< minrep 10)", dict(minrep=8)),
    ("arm length 6-7", dict(minrep=6)),
    ("mirror mismatch 10-20%", dict(mismatch=0.20)),
    ("purine/pyrimidine purity 80-90%", dict(purity=0.80)),
]


def relaxes(base: dict, step: dict) -> bool:
    """True if *step* loosens *base* (a filter that is on, or a stricter threshold)."""
    for k, v in step.items():
        cur = base.get(k)
        if k == "at_threshold" and cur is not None:
            return True
        if k == "filter_homopolymers" and cur:
            return True
        if k == "minrep" and cur > v:
            return True
        if k in ("mismatch",) and cur < v:
            return True
        if k == "purity" and cur > v:
            return True
    return False


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
    ap.add_argument("--config", choices=rk.CONFIGS, default="cli")
    args = ap.parse_args()
    g, config = args.genome, args.config
    base = rk.CONFIGS[config]
    ladder = [(name, relax) for name, relax in LADDER if relaxes(base, relax)]
    cfg = rp.SETS[g]
    loci, hidx = rp.load_hseeker(g, config)
    primary = set(rk.GENOMES[g]["primary"])
    called = sorted(m for m in rp.read_bed(rp.AUTH / "nonB_DNA_ssDNA_enriched" / cfg["species"]
                                           / f"{cfg['called']}.bed.gz") if m[0] in primary)
    genome = dict(rk.read_fasta(rk.KOUZINE / "work" / g / f"{g}.primary.fa", primary))
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
            params = dict(base)
            reason = "not recovered by any relaxation"
            if overlaps(hseeker.scan_sequence(window, **params), ms, me):
                # with +-50 bp context the motif's candidate survives; in the genome scan greedy
                # overlap removal keeps a longer neighbouring candidate that does not cover the motif
                reason = "overlap removal (a longer neighbouring locus was kept)"
            else:
                for name, relax in ladder:
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

    # AT stratification on the authors' analysed universe (SINE filter as in P2); repeated on the
    # motifs 36 bp reads can map to, since unmappable motifs cannot be ssDNA+ (sensitivity)
    called_set = set(called)
    predicted = [m for m in rp.read_bed(auth / "nonB_DNA_predicted" / cfg["species"] / "H-DNA.bed.gz")
                 if m[0] in primary]
    analysed = rp.analysed_universe(g, predicted)
    mp = rp.mappability(g, analysed)
    rows = [(composition_group(genome[c][s:e].decode()), sidd_pred.any(c, s, e), (c, s, e) in called_set,
             hidx.any(c, s, e), mp[(c, s, e)]) for c, s, e in analysed]

    def stratify(subset) -> dict:
        strata = collections.defaultdict(lambda: [0, 0, 0])  # n, ssDNA+, HSeeker-called
        for grp, in_sidd, plus, hcall, _ in subset:
            row = strata[f"{grp} | {'SIDD' if in_sidd else 'no SIDD'}"]
            row[0] += 1
            row[1] += plus
            row[2] += hcall
        tot_n = sum(r[0] for r in strata.values())
        tot_p = sum(r[1] for r in strata.values())

        def pooled(grp):
            n = sum(strata[f"{grp} | {x}"][0] for x in ("SIDD", "no SIDD"))
            p = sum(strata[f"{grp} | {x}"][1] for x in ("SIDD", "no SIDD"))
            return p, n - p

        ref = pooled("A/T<0.8")
        by_composition = {}
        for grp in ("A/T<0.8", "A/T>=0.8", "homopolymer"):
            pos, neg = pooled(grp)
            sd, nsd = strata[f"{grp} | SIDD"], strata[f"{grp} | no SIDD"]
            o, lo, hi = rk.odds_ratio(sd[1], sd[0] - sd[1], nsd[1], nsd[0] - nsd[1])
            row = dict(n=pos + neg, ssDNA_plus=pos, ssDNA_plus_rate=round(pos / (pos + neg), 4),
                       share_of_ssDNA_plus=round(pos / tot_p, 3), sidd_share_of_motifs=round(sd[0] / (pos + neg), 3),
                       or_sidd_vs_no_sidd=[round(o, 3), round(lo, 3), round(hi, 3)])
            if grp != "A/T<0.8":
                o, lo, hi = rk.odds_ratio(pos, neg, *ref)
                row["or_vs_at_below_0.8"] = [round(o, 3), round(lo, 3), round(hi, 3)]
                row["fisher_p"] = float(stats.fisher_exact([[pos, neg], list(ref)]).pvalue)
            by_composition[grp] = row
        return dict(motifs=tot_n, ssDNA_plus=tot_p, overall_ssDNA_plus_rate=round(tot_p / tot_n, 4),
                    by_composition=by_composition,
                    strata={k: dict(n=n, ssDNA_plus=p, ssDNA_plus_rate=round(p / n, 4),
                                    hseeker_called_rate=round(h / n, 4))
                            for k, (n, p, h) in sorted(strata.items())})

    at_strat = dict(note="authors' analysed H-DNA motifs (SINE filter re-applied); SIDD = overlaps a "
                         "predicted SIDD motif; mappable_* = mean CRG 36-mer mappability over read starts "
                         "overlapping the motif at or above the threshold",
                    all=stratify(rows))
    for k, t in rp.MAPPABLE.items():
        at_strat[f"mappable_{k}"] = dict(threshold=t, **stratify([r for r in rows if r[4] >= t]))

    res = dict(genome=g, hseeker_config=config, hseeker_params=base, hseeker_loci=len(loci),
               called=len(called), caught=len(caught), missed=len(missed),
               missed_fraction=len(missed) / len(called),
               missed_by_reason={k: dict(n=v, share=round(v / len(missed), 3),
                                         hseeker_locus_within_100bp=round(near[k] / v, 3))
                                 for k, v in reasons.most_common()},
               sidd_overlap={k: dict(n=n, sidd_motif=round(p / n, 3), ssDNA_plus_sidd=round(q / n, 3))
                             for k, (n, p, q) in sorted(sidd.items())},
               at_stratification=at_strat,
               examples=examples, features=dict(caught=summary(feats["caught"]), missed=summary(feats["missed"])))
    out = rk.ROOT / "analysis" / "results" / "kouzine_v1" / "published"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{g}.{config}.missed.json").write_text(json.dumps(res, indent=1))
    print(json.dumps({k: v for k, v in res.items() if k != "examples"}, indent=1))
    for k, v in examples.items():
        print(f"\n[{k}]")
        for x in v[:4]:
            print("  ", x)


if __name__ == "__main__":
    main()
