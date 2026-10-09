#!/usr/bin/env python3
"""HSeeker vs the authors' published Kouzine et al. 2017 calls (X8 primary, deviation 4).

Data: https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/nonbdna.html
  nonB_DNA_predicted/      predicted non-B motifs (SMnB), mm9 / hg19
  nonB_DNA_ssDNA_enriched/ motifs the authors called ssDNA+ (formed non-B structures)
  ssDNA_wiggle/            authors' ssDNA-seq coverage (100 bp tag extension, 10 bp bins)

Tests (pre-specified in RTR/evidence/X8_kouzine_analysis_design_2026-10-09.md):
  P1 recall       fraction of the authors' ssDNA+ H-DNA overlapped by an HSeeker locus
  P2 selectivity  among the authors' analysed H-DNA motifs (their SINE filter re-applied:
                  no SINE within a 500 bp window centred on the motif), odds that a motif is
                  ssDNA+ when HSeeker calls it vs when it does not; score of ssDNA+ vs ssDNA-
  P3 enrichment   authors' ssDNA signal at HSeeker loci vs composition-matched controls
                  (controls from run_kouzine_comparison.py), activated vs resting

Usage:  python3 analysis/scripts/run_kouzine_published.py --genome mm9|hg19
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_kouzine_comparison as rk  # noqa: E402

AUTH = rk.KOUZINE / "authors"
SINE_HALF_WINDOW = 250
BIN = 10
SETS = {
    "mm9": dict(species="mouse_mm9", called="actB_ssDNA_enriched_H-DNA",
                wig=dict(activated="actB_ssDNA", resting="resB_ssDNA"),
                analysed_in_paper=320_585 + 20_846, called_in_paper=16_876 + 233),
    "hg19": dict(species="human_hg19", called="Raji_ssDNA_enriched_H-DNA",
                 wig=dict(raji="Raji_ssDNA"),
                 analysed_in_paper=144_643 + 13_727, called_in_paper=8_008 + 348),
}


def read_bed(path: Path) -> list[tuple[str, int, int]]:
    with gzip.open(path, "rt") as fh:
        return [(f[0], int(f[1]), int(f[2])) for f in (l.split() for l in fh if not l.startswith("track"))]


class Index:
    """Overlap queries against possibly overlapping intervals (sorted starts + running max end)."""

    def __init__(self, intervals, values=None):
        self.by = {}
        for k, (c, s, e) in enumerate(intervals):
            self.by.setdefault(c, []).append((s, e, k))
        self.values = values
        self.arr = {}
        for c, v in self.by.items():
            v.sort()
            st = np.array([x[0] for x in v])
            en = np.array([x[1] for x in v])
            self.arr[c] = (st, en, np.maximum.accumulate(en), np.array([x[2] for x in v]))

    def any(self, c, s, e) -> bool:
        if c not in self.arr:
            return False
        st, _, rmax, _ = self.arr[c]
        i = np.searchsorted(st, e) - 1
        return bool(i >= 0 and rmax[i] > s)

    def hits(self, c, s, e) -> list[int]:
        """Indices of intervals overlapping [s, e) (linear scan back while possible)."""
        if c not in self.arr:
            return []
        st, en, rmax, ix = self.arr[c]
        i = np.searchsorted(st, e) - 1
        out = []
        while i >= 0 and rmax[i] > s:
            if en[i] > s:
                out.append(int(ix[i]))
            i -= 1
        return out


def sine_index(genome: str) -> Index:
    rows = []
    files = sorted((AUTH / "mm9_rmsk").glob("chr*_rmsk.txt.gz")) if genome == "mm9" else [AUTH / "hg19_rmsk.txt.gz"]
    for f in files:
        with gzip.open(f, "rt") as fh:
            for l in fh:
                x = l.rstrip("\n").split("\t")
                if x[11] == "SINE":
                    rows.append((x[5], int(x[6]), int(x[7])))
    return Index(rows)


def load_hseeker(genome: str):
    loci = rk.load_loci(rk.KOUZINE / "work" / genome / "cli" / "loci.tsv")
    return loci, Index([(l["chrom"], l["start"], l["end"]) for l in loci])


def wig_track(path: Path, lengths: dict[str, int]) -> dict[str, np.ndarray]:
    """Prefix sums of the wiggle coverage per chromosome at 10 bp bins."""
    arrs = {c: np.zeros(L // BIN + 2, dtype=np.float64) for c, L in lengths.items()}
    cur, step = None, BIN
    with gzip.open(path, "rt") as fh:
        for l in fh:
            if l[0].isdigit():
                if cur is not None:
                    p, v = l.split()
                    cur[(int(p) - 1) // BIN] += float(v)
            elif l.startswith("variableStep"):
                f = dict(kv.split("=") for kv in l.split()[1:])
                cur = arrs.get(f["chrom"])
                step = int(f.get("span", BIN))
                assert step == BIN, step
    return {c: np.concatenate(([0.0], np.cumsum(a))) for c, a in arrs.items()}


def wig_signal(cs: dict[str, np.ndarray], c: str, s: int, e: int, total_m: float) -> float:
    a, b = s // BIN, (e - 1) // BIN + 1
    return (cs[c][b] - cs[c][a]) / (b - a) / total_m


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", choices=SETS, required=True)
    args = ap.parse_args()
    g = args.genome
    cfg = SETS[g]
    res: dict = dict(genome=g, hseeker_config="cli", hseeker_params=rk.CONFIGS["cli"])
    loci, hidx = load_hseeker(g)
    res["hseeker_loci"] = len(loci)

    predicted = read_bed(AUTH / "nonB_DNA_predicted" / cfg["species"] / "H-DNA.bed.gz")
    called = set(read_bed(AUTH / "nonB_DNA_ssDNA_enriched" / cfg["species"] / f"{cfg['called']}.bed.gz"))
    primary = set(rk.GENOMES[g]["primary"])
    predicted = [m for m in predicted if m[0] in primary]
    called = {m for m in called if m[0] in primary}

    # --- P1 recall
    rec = sum(hidx.any(*m) for m in called)
    res["P1_recall"] = dict(called_ssDNA_plus=len(called), called_in_paper=cfg["called_in_paper"],
                            overlapped_by_hseeker=rec, fraction=rec / len(called))

    # --- P2 selectivity on the authors' analysed universe (SINE filter re-applied)
    sidx = sine_index(g)
    analysed = [m for m in predicted
                if not sidx.any(m[0], (m[1] + m[2]) // 2 - SINE_HALF_WINDOW, (m[1] + m[2]) // 2 + SINE_HALF_WINDOW)]
    a = b = c = d = 0
    pos_scores, neg_scores = [], []
    for m in analysed:
        hs = hidx.hits(*m)
        plus = m in called
        if hs:
            sc = max((loci[i]["score"] for i in hs if not math.isnan(loci[i]["score"])), default=math.nan)
            (pos_scores if plus else neg_scores).append(sc)
            a += plus
            b += not plus
        else:
            c += plus
            d += not plus
    o, lo, hi = rk.odds_ratio(a, b, c, d)
    fisher_p = stats.fisher_exact([[a, b], [c, d]], alternative="greater").pvalue
    ps, ns = np.array([x for x in pos_scores if not math.isnan(x)]), np.array([x for x in neg_scores if not math.isnan(x)])
    res["P2_selectivity"] = dict(
        analysed_motifs=len(analysed), analysed_in_paper=cfg["analysed_in_paper"],
        called_in_analysed=sum(m in called for m in analysed),
        hseeker_called=dict(ssDNA_plus=a, ssDNA_minus=b, rate=a / (a + b) if a + b else None),
        not_hseeker_called=dict(ssDNA_plus=c, ssDNA_minus=d, rate=c / (c + d) if c + d else None),
        odds_ratio=o, or_ci95=[lo, hi], fisher_p_greater=float(fisher_p),
        score_ssDNA_plus_median=float(np.median(ps)) if len(ps) else None,
        score_ssDNA_minus_median=float(np.median(ns)) if len(ns) else None,
        score_mannwhitney_p_greater=float(stats.mannwhitneyu(ps, ns, alternative="greater").pvalue)
        if len(ps) and len(ns) else None)

    # --- P3 enrichment with the authors' ssDNA coverage (needs the matched controls)
    work = rk.KOUZINE / "work" / g / "cli"
    if not (work / "controls_meta.json").exists():
        res["P3_enrichment_authors_signal"] = "pending: controls not yet sampled"
        print(json.dumps(res, indent=1))
        return
    lengths = {c_: len(s_) for c_, s_ in rk.read_fasta(work.parent / f"{g}.primary.fa", primary)}
    intervals = {}
    for line in open(work / "intervals.bed"):
        c_, s_, e_, name = line.split()
        intervals[name] = (c_, int(s_), int(e_))
    rng = np.random.default_rng(rk.SEED)
    p3 = {}
    for label, wig in cfg["wig"].items():
        cs = wig_track(AUTH / "ssDNA_wiggle" / cfg["species"] / f"{wig}.wig.gz", lengths)
        total_m = sum(v[-1] for v in cs.values()) / 1e6
        sig = {n: wig_signal(cs, *iv, total_m) for n, iv in intervals.items()}
        del cs
        p3[label] = {}
        for cset, prefix in (("asymmetry_matched", "C"), ("gc_matched", "G")):
            pairs = {}
            for n in sorted(intervals):
                if n[0] not in ("L", prefix):
                    continue
                i = int(n[1:].split("_")[0])
                pairs.setdefault(i, {"L": None, "C": []})
                if n[0] == "L":
                    pairs[i]["L"] = n
                else:
                    pairs[i]["C"].append(n)
            use = [i for i, p in pairs.items() if p["L"] and p["C"]]
            lv = np.array([sig[pairs[i]["L"]] for i in use])
            cv = [np.array([sig[n] for n in pairs[i]["C"]]) for i in use]
            est, blo, bhi = rk.boot_ratio(lv, cv, rng)
            flat = np.concatenate(cv)
            p3[label][cset] = dict(loci_used=len(use), locus_mean=float(lv.mean()), control_mean=float(flat.mean()),
                                   ratio_of_means=est, ratio_ci95=[blo, bhi],
                                   mannwhitney_p_greater=float(stats.mannwhitneyu(lv, flat, alternative="greater").pvalue))
    if "activated" in p3 and "resting" in p3:
        for cset in ("asymmetry_matched", "gc_matched"):
            p3[f"activated_over_resting_{cset}"] = (p3["activated"][cset]["ratio_of_means"]
                                                    / p3["resting"][cset]["ratio_of_means"])
    res["P3_enrichment_authors_signal"] = p3
    res["controls_meta"] = json.load(open(work / "controls_meta.json"))
    out = rk.ROOT / "analysis" / "results" / "kouzine_v1" / "published"
    out.mkdir(parents=True, exist_ok=True)
    (out / f"{g}.json").write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
