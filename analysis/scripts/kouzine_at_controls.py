#!/usr/bin/env python3
"""Do A/T-rich mirror repeats carry more ssDNA signal than A/T-rich DNA without one? (X8, exploratory)

HSeeker's AT filter drops motifs >= 80 % A/T. The authors of Kouzine et al. 2017 call some such
motifs ssDNA+. If that signal comes from composition (duplex melting of A/T-rich DNA, KMnO4's
preference for unpaired thymine), then A/T-rich DNA *without* a mirror repeat should carry the same
signal, and the data do not conflict with the filter. If A/T-rich mirror repeats carry more signal
than composition-matched DNA, the filter discards structure-associated signal.

For a seeded sample of the authors' analysed, mappable H-DNA motifs in each composition group
(A/T < 0.8, A/T >= 0.8, single-base), control windows are drawn that have
  - the same chromosome and the same length,
  - exactly the same A/T count,
  - "C" set: purine-strand asymmetry max(A+G, C+T) within +-1 count; "G" set: A/T count only,
  - no N, no overlap with any of the authors' predicted H-DNA motifs (+-100 bp, the tag extension),
  - mean 36-mer mappability >= 0.5, as for the motifs.
Two neighbourhoods: "local" (all positions within +-50 kb, so the same chromatin and
transcriptional context) and "chromosome" (random positions on the same chromosome).
Signal: the authors' ssDNA wiggle; ratio of means motif/controls with a bootstrap CI over motifs.

Usage:  python3 analysis/scripts/kouzine_at_controls.py --genome mm9|hg19
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
from kouzine_missed_hdna import composition_group  # noqa: E402

N_PER_GROUP = 3000
N_CONTROLS = 10
CANDIDATES = 30          # drawn before the mappability check
PAD = 100                # exclusion around predicted H-DNA motifs (authors' tag extension)
LOCAL = 50_000
DRAWS = 20_000
ROUNDS = 5
MAP_MIN = 0.5


def exclusion(predicted):
    by = collections.defaultdict(list)
    for c, s, e in predicted:
        by[c].append((s - PAD, e + PAD))
    out = {}
    for c, v in by.items():
        v.sort()
        st = np.array([x[0] for x in v], dtype=np.int64)
        out[c] = (st, np.maximum.accumulate(np.array([x[1] for x in v], dtype=np.int64)))
    return out


def candidates(p, L, at_n, asym, pref, ex, use_asym):
    at_pref, pu_pref, n_pref = pref
    at_c = at_pref[p + L] - at_pref[p]
    ok = (n_pref[p + L] - n_pref[p] == 0) & (at_c == at_n)
    if use_asym:
        pu = pu_pref[p + L] - pu_pref[p]
        ok &= np.abs(np.maximum(pu, L - pu) - asym) <= 1
    st, rmax = ex
    i = np.searchsorted(st, p + L) - 1
    ok &= ~((i >= 0) & (rmax[np.maximum(i, 0)] > p))
    return p[ok]


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", choices=rp.SETS, default="mm9")
    g = ap.parse_args().genome
    cfg = rp.SETS[g]
    rng = np.random.default_rng(rk.SEED)
    primary = set(rk.GENOMES[g]["primary"])
    predicted = [m for m in rp.read_bed(rp.AUTH / "nonB_DNA_predicted" / cfg["species"] / "H-DNA.bed.gz")
                 if m[0] in primary]
    called = {m for m in rp.read_bed(rp.AUTH / "nonB_DNA_ssDNA_enriched" / cfg["species"]
                                     / f"{cfg['called']}.bed.gz") if m[0] in primary}
    analysed = rp.analysed_universe(g, predicted)
    mp = rp.mappability(g, analysed)
    genome = dict(rk.read_fasta(rk.KOUZINE / "work" / g / f"{g}.primary.fa", primary))

    groups = collections.defaultdict(list)
    for m in sorted(set(analysed)):
        if mp[m] >= MAP_MIN:
            c, s, e = m
            groups[composition_group(genome[c][s:e].decode())].append(m)
    sample = {}
    for grp, v in sorted(groups.items()):
        ix = rng.choice(len(v), size=min(N_PER_GROUP, len(v)), replace=False)
        sample[grp] = sorted(v[i] for i in ix)

    ex = exclusion(predicted)
    # controls[(grp, scope, set)][motif] -> candidate starts
    cand = collections.defaultdict(dict)
    by_chrom = collections.defaultdict(list)
    for grp, v in sample.items():
        for m in v:
            by_chrom[m[0]].append((grp, m))
    for c in sorted(by_chrom):
        arr = np.frombuffer(genome[c], dtype=np.uint8)
        pref = tuple(np.concatenate(([0], np.cumsum(x, dtype=np.int32))) for x in
                     ((arr == 65) | (arr == 84), (arr == 65) | (arr == 71), arr == 78))
        n = len(arr)
        for grp, (_, s, e) in by_chrom[c]:
            L = e - s
            seq = arr[s:e]
            at_n = int(((seq == 65) | (seq == 84)).sum())
            pu = int(((seq == 65) | (seq == 71)).sum())
            asym = max(pu, L - pu)
            for use_asym, cset in ((True, "C"), (False, "G")):
                p = np.arange(max(0, s - LOCAL), min(n - L, e + LOCAL), dtype=np.int64)
                hit = candidates(p, L, at_n, asym, pref, ex[c], use_asym)
                cand[(grp, "local", cset)][(c, s, e)] = rng.permutation(hit)[:CANDIDATES]
                found = []
                for _ in range(ROUNDS):
                    p = rng.integers(0, n - L, DRAWS)
                    found.extend(candidates(p, L, at_n, asym, pref, ex[c], use_asym).tolist())
                    if len(found) >= CANDIDATES:
                        break
                cand[(grp, "chromosome", cset)][(c, s, e)] = np.array(found[:CANDIDATES], dtype=np.int64)
        del pref, arr

    windows = sorted({(m[0], int(p), int(p) + m[2] - m[1]) for d in cand.values()
                      for m, ps in d.items() for p in ps})
    cmp = rp.mappability(g, windows, tag="at_controls")
    controls = {key: {m: [(m[0], int(p), int(p) + m[2] - m[1]) for p in ps
                          if cmp[(m[0], int(p), int(p) + m[2] - m[1])] >= MAP_MIN][:N_CONTROLS]
                      for m, ps in d.items()} for key, d in cand.items()}

    lengths = {c: len(s) for c, s in genome.items()}
    del genome
    res = dict(genome=g, n_per_group=N_PER_GROUP, n_controls=N_CONTROLS, pad=PAD, local_window=LOCAL,
               mappability_min=MAP_MIN,
               groups={grp: dict(motifs_mappable=len(groups[grp]), sampled=len(v),
                                 ssDNA_plus_rate_in_sample=round(sum(m in called for m in v) / len(v), 4))
                       for grp, v in sample.items()},
               signal={})
    for label, wig in cfg["wig"].items():
        cs = rp.wig_track(rp.AUTH / "ssDNA_wiggle" / cfg["species"] / f"{wig}.wig.gz", lengths)
        total_m = sum(v[-1] for v in cs.values()) / 1e6
        out = {}
        for (grp, scope, cset), d in sorted(controls.items()):
            use = [m for m, cl in d.items() if cl]
            if not use:
                out[f"{grp} | {scope} | {cset}"] = dict(motifs_with_controls=0)
                continue
            lv = np.array([rp.wig_signal(cs, *m, total_m) for m in use])
            cv = [np.array([rp.wig_signal(cs, *w, total_m) for w in d[m]]) for m in use]
            est, lo, hi = rk.boot_ratio(lv, cv, rng)
            cm = np.array([v.mean() for v in cv])
            out[f"{grp} | {scope} | {cset}"] = dict(
                motifs_with_controls=len(use), controls_per_motif=round(float(np.mean([len(v) for v in cv])), 2),
                motif_mean=float(lv.mean()), control_mean=float(cm.mean()),
                ratio_of_means=round(float(est), 3), ratio_ci95=[round(lo, 3), round(hi, 3)],
                share_motif_above_its_controls=round(float((lv > cm).mean()), 3))
        res["signal"][label] = out
        del cs
    path = rk.ROOT / "analysis" / "results" / "kouzine_v1" / "published" / f"{g}.at_controls.json"
    path.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
