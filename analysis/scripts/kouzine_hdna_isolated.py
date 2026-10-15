#!/usr/bin/env python3
"""P4: ssDNA signal at HSeeker loci isolated from other sources of ssDNA (X8, design deviation 5).

ssDNA-seq detects unpaired DNA from any cause. P4 removes the other recognised sources from both
the HSeeker loci and their controls, then asks whether the loci still carry more of the authors'
ssDNA signal than composition-matched DNA without a mirror repeat.

Tiers (applied identically to loci and controls):
  A  mappable: mean CRG 36-mer mappability over the starts of every read overlapping the
     interval >= 0.5
  B  isolated (primary): A, plus no predicted Quadruplex / Z-DNA / Cruciform / SIDD motif of the
     authors within +-100 bp, plus no RefSeq TSS within +-2 kb

Controls: same chromosome, same length, same A/T count, purine asymmetry max(A+G, C+T) within
+-1 ("C"; "G" = A/T count only), no N, no HSeeker locus or predicted H-DNA within +-100 bp, the
locus's tier. Neighbourhoods: local (+-50 kb, primary) and chromosome.
Checks: placebo (first control as pseudo-locus vs the rest) and a positive control (the authors'
ssDNA+ H-DNA that pass tier B).

Data root: $KOUZINE_DIR (see run_kouzine_comparison.py).

Usage:  python3 analysis/scripts/kouzine_hdna_isolated.py --genome mm9|hg19 --config cli [--scan-only]
        [--chroms chr19]  (smoke test: the data root's genome.fa must hold only these chromosomes)
"""

from __future__ import annotations

import argparse
import collections
import csv
import gzip
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy import stats

sys.path.insert(0, str(Path(__file__).resolve().parent))
import run_kouzine_comparison as rk  # noqa: E402
import run_kouzine_published as rp  # noqa: E402

OTHER_CLASSES = ("Quadruplex", "Z-DNA", "Cruciform", "SIDD")
PAD = 100                 # authors' tag extension
TSS_PAD = 2000
LOCAL = 50_000
N_SAMPLE = 20_000
N_CONTROLS = 10
CANDIDATES = 30
DRAWS = 20_000
ROUNDS = 5
MAP_MIN = 0.5
ARM_CLASSES = ((8, 9), (10, 11), (12, 14), (15, 19), (20, 10**9))


# ---------------------------------------------------------------- intervals

def interval_index(intervals, pad: int) -> dict:
    """Per chromosome: sorted padded starts and the running maximum of padded ends."""
    by = collections.defaultdict(list)
    for c, s, e in intervals:
        by[c].append((s - pad, e + pad))
    out = {}
    for c, v in by.items():
        v.sort()
        out[c] = (np.array([x[0] for x in v], dtype=np.int64),
                  np.maximum.accumulate(np.array([x[1] for x in v], dtype=np.int64)))
    return out


def blocked(ix: dict, c: str, p: np.ndarray, L) -> np.ndarray:
    """True where [p, p + L) overlaps a padded interval of *ix*."""
    if c not in ix:
        return np.zeros(len(p), dtype=bool)
    st, rmax = ix[c]
    i = np.searchsorted(st, p + L) - 1
    return (i >= 0) & (rmax[np.maximum(i, 0)] > p)


def tss_sites(genome: str, primary: set[str]):
    rows = set()
    with gzip.open(rp.AUTH / "refgene" / f"{genome}_refGene.txt.gz", "rt") as fh:
        for line in fh:
            x = line.split("\t")
            if x[2] in primary:
                t = int(x[4]) if x[3] == "+" else int(x[5]) - 1
                rows.add((x[2], t, t + 1))
    return sorted(rows)


def load_loci(path: Path) -> list[dict]:
    with open(path) as fh:
        return [dict(chrom=r["chrom"], start=int(r["start0"]), end=int(r["end"]), arm=int(r["arm_length"]),
                     score=float(r["total_score"]) if r["total_score"] not in ("", "None") else math.nan)
                for r in csv.DictReader(fh, delimiter="\t")]


# ---------------------------------------------------------------- controls

def sequence_ok(p, L, at_n, asym, pref, use_asym):
    at_pref, pu_pref, n_pref = pref
    ok = (n_pref[p + L] - n_pref[p] == 0) & (at_pref[p + L] - at_pref[p] == at_n)
    if use_asym:
        pu = pu_pref[p + L] - pu_pref[p]
        ok &= np.abs(np.maximum(pu, L - pu) - asym) <= 1
    return ok


def draw_candidates(targets, genome, blocks, rng) -> dict:
    """targets: {key: (c, s, e, isolated)} -> {(key, scope, set): candidate starts}."""
    by_chrom = collections.defaultdict(list)
    for key, t in targets.items():
        by_chrom[t[0]].append((key, t))
    cand = {}
    for c in sorted(by_chrom):
        arr = np.frombuffer(genome[c], dtype=np.uint8)
        pref = tuple(np.concatenate(([0], np.cumsum(x, dtype=np.int32))) for x in
                     ((arr == 65) | (arr == 84), (arr == 65) | (arr == 71), arr == 78))
        n = len(arr)
        for key, (_, s, e, isolated) in by_chrom[c]:
            L = e - s
            seq = arr[s:e]
            at_n = int(((seq == 65) | (seq == 84)).sum())
            pu = int(((seq == 65) | (seq == 71)).sum())
            asym = max(pu, L - pu)
            use = blocks["isolated"] if isolated else blocks["mappable"]

            def keep(p, use_asym):
                ok = sequence_ok(p, L, at_n, asym, pref, use_asym)
                for ix in use:
                    ok &= ~blocked(ix, c, p, L)
                return p[ok]

            for use_asym, cset in ((True, "C"), (False, "G")):
                p = np.arange(max(0, s - LOCAL), min(n - L, e + LOCAL), dtype=np.int64)
                cand[(key, "local", cset)] = rng.permutation(keep(p, use_asym))[:CANDIDATES]
                found = []
                for _ in range(ROUNDS):
                    found.extend(keep(rng.integers(0, n - L, DRAWS), use_asym).tolist())
                    if len(found) >= CANDIDATES:
                        break
                cand[(key, "chromosome", cset)] = np.array(found[:CANDIDATES], dtype=np.int64)
        del pref, arr
    return cand


# ---------------------------------------------------------------- statistics

def compare(lv: np.ndarray, cv: list[np.ndarray], rng) -> dict:
    if len(lv) == 0:
        return dict(n=0)
    est, lo, hi = rk.boot_ratio(lv, cv, rng)
    cm = np.array([v.mean() for v in cv])
    # Probability that the locus beats one of its own controls (ties count 1/2). Exactly 0.5 in
    # expectation when locus and controls are exchangeable, also for zero-inflated signal, where
    # locus - mean(controls) is not (deviation 5, correction).
    ps = np.array([((v < x).sum() + 0.5 * (v == x).sum()) / len(v) for x, v in zip(lv, cv)])
    boots = [ps[ix].mean() for ix in (rng.integers(0, len(ps), len(ps)) for _ in range(rk.N_BOOT))]
    tp = float(stats.ttest_1samp(ps, 0.5, alternative="greater").pvalue) if np.ptp(ps) > 0 else None
    return dict(n=len(lv), controls_per_locus=round(float(np.mean([len(v) for v in cv])), 2),
                locus_mean=float(lv.mean()), control_mean=float(cm.mean()),
                ratio_of_means=round(float(est), 4), ratio_ci95=[round(lo, 4), round(hi, 4)],
                prob_superiority=round(float(ps.mean()), 4),
                prob_superiority_ci95=[round(float(np.percentile(boots, 2.5)), 4),
                                       round(float(np.percentile(boots, 97.5)), 4)],
                prob_superiority_p_greater=tp,
                share_locus_signal_nonzero=round(float((lv > 0).mean()), 4),
                share_control_signal_nonzero=round(float(np.mean([(v > 0).mean() for v in cv])), 4))


def placebo(cv: list[np.ndarray], rng) -> dict:
    keep = [v for v in cv if len(v) >= 3]
    return compare(np.array([v[0] for v in keep]), [v[1:] for v in keep], rng)


# ---------------------------------------------------------------- main

def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", choices=rp.SETS, required=True)
    ap.add_argument("--config", choices=("cli", "cli_minrep8"), default="cli")
    ap.add_argument("--scan-only", action="store_true")
    ap.add_argument("--chroms", nargs="+", help="restrict to these chromosomes (smoke tests only)")
    args = ap.parse_args()
    g, cfg = args.genome, rp.SETS[args.genome]
    work = rk.KOUZINE / "work" / g / args.config
    work.mkdir(parents=True, exist_ok=True)
    rk.scan(g, args.config, work)
    if args.scan_only:
        return

    rng = np.random.default_rng(rk.SEED)
    primary = set(args.chroms or rk.GENOMES[g]["primary"])
    loci = load_loci(work / "loci.tsv")
    rk.log(f"{g}/{args.config}: {len(loci)} loci")
    predicted = [m for m in rp.read_bed(rp.AUTH / "nonB_DNA_predicted" / cfg["species"] / "H-DNA.bed.gz")
                 if m[0] in primary]
    called = sorted({m for m in rp.read_bed(rp.AUTH / "nonB_DNA_ssDNA_enriched" / cfg["species"]
                                            / f"{cfg['called']}.bed.gz") if m[0] in primary})
    other = [m for k in OTHER_CLASSES
             for m in rp.read_bed(rp.AUTH / "nonB_DNA_predicted" / cfg["species"] / f"{k}.bed.gz")
             if m[0] in primary]
    tss = tss_sites(g, primary)
    ix_other, ix_tss = interval_index(other, PAD), interval_index(tss, TSS_PAD)
    ix_struct = interval_index([(l["chrom"], l["start"], l["end"]) for l in loci] + predicted, PAD)
    blocks = dict(mappable=[ix_struct], isolated=[ix_struct, ix_other, ix_tss])

    def isolated(c, s, e) -> bool:
        p = np.array([s], dtype=np.int64)
        return not (blocked(ix_other, c, p, e - s)[0] or blocked(ix_tss, c, p, e - s)[0])

    # --- tiers and samples
    iv = [(l["chrom"], l["start"], l["end"]) for l in loci]
    mp = rp.mappability(g, iv, tag=f"p4_{args.config}_loci")
    mapp = [i for i, t in enumerate(iv) if mp[t] >= MAP_MIN]
    iso = [i for i in mapp if isolated(*iv[i])]
    mp_pos = rp.mappability(g, called, tag="p4_called")
    pos = [m for m in called if mp_pos[m] >= MAP_MIN and isolated(*m)]
    funnel = dict(loci=len(loci), tier_A_mappable=len(mapp), tier_B_isolated=len(iso),
                  ssDNA_plus_hdna=len(called), ssDNA_plus_hdna_tier_B=len(pos))
    rk.log(f"funnel {funnel}")
    samples = {}
    for tier, pool in (("A", mapp), ("B", iso)):
        ix = sorted(rng.choice(len(pool), size=min(N_SAMPLE, len(pool)), replace=False))
        samples[tier] = [pool[i] for i in ix]
    pos_sample = [pos[i] for i in sorted(rng.choice(len(pos), size=min(N_SAMPLE, len(pos)), replace=False))]

    targets = {}
    for tier, v in samples.items():
        for i in v:
            targets[(tier, i)] = (*iv[i], tier == "B")
    for k, m in enumerate(pos_sample):
        targets[("pos", k)] = (*m, True)

    genome = dict(rk.read_fasta(rk.KOUZINE / "work" / g / f"{g}.primary.fa", primary))
    cand = draw_candidates(targets, genome, blocks, rng)
    lengths = {c: len(s) for c, s in genome.items()}
    del genome
    windows = sorted({(targets[key][0], int(p), int(p) + targets[key][2] - targets[key][1])
                      for (key, _, _), ps in cand.items() for p in ps})
    rk.log(f"{len(windows)} candidate control windows")
    cmp = rp.mappability(g, windows, tag=f"p4_{args.config}_controls")
    controls = {}
    for (key, scope, cset), ps in cand.items():
        c, s, e, _ = targets[key]
        w = [(c, int(p), int(p) + e - s) for p in ps]
        controls[(key, scope, cset)] = [x for x in w if cmp[x] >= MAP_MIN][:N_CONTROLS]

    # --- signal and statistics
    res = dict(genome=g, hseeker_config=args.config, hseeker_params=rk.CONFIGS[args.config],
               design="RTR/evidence/X8_kouzine_analysis_design_2026-10-09.md, deviation 5",
               parameters=dict(other_classes=OTHER_CLASSES, pad=PAD, tss_pad=TSS_PAD, local=LOCAL,
                               n_sample=N_SAMPLE, n_controls=N_CONTROLS, candidates=CANDIDATES,
                               mappability_min=MAP_MIN, seed=rk.SEED),
               funnel=funnel, signal={})
    for label, wig in cfg["wig"].items():
        cs = rp.wig_track(rp.AUTH / "ssDNA_wiggle" / cfg["species"] / f"{wig}.wig.gz", lengths)
        total_m = sum(v[-1] for v in cs.values()) / 1e6
        sig = lambda t: rp.wig_signal(cs, t[0], t[1], t[2], total_m)  # noqa: E731
        out = {}
        groups = [("tier_A", [("A", i) for i in samples["A"]]), ("tier_B", [("B", i) for i in samples["B"]]),
                  ("positive_control_tier_B", [("pos", k) for k in range(len(pos_sample))])]
        for name, keys in groups:
            for scope in ("local", "chromosome"):
                for cset in ("C", "G"):
                    use = [k for k in keys if controls[(k, scope, cset)]]
                    lv = np.array([sig(targets[k]) for k in use])
                    cv = [np.array([sig(w) for w in controls[(k, scope, cset)]]) for k in use]
                    r = compare(lv, cv, rng)
                    r["targets_without_controls"] = len(keys) - len(use)
                    r["placebo"] = placebo(cv, rng)
                    if name != "positive_control_tier_B":
                        arm = np.array([loci[k[1]]["arm"] for k in use])
                        score = np.array([loci[k[1]]["score"] for k in use])
                        d = np.array([((v < x).sum() + 0.5 * (v == x).sum()) / len(v) for x, v in zip(lv, cv)])
                        r["by_arm_length"] = {
                            f"{lo}-{hi}" if hi < 10**9 else f">={lo}":
                                compare(lv[m], [cv[j] for j in np.flatnonzero(m)], rng)
                            for lo, hi in ARM_CLASSES for m in [(arm >= lo) & (arm <= hi)]}
                        rho = stats.spearmanr(arm, d)
                        r["arm_length_vs_prob_superiority_spearman"] = dict(rho=float(rho.statistic),
                                                                      p=float(rho.pvalue))
                        ok = ~np.isnan(score)
                        if ok.sum() >= 4:
                            q = np.quantile(score[ok], [0.25, 0.5, 0.75])
                            qi = np.digitize(score, q)
                            r["by_score_quartile"] = {
                                f"Q{j + 1}": compare(lv[m], [cv[i] for i in np.flatnonzero(m)], rng)
                                for j in range(4) for m in [ok & (qi == j)]}
                    out[f"{name} | {scope} | {cset}"] = r
                    rk.log(f"{label} {name} {scope} {cset}: ratio {r.get('ratio_of_means')} {r.get('ratio_ci95')} "
                           f"P(sup) {r.get('prob_superiority')} | placebo ratio {r['placebo'].get('ratio_of_means')} "
                           f"P(sup) {r['placebo'].get('prob_superiority')}")
        res["signal"][label] = out
        del cs
    if "activated" in res["signal"] and "resting" in res["signal"]:
        res["activated_over_resting"] = {
            k: res["signal"]["activated"][k]["ratio_of_means"] / res["signal"]["resting"][k]["ratio_of_means"]
            for k in res["signal"]["activated"] if res["signal"]["resting"][k].get("ratio_of_means")}
    res["chromosomes"] = sorted(primary)
    path = (rk.ROOT / "analysis" / "results" / "kouzine_v1" / "published" / f"{g}.{args.config}.p4.json"
            if not args.chroms else work / "p4_smoke.json")
    path.write_text(json.dumps(res, indent=1))
    rk.log(f"wrote {path}")


if __name__ == "__main__":
    main()
