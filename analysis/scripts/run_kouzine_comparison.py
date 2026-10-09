#!/usr/bin/env python3
"""Compare HSeeker predictions with Kouzine et al. 2017 ssDNA-seq (X8, DEC-14).

Implements the pre-specified design in
RTR/evidence/X8_kouzine_analysis_design_2026-10-09.md:

  1. scan the primary chromosomes with HSeeker (named configuration)
  2. sample 10 composition-matched control intervals per locus
     (same chromosome and length, GC and purine asymmetry within +-0.02,
     no N, no HSeeker locus within 1 kb)
  3. count MAPQ>=10 reads per interval (samtools bedcov -c) for each group
  4. enrichment (loci vs controls), subset test against ssDNA-enriched
     regions (200 bp windows, Poisson vs background, BH FDR 1 %), score
     correlation, activated vs resting

Every step caches its output in --outdir, so the script can be re-run.

Usage:
  python3 analysis/scripts/run_kouzine_comparison.py --genome mm10 --config cli
  python3 analysis/scripts/run_kouzine_comparison.py --genome hg38 --config cli
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
from scipy import stats

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))
import hseeker  # noqa: E402

KOUZINE = Path("/workspaces/.hseeker-pr/kouzine")
SAMTOOLS = "/workspaces/.hseeker-pr/env/bin/samtools"
SEED = 20261009
N_CONTROLS = 10
MAX_DRAWS = 1000
TOL = 0.02
BLOCK_PAD = 1000
WINDOW = 200
FDR = 0.01
N_BOOT = 1000

GENOMES = {
    "mm10": dict(primary=[f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"],
                 groups=dict(activated="mouse_activated_ssDNA", resting="mouse_resting_ssDNA",
                             input="mouse_activated_input")),
    "hg38": dict(primary=[f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"],
                 groups=dict(raji="human_raji_ssDNA")),
    # builds of the authors' published calls (deviation 4); used for scan + controls only
    "mm9": dict(primary=[f"chr{i}" for i in range(1, 20)] + ["chrX", "chrY"], groups={}),
    "hg19": dict(primary=[f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"], groups={}),
}
CONFIGS = {
    # shipped command-line defaults
    "cli": dict(minrep=10, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10,
                at_threshold=0.8, filter_homopolymers=True, overlap_strategy="greedy"),
    # benchmark configuration (analysis scripts)
    "benchmark": dict(minrep=8, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10,
                      overlap_strategy="greedy"),
}


def log(msg: str) -> None:
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# ---------------------------------------------------------------- genome

def read_fasta(path: Path, keep: set[str]):
    """Yield (name, uppercase bytes) for records in *keep*."""
    name, parts = None, []
    with open(path, "rb") as fh:
        for line in fh:
            if line.startswith(b">"):
                if name in keep:
                    yield name, b"".join(parts).upper()
                name, parts = line[1:].split()[0].decode(), []
            elif name in keep:
                parts.append(line.strip())
    if name in keep:
        yield name, b"".join(parts).upper()


def primary_fasta(genome: str, out: Path) -> Path:
    path = out / f"{genome}.primary.fa"
    if not path.exists():
        src = KOUZINE / "index" / genome / "genome.fa"
        with open(path, "wb") as fo:
            for name, seq in read_fasta(src, set(GENOMES[genome]["primary"])):
                fo.write(b">" + name.encode() + b"\n")
                for i in range(0, len(seq), 80):
                    fo.write(seq[i:i + 80] + b"\n")
        log(f"wrote {path}")
    return path


# ---------------------------------------------------------------- 1. scan

def scan(genome: str, config: str, out: Path) -> Path:
    path = out / "loci.tsv"
    if path.exists():
        return path
    fa = primary_fasta(genome, out.parent)
    t = time.perf_counter()
    hits = hseeker.scan_fasta_parallel(str(fa), workers=16, score=True, purity_rmq=True,
                                       **CONFIGS[config])
    with open(path, "w") as fo:
        fo.write("chrom\tstart0\tend\tarm_length\tspacer_length\ttotal_score\n")
        for h in hits:
            fo.write(f"{h['seq_id']}\t{h['start'] - 1}\t{h['end']}\t{h['arm_length']}\t"
                     f"{h['spacer_length']}\t{h['total_score']}\n")
    log(f"scan {genome}/{config}: {len(hits)} loci in {time.perf_counter() - t:.0f} s")
    return path


def load_loci(path: Path) -> list[dict]:
    with open(path) as fh:
        return [dict(chrom=r["chrom"], start=int(r["start0"]), end=int(r["end"]),
                     score=float(r["total_score"]) if r["total_score"] not in ("", "None") else math.nan)
                for r in csv.DictReader(fh, delimiter="\t")]


# ---------------------------------------------------------------- 2. controls

EXACT_MAX = 60            # loci up to this length use windows of their exact length
LENGTH_CLASSES = [(60, 100), (100, 200), (200, 10**9)]   # longer loci: class pools
POOL_STRIDE = 10
POOL_CAP = 200            # windows kept per cell
WIDEN_STEPS = (0, 1, 2)   # asymmetry widened only if a cell is empty (counts, or 0.02 bins)


def _bin(x):
    return np.minimum((np.asarray(x) / TOL).astype(int), int(1 / TOL))


def sample_controls(genome: str, loci: list[dict], out: Path) -> Path:
    """Composition-matched controls from a genome-wide window pool (design deviations 1-3).

    Windows are enumerated at a 10 bp stride on all primary chromosomes, excluding N and anything
    within 1 kb of an HSeeker locus. A locus of length <= 60 bp is matched to windows of its exact
    length with the same GC count and purine-asymmetry count max(#A+G, #C+T); longer loci are
    matched within length classes (window length = class median) on 0.02 bins of GC and asymmetry.
    Two control sets per locus, 10 each, drawn with replacement:
      C: GC and asymmetry matched (primary); asymmetry widened by 1-2 units only if a cell is empty
      G: GC matched only (deviation 2), covering every locus
    """
    path = out / "intervals.bed"
    if path.exists():
        return path
    rng = np.random.default_rng(SEED)
    fa = out.parent / f"{genome}.primary.fa"
    lengths = np.array([l["end"] - l["start"] for l in loci])
    exact_lengths = sorted({int(x) for x in lengths if x <= EXACT_MAX})
    long_cls = np.searchsorted([b for _, b in LENGTH_CLASSES], lengths, side="right")
    cls_len = {int(k): int(np.median(lengths[(lengths > EXACT_MAX) & (long_cls == k)]))
               for k in np.unique(long_cls[lengths > EXACT_MAX])}
    # window plans: key -> window length
    plans = {("L", L): L for L in exact_lengths}
    plans.update({("K", k): W for k, W in cls_len.items()})

    def cell_of(kind, a_gc, a_pur, W):
        if kind == "L":
            return a_gc, np.maximum(a_pur, W - a_pur)
        p = a_pur / W
        return _bin(a_gc / W), _bin(np.maximum(p, 1 - p))

    by_chrom: dict[str, list[int]] = {}
    for i, l in enumerate(loci):
        by_chrom.setdefault(l["chrom"], []).append(i)
    pool: dict[tuple, list[tuple[str, int]]] = {}
    locus_key = {}
    for chrom, seq in read_fasta(fa, set(GENOMES[genome]["primary"])):
        a = np.frombuffer(seq, dtype=np.uint8)
        n = len(a)
        cum = lambda m: np.concatenate(([0], np.cumsum(m, dtype=np.int64)))  # noqa: E731
        gc = cum((a == ord("G")) | (a == ord("C")))
        pur = cum((a == ord("A")) | (a == ord("G")))
        bad = cum(~np.isin(a, np.frombuffer(b"ACGT", dtype=np.uint8)))
        blocked = np.zeros(n, dtype=bool)
        for i in by_chrom.get(chrom, []):
            blocked[max(0, loci[i]["start"] - BLOCK_PAD):min(n, loci[i]["end"] + BLOCK_PAD)] = True
        blk = cum(blocked)
        del a, blocked
        for i in by_chrom.get(chrom, []):
            s, e = loci[i]["start"], loci[i]["end"]
            L = e - s
            if bad[e] - bad[s]:
                continue  # loci containing N are not matched
            if L <= EXACT_MAX:
                kind, W, plan = "L", L, ("L", L)
            else:
                kind, W, plan = "K", L, ("K", int(long_cls[i]))
            g, asym = cell_of(kind, int(gc[e] - gc[s]), int(pur[e] - pur[s]), W)
            locus_key[i] = (plan, int(g), int(asym))
        for plan, W in plans.items():
            st = np.arange(0, n - W, POOL_STRIDE)
            en = st + W
            ok = (bad[en] - bad[st] == 0) & (blk[en] - blk[st] == 0)
            st, en = st[ok], en[ok]
            g, asym = cell_of(plan[0], gc[en] - gc[st], pur[en] - pur[st], W)
            cells = g.astype(np.int64) * 10000 + asym
            keyed = np.lexsort((rng.random(len(cells)), cells))
            cells, st = cells[keyed], st[keyed]
            first = np.r_[0, np.flatnonzero(np.diff(cells)) + 1]
            for f, nxt in zip(first, np.r_[first[1:], len(cells)]):
                c = int(cells[f])
                pool.setdefault((plan, c // 10000, c % 10000), []).extend(
                    (chrom, int(x)) for x in st[f:min(nxt, f + POOL_CAP)])
        log(f"pool {chrom}: done")
    gc_pool: dict[tuple, list[tuple[str, int]]] = {}
    for (plan, g, _), v in pool.items():
        gc_pool.setdefault((plan, g), []).extend(v)
    rows, widened, unmatched, gc_unmatched, n_bad = [], {1: 0, 2: 0}, 0, 0, 0
    for i in range(len(loci)):
        if i not in locus_key:
            n_bad += 1
            continue
        plan, g, asym = locus_key[i]
        W = plans[plan]
        Lc = loci[i]
        rows.append((Lc["chrom"], Lc["start"], Lc["end"], f"L{i}"))
        gcand = gc_pool.get((plan, g), [])
        if gcand:
            for j, ix in enumerate(rng.integers(0, len(gcand), N_CONTROLS)):
                c, s0 = gcand[ix]
                rows.append((c, s0, s0 + W, f"G{i}_{j}"))
        else:
            gc_unmatched += 1
        cand = None
        for w in WIDEN_STEPS:
            cand = [x for d in range(-w, w + 1) for x in pool.get((plan, g, asym + d), [])]
            if cand:
                if w:
                    widened[w] += 1
                break
        if not cand:
            unmatched += 1
            continue
        for j, ix in enumerate(rng.integers(0, len(cand), N_CONTROLS)):
            c, s0 = cand[ix]
            rows.append((c, s0, s0 + W, f"C{i}_{j}"))
    with open(path, "w") as fo:
        for c, s, e, name in rows:
            fo.write(f"{c}\t{s}\t{e}\t{name}\n")
    meta = dict(method="genome-wide pool; exact length + GC/purine counts for loci <= 60 bp, "
                       "length classes + 0.02 bins above",
                exact_max=EXACT_MAX, long_classes={str(LENGTH_CLASSES[k]): v for k, v in cls_len.items()},
                loci=len(loci), loci_with_N_excluded=n_bad,
                asymmetry_matched=len(loci) - n_bad - unmatched, asymmetry_unmatched=unmatched,
                widened_by_1=widened[1], widened_by_2=widened[2],
                gc_matched=len(loci) - n_bad - gc_unmatched, gc_unmatched=gc_unmatched,
                pool_cells=len(pool), pool_windows=sum(len(v) for v in pool.values()))
    json.dump(meta, open(out / "controls_meta.json", "w"), indent=1)
    log(f"controls: {meta}")
    return path


# ---------------------------------------------------------------- 3. counts

def bam_for(genome: str, group: str) -> Path:
    return KOUZINE / "bam" / f"{group}.{genome}.q10.bam"


def mapped_millions(bam: Path, primary: list[str]) -> float:
    out = subprocess.run([SAMTOOLS, "idxstats", str(bam)], capture_output=True, text=True, check=True).stdout
    return sum(int(l.split("\t")[2]) for l in out.splitlines() if l.split("\t")[0] in primary) / 1e6


def count(genome: str, group: str, bed: Path, out: Path) -> Path:
    path = out / f"counts.{group}.tsv.gz"
    if not path.exists():
        res = subprocess.run([SAMTOOLS, "bedcov", "-c", "-Q", "10", str(bed), str(bam_for(genome, group))],
                             capture_output=True, text=True, check=True).stdout
        with gzip.open(path, "wt") as fo:
            fo.write(res)
        log(f"counted {group}")
    return path


def load_counts(path: Path) -> dict[str, int]:
    with gzip.open(path, "rt") as fh:
        return {l.split("\t")[3]: int(l.rstrip("\n").split("\t")[-1]) for l in fh}


# ---------------------------------------------------------------- 4. enriched regions

def read_starts(bam: Path, chrom: str, n_windows: int) -> np.ndarray:
    """Reads per 200 bp window (by leftmost aligned position), streamed."""
    counts = np.zeros(n_windows + 1, dtype=np.int64)
    proc = subprocess.Popen([SAMTOOLS, "view", "-q", "10", str(bam), chrom], stdout=subprocess.PIPE)
    buf = []
    for line in proc.stdout:
        buf.append(int(line.split(b"\t", 4)[3]))
        if len(buf) >= 1_000_000:
            counts += np.bincount((np.array(buf) - 1) // WINDOW, minlength=n_windows + 1)[:n_windows + 1]
            buf = []
    if buf:
        counts += np.bincount((np.array(buf) - 1) // WINDOW, minlength=n_windows + 1)[:n_windows + 1]
    if proc.wait() != 0:
        raise RuntimeError(f"samtools view failed for {bam} {chrom}")
    return counts[:n_windows]


def enriched_regions(genome: str, group: str, control_group: str | None, out: Path) -> Path:
    """200 bp windows, Poisson vs max(input window, input 10 kb mean, genome mean), BH FDR 1 %."""
    path = out.parent / f"enriched.{group}.bed"
    if path.exists():
        return path
    fa = out.parent / f"{genome}.primary.fa"
    lengths = {c: len(s) for c, s in read_fasta(fa, set(GENOMES[genome]["primary"]))}
    bam = bam_for(genome, group)
    per = {c: read_starts(bam, c, math.ceil(L / WINDOW)) for c, L in lengths.items()}
    total = sum(v.sum() for v in per.values())
    nwin = sum(len(v) for v in per.values())
    if control_group:
        cbam = bam_for(genome, control_group)
        ctl = {c: read_starts(cbam, c, len(per[c])) for c in per}
        scale = total / sum(v.sum() for v in ctl.values())
        gmean = sum(v.sum() for v in ctl.values()) / nwin
    pvals, idx = [], []
    for c, k in per.items():
        if control_group:
            x = ctl[c].astype(float)
            local = np.convolve(x, np.ones(50) / 50, mode="same")  # 10 kb = 50 windows
            lam = scale * np.maximum.reduce([x, local, np.full_like(x, gmean)])
        else:
            lam = np.full(len(k), total / nwin)
        p = stats.poisson.sf(k - 1, lam)
        pvals.append(p)
        idx.append(c)
    allp = np.concatenate(pvals)
    order = np.argsort(allp)
    q = np.empty_like(allp)
    ranked = allp[order] * len(allp) / np.arange(1, len(allp) + 1)
    q[order] = np.minimum.accumulate(ranked[::-1])[::-1]
    regions, off = [], 0
    for c, p in zip(idx, pvals):
        sig = q[off:off + len(p)] <= FDR
        off += len(p)
        w = np.flatnonzero(sig)
        if len(w):
            breaks = np.flatnonzero(np.diff(w) > 1)
            for a, b in zip(np.r_[0, breaks + 1], np.r_[breaks, len(w) - 1]):
                regions.append((c, int(w[a]) * WINDOW, (int(w[b]) + 1) * WINDOW))
    with open(path, "w") as fo:
        for c, s, e in regions:
            fo.write(f"{c}\t{s}\t{e}\n")
    log(f"enriched regions {group}: {len(regions)} (windows tested {nwin})")
    return path


def overlaps(intervals: dict[str, tuple[str, int, int]], regions: Path) -> dict[str, bool]:
    reg: dict[str, list[tuple[int, int]]] = {}
    for line in open(regions):
        c, s, e = line.split()
        reg.setdefault(c, []).append((int(s), int(e)))
    starts = {c: np.array([s for s, _ in v]) for c, v in reg.items()}
    ends = {c: np.array([e for _, e in v]) for c, v in reg.items()}
    res = {}
    for name, (c, s, e) in intervals.items():
        if c not in reg:
            res[name] = False
            continue
        i = np.searchsorted(starts[c], e) - 1
        res[name] = bool(i >= 0 and ends[c][i] > s)
    return res


# ---------------------------------------------------------------- 5. statistics

def boot_ratio(locus_vals: np.ndarray, ctrl_vals: list[np.ndarray], rng) -> tuple[float, float, float]:
    """Ratio of means loci/controls; CI by resampling loci together with their controls."""
    cmeans = np.array([v.mean() for v in ctrl_vals])
    est = locus_vals.mean() / cmeans.mean()
    n = len(locus_vals)
    boots = [locus_vals[ix].mean() / cmeans[ix].mean() for ix in (rng.integers(0, n, n) for _ in range(N_BOOT))]
    return est, float(np.percentile(boots, 2.5)), float(np.percentile(boots, 97.5))


def odds_ratio(a: int, b: int, c: int, d: int) -> tuple[float, float, float]:
    a, b, c, d = (x + 0.5 for x in (a, b, c, d))
    o = (a * d) / (b * c)
    se = math.sqrt(1 / a + 1 / b + 1 / c + 1 / d)
    return o, math.exp(math.log(o) - 1.96 * se), math.exp(math.log(o) + 1.96 * se)


def analyse(genome: str, config: str, out: Path) -> dict:
    primary = GENOMES[genome]["primary"]
    groups = GENOMES[genome]["groups"]
    loci = load_loci(out / "loci.tsv")
    intervals = {}
    for line in open(out / "intervals.bed"):
        c, s, e, name = line.split()
        intervals[name] = (c, int(s), int(e))
    counts = {g: load_counts(out / f"counts.{grp}.tsv.gz") for g, grp in groups.items()}
    mm = {g: mapped_millions(bam_for(genome, grp), primary) for g, grp in groups.items()}

    def rpkm(g, name):
        c, s, e = intervals[name]
        return counts[g][name] / ((e - s) / 1000) / mm[g]

    names = list(intervals)
    keep = set(names)
    if "input" in groups:  # drop unmappable intervals symmetrically
        keep = {n for n in names if counts["input"][n] > 0}
    signal_groups = [g for g in groups if g != "input"]
    rng = np.random.default_rng(SEED)
    result = dict(genome=genome, config=config, hseeker_params=CONFIGS[config], n_loci=len(loci),
                  mapped_millions=mm, intervals_kept=len(keep), intervals_total=len(names))
    for cset, prefix in (("asymmetry_matched", "C"), ("gc_matched", "G")):
        result[cset] = _analyse_set(prefix, keep, intervals, loci, counts, groups, signal_groups,
                                    rpkm, mm, genome, out, rng)
    if result["asymmetry_matched"].get("activated") and result["asymmetry_matched"].get("resting"):
        for cset in ("asymmetry_matched", "gc_matched"):
            r = result[cset]
            r["activated_over_resting_ratio"] = r["activated"]["ratio_of_means"] / r["resting"]["ratio_of_means"]
    return result


def _analyse_set(prefix, keep, intervals, loci, counts, groups, signal_groups, rpkm, mm, genome, out, rng):
    result = {}
    pairs = {}
    for name in sorted(keep):  # deterministic order, so bootstrap CIs are reproducible
        if name[0] not in ("L", prefix):
            continue
        i = int(name[1:].split("_")[0])
        pairs.setdefault(i, {"L": None, "C": []})
        if name.startswith("L"):
            pairs[i]["L"] = name
        else:
            pairs[i]["C"].append(name)
    usable = [i for i, p in pairs.items() if p["L"] and p["C"]]
    result["loci_used"] = len(usable)
    for g in signal_groups:
        lv = np.array([rpkm(g, pairs[i]["L"]) for i in usable])
        cv = [np.array([rpkm(g, n) for n in pairs[i]["C"]]) for i in usable]
        est, lo, hi = boot_ratio(lv, cv, rng)
        flat_c = np.concatenate(cv)
        r = dict(locus_mean_rpkm=float(lv.mean()), control_mean_rpkm=float(flat_c.mean()),
                 locus_median_rpkm=float(np.median(lv)), control_median_rpkm=float(np.median(flat_c)),
                 ratio_of_means=est, ratio_ci95=[lo, hi],
                 mannwhitney_p=float(stats.mannwhitneyu(lv, flat_c, alternative="greater").pvalue))
        if "input" in groups:
            enr = lambda n: math.log2((rpkm(g, n) + 0.5) / (counts["input"][n] / ((intervals[n][2] - intervals[n][1]) / 1000) / mm["input"] + 0.5))  # noqa: E731
            le = np.array([enr(pairs[i]["L"]) for i in usable])
            ce = np.array([np.mean([enr(n) for n in pairs[i]["C"]]) for i in usable])
            d = le - ce
            boots = [d[rng.integers(0, len(d), len(d))].mean() for _ in range(N_BOOT)]
            r["log2_enrichment_vs_input"] = dict(loci_mean=float(le.mean()), controls_mean=float(ce.mean()),
                                                 diff=float(d.mean()), diff_ci95=[float(np.percentile(boots, 2.5)),
                                                                                  float(np.percentile(boots, 97.5))])
            sc = np.array([loci[i]["score"] for i in usable])
            ok = ~np.isnan(sc)
            rho, p = stats.spearmanr(sc[ok], le[ok])
            qs = np.quantile(sc[ok], [0.25, 0.5, 0.75])
            quart = np.digitize(sc[ok], qs)
            r["score_vs_enrichment"] = dict(spearman_rho=float(rho), p=float(p),
                                            enrichment_by_score_quartile=[float(le[ok][quart == k].mean()) for k in range(4)])
        reg = enriched_regions(genome, groups[g], groups.get("input"), out)
        ov = overlaps({n: intervals[n] for n in keep}, reg)
        a = sum(ov[pairs[i]["L"]] for i in usable)
        cn = [n for i in usable for n in pairs[i]["C"]]
        c = sum(ov[n] for n in cn)
        o, olo, ohi = odds_ratio(a, len(usable) - a, c, len(cn) - c)
        # fraction of enriched regions that contain (overlap) an HSeeker locus
        regions = {f"R{k}": (c_, int(s_), int(e_)) for k, (c_, s_, e_) in
                   enumerate(line.split() for line in open(reg))}
        locus_bed = out / f"loci_used.{prefix}.{g}.bed"
        with open(locus_bed, "w") as fo:
            for c_, s_, e_ in sorted(intervals[pairs[i]["L"]] for i in usable):
                fo.write(f"{c_}\t{s_}\t{e_}\n")
        reg_hit = overlaps(regions, locus_bed)
        n_reg = len(regions)
        r["subset"] = dict(enriched_regions=n_reg,
                           enriched_regions_with_locus=sum(reg_hit.values()),
                           enriched_regions_with_locus_frac=(sum(reg_hit.values()) / n_reg) if n_reg else None,
                           loci_in_enriched=a, loci_frac=a / len(usable),
                           controls_in_enriched=c, controls_frac=c / len(cn), odds_ratio=o, or_ci95=[olo, ohi])
        result[g] = r
    return result


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome", choices=GENOMES, required=True)
    ap.add_argument("--config", choices=CONFIGS, required=True)
    ap.add_argument("--outdir", type=Path, default=ROOT / "analysis" / "results" / "kouzine_v1")
    args = ap.parse_args()
    # Large intermediates (genome FASTA, intervals, counts) stay outside the repository;
    # only summary.json and controls_meta.json are copied into --outdir.
    out = KOUZINE / "work" / args.genome / args.config
    out.mkdir(parents=True, exist_ok=True)
    results = args.outdir / args.genome / args.config
    results.mkdir(parents=True, exist_ok=True)
    loci = load_loci(scan(args.genome, args.config, out))
    bed = sample_controls(args.genome, loci, out)
    for grp in GENOMES[args.genome]["groups"].values():
        count(args.genome, grp, bed, out)
    res = analyse(args.genome, args.config, out)
    res["hseeker_version"] = hseeker.__version__
    res["git_commit"] = subprocess.run(["git", "-C", str(ROOT), "rev-parse", "HEAD"], capture_output=True, text=True).stdout.strip()
    json.dump(res, open(out / "summary.json", "w"), indent=1)
    for name in ("summary.json", "controls_meta.json"):
        (results / name).write_text((out / name).read_text())
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
