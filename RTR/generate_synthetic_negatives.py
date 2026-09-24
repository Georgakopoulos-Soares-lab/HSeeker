#!/usr/bin/env python3
"""
generate_synthetic_negatives.py — HSeeker benchmark v3 synthetic negative generator.

Generates synthetic non-forming controls to balance the curated experimental
benchmark to exactly 1:1. The number needed is computed from the input CSV
(forming minus non-forming); class counts follow the recipe in target_counts().
Current curated set: 66 forming / 6 non-forming -> 60 synthetic negatives,
132 records total (66:66).

Classes (reviewer R2.1 menu, adapted for DNA; perturbation-focused core):
  A. mirror_disrupted_mutant  (20) — heavy purine-preserving disruption of the
     longest mirror arm of real forming positives. Heavy disruption is required
     because Belotserkovskii 1990 showed single-mismatch mirrors still form H-DNA.
  B. dinucleotide_shuffle     (12) — exact dinucleotide-preserving shuffle
     (random Eulerian path through the fixed dinucleotide multigraph).
  C. random_gc_length_matched (12) — random sequence, same length and GC count as a
     randomly chosen forming positive.
  D. other_structure          (12) — G4 (4), Z-DNA (3), mixed B-DNA (5).
  E. homopolymer_AT           (4)  — poly-A/poly-T; poly-A20 is experimentally
     non-forming (Hanvey 1988, pRW1405 = HDNA0053, already in the experimental
     set); these lengths extend coverage and test the homopolymer filter.
     poly-G/poly-C are NOT used: poly-dG.dC genuinely forms H-DNA (Kohwi 1988;
     benchmark HDNA0054).
  F. perfect_mirror_balanced  (3)  — perfect mirror repeats with ~50% purine:
     mirror symmetry without purine richness must not be called H-DNA.

Mirror detection note: a sequence is H-DNA-competent only if it has a mirror repeat
with arm >= 8 nt at identity >= 0.80 AND a purine-rich (>= 85%) strand. Short
perfect mirrors (6-7 nt) occur by chance and MUST NOT mask longer imperfect ones,
so the strong-mirror search is restricted to arm >= 8 (the tool's own minimum).

Validation is STRUCTURAL ONLY (mirror score, purine fraction, dedup). Sequences are
never filtered by HSeeker's own predictions, to avoid circular benchmarking.
Classes E and F are exempt from the structural rejection check by design:
E rests on experimental evidence (Hanvey 1988) + the homopolymer filter;
F is a deliberate mirror-positive/purity-negative hard control.

Usage: python3 generate_synthetic_negatives.py
Reads:  hdna_benchmark_experimental_v3.csv (same directory)
Writes: hdna_benchmark_balanced_v3.csv (same directory)
Seed:   20260922 (fixed)
"""
import csv
import os
import random
from collections import Counter, defaultdict

SEED = 20260922
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)  # benchmark CSVs live at the repository root
IN_CSV = os.path.join(ROOT, 'hdna_benchmark_experimental_v3.csv')
OUT_CSV = os.path.join(ROOT, 'hdna_benchmark_balanced_v3.csv')

MIN_ARM = 6         # reporting floor for the best-mirror search
MAX_SPACER = 12
STRONG_ARM = 8      # mirror arm length at/above which...
STRONG_ID = 0.80    # ...identity >= this counts as an H-DNA-competent mirror
PURINE_RICH = 0.85  # strand purine fraction at/above which the purity criterion passes

PUR_SWAP = {'G': 'A', 'A': 'G', 'C': 'T', 'T': 'C'}  # purine/pyrimidine-preserving


# ---------------------------------------------------------------- mirror scoring
def find_best_mirror(seq, min_arm=MIN_ARM, max_spacer=MAX_SPACER):
    """Best same-strand mirror repeat: max (identity, arm_len), lexicographic.
    Mirror = left arm == reverse(right arm). Returns (identity, arm_len, i, j, spacer)."""
    best = (0.0, 0, 0, 0, 0)
    n = len(seq)
    for L in range(min_arm, n // 2 + 1):
        for i in range(0, n - 2 * L + 1):
            left = seq[i:i + L]
            for spacer in range(0, max_spacer + 1):
                j = i + L + spacer
                if j + L > n:
                    break
                right_rev = seq[j:j + L][::-1]
                m = sum(1 for a, b in zip(left, right_rev) if a == b)
                fid = m / L
                if (fid, L) > (best[0], best[1]):
                    best = (fid, L, i, j, spacer)
    return best


def best_strong_mirror(seq):
    """Best mirror restricted to arm >= STRONG_ARM — the H-DNA-competent scale.
    Short perfect mirrors cannot mask longer imperfect ones here."""
    return find_best_mirror(seq, min_arm=STRONG_ARM)


def has_strong_mirror(seq):
    fid, L, *_ = best_strong_mirror(seq)
    return L >= STRONG_ARM and fid >= STRONG_ID


def purine_fraction(seq):
    return sum(1 for c in seq if c in 'GA') / len(seq)


def revcomp(seq):
    return seq.translate(str.maketrans('ACGT', 'TGCA'))[::-1]


# ---------------------------------------------------------------- class A: mutants
def disrupt_mirror(seq, rng, max_iter=40):
    """Introduce purine-preserving substitutions in the left arm of the best
    arm>=8 mirror until none at identity >= 0.80 remains. Heavy disruption by
    design (see module docstring)."""
    s = list(seq)
    for _ in range(max_iter):
        fid, L, i, j, _ = best_strong_mirror(''.join(s))
        if not (L >= STRONG_ARM and fid >= STRONG_ID):
            return ''.join(s), fid, L
        left = s[i:i + L]
        right_rev = s[j:j + L][::-1]
        matches = [k for k in range(L) if left[k] == right_rev[k]]
        k_mut = max(3, (len(matches) + 1) // 2)
        for k in rng.sample(matches, min(k_mut, len(matches))):
            s[i + k] = PUR_SWAP[s[i + k]]
    fid, L, *_ = best_strong_mirror(''.join(s))
    return ''.join(s), fid, L


# ------------------------------------------------------- class B: dinuc shuffle
def dinuc_shuffle(seq, rng, max_tries=500):
    """Exact dinucleotide-preserving shuffle: sample a random Eulerian path through
    the FIXED dinucleotide multigraph of the original sequence (random walk over
    unused edges; restart if stuck at the terminal node). The edge multiset — hence
    the dinucleotide composition — is preserved exactly by construction; first/last
    nucleotides are preserved by the degree constraints. A walk can only get stuck
    at the terminal node (degree argument), in which case we retry."""
    edges = [(seq[i], seq[i + 1]) for i in range(len(seq) - 1)]
    for _ in range(max_tries):
        adj = defaultdict(list)
        for idx, (a, b) in enumerate(edges):
            adj[a].append(idx)
        used = [False] * len(edges)
        u = seq[0]
        path = [u]
        n_used = 0
        while n_used < len(edges):
            avail = [idx for idx in adj[u] if not used[idx]]
            if not avail:
                break  # stuck at terminal node with edges unused -> retry
            idx = rng.choice(avail)
            used[idx] = True
            n_used += 1
            u = edges[idx][1]
            path.append(u)
        if n_used == len(edges):
            return ''.join(path)
    raise RuntimeError(f'no Eulerian path found for {seq}')


# ------------------------------------------------------------- class C: random
def random_matched(positive, rng):
    L = len(positive)
    gc = sum(1 for c in positive if c in 'GC')
    g = rng.randint(0, gc)
    at = L - gc
    a = rng.randint(0, at)
    pool = ['G'] * g + ['C'] * (gc - g) + ['A'] * a + ['T'] * (at - a)
    rng.shuffle(pool)
    return ''.join(pool)


def random_bdna(rng, L):
    while True:
        s = ''.join(rng.choice('ACGT') for _ in range(L))
        if 0.4 <= purine_fraction(s) <= 0.6 and not has_strong_mirror(s):
            return s


# ---------------------------------------------------------------- class F: mirror
def perfect_mirror_balanced(rng, arm_len):
    while True:
        left = ''.join(rng.choice('ACGT') for _ in range(arm_len))
        if 0.4 <= purine_fraction(left) <= 0.6:
            break
    return left + 'TATA' + left[::-1]


# ---------------------------------------------------------------- fixed sequences
G4S = {
    'telomeric_(TTAGGG)4': 'TTAGGGTTAGGGTTAGGGTTAGGG',
    'G4_G3-tracts_TT-loops': 'GGGTTGGGTTGGGTTGGG',
    'c-kit21_G4': 'AGGGAGGGCGCTGGGAGGAGGG',
    'thrombin_aptamer_G4': 'GGTTGGTGTGGTTGG',
}
# NOTE: (GGGA)5 was considered and rejected by the built-in safety check: it is
# 100% purine AND contains a strong mirror repeat (periodic GGA mirrors), i.e. it
# is H-DNA-competent (cf. forming GGAA/GAAA repeats HDNA0073/0074) — not a valid
# negative. The G4s above all fail at least one H-DNA criterion (verified below).
ZDNAS = {'Z-DNA_(CG)8': 'CG' * 8, 'Z-DNA_(CG)10': 'CG' * 10, 'Z-DNA_(CG)12': 'CG' * 12}
# NOTE: poly-A20 is NOT generated — it already exists in the experimental set as
# HDNA0053 (pRW1405, Hanvey 1988, non-forming). These lengths add coverage without
# duplicating it.
HOMOPOLYMERS = {'poly-A15': 'A' * 15, 'poly-A30': 'A' * 30,
                'poly-T25': 'T' * 25, 'poly-T35': 'T' * 35}

# Class-count recipe. Fixed classes reflect the reviewer's explicit asks; the
# perturbation/random classes fill the remainder needed for a 1:1 balance, split
# ~50/25/25 (perturbation-focused core). For the current curated set
# (66 forming / 6 non-forming) this yields 60 synthetic negatives:
# 20 mutants, 10 shuffles, 11 random, 12 other-structure, 4 homopolymer, 3 mirror.
FIXED_COUNTS = {'other_structure': 12, 'homopolymer_AT': 4, 'perfect_mirror_balanced': 3}


def target_counts(n_needed):
    fixed = sum(FIXED_COUNTS.values())
    rem = n_needed - fixed
    assert rem >= 6, f'need at least {fixed + 6} synthetic negatives, got {n_needed}'
    n_mut = round(rem * 0.5)
    n_shuf = round(rem * 0.25)
    return {'mirror_disrupted_mutant': n_mut,
            'dinucleotide_shuffle': n_shuf,
            'random_gc_length_matched': rem - n_mut - n_shuf,
            **FIXED_COUNTS}


def main():
    rng = random.Random(SEED)
    rows = list(csv.DictReader(open(IN_CSV)))
    kept = [r for r in rows if r['curation_decision'] == 'kept']
    forming = [r for r in kept if r['label'] == 'forming']
    n_nonforming = len(kept) - len(forming)
    n_needed = len(forming) - n_nonforming
    COUNTS = target_counts(n_needed)
    print(f'experimental: {len(kept)} kept = {len(forming)} forming / '
          f'{n_nonforming} non-forming -> generating {n_needed} synthetic negatives')

    # dedup registry: exact + reverse-complement, experimental first
    taken = set()
    for r in kept:
        taken.add(r['sequence_5to3'])
        taken.add(revcomp(r['sequence_5to3']))

    # eligible sources for perturbation classes: forming, len>=20, >=3 distinct nt
    sources = [r for r in forming
               if len(r['sequence_5to3']) >= 20 and len(set(r['sequence_5to3'])) >= 3]
    rng.shuffle(sources)
    print(f'eligible perturbation sources: {len(sources)}')

    synth = []

    def register(seq, name, method, source, rationale, family):
        fid_any, L_any, *_ = find_best_mirror(seq)          # best mirror, any arm >= 6
        fid_s, L_s, *_ = best_strong_mirror(seq)            # best mirror, arm >= 8
        rec = {k: '' for k in rows[0].keys()}
        rec.update({
            'record_id': f'SYN{len(synth)+1:04d}',
            'sequence_name': name,
            'sequence_5to3': seq,
            'sequence_length': str(len(seq)),
            'label': 'non_forming',
            'label_scope': 'synthetic_design',
            'confidence_tier': 'synthetic',
            'study_id': 'synthetic_v3',
            'family_id': family,
            'flags': 'synthetic_negative',
            'relationships': f'derived_from:{source}' if source else '',
            'verification_status': 'not_applicable_synthetic',
            'resolved_origin': f'synthetic ({method}, seed={SEED})',
            'disposition': 'synthetic',
            'curation_decision': 'kept',
            'label_conflict': 'false',
            'record_type': 'synthetic',
            'generation_method': method,
            'source_record': source or '',
            'seed': str(SEED),
            'expected_label_rationale': rationale,
            'validation_mirror_identity': f'{fid_any:.3f}',
            'validation_mirror_arm': str(L_any),
            'validation_strong_mirror_identity': f'{fid_s:.3f}',
            'validation_strong_mirror_arm': str(L_s),
            'validation_purine_fraction': f'{purine_fraction(seq):.3f}',
        })
        synth.append(rec)
        taken.add(seq)
        taken.add(revcomp(seq))

    def novel(seq):
        return seq not in taken and revcomp(seq) not in taken

    # ---- class A: mirror-disrupted mutants (20)
    src_iter = iter(sources)
    made = 0
    while made < COUNTS['mirror_disrupted_mutant']:
        src = next(src_iter)
        mut, fid, L = disrupt_mirror(src['sequence_5to3'], rng)
        if L >= STRONG_ARM and fid >= STRONG_ID:
            continue  # disruption failed; next source
        if not novel(mut):
            continue  # e.g. a mutant that accidentally equals another pXY32 variant
        register(mut, f"{src['sequence_name']}_mirrordisrupted", 'mirror_disrupted_mutant',
                 src['record_id'],
                 f"Heavy purine-preserving disruption of the longest mirror arm of "
                 f"{src['record_id']} (residual arm>=8 mirror: arm={L}, identity={fid:.2f} "
                 f"< {STRONG_ID}). Heavy disruption required: single-mismatch mirrors still "
                 f"form H-DNA (Belotserkovskii 1990). Expected non-forming by design.",
                 f"SYN:{src['family_id']}")
        made += 1

    # ---- class B: dinucleotide shuffles (12)
    made = 0
    while made < COUNTS['dinucleotide_shuffle']:
        src = next(src_iter)
        for _ in range(25):
            sh = dinuc_shuffle(src['sequence_5to3'], rng)
            if sh != src['sequence_5to3'] and not has_strong_mirror(sh) and novel(sh):
                break
        else:
            continue
        assert Counter(sh[i:i+2] for i in range(len(sh)-1)) == \
               Counter(src['sequence_5to3'][i:i+2] for i in range(len(src['sequence_5to3'])-1))
        fid, L, *_ = best_strong_mirror(sh)
        register(sh, f"{src['sequence_name']}_dinucshuffle", 'dinucleotide_shuffle',
                 src['record_id'],
                 f"Dinucleotide-preserving shuffle of {src['record_id']} (exact dinucleotide "
                 f"composition, same length); arm>=8 mirror destroyed (residual arm={L}, "
                 f"identity={fid:.2f}). Expected non-forming by design.",
                 f"SYN:{src['family_id']}")
        made += 1

    # ---- class C: random GC/length-matched (12)
    made = 0
    while made < COUNTS['random_gc_length_matched']:
        src = forming[rng.randrange(len(forming))]
        s = random_matched(src['sequence_5to3'], rng)
        if has_strong_mirror(s) or not novel(s):
            continue
        fid, L, *_ = best_strong_mirror(s)
        register(s, f"random_matched_to_{src['record_id']}", 'random_gc_length_matched',
                 src['record_id'],
                 f"Random sequence matched to {src['record_id']} for length ({len(s)} nt) and GC "
                 f"count; no arm>=8 mirror at identity>={STRONG_ID} (best arm={L}, "
                 f"identity={fid:.2f}). Reviewer-requested random control. "
                 f"Expected non-forming by design.",
                 'SYN:random')
        made += 1

    # ---- class D: other structures (12 = 4 G4 + 3 Z-DNA + 5 B-DNA)
    for name, seq in G4S.items():
        fid, L, *_ = best_strong_mirror(seq)
        assert novel(seq), name
        assert not (L >= STRONG_ARM and fid >= STRONG_ID
                    and purine_fraction(seq) >= PURINE_RICH), name
        register(seq, name, 'other_structure:G4', '',
                 f"Known G-quadruplex-forming sequence; fails H-DNA criteria (arm>=8 mirror: "
                 f"arm={L}, identity={fid:.2f}; purine fraction={purine_fraction(seq):.2f}). "
                 f"Reviewer-requested other-structure control.", 'SYN:G4')
    for name, seq in ZDNAS.items():
        assert novel(seq), name
        register(seq, name, 'other_structure:Z-DNA', '',
                 f"Alternating CG (Z-DNA-forming); purine fraction "
                 f"{purine_fraction(seq):.2f} fails the H-DNA purity criterion. "
                 f"Reviewer-requested other-structure control.", 'SYN:ZDNA')
    made = 0
    while made < 5:
        s = random_bdna(rng, rng.choice([25, 30, 35, 40, 45]))
        if not novel(s):
            continue
        register(s, f'B-DNA_random_{made+1}', 'other_structure:B-DNA', '',
                 f"Random mixed-sequence B-DNA (~50% GC, purine fraction "
                 f"{purine_fraction(s):.2f}, no arm>=8 mirror). Reviewer-requested "
                 f"B-DNA control.", 'SYN:BDNA')
        made += 1

    # ---- class E: poly-A/poly-T homopolymers (4)
    for name, seq in HOMOPOLYMERS.items():
        assert novel(seq), name
        register(seq, name, 'homopolymer_AT', '',
                 "Homopolymer A/T. Poly-A20 experimentally non-forming (Hanvey 1988, pRW1405 = "
                 "HDNA0053, already in the experimental set); these lengths extend homopolymer "
                 "coverage and test the homopolymer filter. poly-G/poly-C deliberately excluded: "
                 "poly-dG.dC forms H-DNA (Kohwi 1988).", 'SYN:homopolymer')

    # ---- class F: perfect mirror, balanced composition (3)
    made = 0
    while made < COUNTS['perfect_mirror_balanced']:
        s = perfect_mirror_balanced(rng, rng.choice([10, 11, 12]))
        if not novel(s):
            continue
        assert has_strong_mirror(s)
        assert purine_fraction(s) < 0.6
        register(s, f'perfect_mirror_balancedGC_{made+1}', 'perfect_mirror_balanced', '',
                 f"Perfect mirror repeat (arm={(len(s)-4)//2} nt) but ~50% purine "
                 f"({purine_fraction(s):.2f}): mirror symmetry WITHOUT purine richness must not "
                 f"be called H-DNA. Reviewer-requested mirror control.", 'SYN:mirror')
        made += 1

    # ---------------------------------------------------------------- write
    assert len(synth) == n_needed, (len(synth), n_needed)
    extra_cols = ['record_type', 'generation_method', 'source_record', 'seed',
                  'expected_label_rationale', 'validation_mirror_identity',
                  'validation_mirror_arm', 'validation_strong_mirror_identity',
                  'validation_strong_mirror_arm', 'validation_purine_fraction']
    fieldnames = list(rows[0].keys()) + [c for c in extra_cols if c not in rows[0]]
    for r in kept:
        r['record_type'] = 'experimental'
    out = kept + synth
    try:
        with open(OUT_CSV, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(out)
    except PermissionError:
        # S3-backed mounts may refuse truncate-in-place but allow unlink+recreate.
        import shutil
        staging = '/workspace/' + os.path.basename(OUT_CSV)
        with open(staging, 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=fieldnames)
            w.writeheader()
            w.writerows(out)
        if os.path.exists(OUT_CSV):
            os.remove(OUT_CSV)
        shutil.copy(staging, OUT_CSV)
        print(f'(staged via {staging} — mount blocked direct overwrite)')

    # ---------------------------------------------------------------- summary
    print(f'\nwrote {OUT_CSV}: {len(out)} rows '
          f'({len(kept)} experimental + {len(synth)} synthetic)')
    print('label balance:', dict(Counter(r['label'] for r in out)))
    print('\nper-class summary (strong mirror = arm>=8):')
    for method in COUNTS:
        rs = [r for r in synth if r['generation_method'].startswith(method.split(':')[0])]
        sfids = [float(r['validation_strong_mirror_identity']) for r in rs]
        purs = [float(r['validation_purine_fraction']) for r in rs]
        print(f'  {method:28s} n={len(rs):2d}  strong_mirror_id max={max(sfids):.2f}  '
              f'purine range=[{min(purs):.2f},{max(purs):.2f}]')
    # A synthetic negative is INVALID only if it is H-DNA-competent:
    # strong mirror AND purine-rich. Classes E and F are exempt by design (see docstring).
    bad = [r['record_id'] for r in synth
           if r['generation_method'] not in ('perfect_mirror_balanced', 'homopolymer_AT')
           and float(r['validation_strong_mirror_identity']) >= STRONG_ID
           and int(r['validation_strong_mirror_arm']) >= STRONG_ARM
           and float(r['validation_purine_fraction']) >= PURINE_RICH]
    print('\nH-DNA-competent synthetic rows (must be []):', bad)


if __name__ == '__main__':
    main()
