# WS6-A4 / WS4-A4: chr1 runtime, memory and overlap-rule comparison on the merged code

**Status:** measured 2026-10-07; PROVISIONAL (development machine; see "Hardware").
**Code:** branch `RTR-merge-main` = origin/main `d237bca` (Nikol, PR #18) + DEC-9 score-informed
default + score-selection optimisation + PR #18 memory-regression fix (§3). Python 3.14.
**Default since 2026-10-08 (DEC-10):** greedy again, so the "greedy" rows below are the default-mode numbers and the "score" rows describe the optional score mode.
**Input:** hg38 chr1, UCSC `goldenPath/hg38/chromosomes/chr1.fa.gz` (248,956,422 bp).

## Hardware

AMD EPYC 7763, 16 vCPUs (lscpu), 62 GB RAM — the same CPU model and vCPU count as the
manuscript's benchmark node (`benchmarks/benchmark_comparison.md`: "AMD EPYC 7763, 32 GB RAM,
16 vCPUs"). Timings are wall-clock `time.perf_counter()` around `scan_fasta_parallel`; peak
memory is `resource.getrusage(RUSAGE_SELF).ru_maxrss` of a fresh process per run. The published
RSS used `/usr/bin/time -v`; both report the process's maximum resident set size.

## 1. Runtime and hits

Parameters as in the deposited chr1 run (`analysis/results/chr1/chr1_run_summary.json`): minrep 10,
maxrep 1000, maxspacer 10, purity 0.90, mismatch 0.10, `purity_rmq=True`, scoring on.
"Library" = no composition filters (as the deposited run); "CLI" = `at_threshold=0.8`,
`filter_homopolymers=True` (the command-line defaults).

| Mode | Workers | Wall (s) | Hits | Peak RSS (MB) |
|---|---|---|---|---|
| greedy, library | 16 | 26.3 / 26.5 | **96,729** | 1,979 (after fix; 2,704 before) |
| score, library (**new default**) | 16 | 58.6 / 60.1 / 61.0 | 98,223 | 1,969 (after fix) |
| greedy, library | 1 | 172.1 | 96,729 | 2,379 (before fix) |
| score, library | 1 | 205.2 | 98,223 | 2,672 (before fix) |
| greedy, CLI | 16 | 23.6 | 41,937 | 1,979 (before fix) |
| score, CLI | 16 | 40.4 | 43,150 | 1,986 (before fix) |

Published (manuscript benchmark, hseeker 0.1.0, greedy, library parameters): 16 workers 25.2 s,
1 worker 187.7 s, RSS 1,965 MB, 96,729 hits; triplex 57.7 s, Triplexator 67.5 s (16 cores).

**Reading.**
- Greedy reproduces the published result exactly (96,729 hits) and its runtime (26.3 s vs 25.2 s).
- Under the score-informed default (DEC-9) chr1 takes ~60 s at 16 workers — about the published
  triplex time (57.7 s). The "2.3× faster than triplex" statement holds only for greedy mode.
- The score-mode overhead (~33 s) is the same at 1 and 16 workers: it is serial Python.
  Breakdown at 16 workers: C scan ≈ 19–21 s in both modes (score mode returns 1,393,237 raw
  candidates instead of 96,729 hits); greedy then spends 1.3 s on overlap removal and 5.3 s scoring
  96,729 hits; score mode spends 41.4 s in lazy selection, of which 18.4 s is exact scoring and the
  rest is per-candidate work (bounds, heap, overlap checks) on 1.39 M candidates.
- A further pure-Python attempt (staged length → matched-pair → Kadane window bounds) returned
  identical output but no speedup (3.56 s vs 3.58 s on a 20 Mb chr1 slice) and was reverted.
  Closing the gap requires moving bound computation and the selection loop into the C extension
  (bounds computed in the parallel scan threads), estimated ~60 s → ~35 s; not done.

## 2. Score vs greedy: what changes in the output (chr1)

| | library | CLI |
|---|---|---|
| greedy hits | 96,729 | 41,937 |
| score hits | 98,223 (+1.5 %) | 43,150 (+2.9 %) |
| greedy hits with identical coordinates in score output | 80,932 (83.7 %) | 38,385 (91.5 %) |
| greedy hits overlapping some score hit | 96,729 (100 %) | 41,937 (100 %) |
| score-only loci (no greedy hit overlaps) | 1,076 | 936 |
| of these, ≥ 90 % one base | 32 | 5 |

Score selection loses no greedy locus; it reports a different (higher-scoring) hit in ~16 % (library)
of loci and adds ~1 k loci, mostly GA-rich mirror repeats (e.g. `AGAGAAAGAAAGAAAGAAGAAAGAGAAAG…`)
that greedy's longest-arm choice suppressed. Near-homopolymer "edge hits" (DEC-9 evidence §1) are
negligible on a real chromosome: the longest homopolymer runs on chr1 are G 27, C 21, A 69, T 57 bp.

## 3. Memory regression in PR #18, and its fix

Peak RSS, greedy, 16 workers, chr1, measured per commit:

| Code | RSS (MB) | Wall (s) |
|---|---|---|
| d1aca31 (pre-September) | 1,946 | 34.9 |
| 0a35aa6 | 1,954 | 33.4 |
| 063e3b5 | 1,974 | 29.5 |
| 5386b9f | 1,969 | 29.5 |
| 4784c4b | 1,961 | 29.9 |
| **d237bca (PR #18)** | **2,704** | 29.9 |
| merged + fix | 1,979 | 26.5 |

PR #18 scans each batch while the FASTA parser generator is suspended and still references the
current record's full sequence (previously the generator was exhausted before scanning). Fix:
flush a full batch when the next record arrives (or at end of file) instead of immediately —
batch contents are unchanged, so output is unchanged; Nikol's batching tests
(`tests/test_parallel_batching.py`, incl. golden outputs) pass. Not caused by the FASTA parser
(in-house 2,700 MB vs biopython 2,696 MB before the fix).

## Reproduce

```
python3 - <<'EOF'
import sys, time, resource; import hseeker
t = time.perf_counter()
h = hseeker.scan_fasta_parallel("chr1.fa", workers=16, minrep=10, maxrep=1000, maxspacer=10,
        purity=0.9, mismatch=0.1, score=True, purity_rmq=True, overlap_strategy="score")  # or "greedy"
print(time.perf_counter() - t, len(h), resource.getrusage(resource.RUSAGE_SELF).ru_maxrss // 1024)
EOF
```

## 4. Upstream baseline and PRs (2026-10-09)

Measured on Nikol's unmodified main (`d237bca`) and on the two PR branches, each tested against
its own source tree (`PYTHONPATH=src`), fresh process per run, chr1, 16 workers, minrep 10:

| Code | greedy | score | output sha256 (greedy / score) |
|---|---|---|---|
| main `d237bca` | 29.4 s, 2,708 MB | 193.0 s, 2,779 MB | `da869cb2…` / `8c5d8f12…` |
| PR #19 memory fix | 29.3 s, 1,978 MB | 193.3 s, 2,043 MB | identical |
| PR #20 score speed-up | 26.5 s, 2,702 MB | 60.1 s, 2,712 MB | identical |

Unmodified main's score mode takes 193 s on chr1; the RTR branch (both changes) takes ~60 s.
Earlier "origin/main" memory figures in §3 were measured the same way and agree (2,704 MB).
Caution for reruns: inside a git worktree, `import hseeker` resolves to the editable install of
the main checkout unless `PYTHONPATH=src` is set.
