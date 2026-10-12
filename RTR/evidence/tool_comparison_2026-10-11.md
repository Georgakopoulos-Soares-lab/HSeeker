# Comparison with triplex, Triplexator, NeSSie and nBMST — results (2026-10-11)

**Design:** `RTR/TOOL_COMPARISON_PLAN_2026-10-11.md` (settings, hit rules, statistics and the wording rule
committed in `caa7abb` and `ca63f0b` before any comparator was run on the benchmark).
**Script:** `analysis/scripts/run_tool_comparison.py` (+ `triplex_search_params.R`); builds
`benchmarks/tools/build_comparators.sh`. **Outputs:** `analysis/results/tool_comparison/`
(`summary_metrics_v2.csv`, `paired_differences.csv`, `tuned_out_of_fold.csv`, `by_negative_class.csv`,
`homopolymer_panel_tools.csv`, `per_record_calls.csv`, `verdict.csv`, `raw/`, `metadata.json`).
No tool run failed.

**Benchmark:** `hdna_benchmark_balanced_v3.csv` (sha256 `c3f3b842…`), 130 records (65 forming / 65
non-forming); also the 70 experimental and 58 primary records. Call rule for every tool: at least one
qualifying hit in the record. HSeeker at its shipped defaults.

## 1. Main comparison (each tool at its documented defaults, setting A), balanced benchmark

| Tool (hit type) | TP | FN | TN | FP | MCC (95 % CI, family bootstrap) | HSeeker − tool, MCC (95 % CI, family) | McNemar (HSeeker right / tool right) |
|---|---|---|---|---|---|---|---|
| **HSeeker** | 63 | 2 | 65 | 0 | **0.970** (0.909–1.000) | — | — |
| triplex (any hit) | 62 | 3 | 56 | 9 | 0.819 (0.743–0.935) | +0.151 (0.027–0.224) | 10 / 0, p = 0.002 |
| NeSSie (`-T` motif) | 33 | 32 | 65 | 0 | 0.583 (0.421–0.779) | +0.386 (0.178–0.559) | 30 / 0, p < 10⁻⁶ |
| nBMST (triplex subset) | 24 | 41 | 64 | 1 | 0.449 (0.351–0.561) | +0.521 (0.400–0.626) | 40 / 0, p < 10⁻⁶ |
| Triplexator (`-ds` target site) | 8 | 57 | 58 | 7 | 0.024 (−0.132–0.250) | +0.946 (0.709–1.108) | 62 / 0, p < 10⁻⁶ |

## 2. Matched geometry (setting B), manuscript settings, secondary hit types

| Tool, setting | TP | FN | TN | FP | MCC |
|---|---|---|---|---|---|
| triplex, matched (min_len 8, max_len 50) | 62 | 3 | 56 | 9 | 0.819 |
| triplex, manuscript (… prokaryotic) | 62 | 3 | 56 | 9 | 0.819 |
| Triplexator, matched (`-l 16 -e 10`) | 25 | 40 | 42 | 23 | 0.032 |
| Triplexator, manuscript (`-l 10 -L -1 -e 15 -fr off`) | 64 | 1 | 23 | 42 | 0.436 |
| NeSSie `-T`, matched (`-p 10 -m 10`) | 39 | 26 | 65 | 0 | 0.655 |
| NeSSie `-M` (mirror), defaults / matched | 34 / 65 | 31 / 0 | 62 / 25 | 3 / 40 | 0.529 / 0.488 |
| nBMST triplex subset, matched (rep 8, spacer ≤ 10) | 30 | 35 | 62 | 3 | 0.477 |
| nBMST any mirror repeat, defaults / matched | 41 / 47 | 24 / 18 | 61 / 58 | 4 / 7 | 0.598 / 0.624 |

## 3. Tuned out of fold (setting C; selection inside training folds only)

| Tool (grid size) | 5-fold × 20 | Leave-one-family-out (95 % CI) | Leave-one-study-out |
|---|---|---|---|
| HSeeker (576) | 0.941 | 0.939 (0.848–0.988) | 0.909 |
| triplex (16) | 0.819 | 0.819 (0.744–0.937) | 0.819 |
| nBMST (12) | 0.598 | 0.597 (0.446–0.812) | 0.556 |
| NeSSie (21) | 0.573 | 0.547 (0.391–0.732) | 0.524 |
| Triplexator (18) | 0.465 | 0.520 (0.376–0.746) | 0.520 |

Paired, leave-one-family-out, HSeeker − tool MCC: triplex +0.120 (−0.024–0.210), Triplexator +0.419
(0.170–0.593), NeSSie +0.392 (0.193–0.539), nBMST +0.342 (0.106–0.507).

## 4. Where the errors come from (setting A)

- **triplex**: false positives are 7 of 22 mirror-disrupted mutants and 2 of 11 dinucleotide shuffles (none
  among random, G4, Z-DNA, B-DNA, balanced perfect mirrors or the 5 experimental negatives); misses
  HDNA0042 (pRW1412, an interrupted (GAA)9 tract) in addition to HSeeker's two misses.
- **Triplexator** (defaults) finds 8 of 65 forming records; at the manuscript settings it finds 64 but calls
  42 negatives, including all 22 mirror-disrupted mutants and 9 of 11 shuffles. It reports purine target
  sites for a third strand and does not require a mirror repeat.
- **NeSSie `-T`** makes no false call but finds 33 of 65 (e.g. none of the pXY32 series at the default
  lengths); `-M` with mismatches and gaps finds all 65 but calls 40 negatives.
- **nBMST** triplex subset finds 4 of the 16 pXY32 records and **none of the 19 forming records written on
  the pyrimidine-rich strand** (its any-mirror-repeat output finds 15 of those 19). This is the
  integer-division rule of plan §11, which only accepts mirror repeats with more purines than pyrimidines
  on the + strand. Its any-mirror-repeat output calls all three balanced perfect mirrors.
- **NeSSie `-T`** finds none of the 16 pXY32 records; triplex and HSeeker find all 16.
- **HSeeker** makes no false call in any negative class.

## 5. Experimental and primary records only (setting A)

| Tool | Experimental (70): TP/FN/TN/FP, MCC | Primary (58): MCC | HSeeker − tool, experimental (family CI) |
|---|---|---|---|
| HSeeker | 63/2/5/0, 0.832 | 0.829 | — |
| triplex | 62/3/5/0, 0.772 | 0.768 | +0.060 (0.000–0.223); McNemar 1 / 0 |
| NeSSie | 33/32/5/0, 0.262 | 0.231 | +0.570 (0.285–0.802) |
| nBMST | 24/41/5/0, 0.200 | 0.198 | +0.632 (0.395–0.844) |
| Triplexator | 8/57/5/0, 0.100 | 0.094 | +0.732 (0.478–0.935) |

## 6. Homopolymer panel (setting A)

HSeeker: no hit in any of the six (correct for the five A/T homopolymers, wrong for poly-(dC)30). triplex
and NeSSie: hit in A20, A30, T25, T35 and poly-(dC)30 (4 of 5 non-forming called forming). Triplexator: no
hit. nBMST triplex subset: hit in A30 only.

## 7. Verdict under the pre-specified wording rule (plan §2.4)

"HSeeker better" for all four tools: the family-cluster CI of the paired MCC difference under setting A
excludes 0 and the difference is positive under B and C. **Qualification that the text must carry:**
against triplex the difference arises entirely from the designed negatives (mirror-disrupted mutants and
shuffles); on the experimental records alone the two differ by one record and are comparable, and the
out-of-fold (setting C) interval for triplex includes 0. So the accurate wording is: HSeeker outperforms
Triplexator, NeSSie and nBMST on every set, and outperforms triplex on the balanced benchmark through
higher specificity on mirror-disrupted negatives, with comparable performance on experimental records.

## 8. Speed and memory on hg38 chr1 (same node and day as HSeeker)

| Tool (version) | Cores | Wall time, s (mean of 3; range) | Peak RSS (MB) | Output on chr1 |
|---|---|---|---|---|
| HSeeker (final code, defaults) | 16 threads | **16.3** (16.3–16.4) | 1,334 | 45,927 loci |
| Triplexator 1.3.2, defaults (`-ds -rm 1 -p 16`) | 16 threads (OpenMP) | 45.0 (43.5–47.8) | 1,714 | 275,476 target sites |
| triplex 1.50.0, defaults (`benchmarks/triplex_benchmark.R`, 1-Mb windows) | 16 forked processes | 64.7 (63.2–65.7) | 660 per process | 81,906 hits |
| nBMST / non-B_gfa `a891b59`, mirror-repeat module only | 1 (no parallel mode) | 96.5 (96.2–96.9) | 866 | 147,222 mirror repeats (33,519 classified triplex) |
| NeSSie `dbe6cb2`, `-T`, k = 16–20 | 1 (no parallel mode) | 295.7 (one run) | 1,056 | 720,811 motif occurrences |

Same node (AMD EPYC 7763, 16 vCPUs, Rocky Linux 8.10), same day, nothing else running; `/usr/bin/time -v`.
Every tool gave identical output in each repeat; triplex and Triplexator reproduce the hit counts of the
June 2026 measurement exactly (81,906 and 275,476). For triplex, `/usr/bin/time` reports the largest single
forked process, not the sum over the 16 workers. Raw output: `benchmarks/results_2026-10-11/comparators/`.
