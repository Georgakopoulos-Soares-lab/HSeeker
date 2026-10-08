# WS2 pre-work — first HSeeker baseline on benchmark v3

**Regenerated 2026-10-07 on 128 records (64:64) after DEC-8 rebalance and DEC-9
score-informed overlap default; 124-record values in this file's 2026-10-02 version.**
DEC-8 added four synthetic non-forming records (SYN0060–SYN0063: two mirror-disrupted
mutants, one dinucleotide shuffle, one length/GC-matched random control) in place of
the excluded homopolymer negatives, restoring exact 1:1 balance; the 124 earlier rows
are byte-identical. DEC-9 made `overlap_strategy="score"` the default whenever scoring
was on; no detection call changed on the benchmark (0 of 73,728 grid calls), although
in general it can with the homopolymer filter on (see "Overlap rule" below). DEC-10
(2026-10-08) made greedy the default again; the numbers below are the same under either
rule.

**Label corrected 2026-10-08.** Earlier versions of this file called the minrep-8
configuration "library defaults". It is the **benchmark/script configuration** (minrep 8,
mismatch 0.10, no filters), i.e. what the analysis scripts use. The software defaults use
minrep 10 (library: no filters; CLI: both filters) and give 56/8/64/0, MCC 0.882; both are
shown in the tables below.

**Regenerated 2026-10-02** after merging main (Nikol's scanner) and removing the 6
pure homopolymers; previous values in git history at `007288b`. Scanner changes that
matter here: the AT filter now tests the full motif (both arms + spacer) rather than
the left arm, the scoring boundary search no longer shrinks an arm below `min_al`,
and the homopolymer filter has PRE/POST modes. Benchmark: 130 → 124 (2026-10-02) →
**128 records (64 forming / 64 non-forming — 1:1 again)** (2026-10-07).

**Status:** PROVISIONAL. This is **not** WS2-A2 and must not be cited in the
response letter. See "Why provisional" below.
**Date:** 2026-09-23 (first run); regenerated 2026-10-02 and 2026-10-07
**Input:** `hdna_benchmark_balanced_v3.csv`, sha256 `7bb495444e91709f…` (was
`7f78555dd7c788a7…` at 124 records, `f4602d6453d0926c…` before the homopolymer removal)

## Why this exists

The deposited results in `analysis/results/` were run on **2026-06-24** against the
**old** CSV (`hdna_experimental_sequences_final.csv`, 80 rows, 71 forming / 9
non-forming), on a different machine, with code predating the AT/homopolymer filter
commit (`6ecaa84`, 2026-07-02). **No HSeeker result had ever been produced on the v3
benchmark.** Establishing that baseline is a prerequisite for every comparator
comparison (WS3-A3 depends on WS2-A2), so it was run before going further with
NeSSie/nBMST.

## Result — HSeeker detection call, current script parameters

`minrep=8, maxspacer=10, purity=0.90, mismatch=0.10, purity_rmq=True` — the values
the direct script uses today, **not** the manuscript's (DEC-5 unresolved), and **not** the
software defaults (minrep 10; see the configuration tables below).

| subset | n | pos | neg | TP | FN | TN | FP | sens | spec | prec | MCC | acc | AUC |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| balanced (all) | 128 | 64 | 64 | 62 | 2 | 64 | 0 | 0.969 | 1.000 | 1.000 | 0.969 | 0.984 | 0.984 |
| experimental | 69 | 64 | 5 | 62 | 2 | 5 | 0 | 0.969 | 1.000 | 1.000 | 0.832 | 0.971 | 0.984 |
| synthetic | 59 | 0 | 59 | 0 | 0 | 59 | 0 | — | 1.000 | — | — | 1.000 | n/a |

AUC is computed as in `run_direct_sensitivity_analysis.py`
(`run_hseeker_direct` → best hit's `total_score`, records with no hit scored −∞,
`auc_roc_from_scores`). Every hit now carries a score (none is `None` after the
scoring fix). **Read the AUC with care:** no negative has any hit, so every negative
sits at −∞ and the AUC reduces to (1 + sensitivity)/2 = 0.984 — the two unhit
positives tie with all negatives and score 0.5 each. It says nothing about how well
`total_score` ranks; on this benchmark there is nothing for it to rank. (Old 130-record
values: 0.963 balanced, 0.977 experimental; 124-record balanced MCC was 0.968 on
62/2/60/0.)

The balanced set does what it was built to do: specificity rests on **64**
negatives instead of the 5 experimental ones, so it is no longer hostage to one or
two records — the exact objection in R3.1. The homopolymer removal (five negatives,
one positive) had left the classes at 64:60; DEC-8's four replacement negatives
restore exactly 1:1 (64:64). **All four new negatives are true negatives** under every
configuration below (no hit at all), so the only metric changes from the 124-record
round are the TN count and the MCC that depends on it (0.968 → 0.969).

## What happened to the false positives — the benchmark no longer contains any

On the original 130-record v3 benchmark the direct-script configuration produced four
false positives, and **every one was a pure homopolymer**:

| Record | Type | Sequence | status now |
|---|---|---|---|
| HDNA0053 `pRW1405` | experimental | `A`×20 | removed |
| SYN0054 `poly-A30` | synthetic | `A`×30 | removed |
| SYN0055 `poly-T25` | synthetic | `T`×25 | removed |
| SYN0056 `poly-T35` | synthetic | `T`×35 | removed |

By author decision (2026-10-02) all six pure single-base homopolymers were removed
from the benchmark — these four, plus SYN0053 `poly-A15` (which was already a true
negative) and the positive HDNA0054 `poly-(dC)30`. **With them gone there are no
false positives left**: the script configuration (and equally the software defaults)
now scores 64/64 negatives correctly (the 60 that survived the removal plus the four
DEC-8 replacements).

The consequence for the filters is direct, and must be stated plainly rather than
read as a success: the AT-content and homopolymer filters target exactly the negatives
that were removed, so **on this benchmark they can no longer be evaluated** — there is
nothing left for them to catch. Switching them on changes no call:

### Balanced set (n=128)

| configuration | TP | FN | TN | FP | sens | spec | MCC |
|---|---|---|---|---|---|---|---|
| benchmark/script configuration (minrep 8, no filters; **what the direct script uses**) | 62 | 2 | 64 | 0 | 0.969 | 1.000 | 0.969 |
| + AT filter only | 62 | 2 | 64 | 0 | 0.969 | 1.000 | 0.969 |
| + homopolymer filter only | 62 | 2 | 64 | 0 | 0.969 | 1.000 | 0.969 |
| software defaults, library (minrep 10, no filters) | 56 | 8 | 64 | 0 | 0.875 | 1.000 | 0.882 |
| **software defaults, production CLI (minrep 10, both filters)** | 56 | 8 | 64 | 0 | 0.875 | 1.000 | 0.882 |

### Experimental only (n=69)

| configuration | TP | FN | TN | FP | sens | spec | MCC |
|---|---|---|---|---|---|---|---|
| benchmark/script configuration (minrep 8, no filters) | 62 | 2 | 5 | 0 | 0.969 | 1.000 | 0.832 |
| software defaults, library (minrep 10, no filters) | 56 | 8 | 5 | 0 | 0.875 | 1.000 | 0.580 |
| software defaults, production CLI (minrep 10, both filters) | 56 | 8 | 5 | 0 | 0.875 | 1.000 | 0.580 |

Synthetic only (n=59): 59 TN / 0 FP under every configuration (including the four
DEC-8 records SYN0060–SYN0063).

**Reading:** the filters are no-ops on the current benchmark. The lower sensitivity of
the software defaults (production CLI and library alike; 0.969 → 0.875, MCC 0.969 →
0.882) comes **only from `minrep=10`**: re-running the production configuration with both filters switched off
gives the identical 56/8/64/0, and the six extra misses (HDNA0022, HDNA0023,
HDNA0024, HDNA0025, HDNA0027, HDNA0030) are lost to the longer minimum arm, not to a
filter. On the experimental-only set MCC falls 0.832 → 0.580 for the same reason;
with five negatives MCC is dominated by the positive class there, which is itself an
argument for dual reporting.

(124-record values, 2026-10-02: script configuration 62/2/60/0 MCC 0.968; production CLI
56/8/60/0 MCC 0.879. Old 130-record result, for the record: script configuration
63/2/61/4 MCC 0.908; either filter alone 62/3/65/0 MCC 0.955; production CLI
61/4/65/0 MCC 0.940. There the filters removed the four homopolymer FPs at a cost of 1–2 TPs.)

**Consequence:** register **D-4** (direct script without filters vs production filters)
no longer moves the headline on this benchmark — both arms give specificity 1.000. The two arms
still have to be pinned to one configuration before any number is quoted (WS2-A1),
and Methods must state which, but the difference between "direct script" and
"production CLI" is now entirely the `minrep` setting (8 vs 10). Which configuration the
manuscript reports, and whether the shipped default minrep should change from 10 to 8, is
an open `TODO(authors)` (plan, under WS2; added 2026-10-08). Any claim about what
the AT/homopolymer filters contribute needs evidence from somewhere other than this
benchmark (e.g. genome-scale scans), and Methods should say that the benchmark
excludes pure homopolymers by design.

## Overlap rule (DEC-9, DEC-10) — no effect on these numbers

These numbers were generated on 2026-10-07, when the library and CLI resolved
overlapping hits by score (`overlap_strategy="score"`) whenever scoring was on (DEC-9).
Since 2026-10-08 (DEC-10) the default is greedy (longest arm) again, with score mode as an
option. The rule changes **which** hit is reported inside an overlapping cluster. Without
the homopolymer filter it cannot change **whether** a record has a hit, because score
selection always keeps the top-scoring candidate. With the filter on it can in general:
greedy's longest-arm pick may be a pure homopolymer that the filter drops, leaving nothing,
whereas score selection can report an overlapping near-homopolymer candidate (see the
DEC-9 evidence for an example). On the benchmark this never happens: 0 binary-call
differences between greedy and score on all 128 records, and 0 of 73,728 calls across the
full WS4-A3 grid (576 configurations, all four filter states). Every number in this file
(run with scoring on) is therefore unaffected by DEC-9 and by DEC-10.

## The two false negatives are the same motif class

| Record | Sequence | Note |
|---|---|---|
| HDNA0058 `cMYC_NSE_ACCCTCCCC4` | `ACCCTCCCC` ×4 | pyrimidine-rich tandem **direct** repeat |
| HDNA0078 `gamma_globin_228_189` | `CTCCCAGC…` ×4 | pyrimidine-rich tandem **direct** repeat |

Both are ~36–37 nt tandem repeats of a ~9-nt unit, i.e. direct-repeat architecture
rather than the arm–spacer–arm mirror HSeeker searches for. Worth an explicit
sentence in the FN discussion (WS2-A5) rather than leaving them unexplained; a
reviewer who checks will notice they share a structure.

## Why provisional — four reasons it is not WS2-A2

1. **DEC-5 and the reporting configuration unresolved.** Run with the script's current
   `mismatch=0.10, minrep=8`, not the manuscript's `mismatch=0.15` and not the software
   defaults' `minrep=10`. Numbers will move (on the 128-record set,
   mismatch 0.15 with minrep 8 gives 63/1/63/1, MCC 0.969 with or without the
   filters: it recovers HDNA0058 but admits HDNA0028 `B-33` as a false positive —
   see WS4-A3).
2. **Benchmark not frozen (WS0-A1).** Blocked on DEC-2/DEC-3. If adopted, DEC-3
   would add a positive (→129 rows, 65:64 on the current 128-record base); the pCON
   terminal-C correction (D-9) changes one negative's sequence. The homopolymer
   removal (2026-10-02) and the DEC-8 rebalance (2026-10-07) have already changed
   every number twice.
3. **HSeeker only — no Triplex comparator.** R/Rscript is not installed in this
   environment, so the direct script's Triplex arm could not run. Every
   tool-vs-tool claim is still unsupported.
4. **No confidence intervals, no family-aware evaluation.** Point estimates only;
   WS2-A2/WS2-A3 remain outstanding here, and R2.7/R3.1/R3.2/R3.3 require them
   (CIs and family-aware evaluation are in `WS4-A3_parameter_robustness.md` and
   `WS2-A3_family_aware_validation.md`).

## Reproduce

```
python3 -c "…"   # uses load_rows/run_hseeker_direct/metrics_from_binary
                  # from analysis/scripts/run_direct_sensitivity_analysis.py
                  # with --subset all|experimental|synthetic; the filter rows call
                  # hseeker.scan_sequence with at_threshold=0.8 /
                  # filter_homopolymers=True; AUC via auc_roc_from_scores on the
                  # run_hseeker_direct hseeker_score column
```
