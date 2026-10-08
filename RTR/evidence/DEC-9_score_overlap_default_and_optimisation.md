# DEC-9 — score-informed overlap removal becomes the default; its cost is removed without changing results

**Status:** done, verified (2026-10-07). **Default superseded 2026-10-08 by DEC-10:** the longest-arm (greedy) rule is the default again in the library, CLI and webapp (author instruction "make the greedy the default"), because score mode costs ~2.3x greedy on a genome (chr1, 16 workers: ~60 s vs 26.3 s; `WS6-A4_chr1_runtime_memory_2026-10-07.md`). Score mode remains available (`overlap_strategy="score"`, CLI `-overlap-strategy score`) and keeps the result-preserving optimisation in §3-4. §1 below describes the 2026-10-07 state.
**Decision:** author instruction, 2026-10-07 — "In overlap remove pick the default (scoring
informed)", then "optimize the score based overlapping while keeping the logic and results the
same."
**Code:** `src/hseeker/__init__.py`, `src/hseeker/_scoring.py`, `src/hseeker/__main__.py`,
`webapp/main.py`, `webapp/templates/index.html`, README; tests in `tests/test_score_overlap.py`
and `tests/test_hdna.py`.

## 1. What changed for users

Nikol's main (PR #16) added an optional score-based overlap strategy alongside the historical
longest-arm rule. DEC-9 makes it the default:

| Entry point | Before | Now |
|---|---|---|
| `scan_sequence` / `scan_fasta` / `scan_fasta_iter` / `scan_fasta_parallel` | `overlap_strategy="greedy"` | `overlap_strategy=None` → `"score"` when `score=True`, `"greedy"` when `score=False` |
| CLI `-overlap-strategy` | default `greedy` | default `score`; with `-no-score`, `greedy`; startup banner prints the mode |
| Webapp form | "Longest arm (greedy)" preselected; scoring off + score → HTTP 400 | "Highest scores first" preselected; scoring off → falls back to greedy |
| `hseeker.ablation` | explicit per-model handling | unchanged (ablation models without scores cannot rank by score) |

Explicitly requesting `overlap_strategy="score"` with `score=False` still raises `ValueError`, as
before; only the *default* adapts, because without scores there is nothing to rank by.

**Score-informed selection** (Nikol's definition, unchanged): every candidate is scored, then
hits are accepted in descending `total_score` order when they overlap no accepted hit; ties go to
the longer arm, the shorter spacer, then the earlier coordinate.

### Effect on reported results

- **No detection call changed on the benchmark.** Greedy and score were compared on every
  record under every configuration of the WS4-A3 grid (576 configs × 128 records = 73,728 calls,
  covering all four filter states): **0** differences in the binary call. Every metric in WS2-PRE,
  WS4-A3 and WS2-A3 is therefore unaffected by DEC-9 (they were nevertheless regenerated with it).
- **In general a call *can* differ when the homopolymer filter is on.** Without that filter,
  score selection always keeps the top-scoring candidate of an overlapping cluster, so "is there a
  hit" cannot change. With it, greedy removes overlaps *before* the filter: if its longest-arm pick
  is a pure homopolymer, the filter drops it and nothing is left, whereas score selection can still
  report an overlapping non-homopolymer candidate. Example: `ACGT×10 + G×100 + C×100 + ACGT×10`
  with the CLI filters gives no hit under greedy and 3 under score — the G run extended by one
  flanking base (`gtggggggggg[gggggggggg]ggggggggggg`, `gggggggggg[gggggggggg]gggggggggc`, and
  the C run with `…ac`, reported on the purine strand). The same happens for a G run flanked by
  random sequence. These are near-homopolymer edge hits the homopolymer filter is not designed to
  catch; the authors may want to say whether they should be reported.
- **Which hit is reported can change** inside overlapping clusters. Examples: on
  `GAAGGAAAGAAAGAAAGGG` score selection reports (3, 19) with a higher `total_score` than greedy's
  (1, 17); on the GAMIR test motif it reports the higher-scoring arm-6 hit where greedy kept the
  arm-7 hit. This is the purpose of the change.
- The WS3 smoke-test values for HSeeker are unchanged under the new default (pGG32 arm 15 /
  spacer 2 / score 152.55; pcMyc arm 8 / spacer 7 / 92.03; poly-A20 83.30; SYN0001 no hit;
  FXN_GAA66 arm 98 / 984.75).

## 2. The problem the default exposed

Score mode has to score every raw candidate before it can choose, and in long low-complexity
runs almost every sub-window is a candidate. Each candidate's scoring was itself O(arm³) (the
boundary search rebuilt the stacking sum for every window). Measured on this machine before
optimisation:

| Input | greedy | score (before) |
|---|---|---|
| poly-G 100 bp | 0.002 s | 0.43 s |
| poly-G 200 bp | 0.014 s | 6.5 s |
| poly-G 400 bp | 0.1 s | **96 s** |
| `G×3000 + C×2000 + A×5000` | 1.6 s | did not finish |
| 2 Mb random sequence with GA/CT repeat islands | 8.9 s | 39 s |

With score as the default, 11 parallel-scan tests timed out (> 60 s each) and the suite went
from ~11 s to > 11 min. Real genomes contain long homopolymer runs and microsatellites, so this
would have made default genome-wide scans impractical.

## 3. The optimisation (same logic, same results)

1. **Incremental window search** (`_Scorer._score_impl`). For each left boundary L, the window
   [L, R) is extended one base at a time and its pairing and stacking sums are updated with the
   next term. The terms are added in the order `sum(pairing_scores[L:R])` and
   `sum(_calc_stacking(scoring_array[L:R]))` add them, and windows are compared after
   `round(..., 3)`, visited in the same order with the same strict `>` tie rule. On Python ≥ 3.12
   `sum()` uses compensated summation, so a raw running sum can differ from `sum()` in the last
   bits; that cannot change a comparison unless a window score lies within ~1e-12 of a rounding
   boundary, and the reported `pairing_score`/`stacking_score` are still computed with `sum()`.
   The equivalence harness (§4, run on Python 3.14) and Nikol's golden-output tests
   (`tests/test_parallel_batching.py`) found no difference.
   Cost per hit: O(arm³) → O(arm²).
2. **A cheap upper bound** (`_Scorer.upper_bound`). No window can score more than the sum of the
   positive per-position pairing scores plus +5 for each adjacent pair of matched positions
   (every other stacking term is a penalty). Computed in O(arm) with string-translate bitmasks and
   rounded to 3 decimals like the scores themselves (rounding is monotone, so it stays a valid
   bound, and a candidate that attains it ties exactly).
3. **Lazy score-ordered selection** (`_score_and_select`). Candidates enter a heap keyed by the
   selection priority with the score replaced by its bound. A candidate popped with its exact
   score known is provably the next one in the full sort, so it is accepted or rejected exactly
   as before. A candidate popped with only its bound is **discarded unscored** if it already
   overlaps an accepted hit — accepted hits only accumulate, so it would be rejected at its true
   turn as well; otherwise it is scored and pushed back. The POST homopolymer filter is applied
   per hit, which is equivalent because it depends only on the hit itself. This replaces
   "score everything → filter → select" in `scan_sequence`, `scan_fasta` and
   `scan_fasta_parallel`; `_filter_overlapping_hits` itself is unchanged.

### Timings after optimisation

| Input | greedy | score (before) | score (now) |
|---|---|---|---|
| poly-G 400 bp | 0.01 s | 96 s | 0.06 s |
| poly-G 1000 bp | 0.08 s | — | 0.27 s |
| poly-G 3000 bp | 0.44 s | — | 1.5 s |
| `G×3000 + C×2000 + A×5000` (no filters) | 1.6 s | did not finish | 5.4 s |
| 2 Mb with repeat islands | 7.8 s | 39 s | 8.8 s |
| Test suite | — | > 11 min, 11 timeouts | **225 passed, 3 skipped, ~21 s** |

Two further changes, also result-preserving, address the homopolymer filter (always on in the
CLI): candidates whose detected `full_sequence` is a single base are dropped *before* scoring when
the filter is on (the POST filter removes them whatever their score), and the filter's
single-base test is now `len(set(seq.lower())) == 1` instead of a per-base Python loop (identical
on 791 edge cases, including the empty string and mixed case). Poly-G 300 bp with filters went from
6.4 s to 0.07 s.

**Remaining worst case.** A long G or C homopolymer run *with the homopolymer filter on* still
costs far more under score than under greedy. Nearly every candidate near the run's edge is
"run + one flanking base"; its best window trims that base, so its triplex is a homopolymer and the
filter drops it — but only after it has been scored, and since none is ever accepted, nothing is
pruned. Cost grows roughly with the cube of the run length:

| Input (CLI filters on) | greedy | score |
|---|---|---|
| G100 + C100 | < 0.01 s | 0.57 s |
| G200 + C200 | 0.01 s | 4.8 s |
| G400 + C400 | 0.02 s | 40 s |
| G800 + C800 | 0.08 s | 330 s |
| random flank + G200 + random flank | 0.01 s | 2.4 s |
| G3000 alone | 0.45 s | 3.7 s |

The original implementation scored every candidate as well, so this is not a regression of the
optimisation, but it is a cost of making score the default. G/C homopolymer runs of hundreds of
bases are rare in natural genomes; for inputs that contain them, `-overlap-strategy greedy` avoids
it. Removing it entirely would require replacing the per-candidate O(arm²) window search with an
O(arm) algorithm, which was not done here.

Greedy also became slightly faster (8.9 → 7.8 s on 2 Mb) because every surviving hit is scored
with the O(arm²) search.

## 4. Verification that results are unchanged

**In-repo tests (permanent guards, all passing):**

- `test_lazy_score_selection_matches_scoring_every_candidate` — 12 cases (poly-G 120 bp,
  (GA)60, (GAA)40, a CT/G composite, two real motifs; with and without the homopolymer filter):
  lazy selection returns exactly what "score every candidate → POST filter →
  `_filter_overlapping_hits(score)`" returns.
- `test_incremental_window_search_matches_naive_recomputation` — 300 random arms: the chosen
  window, pairing, stacking and total scores equal a from-scratch per-window recomputation, and
  the bound is never below the score.
- `test_default_is_score_with_scoring_and_greedy_without` and
  `test_cli_default_is_score_and_no_score_falls_back_to_greedy` — the new default semantics in the
  API and CLI.
- 24 detector-geometry tests in `tests/test_hdna.py` (GAMIR/CTMIR arm-7 hit, N-flanked motif,
  dual motif, embedded coordinates) now request `overlap_strategy="greedy"` explicitly, because
  they assert the longest-arm hit; their expectations are unchanged.

**Equivalence against the pre-optimisation code** (one-off harness, 2026-10-07): the
pre-optimisation `__init__.py` and `_scoring.py` were loaded side by side with the optimised
modules and compared dict-for-dict.

Final run on the final code, all **passed**:

| Comparison | Scope | Result |
|---|---|---|
| Scorer (`_Scorer.score`) | every raw candidate (minrep 6/8/10 × maxspacer 15/10/5) of 438 sequences — random, purine-rich, pyrimidine-rich, homopolymers to 150 bp, 8 repeat units × 3 lengths, N-containing repeats, all 128 benchmark and all 80 experimental-CSV sequences — **70,984 candidates × 4 flag combinations** (stacking on/off × boundary optimisation on/off) | identical dicts; `upper_bound` ≥ `total_score` for every candidate |
| `scan_sequence` | the same 438 sequences × 5 parameter sets × 3 strategies (`score`, default, `greedy`) × 3 filter settings, alternating `purity_rmq` — **19,710 calls** | identical |
| `scan_fasta`, `scan_fasta_iter`, `scan_fasta_parallel` | 41-record FASTA × 3 parameter sets × 2 strategies × homopolymer filter on/off; parallel at chunk sizes 60, 150 and 10,000 with 3 workers | identical |
| Long runs | G×300, (GA)×150, (GAA)×100, (CT)×120 + ACGT×10 + G×200, filter on and off | identical; e.g. G×300: 0.07 s vs 32.8 s; (GA)×150: 0.04 s vs 14.7 s |

The harness (`equiv.py`) was a one-off script outside the repository, because the
pre-optimisation code is not committed; the in-repo tests above are the permanent guards.

## Reproduce

```
python -m pytest tests -q                       # 225 passed, 3 skipped
python -m pytest tests/test_score_overlap.py -q # default semantics + equivalence guards
```
