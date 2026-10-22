# WS1-A2 — profiling counters must not be shared across threads (reviewer R2.9)

**Status (2026-10-02):** resolved — RTR implementation superseded by Nikol's fix on main (merged)
**Date:** 2026-09-23 (RTR fix); updated 2026-10-02 after merging main
**Reviewer item R2.9 (verbatim):** "There is no synchronization in mutating
prof_inner_iters and prof_ctr_sp_pairs in src/hseeker/_hdna.c. It may lead to
losing the counting."

## Reviewer's reading confirmed — and the defect is larger than "may"

Before either fix, `_hdna.c:54-61` declared both counters as plain static globals.
They are incremented inside the two scan cores (`findHDNA_core`,
`findHDNA_purity_rmq_core`) which run **with the GIL released**
(`Py_BEGIN_ALLOW_THREADS` in `py_scan_sequence`), and `scan_fasta_parallel()`
dispatches chunk scans across a `ThreadPoolExecutor`. So several OS threads
incremented the same two words concurrently, with no synchronization.

`prof_reset()` made it worse: each new scan zeroed the shared counters, so a scan
starting mid-flight reset a concurrent scan's totals.

**Measured, not assumed.** Eight concurrent scans of *identical* input, compared
against the single-threaded total for that same input (original, pre-fix build):

```
single-thread baseline : inner_iters = 1,823,346   ctr_sp_pairs = 659,736

8 concurrent scans:
  distinct ctr_sp_pairs values: 695711, 707816, 713698, 720248, 743261, 744255, ...  (8 distinct)
  distinct inner_iters  values: 6526379, 6646285, 6723703, 6757692, 6905192, ...     (8 distinct)
  max inner_iters 6,993,228 vs baseline 1,823,346  ->  3.84x inflated
```

Every thread reported a different number, none correct. The counts were not merely
"lost" — they were cross-contaminated, each thread reading a total that included
other threads' work.

## Resolution (current tree) — per-scan counters, published under the GIL

The shared static counters are gone. In the merged tree (`src/hseeker/_hdna.c`):

- `ScanProfile` (`_hdna.c:54-58`, `{inner_iters, ctr_sp_pairs}`) is a local in each
  `py_scan_sequence` call (`_hdna.c:599`) and is passed by pointer to the scan cores
  (`findHDNA_core`, `findHDNA_purity_rmq_core`; `_hdna.c:266`, `:383`). While the
  GIL is released (`Py_BEGIN_ALLOW_THREADS`, `_hdna.c:680-699`) each scan writes only
  its own struct — no shared mutation, no atomics in the hot loop.
- Module state `ProfileState` (`_hdna.c:61-66`: `last`, `total`, `scans_completed`)
  is updated by `profile_publish()` (`_hdna.c:68-75`) only **after the GIL is
  reacquired** (`_hdna.c:710`; empty input publishes at `:641`).
- Python API: `_hdna.profiling_info()` returns `inner_iters`, `ctr_sp_pairs` (last
  completed scan), `total_inner_iters`, `total_ctr_sp_pairs`, `scans_completed`;
  `_hdna.reset_profiling()` zeroes the module state (`_hdna.c:838-862`).

So a concurrent scan can neither reset nor contaminate another's counts, and the
totals across threads are exact. Still diagnostic-only: no consumer outside the test.

## Verification

**1. Totals across threads are exact** (merged tree, 2026-10-02): 300 kb random
sequence, default parameters; 8 concurrent scans vs one scan:

```
single scan : inner_iters = 2,524,852,128   ctr_sp_pairs = 6,299,475
8 concurrent: total_inner_iters = 20,198,817,024 (= 8x)   total_ctr_sp_pairs = 50,395,800 (= 8x)
scans_completed = 8; last-scan counts equal the single-scan counts
```

**2. Regression tests** — `tests/test_profiling_threadsafe.py` (2 tests, both pass):

- `test_parallel_scans_preserve_every_profile_count`: 40 scans (both scan cores,
  8 workers); `scans_completed` and both totals equal the sum of the per-scan
  single-threaded counts.
- `test_empty_scan_has_zero_last_counts_and_keeps_totals`: an empty scan reports
  zero last-scan counts and leaves totals unchanged.

**3. Full suite on the merged tree: 211 passed, 3 skipped** (chr1 tests skip without
`benchmarks/data/chr1.fa`).

## History (superseded RTR implementation)

RTR (2026-09-23) declared both static counters thread-local (`HDNA_THREAD_LOCAL`
macro: `__declspec(thread)` / `__thread` / `_Thread_local`), leaving the module
surface unchanged. Verified then: 8 concurrent identical scans each reported exactly
the single-thread baseline (1,823,346 / 659,736); parallel `scan_fasta` output
identical to serial for 2–16 workers; regression test
`test_profiling_counters_are_thread_local` (tests/test_hdna.py) failed on the pre-fix
build and passed on the fix; suite 150 passed, 3 skipped. Nikol fixed the same defect
independently on main; on merging, the authors chose Nikol's version (2026-10-02).
It additionally gives exact cross-thread totals (thread-local copies could not be
summed from Python). The macro and `test_profiling_counters_are_thread_local` were
removed (the test was redundant with the new tests).

## Note for the response letter

The honest framing is that the reviewer under-stated the problem: it is not only
that counts may be lost, but that any reported count from a parallel run was wrong
by an unbounded factor. Since the counters are diagnostic and unused, **no published
number in the manuscript is affected** — the benchmark timings come from
`/usr/bin/time -v` and Python `perf_counter`, not from these counters.

## Files changed (current tree, Nikol's fix)

- `src/hseeker/_hdna.c` (`ScanProfile` / `ProfileState`, `profile_publish`,
  `reset_profiling`, extended `profiling_info`)
- `tests/test_profiling_threadsafe.py` (new, 2 tests)

## Reproduce

```
python3 -m pytest -q tests/test_profiling_threadsafe.py
```
