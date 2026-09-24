# WS1-A2 — profiling counters must not be shared across threads (reviewer R2.9)

**Status:** done, verified
**Date:** 2026-09-23
**Reviewer item R2.9 (verbatim):** "There is no synchronization in mutating
prof_inner_iters and prof_ctr_sp_pairs in src/hseeker/_hdna.c. It may lead to
losing the counting."

## Reviewer's reading confirmed — and the defect is larger than "may"

Before the fix, `_hdna.c:54-61` declared both counters as plain static globals.
They are incremented inside the two scan cores (`findHDNA_core`,
`findHDNA_purity_rmq_core`) which run **with the GIL released**
(`Py_BEGIN_ALLOW_THREADS` in `py_scan_sequence`), and `scan_fasta_parallel()`
dispatches chunk scans across a `ThreadPoolExecutor`. So several OS threads
incremented the same two words concurrently, with no synchronization.

`prof_reset()` made it worse: each new scan zeroed the shared counters, so a scan
starting mid-flight reset a concurrent scan's totals.

**Measured, not assumed.** Eight concurrent scans of *identical* input, compared
against the single-threaded total for that same input (pre-fix build):

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

## Fix — thread-local storage, not atomics

Both counters are now declared `HDNA_THREAD_LOCAL`, via a portable macro:

```c
#if defined(_MSC_VER)
#  define HDNA_THREAD_LOCAL __declspec(thread)
#elif defined(__GNUC__) || defined(__clang__)
#  define HDNA_THREAD_LOCAL __thread
#elif defined(__STDC_VERSION__) && __STDC_VERSION__ >= 201112L
#  define HDNA_THREAD_LOCAL _Thread_local
#else
#  define HDNA_THREAD_LOCAL   /* single-threaded fallback */
#endif
```

Each thread accumulates privately; `prof_reset()` and `profiling_info()` act on the
calling thread's own copy.

**Why thread-local rather than atomics:** these counters are incremented in the
innermost scan loop (millions of times per scan — 1.8M in the measurement above).
An atomic read-modify-write per iteration would tax the exact hot path whose speed
is this tool's central claim. Thread-local storage costs essentially nothing and
removes the race completely.

**Why not delete them:** they are diagnostic-only and have **no consumer anywhere**
in the repository (verified: no reference in any `.py`, `.ipynb`, `.md` outside the
revision plan). Deletion is therefore also defensible and remains open to the
authors; thread-local was chosen as the smaller change, since it fixes the defect
without altering the `_hdna` module's surface.

The macro covers the full wheel matrix (MSVC / clang / gcc), so no platform loses
the counters.

## Verification

**1. Post-fix, every thread reports exactly the single-threaded total:**

```
baseline: inner_iters = 1,823,346   ctr_sp_pairs = 659,736
8 concurrent -> distinct ctr_sp_pairs: [659736]   distinct inner_iters: [1823346]
ALL threads exactly match baseline: True
```

**2. Parallel scan results unchanged and deterministic.** 6-record FASTA (~720 kb),
digest over `(seq_id, start, end, arm_length, full_sequence, total_score)`:

```
serial hits: 9   digest b6582aeefa925b9a
workers= 2  3 repeats identical to serial: True
workers= 4  3 repeats identical to serial: True
workers= 8  3 repeats identical to serial: True
workers=16  3 repeats identical to serial: True
```

**3. Regression test added and proven to bite.**
`test_profiling_counters_are_thread_local` (tests/test_hdna.py, section 24):

- against the **pre-fix** build: **FAILED** (as required)
- against the **fixed** build: **passed**

A regression test that passed on the broken code would guard nothing; this one was
explicitly checked both ways.

**4. Full suite: 150 passed, 3 skipped** (chr1 tests skip without
`benchmarks/data/chr1.fa`). Was 149/3 before this action.

## Note for the response letter

The honest framing is that the reviewer under-stated the problem: it is not only
that counts may be lost, but that any reported count from a parallel run was wrong
by an unbounded factor. Since the counters are diagnostic and unused, **no published
number in the manuscript is affected** — the benchmark timings come from
`/usr/bin/time -v` and Python `perf_counter`, not from these counters.

## Files changed

- `src/hseeker/_hdna.c` (counters → thread-local + rationale comment; +26/-3)
- `tests/test_hdna.py` (+1 regression test)

## Reproduce

```
python3 -m pytest -q tests -k thread_local
```
