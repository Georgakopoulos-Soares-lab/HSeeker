# Implementation Directive — Bounded DP / LCE Optimization of `findHDNA_core`

**Audience:** GitHub Copilot (agent mode) working inside the `HSeeker` repository.
**Target file:** `src/hseeker/_hdna.c` (plus new supporting files as specified below).
**Read `.github/copilot-instructions.md` in full before touching anything.** Every invariant listed there (struct layout, `is_perfect` float-equality semantics, 1-based coordinates, N-base handling, `cleanup:` error-path discipline) still applies. This document does not override that file — it adds a performance work-stream on top of it.

This is a **performance-only change**. There is no new feature, no new CLI flag behavior, and no change to any returned field. If you cannot prove bit-for-bit output parity with the current implementation, do not merge.

---

## 0. Non-negotiable constraints (read this twice)

1. **Byte-for-byte output parity.** For every `(seq, minrep, maxrep, maxspacer, purity, mismatch, seq_offset)` combination, the new code must produce exactly the same hit list (same count, same order, same field values, same float rounding) as the current `findHDNA_core` before overlap removal. Overlap removal and scoring are untouched — do not modify `remove_overlaps_core`, `_scoring.py`, or anything in `__init__.py`'s chunking logic as part of this work, except where Section 6 explicitly says to add a feature flag.
2. **No regression on the existing 120+ tests in `tests/test_hdna.py`.** All must pass unmodified (the test file itself is the oracle; do not edit expected values to make tests pass).
3. **Memory must scale with `chunk_size`, not with total genome length.** Do not build a single global suffix array over an entire chromosome/genome. See Section 4.2 and Section 5 for why and how.
4. **Land this behind a feature flag** (`HSEEKER_FAST_HDNA` env var or a private `_use_dp` kwarg threaded through `py_scan_sequence` — your choice, but it must be possible to flip back to the legacy loop with zero code changes, for bisection if a correctness bug surfaces post-merge).
5. **Do this in the phases in Section 6, as separate, independently testable, independently revertible commits/PRs.** Do not attempt a single monolithic rewrite of `findHDNA_core`. If you cannot fit a phase in one self-contained, testable diff, the phase is too big — split it further.
6. **When in doubt about a corner case, fall back to the legacy character-by-character loop for that specific `(ctr, sp)` pair rather than guessing.** Correctness > speed for this codebase (it's used for genomic motif calling; silent false negatives/positives are worse than a slow scan).

---

## 1. Formal specification of current behavior (the "oracle")

Before writing any new code, write this specification down as a comment block at the top of the new module and as a set of assertions in a new test file (`tests/test_hdna_dp_equivalence.py`). This is the contract you are optimizing — not your mental model of what the algorithm "should" do.

For a fixed `(ctr, sp)`, define `right0 = ctr + sp + 1`. At trial offset `t = 1, 2, 3, ...` the loop compares:

```
left_char(t)  = dna[ctr - t + 1]
right_char(t) = dna[right0 + t - 1]
```

`k` after `t` iterations equals `t`. The loop's *reachable range* for `k`, ignoring the mismatch-budget check, is bounded by three independent ceilings — compute and unit-test each one separately:

- **Left ceiling:** `K_left = ctr + 1` (array bound only — there is **no N/ambiguity check on the left character at all**. A non-ACGT base on the left side does not stop anything; it simply can never satisfy the equality test against a right character that is guaranteed ACGT by construction — see next point — so it always counts as a mismatch. Do not special-case it; do not add a "left is N → break" rule. That would NOT match current behavior.)
- **Right hard ceiling:** the loop calls `BT[dna[right_char(t)]]` and **breaks immediately, before any comparison or counting**, the first time this is `0` (non-ACGT, including lowercase `n` and any IUPAC ambiguity letter). This position is **not counted** — `k` does not advance to include it. Precompute `next_bad[i]` = smallest `idx >= i` with `BT[dna[idx]] == 0`, sentinel `= total_bases` if none exists. Then `K_right_hard = next_bad[right0] - right0`.
- **`maxrep` ceiling:** `K_maxrep = maxrep`.

`K_total = min(K_left, K_right_hard, K_maxrep)`.

**Mismatch counting:** `mismatches(k) = #{t in [1,k] : left_char(t) != right_char(t)}`. Because the right character is guaranteed ACGT for all `t <= K_total` (that's what `K_right_hard` guarantees), `mismatches(k)` is a simple, well-defined non-decreasing step function on `[1, K_total]`.

**Early-exit semantics — read carefully, this is asymmetric and easy to get wrong:** the check `if (mismatches > mismatch_budget) break;` only executes **inside** the `if (k >= minrep)` block. `mismatch_budget = (int)(mismatch_tol * maxrep)` is computed **once, from `maxrep`, not from the running `k`** — it is an absolute cap, not a per-position ratio. Consequently:
- For `k < minrep`, no early exit check happens at all (other than the three ceilings above).
- `K_mismatch_stop` = smallest `k >= minrep` with `mismatches(k) > mismatch_budget`, if such a `k <= K_total` exists; otherwise `K_mismatch_stop = K_total`.

`K_final = min(K_total, K_mismatch_stop)`.

**Purity test** at a given `k` (only meaningful for `k >= minrep`): let `ga(k)`/`ct(k)` be the count of GA/CT bases in `dna[right0 : right0+k)` (note: **right arm composition only**, never the left arm — this is explicit in the existing `.github/copilot-instructions.md` invariant list, do not violate it). Using the integer-scaled threshold exactly as today:

```
ga_ok(k) = ga(k) * 100 >= purity_thresh_int * k
ct_ok(k) = ct(k) * 100 >= purity_thresh_int * k
purity(k) = ga_ok(k) || ct_ok(k)
ratio_ok(k) = mismatches(k) * 100 <= mismatch_tol_int * k
valid(k) = purity(k) && ratio_ok(k)
```

**`best_k` definition — this is the subtle part:** the legacy loop does **not** stop scanning at the first `k` where `valid(k)` becomes false and resume only at run boundaries — it evaluates `valid(k)` at *every* `k` from `minrep` to `K_final`, and **overwrites** `best_k = k` every time `valid(k)` is true. Because `k` strictly increases through the loop, this means:

> `best_k` = the **largest** `k` in `[minrep, K_final]` such that `valid(k)` is true, or `0` if no such `k` exists.

It is **not** "the first valid `k`", and it is **not** "the longest contiguous valid prefix" — `valid(k)` can be false, then true again at a larger `k`, and `best_k` will track that larger `k`. Any optimization that assumes monotonicity of `valid(k)` or stops at the first failure is **wrong**. Write a brute-force Python or C reference function that reproduces exactly this definition (scanning every `k`) and use it as the ground truth in property tests — do not derive the ground truth from your optimized code.

When `best_k >= minrep`, the hit's `ga_pct`/`ct_pct`/`mirror_id`/`is_perfect` are computed from the **floats at that specific `best_k`**, recomputed via `ga(best_k)`, `ct(best_k)`, `mismatches(best_k)` — not accumulated incrementally. `is_perfect` uses **exact float equality** against `1.0f` (`best_ga == 1.0f || best_ct == 1.0f`) **and** `best_mir == 1.0f` — preserve this exactly, including the fact that it is computed from `(float)x * (1.0f/k)`, not `(float)x/k`, since the reciprocal-multiplication path can produce a different rounding than direct division in rare cases. Use the **exact same arithmetic expression**, not an equivalent-looking one, or `is_perfect` can flip on edge cases.

---

## 2. Performance characterization — set correct expectations before coding

Do not oversell this to the team. Be precise about what changes and what doesn't:

- `mismatch_budget = mismatch_tol * maxrep` scales **linearly with `maxrep`**. The optimization in Section 4.2 turns the per-`(ctr,sp)` matching cost from `O(maxrep)` character comparisons into `O(mismatch_budget / "useful jump length")` LCE queries. In the **worst case** (long exact-match runs — satellite DNA, large `maxrep`, low `mismatch_tol`), a single LCE query collapses an entire run that would otherwise cost up to `maxrep` character comparisons, so the speedup there can be on the order of `maxrep` itself (i.e., dramatic — this is the whole point, since this is exactly the pathological case that currently makes the tool slow on repeat-rich regions).
- On **ordinary, low-repetitivity background sequence**, runs between mismatches are short (a handful of bases for random ACGT), so the new code issues roughly as many LCE queries as the old code does character comparisons. Each LCE query costs more in absolute terms (array-indexed RMQ lookup, less cache-friendly) than a single sequential byte comparison. **Do not assume this case gets faster — benchmark it. It may be a wash or mildly negative**, and that is an acceptable, expected outcome, not a bug, as long as the worst case improves substantially and the feature flag lets users opt out if their workload is dominated by background sequence.
- Section 4.1 (purity memoization) is the separate, lower-risk piece. It removes a genuine `O(maxspacer)` multiplicative factor from the purity-evaluation cost (not from the matching/mismatch cost) by exploiting that `right0 = ctr + sp + 1` is shared by up to `maxspacer + 1` different `(ctr, sp)` pairs, and purity at a given `right0` is independent of which pair produced it. Implement and benchmark this **independently** of Section 4.2 — it has a much smaller blast radius and should land first.

Report both regimes (synthetic random ACGT sequence, and synthetic/real tandem-repeat-rich sequence such as a microsatellite or ALU-dense region) in every benchmark you run. A PR that only reports one regime is incomplete.

---

## 3. New file layout

Do **not** grow `_hdna.c` into an unreadable single file. Add:

```
src/hseeker/
├── _hdna.c                 ← existing; only the call site inside findHDNA_core
│                              changes (gated by feature flag), plus next_bad[]
│                              and prefix-sum array construction calls.
├── _hdna_dp.h               ← new; struct/function declarations shared below
├── _hdna_purity.c           ← new; Optimization A — prefix sums + per-right0
│                              memoized "largest valid-purity k at or before B"
│                              table, with an LRU/ring cache (see 4.1.3).
├── _hdna_sa.c                ← new; Optimization B — suffix array + Kasai LCP
│                              + sparse-table RMQ, scoped to one chunk buffer.
└── _hdna_lce.c                ← new; the run-boundary enumerator that walks
                                  mismatch runs for a given (ctr, sp) using
                                  _hdna_sa.c's O(1) LCE query, applies the
                                  closed-form ratio threshold, and looks up
                                  _hdna_purity.c's cache. This is the function
                                  findHDNA_core calls instead of the inner
                                  while-loop when the feature flag is on.
```

Update `setup.py`'s `sources=[...]` list to include the new `.c` files, and add the new `.c`/`.h` files to `MANIFEST.in` / `pyproject.toml`'s `package-data` so they ship in the sdist (check both — the existing config only lists `_hdna.c` explicitly in `[tool.setuptools.package-data]`; this must be extended or the sdist will be broken for source builds).

---

## 4. Design

### 4.1 Optimization A — purity memoization (implement and merge first)

**4.1.1 Prefix sums.** Build once, per chunk, two `int32_t[total_bases+1]` arrays `prefix_ga`, `prefix_ct` (or pack into one `int32_t` per base with two 8-bit-ish counters if you want to halve memory — measure first, only micro-optimize if profiling says so). `ga(right0, k) = prefix_ga[right0+k] - prefix_ga[right0]`, same for `ct`. This replaces the incremental `ga_count++`/`ct_count++` in the legacy loop. Memory: `2 * 4 * total_bases` bytes ≈ `8N` bytes — trivial, always allocate this regardless of which other optimization is enabled.

**4.1.2 Per-`right0` purity table.** For a given `right0`, build `int32_t best_purity_idx[maxrep+1]` where `best_purity_idx[B]` = largest `k <= B` with `purity(k)` true, or `0`/sentinel if none. Build it with a single backward pass over `k = K (=min(maxrep, total_bases-right0)) down to 1`, carrying forward the running max. Cost: `O(maxrep)` time, `O(maxrep)` memory, **per distinct `right0` actually visited** (don't precompute for all `right0` up front — that's `O(N * maxrep)` memory, infeasible).

**4.1.3 Cache, do not recompute, across `sp`.** Note that `right0 = ctr + sp + 1`. If you keep the **existing loop order** (`for ctr { for sp { ... } }`), a given `right0` is visited for `sp = 0..maxspacer` as `ctr` decreases — i.e., `right0` values are **not** contiguous in a simple way across the `ctr` loop in the existing nesting. You have two options, evaluate both and pick based on measured cache hit rate:
   - **(a) Keep loop order, add a small LRU cache** (size ~ `2 * maxspacer`, since within a sliding window of `ctr` values, the same `right0` can recur) keyed by `right0`, storing the `best_purity_idx` table. Cheaper to implement, may have a lower hit rate.
   - **(b) Restructure the outer loop to iterate over `right0` directly** (`for right0 in [minrep, total_bases-minrep]`, then derive the `(ctr, sp)` pairs that map to it: `ctr = right0 - sp - 1` for `sp = 0..min(maxspacer, right0-minrep-1)`), build the purity table once per `right0`, use it for all `(ctr,sp)` pairs that map to it, then discard. This guarantees a 100% reuse rate (each table built exactly once) at the cost of restructuring `findHDNA_core`'s loop nest — **this changes iteration order, which changes hit *discovery* order but must not change the final hit *set/values* after the existing overlap-removal sort.** Verify this explicitly: `remove_overlaps_core` sorts by `start` before sweeping, so output order is independent of discovery order — confirm this with a test, don't just assume it.
   
   Recommendation: prototype (a) first since it's a smaller diff; only do (b) if profiling shows the LRU hit rate is poor.

**4.1.4 Wire it in.** Inside the legacy-loop replacement, when you need "does any `k` in `[lo, hi]` satisfy `purity(k)`, and what's the largest one", look up `best_purity_idx[hi]` (clip `hi` to the table's built range) and check `>= lo`. This is the O(1) replacement for re-scanning purity at every `k`.

### 4.2 Optimization B — per-chunk suffix array + LCP + sparse-table RMQ for O(1) LCE

**4.2.1 Construction, scoped to one chunk.** For the chunk buffer `dna[0..C)` (the same buffer/size already produced by `scan_fasta_parallel`'s `chunk_size` + `overlap` logic in `__init__.py` — **do not introduce a second, different chunking scheme**; reuse the existing one so the memory bound below is the one the rest of the codebase already assumes), build the generalized string:

```
S = dna[0..C) + SEP + reverse(dna[0..C))
```

`SEP` must be a byte value that cannot appear in `dna` after lowercasing (the input is restricted to `ACGTNacgtn` plus possibly other IUPAC letters via `tolower`; pick a byte outside the lowercase-alpha range, e.g. `0x00` or `0x01`, and assert this at construction time rather than assuming it). Build:
   - Suffix array `SA` of `S` (length `M = 2C+1`). Use a standard `O(M log M)` doubling algorithm or pull in a small, dependency-free `O(M)` SA-IS implementation — your call, but whichever you choose must have **no new third-party dependency** added to `pyproject.toml`'s core `dependencies` (the project's stated invariant is the C extension has no runtime deps; check `_scoring.py`'s "zero external dependencies" constraint and apply the same discipline here).
   - Rank array (inverse of `SA`).
   - LCP array via Kasai's algorithm, `O(M)`.
   - Sparse table over the LCP array for `O(1)` range-minimum queries: `sparse[j][i] = min(sparse[j-1][i], sparse[j-1][i + 2^(j-1)])`, `j = 0..floor(log2(M))`. **This is the dominant memory cost — see Section 5 for the exact budget and a hard cap.**

**4.2.2 LCE query.** `LCE(a, b)` for two starting positions in `S` = RMQ over `LCP[rank[a]+1 .. rank[b]]` (standard formula; if `rank[a] == rank[b]` they're the same suffix, handle as a degenerate case — should not occur here since `a` and `b` always come from different halves of `S`, but assert it rather than silently mishandling it).

To get the legacy loop's match-run length starting at trial offset `t0+1` (i.e., how many further consecutive *exact* matches occur starting right after a known mismatch at `t0`), map `left_char(t0+1) = dna[ctr-t0]` and `right_char(t0+1) = dna[right0+t0]` to their positions in `S`: the right side maps directly (`b = right0+t0`), the left side maps into the **reversed half** of `S` (`a = C + 1 + (C - 1 - (ctr - t0))`, i.e., the position in `reverse(dna)` corresponding to `dna[ctr-t0]` — **derive and unit-test this index arithmetic in isolation before wiring it into anything else; an off-by-one here silently corrupts every result and is the single most likely bug in this whole project**). `LCE(a, b)` then gives the number of further matching characters, **before applying the hard right ceiling from `next_bad[]` — you must clip the LCE result by `next_bad[right0+t0] - (right0+t0)` regardless of what the LCE query returns**, because `S`'s suffix-array machinery has no concept of "this byte is N and must hard-stop" — it will happily report exact matches across two `n` bytes (or two identical IUPAC ambiguity bytes) that the legacy loop would never have reached (since it breaks the instant it sees `BT==0` on the right, before any comparison). **This clipping step is not optional — write a specific unit test with two `N`s at the same relative offset on both arms and confirm the new code produces zero match length there, not whatever the raw LCP value says.**

**4.2.3 Run enumeration driver (`_hdna_lce.c`).** For a given `(ctr, sp)`:
1. Compute `K_total` (Section 1) using `next_bad[]` and the array-bound ceilings.
2. Starting at `t0 = 0`, repeatedly: query `LCE` from the current position, clip by the right-hard-ceiling, clip by `K_total - t0` (don't overrun); this gives the length of the current exact-match run. Record the run boundary `[run_start, run_end]`. If `run_end == K_total`, this is the last run — stop. Otherwise, the character immediately after the run is, by definition of LCE, a mismatch (or you've hit `K_total` exactly) — increment a mismatch counter `i`, advance `t0 = run_end + 1`, and continue.
3. Stop enumerating runs as soon as `mismatches(k) > mismatch_budget` becomes true for `k >= minrep` (matches the legacy early-exit) **or** `i` (the running mismatch count) exceeds `mismatch_budget` regardless of `k >= minrep`, mirroring the exact nested-`if` ordering in Section 1 — re-read that section, the check is conditional on `k >= minrep`, do not apply it unconditionally or you will diverge from the oracle for short, dense-mismatch regions below `minrep`.
4. For each enumerated run (most recent first, since `best_k` wants the largest valid `k` — see 4.1.4 on why scanning from the end first lets you stop at the first hit), compute `lo' = max(run_start, minrep, ceil(100*i / mismatch_tol_int))` (closed-form solution of `ratio_ok(k)` being constant-`i`-mismatch and monotonically easier as `k` grows within a run — derive and **prove this formula in a code comment with the algebra**, don't just assert it) and `hi = run_end`. Look up `best_purity_idx_for_right0[hi]` (Optimization A) and check `>= lo'`. If yes, that value is `best_k` for this `(ctr, sp)` — stop. If no run yields a hit, `best_k = 0`.
5. **Special case `mismatch_tol_int == 0`:** the closed-form `lo'` formula divides by `mismatch_tol_int`. When mismatch tolerance is `0.0` (exact-mirror mode, explicitly supported per the README: `mismatch=0.0` for exact mirrors only), `ratio_ok(k)` requires `i == 0` exactly — i.e., **only the first run (i=0) can ever be valid**; do not enter the loop with a division by zero, special-case `mismatch_tol_int == 0` to only ever consider the first run and return immediately after.

### 4.3 What `findHDNA_core` looks like after this lands

The inner `while` loop's body is replaced by a single call to the Section 4.2.3 driver when the feature flag is on, falling back to the existing character-by-character loop verbatim when it's off (keep the old code path **in the same file, unconditionally compiled**, not deleted — this is your regression safety net and your fallback for OOM/construction-failure cases in Section 5).

---

## 5. Memory budget — hard requirement, not a suggestion

A **global** suffix array + sparse table over an entire chromosome is infeasible: for `N = 100` Mbp, `M ≈ 2×10^8`, `log2(M) ≈ 28`, sparse table alone is `M × 28 × 4` bytes ≈ **22 GB**. For the STANDALONE binary's documented ceiling (`MAX_DNA_STANDALONE = 600,000,000`), this would be well over 100 GB. **Do not build this globally.**

Instead, scope construction to the **existing chunk** (`chunk_size`, default `1,000,000`, already used by `scan_fasta_parallel` in `__init__.py`). For `C = 1,000,000`: `M ≈ 2.001×10^6`, `log2(M) ≈ 21`, sparse table ≈ `2×10^6 × 21 × 4` bytes ≈ **168 MB** per chunk. With `n_workers` chunks in flight concurrently (the existing `ThreadPoolExecutor`), peak RAM for this structure alone is `n_workers × ~170 MB`. **Compute and log this number at runtime** (e.g., a debug log line: `"DP scan: chunk=%d bases, SA+RMQ memory=%zu MB"`), and:

- Add a **hard cap**: if the projected sparse-table memory for a given chunk size exceeds a configurable threshold (default suggestion: 512 MB per chunk — make this a named constant, not a magic number), **do not build the structure for that chunk; fall back to the legacy loop for that chunk only.** This must be a per-chunk decision, not a global abort — a few oversized chunks falling back should not crash or disable the optimization for the rest of the run.
- This cap interacts with `chunk_size`: if a user calls `scan_fasta_parallel(chunk_size=50_000_000, ...)`, the sparse table would be ~8.4 GB for that chunk and should trip the fallback. Write a test that constructs a chunk sized to deliberately exceed the cap and asserts the fallback path is taken (check via the profiling counters in Section 9.3, not by parsing log strings).
- Also account for: `SA` (`M × 4` bytes, or `8` bytes if `M > 2^31`, which it will not be for any reasonable chunk size — assert `M < INT32_MAX` and fail loudly rather than silently truncating if someone sets an enormous `chunk_size`), `rank` array (same size), `LCP` array (same size), and the prefix-sum arrays from 4.1.1 (`8 × C` bytes — negligible by comparison). Sum all of these in the logged number, not just the sparse table.
- Document the final memory formula in `.github/copilot-instructions.md` once this lands (Section 9.5 covers this).

---

## 6. Phased implementation plan

Each phase is a separate, independently mergeable, independently revertible unit of work. Do not start phase *n+1* until phase *n*'s tests pass and its benchmark numbers (Section 2) are recorded.

**Phase 0 — Baseline harness (no algorithm changes).**
Generate golden-output fixtures: run the current (unmodified) `scan_sequence`/`scan_fasta` on (a) every sequence already used in `tests/test_hdna.py`, (b) a battery of randomly generated sequences spanning the edge cases in Section 7, (c) at least one real repeat-rich region (e.g., a known microsatellite locus or the validation FASTA already in `sensitivity_analysis/hseeker_validation_motifs.fasta`), at multiple `(minrep, maxrep, maxspacer, purity, mismatch)` parameter combinations including the boundary values `mismatch=0.0` and `purity=1.0`. Serialize full hit lists (all fields) as JSON fixtures checked into `tests/fixtures/`. Write `tests/test_hdna_dp_equivalence.py` that will, in later phases, re-run the same inputs through the new code path and diff against these fixtures field-by-field.

**Phase 1 — Optimization A only (Section 4.1).** Implement prefix sums + per-`right0` purity memoization. Wire it in so it changes *only* how purity is evaluated — the mismatch-counting and run-discovery logic stays the legacy char-by-char loop. This phase has a much smaller surface area than Optimization B and should be fully landed, tested, and benchmarked before touching suffix arrays. Run the Phase 0 fixtures through it; any diff is a stop-the-line bug.

**Phase 2 — Suffix array / LCP / sparse table as a standalone, generically-testable component.** Implement `_hdna_sa.c` with no dependency on `findHDNA_core` at all — just `build_sa_lcp_rmq(const char *s, int m) -> handle` and `lce_query(handle, int a, int b) -> int`. Unit-test this in isolation against a brute-force `O(n)` LCE reference (direct character comparison) on hundreds of random strings of varying length and alphabet (including strings with repeated substrings, palindromes, and all-same-character strings — the degenerate case that stresses LCP-array correctness most). Do not let any HSeeker-specific concept (arms, spacers, purity) leak into this file — it should be reusable, generic string-algorithm code.

**Phase 3 — Run-boundary enumerator (`_hdna_lce.c`), tested against the Phase 0 oracle in isolation.** Implement the Section 4.2.3 driver. Before wiring it into `findHDNA_core`, write a standalone test harness that, for thousands of random small sequences (length 50–500, so brute force is cheap) and random `(ctr, sp, minrep, maxrep, mismatch_tol, purity_thresh)` combinations, computes `best_k` via (a) the brute-force oracle from Section 1 and (b) the new enumerator, and asserts equality. This is where the off-by-one risk flagged in 4.2.2 will surface — **do not proceed to Phase 4 until this fuzz test has run at least 100,000 random cases with zero mismatches.**

**Phase 4 — Integration into `findHDNA_core` behind the feature flag.** Wire Phases 1–3 together, gated by the flag from Constraint 4. Re-run all Phase 0 fixtures through `py_scan_sequence` with the flag on; diff every field. Re-run the full existing `tests/test_hdna.py` suite with the flag on (temporarily monkey-patch the env var in a fixture, or parametrize the test module to run twice — once per flag state — so both code paths get full coverage from the existing suite, not just the new equivalence test).

**Phase 5 — Chunk-scoped construction wired into `scan_fasta_parallel`'s existing chunking, with the memory cap and fallback from Section 5.** Test specifically: (a) a chunk that fits comfortably under the cap, (b) a chunk deliberately sized to exceed it (assert fallback engaged and output is still correct), (c) the existing chunk-boundary/overlap correctness tests in `tests/test_parallel_edges.py` (note this file existed in a prior repo snapshot and may need to be re-added/ported if missing in the current `tests/` — check) still pass unchanged, since the chunking/overlap math itself is not something this work should touch.

**Phase 6 — Benchmarking and documentation.** Extend `benchmarks/benchmark.py` (don't fork a new script) with a flag-state comparison mode: run the same chromosome/region through both code paths, report wall time and peak RSS for both, on at least one low-repetitivity chromosome and one region known to be repeat-dense. Update `.github/copilot-instructions.md` with: the new file layout, the new feature flag, the memory formula and cap, and the calibrated performance expectations from Section 2 (do not let documentation overpromise a uniform speedup — state the repeat-rich-vs-background distinction explicitly, as future contributors will otherwise "fix" the background-case overhead by reintroducing correctness bugs).

---

## 7. Edge-case checklist — every item needs an explicit test

- [ ] Sequence shorter than `2*minrep` (no valid center exists at all).
- [ ] Sequence consisting entirely of `N`/non-ACGT — `next_bad[]` degenerates correctly, no hits, no crash.
- [ ] A single `N` immediately at `right0` (zero-length right arm, `K_right_hard = 0`).
- [ ] `N` appearing only on the **left** side, never the right — must still produce a hit if the right side alone would otherwise qualify and the left-side `N` doesn't push mismatches over budget; must **not** be treated as a hard stop.
- [ ] `ctr` at position `0` and at `total_bases - 1` (left/right array-bound ceilings binding exactly).
- [ ] `maxspacer = 0`.
- [ ] `mismatch = 0.0` (exact mirror mode) — confirm the `mismatch_tol_int == 0` special case (4.2.3 step 5) and that no run after the first is ever considered.
- [ ] `purity = 1.0` combined with `mismatch = 0.0` (strictest legacy-equivalent mode, matches the README's "Set to 1.0 to reproduce strict non-B_gfa mirror-repeat results").
- [ ] `purity = 0.0` (always satisfied — confirm the memoized table still degenerates correctly rather than dividing by zero or similar).
- [ ] An arm that is exactly `minrep` long and exactly satisfies thresholds at that boundary (off-by-one in the `lo'` formula would show up here first).
- [ ] An arm where `valid(k)` is true, then false, then true again at a larger `k` within `K_final` — confirm `best_k` picks the **larger**, later one (this directly tests the non-monotonicity warning in Section 1).
- [ ] Lowercase, uppercase, and mixed-case input (the existing code lowercases before calling the core; confirm the new arrays are built from the lowercased buffer, not the raw input).
- [ ] IUPAC ambiguity codes other than `N` (e.g., `R`, `Y`, `W`, `S`) on both left and right sides — same `BT == 0` treatment as `N`, confirm no special-casing was accidentally added that distinguishes `N` from other ambiguity letters (the legacy `BT[256]` table treats them identically; your new code must too).
- [ ] `total_bases == maxrep` and `total_bases` slightly less than `maxrep` (boundary interaction between the sequence-end ceiling and `maxrep`).
- [ ] A run that ends exactly at `K_total` with zero mismatches (the "last run, no trailing mismatch" branch in 4.2.3 step 2).
- [ ] Two adjacent `(ctr, sp)` pairs sharing the same `right0` where one is valid and the other is not, to confirm the memoized purity table is correctly reused without leaking the wrong `lo'`/mismatch context between them (purity table is shared; the run/mismatch context is not — make sure your code never accidentally caches the latter).
- [ ] Chunk-boundary cases once Phase 5 lands: a hit whose left arm sits in the overlap region from the previous chunk and whose right arm sits in this chunk's owned region (this is exactly what the existing `overlap = 2*maxrep + maxspacer + 1` formula in `__init__.py` is designed to cover — confirm the new per-chunk SA/LCP/RMQ scope doesn't accidentally exclude positions the legacy per-chunk scan would have included).

---

## 8. Testing requirements summary

- `tests/test_hdna_dp_equivalence.py`: fixture-replay equivalence test (Phase 0/4) — exact field-by-field diff, both flag states.
- A standalone, HSeeker-agnostic LCE correctness fuzz test for `_hdna_sa.c` (Phase 2) — brute force vs. SA/LCP/RMQ, hundreds of random strings, must include degenerate all-same-character and highly repetitive inputs.
- A standalone `best_k` fuzz test for `_hdna_lce.c` (Phase 3) — brute-force oracle (Section 1, written independently from the optimized code) vs. enumerator, ≥100,000 random `(seq, ctr, sp, params)` tuples, zero mismatches required before integration.
- Every item in Section 7 as a named, individually-readable test case (not folded into one giant parametrized blob with unclear failure messages — if a future change regresses "N appearing only on the left side", the test name should say exactly that).
- Re-run of the full existing `tests/test_hdna.py` suite under both flag states (Phase 4).
- Memory-cap fallback test (Phase 5).
- Updated `benchmarks/benchmark.py` run, both regimes, both flag states (Phase 6) — numbers go in the PR description, not just in a local notebook.

---

## 9. Engineering process requirements

1. **Compiler flags unchanged** (`-O3 -march=native -Wall`, per `setup.py`) — do not add `-ffast-math` or similar to "help" the new code; the existing `is_perfect` exact-float-equality contract is fragile enough already.
2. **Error-path discipline.** Every new `malloc` (SA, rank, LCP, sparse table, prefix sums, purity tables) needs a checked allocation with a `cleanup:`-style goto/early-return on failure that frees everything allocated so far, matching the existing convention in `py_scan_sequence`. An OOM building the chunk-local SA/RMQ structure is **not** a fatal error for the whole scan — per Section 5, it should trigger the legacy-loop fallback for that chunk, not propagate as a Python exception that aborts `scan_fasta_parallel` for the entire genome.
3. **Extend the existing profiling counters** (`prof_inner_iters`, `prof_ctr_sp_pairs`) rather than inventing a parallel instrumentation scheme. Add: `prof_lce_queries`, `prof_purity_cache_hits`, `prof_purity_cache_misses`, `prof_chunk_fallback_count` (incremented when Section 5's memory cap trips). Expose these the same way the existing counters are exposed to Python (check how `prof_inner_iters` currently surfaces, if at all, and follow that pattern — if it doesn't currently surface to Python, add a minimal `_hdna.get_profiling_counters()` rather than printing to stderr).
4. **Do not rename any of the canonical parameter names** (`minrep`, `maxrep`, `maxspacer`, `purity`, `mismatch`, `score`) — this is an explicit standing invariant in `.github/copilot-instructions.md`.
5. **Documentation.** Once Phase 6 lands, update `.github/copilot-instructions.md`'s "Performance-critical optimizations in `findHDNA_core`" section to list the new components, the feature flag, and the memory-budget formula, so future Copilot sessions inherit this context instead of rediscovering it.
6. **Rollback plan.** Because everything is behind a flag and the legacy path is never deleted, rollback is "flip the flag back" — confirm this works with a one-line test that sets the flag off and re-runs the full suite, as a literal CI check, not just an assumption.

---

## 10. Definition of done

- All six phases merged, each with passing tests and recorded benchmark numbers for both the background-sequence and repeat-rich regimes.
- Zero diffs against Phase 0 golden fixtures under either flag state.
- Memory for the new structures, at the default `chunk_size=1_000_000`, measured and confirmed under the cap from Section 5, with the cap-and-fallback path exercised by a real test (not just code-read).
- `.github/copilot-instructions.md` updated.
- A short "When this helps and when it doesn't" paragraph added to the README's existing performance/limitations material, calibrated per Section 2 — do not let the README claim a universal speedup.
