/*
 * _hdna_dp.h — shared declarations for the HSeeker DP/LCE optimization.
 *
 * This header is included by _hdna_purity.c, _hdna_sa.c, _hdna_lce.c,
 * and the call-site in _hdna.c.  It does NOT include Python.h so that
 * the SA/LCE modules remain independently testable without a Python build.
 *
 * Feature flag: HSEEKER_FAST_HDNA environment variable.
 *   Unset or "0"  → legacy character-by-character loop (unchanged).
 *   "1"           → new DP path (prefix sums + SA/LCE run-boundary enumerator).
 *
 * Memory formula (per chunk of C bases, log2(2C) levels):
 *   SA:           (2C+1) * 4 bytes
 *   Rank:         (2C+1) * 4 bytes
 *   LCP:          (2C+1) * 4 bytes
 *   Sparse table: (2C+1) * ceil(log2(2C+1)) * 4 bytes  ≈ 21*(2C+1)*4 for C=1e6
 *   next_bad[]:   C * 4 bytes
 *   prefix_ga[]:  (C+1) * 4 bytes
 *   prefix_ct[]:  (C+1) * 4 bytes
 *   Total (C=1e6): ~175 MB per worker thread.
 *
 * Hard cap: DP_MAX_CHUNK_BYTES (default 512 MB).  If a chunk would exceed
 * this, findHDNA_core falls back to the legacy loop for that chunk only.
 */

#ifndef HSEEKER_HDNA_DP_H
#define HSEEKER_HDNA_DP_H

#include <stddef.h>
#include <stdint.h>

/* ── Feature flag ────────────────────────────────────────────────────────── */

/* Returns 1 if the HSEEKER_FAST_HDNA env var is set to "1". */
int hdna_dp_enabled(void);

/* ── Memory cap ──────────────────────────────────────────────────────────── */

/* Hard cap in bytes for the combined SA+RMQ allocation per chunk.
 * If the projected size exceeds this, fall back to the legacy loop. */
#define DP_MAX_CHUNK_BYTES  ((size_t)(512UL * 1024UL * 1024UL))

/* ── Profiling counters (exposed via py_profiling_info) ──────────────────── */
extern long long prof_lce_queries;
extern long long prof_purity_cache_hits;
extern long long prof_purity_cache_misses;
extern long long prof_chunk_fallback_count;

void prof_dp_reset(void);

/* ── Prefix-sum purity arrays (Optimization A) ───────────────────────────── */

/*
 * PurityArrays: prefix sums over GA and CT base counts in dna[0..n).
 * prefix_ga[i] = number of GA bases in dna[0..i-1).
 * prefix_ct[i] = number of CT bases in dna[0..i-1).
 * Length of each array is n+1.
 */
typedef struct {
    int32_t *prefix_ga;   /* length n+1 */
    int32_t *prefix_ct;   /* length n+1 */
    int      n;           /* length of original sequence */
} PurityArrays;

/* Build prefix-sum arrays for dna[0..n).  Returns NULL on OOM. */
PurityArrays *purity_arrays_build(const char *dna, int n);

/* Free a PurityArrays object. */
void purity_arrays_free(PurityArrays *pa);

/* ga(right0, k) = number of GA bases in dna[right0..right0+k). */
static inline int purity_ga(const PurityArrays *pa, int right0, int k)
{
    return pa->prefix_ga[right0 + k] - pa->prefix_ga[right0];
}

/* ct(right0, k) = number of CT bases in dna[right0..right0+k). */
static inline int purity_ct(const PurityArrays *pa, int right0, int k)
{
    return pa->prefix_ct[right0 + k] - pa->prefix_ct[right0];
}

/*
 * PerRight0Table: for a fixed right0 and K = min(maxrep, n-right0),
 * best_purity_k[b] = largest k in [1,b] such that purity(k) is true, or 0.
 * Length: K+1.
 */
typedef struct {
    int32_t *best_k;  /* length K+1; best_k[b] = largest valid-purity k <= b */
    int       K;      /* max index (= min(maxrep, n-right0)) */
    int       right0; /* the right0 this table was built for */
} PerRight0Table;

/* Build the per-right0 purity table using already-built PurityArrays. */
PerRight0Table *per_right0_build(
    const PurityArrays *pa,
    int right0, int maxrep,
    int purity_thresh_int);

/* Free a PerRight0Table. */
void per_right0_free(PerRight0Table *t);

/*
 * LRU cache of PerRight0Table, keyed by right0.
 * Capacity = 2 * maxspacer + 4 (small constant; fits in L1 cache).
 */
typedef struct {
    PerRight0Table **entries; /* array of pointers, length = capacity */
    int              capacity;
    int              count;
    /* round-robin eviction — simple and effective for the access pattern */
    int              next_evict;
} PurityCache;

PurityCache *purity_cache_create(int capacity);
void         purity_cache_free(PurityCache *c);

/*
 * Lookup or build the PerRight0Table for right0.
 * Returns a borrowed pointer (cache owns it); do not free.
 */
const PerRight0Table *purity_cache_get(
    PurityCache *cache,
    const PurityArrays *pa,
    int right0, int maxrep,
    int purity_thresh_int);

/* ── next_bad array ──────────────────────────────────────────────────────── */

/*
 * Build next_bad[i] = smallest idx >= i with BT[dna[idx]]==0, or n if none.
 * Length n+1 (next_bad[n] = n as sentinel).
 * Caller must free the returned pointer.
 */
int32_t *next_bad_build(const char *dna, int n, const unsigned char *BT);

/* ── Suffix array / LCP / sparse-table RMQ (Optimization B) ─────────────── */

/*
 * SaRmqHandle: opaque struct holding the suffix array, rank, LCP array,
 * and sparse table for a combined string S = dna + SEP + rev(dna).
 *
 * All internal arrays are length M = 2*C+1 where C = chunk size.
 */
typedef struct SaRmqHandle SaRmqHandle;

/*
 * Build the SA/LCP/RMQ structure for the chunk buffer dna[0..C).
 * Returns NULL on OOM or if the projected memory exceeds DP_MAX_CHUNK_BYTES
 * (in which case prof_chunk_fallback_count is incremented).
 *
 * On success, *out_fallback = 0.
 * On memory-cap fallback, *out_fallback = 1 and return value is NULL.
 * On true OOM, *out_fallback = 0 and return value is NULL.
 */
SaRmqHandle *sa_rmq_build(const char *dna, int C, int *out_fallback);

/* Free the handle and all internal arrays. */
void sa_rmq_free(SaRmqHandle *h);

/*
 * LCE query: longest common extension of S starting at positions a and b.
 * Returns the number of characters that match before the first mismatch.
 *
 * a and b are positions in the combined string S (not just dna).
 * Positions in the forward half: 0..C-1  → dna[0..C-1].
 * Positions in the reverse half: C+1..2C → rev(dna)[0..C-1] = dna[C-1-i].
 */
int sa_lce(const SaRmqHandle *h, int a, int b);

/* C = the chunk size (length of dna half). */
int sa_chunk_size(const SaRmqHandle *h);

/* ── Run-boundary enumerator (Optimization B driver) ─────────────────────── */

/*
 * dp_best_k: compute best_k for a single (ctr, sp) pair using the
 * run-boundary enumerator (SA/LCE + purity cache).
 *
 * Semantics: exactly matches findHDNA_core's inner while loop for that
 * (ctr, sp) pair — same early-exit conditions, same mismatch-budget check,
 * same best_k definition (largest valid k).
 *
 * Parameters:
 *   dna            : lowercase chunk buffer
 *   total_bases    : C (length of chunk)
 *   ctr, sp        : center and spacer for this iteration
 *   minrep..mismatch_tol_int : scan parameters (already int-scaled where noted)
 *   next_bad       : precomputed next_bad[] array
 *   pa             : precomputed prefix-sum arrays
 *   cache          : per-right0 purity LRU cache (shared across (ctr,sp))
 *   sa             : SA/RMQ handle for LCE queries
 *
 * Returns best_k >= minrep if a hit exists, else 0.
 */
int dp_best_k(
    const char *dna, int total_bases,
    int ctr, int sp,
    int minrep, int maxrep,
    int purity_thresh_int, int mismatch_tol_int,
    int mismatch_budget,
    const int32_t *next_bad,
    const PurityArrays *pa,
    PurityCache *cache,
    const SaRmqHandle *sa);

#endif /* HSEEKER_HDNA_DP_H */
