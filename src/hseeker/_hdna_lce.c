/*
 * _hdna_lce.c — Run-boundary enumerator (Section 4.2.3).
 *
 * Implements dp_best_k(): given a (ctr, sp) pair and the prebuilt SA/LCE +
 * purity-cache structures, enumerate mismatch runs and find the largest valid
 * arm length — exactly matching the semantics of findHDNA_core's inner loop.
 *
 * Design principles (from copilot_dp_optimization_instructions.md):
 *
 * 1. K_total = min(K_left, K_right_hard, maxrep)
 *    K_left       = ctr + 1
 *    K_right_hard = next_bad[right0] - right0   (right hard ceiling)
 *
 * 2. Mismatch budget = int(mismatch_tol * maxrep)  (from maxrep, not running k)
 *
 * 3. Early-exit check only inside k >= minrep block.
 *
 * 4. best_k = LARGEST k in [minrep, K_final] where valid(k) is true.
 *    Scanning runs back-to-front lets us stop at the first hit.
 *
 * 5. Special case: mismatch_tol_int == 0 → only the first run (i=0) can be valid.
 *
 * 6. The closed-form lower bound for the ratio check within a run with i mismatches:
 *      ratio_ok(k) ⟺ i * 100 <= mismatch_tol_int * k
 *                  ⟺ k >= ceil(100 * i / mismatch_tol_int)
 *    So for a run at [run_start, run_end] with i accumulated mismatches before it:
 *      lo' = max(run_start, minrep, ceil(100*i / mismatch_tol_int))
 *      hi  = run_end
 *      If best_purity_k[hi] >= lo' then best_k = best_purity_k[hi].
 *
 *    Proof that ratio_ok is monotonically easier as k grows within a run:
 *      Within a run, mismatches = i (constant, no new mismatches).
 *      ratio_ok(k) = i*100 <= mismatch_tol_int * k.
 *      As k increases, RHS grows, so once true it stays true. QED.
 *    Therefore the largest valid-purity k in [lo', hi] is best_purity_k[hi]
 *    when best_purity_k[hi] >= lo'.
 *
 * 7. LCE position mapping:
 *    The combined string S = dna[0..C) + SEP + rev(dna)[0..C).
 *    At step t (0-indexed, so left_i = ctr-t, right_j = right0+t):
 *      Forward position (right arm) : right0 + t       (= position in S directly)
 *      Reverse position (left arm)  : C + 1 + (C-1-(ctr-t)) = 2C - ctr + t
 *    LCE(2C-ctr+t, right0+t) gives the length of the match run starting
 *    at the CURRENT step t.  We must then clip by:
 *      a) next_bad[right0+t] - (right0+t)  (right hard ceiling from current pos)
 *      b) K_total - t                       (don't overrun K_total)
 *    The result is the run length R; the NEXT position is t + R (the mismatch).
 */

#include <stdlib.h>
#include "_hdna_dp.h"

int dp_best_k(
    const char *dna, int total_bases,
    int ctr, int sp,
    int minrep, int maxrep,
    int purity_thresh_int, int mismatch_tol_int,
    int mismatch_budget,
    const int32_t *next_bad,
    const PurityArrays *pa,
    PurityCache *cache,
    const SaRmqHandle *sa)
{
    int right0  = ctr + sp + 1;
    int C       = sa_chunk_size(sa);

    /* K_total: binding ceiling */
    int K_left       = ctr + 1;
    int K_right_hard = next_bad[right0] - right0;
    int K_total      = K_left;
    if (K_right_hard < K_total) K_total = K_right_hard;
    if (maxrep       < K_total) K_total = maxrep;

    if (K_total < minrep) return 0;

    /* Fetch purity cache entry for this right0 */
    const PerRight0Table *pt = purity_cache_get(
        cache, pa, right0, maxrep, purity_thresh_int);
    if (!pt || pt->K == 0) return 0;  /* OOM or empty right arm */

    /* Clip K_total to the purity table's range */
    if (K_total > pt->K) K_total = pt->K;
    if (K_total < minrep) return 0;

    /* Special case: mismatch_tol_int == 0 — exact mirror only.
     * Only the first run (i=0) is ever valid; it starts at t=0.
     * LCE gives the run length; best_k = min(run_len, K_total) if >= minrep
     * and purity passes. */
    if (mismatch_tol_int == 0) {
        int a_pos = 2 * C - ctr;        /* left arm start in reverse half */
        int b_pos = right0;             /* right arm start in forward half */
        prof_lce_queries++;
        int run_len = sa_lce(sa, a_pos, b_pos);
        /* Clip by right hard ceiling from right0 */
        int right_avail = next_bad[right0] - right0;
        if (run_len > right_avail) run_len = right_avail;
        if (run_len > K_total)     run_len = K_total;
        if (run_len < minrep)      return 0;
        /* Check purity at run_len */
        int32_t best_purity = pt->best_k[run_len];
        return (best_purity >= minrep) ? best_purity : 0;
    }

    /* General case: enumerate runs front-to-back, collect as array,
     * then scan back-to-front to find largest valid best_k quickly. */

    /* Maximum possible runs = (budget+1 match runs) + (budget mismatch-step entries).
     * Each mismatch step is recorded as a single-element entry because the legacy
     * inner loop also records k = (mismatch_position + 1) as a valid arm length.
     * Total entries ≤ 2*budget + 1; allocate 2*budget + 2 for safety. */
    int max_runs = 2 * mismatch_budget + 2;
    /* Allocate on stack for small budgets, heap for large. */
    typedef struct { int run_start; int run_end; int mismatches_before; } RunInfo;

    RunInfo stack_runs[256];
    RunInfo *runs = stack_runs;
    RunInfo *heap_runs = NULL;
    if (max_runs > 256) {
        heap_runs = (RunInfo *)malloc(sizeof(RunInfo) * (size_t)(max_runs));
        if (!heap_runs) return 0;  /* OOM → fall back returns 0; caller uses legacy */
        runs = heap_runs;
    }

    int n_runs   = 0;
    int t        = 0;  /* current trial offset (0-indexed, so k=t at start of step) */
    int i        = 0;  /* accumulated mismatch count */

    while (t < K_total) {
        /* Map positions to combined string S */
        int a_pos = 2 * C - ctr + t;   /* left arm position at step t in rev-half */
        int b_pos = right0 + t;         /* right arm position at step t in fwd-half */

        prof_lce_queries++;
        int raw_lce = sa_lce(sa, a_pos, b_pos);

        /* Clip by right hard ceiling from current right position */
        int right_avail = next_bad[right0 + t] - (right0 + t);
        if (raw_lce > right_avail) raw_lce = right_avail;

        /* Clip by remaining K_total */
        int remaining = K_total - t;
        if (raw_lce > remaining) raw_lce = remaining;

        /* This run covers [t+1, t+run_len] (1-indexed k values) */
        int run_start = t + 1;   /* first k in this run (k = t+1 at trial t=0, first position) */
        int run_end   = t + raw_lce;  /* last k in this run */

        if (run_start <= run_end && run_end >= minrep) {
            if (n_runs < max_runs) {
                runs[n_runs].run_start        = run_start;
                runs[n_runs].run_end          = run_end;
                runs[n_runs].mismatches_before = i;
                n_runs++;
            }
        }

        t += raw_lce;  /* advance past the exact-match run */

        if (t >= K_total) break;  /* hit K_total exactly — done */

        /* t is now at a mismatch position (k = t+1 including this mismatch). */
        i++;  /* accumulate this mismatch */

        /* Record the mismatch step itself as a candidate arm endpoint.
         * The legacy inner loop evaluates k=(t+1) even when the last character
         * is a mismatch, so we must capture it here too. */
        {
            int mis_k = t + 1;  /* arm length ending at this mismatch */
            if (mis_k >= minrep && mis_k <= K_total && mis_k <= pt->K) {
                if (n_runs < max_runs) {
                    runs[n_runs].run_start         = mis_k;
                    runs[n_runs].run_end           = mis_k;
                    runs[n_runs].mismatches_before = i;  /* total, including this step */
                    n_runs++;
                }
            }
        }

        t++;  /* advance past the mismatch */

        /* Early-exit check: mirrors the C code's "if k >= minrep { if mis > budget break }" */
        /* Here k = t after consuming the mismatch */
        if (t >= minrep && i > mismatch_budget) break;
    }

    /* Scan runs back-to-front to find the largest valid best_k */
    int result = 0;
    for (int r = n_runs - 1; r >= 0; r--) {
        int run_start_k = runs[r].run_start;
        int run_end_k   = runs[r].run_end;
        int i_before    = runs[r].mismatches_before;

        /* Clip run_end_k to K_total and purity table range */
        if (run_end_k > K_total)  run_end_k = K_total;
        if (run_end_k > pt->K)    run_end_k = pt->K;

        /* lo' = max(run_start_k, minrep, ceil(100*i_before / mismatch_tol_int)) */
        int lo = run_start_k;
        if (minrep > lo) lo = minrep;

        /* ratio lower bound: k >= ceil(100*i / mismatch_tol_int) */
        if (i_before > 0) {
            int ratio_lo = (100 * i_before + mismatch_tol_int - 1) / mismatch_tol_int;
            if (ratio_lo > lo) lo = ratio_lo;
        }

        if (lo > run_end_k) continue;  /* this run cannot contribute */

        /* best_purity_k[run_end_k] = largest k in [1, run_end_k] with purity true */
        int32_t best_p = pt->best_k[run_end_k];
        if (best_p >= lo) {
            result = best_p;
            break;  /* largest k found — stop */
        }
    }

    if (heap_runs) free(heap_runs);
    return result;
}
