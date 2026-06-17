/*
 * _hdna.c — CPython extension: H-DNA / Triplex Mirror Repeat Detector
 *
 * Core algorithm ported verbatim from findHDNA.c (MirrorHunter / TriplexDetector).
 * Original algorithm inspired by non-B_gfa (Cer et al., Nucleic Acids Res. 2013).
 *
 * When compiled via setup.py / pip install, this becomes the
 * Python C extension module  hseeker._hdna.
 *
 * When compiled with -DSTANDALONE it behaves as the original CLI binary:
 *   gcc -DSTANDALONE -O2 -Wall -o findHDNA _hdna.c -lm
 *
 * Python API exposed by this module:
 *   scan_sequence(seq, *, minrep=8, maxrep=3000, maxspacer=20,
 *                 purity=0.90, mismatch=0.10, remove_overlaps=True,
 *                 seq_offset=1) -> list[dict]
 *
 * Output dict keys per hit:
 *   start, end, arm_length, spacer_length, total_length,
 *   ga_pct, ct_pct, mirror_identity, is_perfect,
 *   left_arm, spacer, right_arm, full_sequence
 */

#define PY_SSIZE_T_CLEAN
#ifndef STANDALONE
#  include <Python.h>
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <math.h>
#include <time.h>

/* ------------------------------------------------------------------ */
/*  Constants                                                          */
/* ------------------------------------------------------------------ */

/*
 * Base-type lookup table — maps any byte to:
 *   1  = GA (A or G, upper or lower)
 *   2  = CT (C or T, upper or lower)
 *   0  = other (N, n, or any non-ACGT character)
 *
 * Using a single table lookup in the inner loop eliminates multiple
 * character comparisons and reduces branch mispredictions.
 */
static const unsigned char BT[256] = {
    ['A']=1, ['a']=1, ['G']=1, ['g']=1,
    ['C']=2, ['c']=2, ['T']=2, ['t']=2,
};

/* ── profiling counters (zeroed per call, exposed via Python) ──────────── */
static long long prof_inner_iters = 0;     /* total inner while-body executions  */
static long long prof_ctr_sp_pairs = 0;    /* how many (ctr, sp) pairs explored */
static long long prof_skipped_invalid_kmax = 0;
static long long prof_skipped_purity_prefilter = 0;
static long long prof_gate_passed = 0;
static long long prof_raw_hits = 0;
static long long prof_final_hits = 0;
static double    prof_runtime_seconds = 0.0;

static void prof_reset(void) {
    prof_inner_iters  = 0;
    prof_ctr_sp_pairs = 0;
    prof_skipped_invalid_kmax = 0;
    prof_skipped_purity_prefilter = 0;
    prof_gate_passed = 0;
    prof_raw_hits = 0;
    prof_final_hits = 0;
    prof_runtime_seconds = 0.0;
}

/*
 * Initial capacity of the per-call hits buffer in py_scan_sequence().
 * The buffer doubles automatically via realloc as hits accumulate, so
 * there is no hard cap — only available heap memory limits the count.
 */
#define INITIAL_HIT_CAPACITY  65536

/* Standalone CLI only: large static buffers to avoid stack overflow */
#ifdef STANDALONE
#  define MAX_DNA_STANDALONE  600000000
#  define MAX_HITS_STANDALONE  5000000
#  define MAX_SEQ_ID           256
#endif

/* ------------------------------------------------------------------ */
/*  Data structures                                                    */
/* ------------------------------------------------------------------ */

/*
 * HDNA_HIT — one H-DNA / triplex mirror repeat candidate.
 *
 * Coordinate convention: 1-based (NCBI / non-B_gfa compatible).
 *   start = first base of left arm
 *   end   = last  base of right arm
 *
 * *_start_idx fields are 0-based offsets into the caller-supplied dna[].
 */
typedef struct {
    long  start;             /* genomic start of left arm (1-based)       */
    long  end;               /* genomic end   of right arm (1-based)      */
    int   arm_len;           /* length of each arm (both arms are equal)  */
    int   spacer_len;        /* length of spacer between arms             */
    float ga_pct;            /* % GA (purine)  in right arm               */
    float ct_pct;            /* % CT (pyrimidine) in right arm            */
    float mirror_id;         /* mirror identity: fraction matched × 100   */
    int   is_perfect;        /* 1 = 100% pure AND exact mirror            */
    int   left_start_idx;    /* 0-based index into dna[] for left arm     */
    int   spacer_start_idx;  /* 0-based index into dna[] for spacer       */
    int   right_start_idx;   /* 0-based index into dna[] for right arm    */
} HDNA_HIT;

/* ------------------------------------------------------------------ */
/*  Core algorithm                                                     */
/* ------------------------------------------------------------------ */

/*
 * findHDNA_core
 *
 * Center-outward mirror repeat scanner with fuzzy matching.
 *
 * For each candidate center position 'ctr', for each spacer length 'sp':
 *   - Virtual left pointer moves left, right pointer moves right.
 *   - Count mirror matches and arm nucleotide composition.
 *   - Record the longest valid arm at each (ctr, sp) if purity and
 *     mismatch thresholds are satisfied.
 *
 * This is identical to the original findHDNA() from findHDNA.c, except:
 *   - dna is a function argument (not a global) — makes it re-entrant.
 *   - hits array and max_hits are also arguments.
 *   - seq_offset is a parameter for genomic coordinate calculation.
 *
 * Parameters
 * ----------
 * dna          : lowercase DNA string (ACGTn)
 * total_bases  : length of dna
 * p_hits       : pointer to a malloc'd HDNA_HIT array; may be realloc'd on overflow
 * p_capacity   : pointer to current buffer capacity; updated on realloc
 * minrep       : minimum arm length
 * maxrep       : maximum arm length
 * maxspacer    : maximum spacer between arms
 * purity_thresh: min fraction of GA or CT in each arm
 * mismatch_tol : max fraction of mismatched mirror positions
 * seq_offset   : 1-based genomic position of dna[0]
 *
 * Returns
 * -------
 * Number of hits written.  Returns -1 on allocation failure (OOM); the
 * buffer at *p_hits is still valid and the caller must free it.
 */
static int findHDNA_core(
    const char *dna, int total_bases,
    HDNA_HIT **p_hits, int *p_capacity,
    int   minrep, int maxrep, int maxspacer,
    float purity_thresh, float mismatch_tol,
    long  seq_offset)
{
    int   ndx             = 0;
    int   mismatch_budget = (int)(mismatch_tol * maxrep);  /* loop-invariant */
    /* integer-scaled thresholds (×100) — eliminate float division from inner loop */
    int   purity_thresh_int = (int)(purity_thresh  * 100.0f + 0.5f);
    int   mismatch_tol_int  = (int)(mismatch_tol   * 100.0f + 0.5f);

    prof_reset();

    for (int ctr = minrep - 1; ctr <= total_bases - minrep - 1; ctr++) {
        /* skip centers that are N or non-ACGT — no sp can succeed */
        if (BT[(unsigned char)dna[ctr]] == 0) continue;

        int max_sp = (ctr + maxspacer < total_bases - minrep)
                      ? maxspacer
                      : (total_bases - minrep - ctr - 1);

        for (int sp = 0; sp <= max_sp; sp++) {
            prof_ctr_sp_pairs++;

            int   left_i    = ctr;
            int   right_j   = ctr + sp + 1;
            int   k         = 0;
            int   mismatches = 0;
            int   ga_count  = 0;
            int   ct_count  = 0;
            int   best_k    = 0;
            float best_mir  = 0.0f;
            float best_ga   = 0.0f;
            float best_ct   = 0.0f;

            /* extend outward one base at a time */
            while (left_i >= 0 && right_j < total_bases && k < maxrep) {
                prof_inner_iters++;

                char lb = dna[left_i];
                unsigned char bt_rb = BT[(unsigned char)dna[right_j]];

                /* N/other bases cannot participate in H-DNA — stop extending */
                if (bt_rb == 0) break;

                if (lb != dna[right_j]) mismatches++;
                k++;

                /* track composition of the right arm (read forward) */
                if (bt_rb == 1) ga_count++; else ct_count++;

                if (k >= minrep) {
                    /* integer threshold checks — no float divisions */
                    int ga_ok = ga_count * 100 >= purity_thresh_int * k;
                    int ct_ok = ct_count * 100 >= purity_thresh_int * k;

                    if ((ga_ok || ct_ok) &&
                        mismatches * 100 <= mismatch_tol_int * k) {
                        best_k = k;
                        /* compute floats only for the winning candidate */
                        float inv_k = 1.0f / k;
                        best_mir = 1.0f - (float)mismatches * inv_k;
                        best_ga  = (float)ga_count * inv_k;
                        best_ct  = (float)ct_count * inv_k;
                    }

                    /* early exit: mismatch budget already fully consumed */
                    if (mismatches > mismatch_budget) break;
                }

                left_i--;
                right_j++;
            }

            if (best_k >= minrep) {
                /* grow buffer on overflow */
                if (ndx >= *p_capacity) {
                    int new_cap = (*p_capacity) * 2;
                    HDNA_HIT *new_buf = (HDNA_HIT *)realloc(
                        *p_hits, sizeof(HDNA_HIT) * (size_t)new_cap);
                    if (!new_buf) return -1;  /* OOM; caller must free *p_hits */
                    *p_hits    = new_buf;
                    *p_capacity = new_cap;
                }

                int left_start  = ctr - best_k + 1;
                int right_start = ctr + sp + 1;

                HDNA_HIT *h     = &(*p_hits)[ndx];
                h->arm_len          = best_k;
                h->spacer_len       = sp;
                h->start            = (long)(left_start + 1) + (seq_offset - 1);
                h->end              = (long)(right_start + best_k) + (seq_offset - 1);
                h->ga_pct           = best_ga * 100.0f;
                h->ct_pct           = best_ct * 100.0f;
                h->mirror_id        = best_mir * 100.0f;
                h->is_perfect       = ((best_ga == 1.0f || best_ct == 1.0f)
                                       && best_mir == 1.0f) ? 1 : 0;
                h->left_start_idx   = left_start;
                h->spacer_start_idx = ctr + 1;
                h->right_start_idx  = right_start;
                ndx++;
            }
        }
    }
    return ndx;
}

static int append_hit(
    HDNA_HIT **p_hits, int *p_capacity, int *p_ndx,
    int ctr, int sp, int best_k,
    float best_mir, float best_ga, float best_ct,
    long seq_offset)
{
    if (*p_ndx >= *p_capacity) {
        int new_cap = (*p_capacity) * 2;
        HDNA_HIT *new_buf = (HDNA_HIT *)realloc(
            *p_hits, sizeof(HDNA_HIT) * (size_t)new_cap);
        if (!new_buf) return -1;
        *p_hits = new_buf;
        *p_capacity = new_cap;
    }

    int left_start  = ctr - best_k + 1;
    int right_start = ctr + sp + 1;

    HDNA_HIT *h = &(*p_hits)[*p_ndx];
    h->arm_len          = best_k;
    h->spacer_len       = sp;
    h->start            = (long)(left_start + 1) + (seq_offset - 1);
    h->end              = (long)(right_start + best_k) + (seq_offset - 1);
    h->ga_pct           = best_ga * 100.0f;
    h->ct_pct           = best_ct * 100.0f;
    h->mirror_id        = best_mir * 100.0f;
    h->is_perfect       = ((best_ga == 1.0f || best_ct == 1.0f)
                           && best_mir == 1.0f) ? 1 : 0;
    h->left_start_idx   = left_start;
    h->spacer_start_idx = ctr + 1;
    h->right_start_idx  = right_start;
    (*p_ndx)++;
    return 0;
}

/*
 * FastGatedHSeeker core
 *
 * This is an additive detector. It gates each (ctr, spacer) pair using
 * local density of mirror equality plus right-arm GA/CT purity, then falls
 * back to the original HSeeker-compatible outward scan once the gate passes.
 */
static int findHDNA_fast_gated_core(
    const char *dna, int total_bases,
    HDNA_HIT **p_hits, int *p_capacity,
    int   minrep, int maxrep, int maxspacer,
    float purity_thresh, float mismatch_tol,
    int   gate_window, int gate_search_limit,
    float gate_mirror_frac, float gate_purity_frac,
    int   use_purity_prefilter,
    long  seq_offset)
{
    int ndx = 0;
    int purity_thresh_int = (int)(purity_thresh * 100.0f + 0.5f);
    int gate_mirror_need = (int)ceilf(gate_mirror_frac * (float)gate_window);
    int gate_purity_need = (int)ceilf(gate_purity_frac * (float)gate_window);
    int *next_invalid_right = NULL;
    unsigned char *bt = NULL;
    unsigned char *right_pure_window_exists = NULL;
    unsigned char *ring_match = NULL;
    unsigned char *ring_bt = NULL;

    prof_reset();
    clock_t start_clock = clock();

    if (gate_window < 1) gate_window = minrep;
    if (gate_search_limit < gate_window) gate_search_limit = gate_window;

    next_invalid_right = (int *)malloc(sizeof(int) * (size_t)(total_bases + 1));
    bt = (unsigned char *)malloc(sizeof(unsigned char) * (size_t)total_bases);
    right_pure_window_exists = (unsigned char *)calloc((size_t)total_bases, 1);
    ring_match = (unsigned char *)malloc((size_t)gate_window);
    ring_bt = (unsigned char *)malloc((size_t)gate_window);
    if (!next_invalid_right || !bt || !right_pure_window_exists ||
        !ring_match || !ring_bt) {
        free(next_invalid_right); free(bt); free(right_pure_window_exists);
        free(ring_match); free(ring_bt);
        return -1;
    }

    next_invalid_right[total_bases] = total_bases;
    for (int i = total_bases - 1; i >= 0; i--) {
        bt[i] = BT[(unsigned char)dna[i]];
        next_invalid_right[i] = (bt[i] == 0) ? i : next_invalid_right[i + 1];
    }

    if (use_purity_prefilter && total_bases >= gate_window) {
        int *ga_prefix = (int *)calloc((size_t)(total_bases + 1), sizeof(int));
        int *ct_prefix = (int *)calloc((size_t)(total_bases + 1), sizeof(int));
        if (!ga_prefix || !ct_prefix) {
            free(ga_prefix); free(ct_prefix);
            free(next_invalid_right); free(bt); free(right_pure_window_exists);
            free(ring_match); free(ring_bt);
            return -1;
        }
        for (int i = 0; i < total_bases; i++) {
            ga_prefix[i + 1] = ga_prefix[i] + (bt[i] == 1);
            ct_prefix[i + 1] = ct_prefix[i] + (bt[i] == 2);
        }
        for (int r = 0; r <= total_bases - gate_window; r++) {
            int stop = r + gate_search_limit;
            if (stop > total_bases - gate_window) stop = total_bases - gate_window;
            for (int w = r; w <= stop; w++) {
                if (next_invalid_right[w] < w + gate_window) continue;
                int ga = ga_prefix[w + gate_window] - ga_prefix[w];
                int ct = ct_prefix[w + gate_window] - ct_prefix[w];
                if (ga >= gate_purity_need || ct >= gate_purity_need) {
                    right_pure_window_exists[r] = 1;
                    break;
                }
            }
        }
        free(ga_prefix); free(ct_prefix);
    }

    for (int ctr = minrep - 1; ctr <= total_bases - minrep - 1; ctr++) {
        if (bt[ctr] == 0) continue;

        int max_sp = (ctr + maxspacer < total_bases - minrep)
                      ? maxspacer
                      : (total_bases - minrep - ctr - 1);

        for (int sp = 0; sp <= max_sp; sp++) {
            prof_ctr_sp_pairs++;

            int right_start = ctr + sp + 1;
            int kmax = maxrep;
            if (ctr + 1 < kmax) kmax = ctr + 1;
            if (total_bases - right_start < kmax) kmax = total_bases - right_start;
            if (next_invalid_right[right_start] - right_start < kmax)
                kmax = next_invalid_right[right_start] - right_start;

            if (kmax < minrep) {
                prof_skipped_invalid_kmax++;
                continue;
            }
            if (use_purity_prefilter && right_pure_window_exists[right_start] == 0) {
                prof_skipped_purity_prefilter++;
                continue;
            }

            int gate_limit = gate_search_limit;
            if (gate_limit > kmax) gate_limit = kmax;
            int rolling_match = 0, rolling_ga = 0, rolling_ct = 0;
            int gate_passed = 0;

            for (int u = 0; u < gate_limit; u++) {
                int slot = u % gate_window;
                if (u >= gate_window) {
                    rolling_match -= ring_match[slot];
                    if (ring_bt[slot] == 1) rolling_ga--;
                    else if (ring_bt[slot] == 2) rolling_ct--;
                }
                int left_i = ctr - u;
                int right_j = right_start + u;
                unsigned char rb = bt[right_j];
                unsigned char m = (dna[left_i] == dna[right_j]) ? 1 : 0;
                ring_match[slot] = m;
                ring_bt[slot] = rb;
                rolling_match += m;
                if (rb == 1) rolling_ga++;
                else if (rb == 2) rolling_ct++;

                prof_inner_iters++;

                if (u + 1 >= gate_window &&
                    rolling_match >= gate_mirror_need &&
                    (rolling_ga >= gate_purity_need || rolling_ct >= gate_purity_need)) {
                    gate_passed = 1;
                    break;
                }
            }

            if (!gate_passed) continue;
            prof_gate_passed++;

            int mismatches = 0, ga_count = 0, ct_count = 0;
            int best_k = 0;
            float best_mir = 0.0f, best_ga = 0.0f, best_ct = 0.0f;

            for (int k = 1; k <= kmax; k++) {
                int left_i = ctr - k + 1;
                int right_j = right_start + k - 1;
                unsigned char rb = bt[right_j];

                if (dna[left_i] != dna[right_j]) mismatches++;
                if (rb == 1) ga_count++; else ct_count++;

                prof_inner_iters++;

                if (k >= minrep) {
                    int ga_ok = ga_count * 100 >= purity_thresh_int * k;
                    int ct_ok = ct_count * 100 >= purity_thresh_int * k;
                    int max_mismatches = (int)floorf(mismatch_tol * (float)k);
                    if ((ga_ok || ct_ok) && mismatches <= max_mismatches) {
                        best_k = k;
                        float inv_k = 1.0f / (float)k;
                        best_mir = 1.0f - (float)mismatches * inv_k;
                        best_ga = (float)ga_count * inv_k;
                        best_ct = (float)ct_count * inv_k;
                    }
                }
            }

            if (best_k >= minrep) {
                if (append_hit(p_hits, p_capacity, &ndx, ctr, sp, best_k,
                               best_mir, best_ga, best_ct, seq_offset) < 0) {
                    free(next_invalid_right); free(bt); free(right_pure_window_exists);
                    free(ring_match); free(ring_bt);
                    return -1;
                }
            }
        }
    }

    prof_raw_hits = ndx;
    prof_final_hits = ndx;
    prof_runtime_seconds = (double)(clock() - start_clock) / (double)CLOCKS_PER_SEC;

    free(next_invalid_right); free(bt); free(right_pure_window_exists);
    free(ring_match); free(ring_bt);
    return ndx;
}

/* ------------------------------------------------------------------ */
/*  Overlap removal                                                    */
/* ------------------------------------------------------------------ */

/*
 * remove_overlaps_core
 *
 * For overlapping hits, keep the one with the longer arm; break ties by
 * keeping the shorter spacer.  Identical semantics to the original algorithm.
 *
 * O(n log n) implementation: sort by genomic start (longest arm first on
 * ties), then a single linear sweep that greedily keeps the best hit in
 * each overlap group.
 *
 * Returns the compacted hit count.
 */
static int hit_cmp_start(const void *a, const void *b)
{
    const HDNA_HIT *ha = (const HDNA_HIT *)a;
    const HDNA_HIT *hb = (const HDNA_HIT *)b;
    if (ha->start < hb->start) return -1;
    if (ha->start > hb->start) return  1;
    /* equal start: longer arm first so the greedy sweep picks it immediately */
    if (ha->arm_len > hb->arm_len) return -1;
    if (ha->arm_len < hb->arm_len) return  1;
    /* equal arm: shorter spacer first */
    return ha->spacer_len - hb->spacer_len;
}

static int remove_overlaps_core(HDNA_HIT *hits, int nhits)
{
    if (nhits <= 1) return nhits;

    /* O(n log n) sort by genomic start */
    qsort(hits, (size_t)nhits, sizeof(HDNA_HIT), hit_cmp_start);

    /*
     * Linear sweep:
     *   - If the current hit does not overlap the last kept hit, always keep it.
     *   - If it overlaps, replace the last kept hit only when the current hit
     *     has a longer arm (or equal arm and shorter spacer).
     * When we replace, we update last->end so the overlap zone shrinks or
     * expands correctly for subsequent hits.
     */
    int n = 0;
    for (int i = 0; i < nhits; i++) {
        if (n == 0) {
            hits[n++] = hits[i];
            continue;
        }
        HDNA_HIT *last = &hits[n - 1];
        if (hits[i].start > last->end) {
            /* no overlap — always keep */
            hits[n++] = hits[i];
        } else if (hits[i].arm_len > last->arm_len ||
                   (hits[i].arm_len == last->arm_len &&
                    hits[i].spacer_len < last->spacer_len)) {
            /* overlaps but is strictly better — replace in place */
            *last = hits[i];
        }
        /* else: overlaps and is not better — discard */
    }
    return n;
}

/* ================================================================== */
/*  CPython extension                                                  */
/* ================================================================== */

#ifndef STANDALONE

/*
 * py_scan_sequence
 *
 * Python signature:
 *   scan_sequence(seq: str, *,
 *                 minrep: int = 8,
 *                 maxrep: int = 3000,
 *                 maxspacer: int = 20,
 *                 purity: float = 0.90,
 *                 mismatch: float = 0.10,
 *                 remove_overlaps: bool = True,
 *                 seq_offset: int = 1) -> list[dict]
 *
 * seq is accepted as a raw byte string (ASCII DNA).  The GIL is released
 * for the compute-intensive C scan so that other Python threads can run.
 */
static PyObject *
py_scan_sequence(PyObject *self, PyObject *args, PyObject *kwargs)
{
    const char *raw_seq   = NULL;
    Py_ssize_t  raw_len   = 0;
    int         minrep    = 8;
    int         maxrep    = 3000;
    int         maxspacer = 20;
    double      purity    = 0.90;
    double      mismatch  = 0.10;
    int         do_overlap = 1;
    long        seq_offset = 1;

    static char *kwlist[] = {
        "seq",
        "minrep", "maxrep", "maxspacer",
        "purity", "mismatch",
        "remove_overlaps", "seq_offset",
        NULL
    };

    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs,
            "s#|iiiddil",
            kwlist,
            &raw_seq, &raw_len,
            &minrep, &maxrep, &maxspacer,
            &purity, &mismatch,
            &do_overlap, &seq_offset))
        return NULL;

    /* parameter validation */
    if (minrep < 1 || maxrep < minrep || maxspacer < 0 ||
        purity  < 0.0 || purity  > 1.0 ||
        mismatch < 0.0 || mismatch > 1.0) {
        PyErr_SetString(PyExc_ValueError,
            "Invalid parameters: minrep>=1, maxrep>=minrep, maxspacer>=0, "
            "0<=purity<=1, 0<=mismatch<=1");
        return NULL;
    }

    if (raw_len == 0)
        return PyList_New(0);

    /* ---- allocate working buffers ---- */

    /* lowercase copy of input */
    char *dna = (char *)malloc((size_t)raw_len + 1);
    if (!dna) { PyErr_NoMemory(); return NULL; }
    for (Py_ssize_t i = 0; i < raw_len; i++)
        dna[i] = (char)tolower((unsigned char)raw_seq[i]);
    dna[raw_len] = '\0';

    /* hits array — starts small and doubles inside findHDNA_core as needed */
    int hit_capacity = INITIAL_HIT_CAPACITY;
    HDNA_HIT *hits = (HDNA_HIT *)malloc(sizeof(HDNA_HIT) * (size_t)hit_capacity);
    if (!hits) { free(dna); PyErr_NoMemory(); return NULL; }

    /* arm / spacer extraction buffers (sized to worst-case arm length) */
    size_t arm_sz    = (size_t)(maxrep + 2);
    size_t spacer_sz = (size_t)(maxspacer + 2);
    size_t full_sz   = arm_sz * 2 + spacer_sz + 2;

    char *left_arm_buf  = (char *)malloc(arm_sz);
    char *spacer_buf    = (char *)malloc(spacer_sz);
    char *right_arm_buf = (char *)malloc(arm_sz);
    char *full_seq_buf  = (char *)malloc(full_sz);

    if (!left_arm_buf || !spacer_buf || !right_arm_buf || !full_seq_buf) {
        free(left_arm_buf); free(spacer_buf);
        free(right_arm_buf); free(full_seq_buf);
        free(hits); free(dna);
        PyErr_NoMemory();
        return NULL;
    }

    /* ---- run C computation (GIL released) ---- */
    int nhits;
    int raw_nhits;
    clock_t start_clock = clock();

    Py_BEGIN_ALLOW_THREADS
        nhits = findHDNA_core(
            dna, (int)raw_len,
            &hits, &hit_capacity,
            minrep, maxrep, maxspacer,
            (float)purity, (float)mismatch,
            seq_offset);

        raw_nhits = nhits;
        if (nhits >= 0 && do_overlap && nhits > 1)
            nhits = remove_overlaps_core(hits, nhits);
    Py_END_ALLOW_THREADS

    prof_raw_hits = raw_nhits >= 0 ? raw_nhits : 0;
    prof_final_hits = nhits >= 0 ? nhits : 0;
    prof_runtime_seconds = (double)(clock() - start_clock) / (double)CLOCKS_PER_SEC;

    if (nhits < 0) {
        /* realloc failed inside findHDNA_core; hits is still the last valid ptr */
        free(hits); free(dna);
        free(left_arm_buf); free(spacer_buf);
        free(right_arm_buf); free(full_seq_buf);
        PyErr_NoMemory();
        return NULL;
    }

    /* ---- build Python result list ---- */
    PyObject *result = PyList_New((Py_ssize_t)nhits);
    if (!result) goto cleanup_err;

    for (int i = 0; i < nhits; i++) {
        HDNA_HIT *h = &hits[i];

        /* extract arm / spacer sequences from the dna buffer */
        strncpy(left_arm_buf,  dna + h->left_start_idx,  (size_t)h->arm_len);
        left_arm_buf[h->arm_len] = '\0';

        if (h->spacer_len > 0) {
            strncpy(spacer_buf, dna + h->spacer_start_idx, (size_t)h->spacer_len);
            spacer_buf[h->spacer_len] = '\0';
        } else {
            spacer_buf[0] = '.';
            spacer_buf[1] = '\0';
        }

        strncpy(right_arm_buf, dna + h->right_start_idx, (size_t)h->arm_len);
        right_arm_buf[h->arm_len] = '\0';

        if (h->spacer_len > 0)
            snprintf(full_seq_buf, full_sz, "%s%s%s",
                     left_arm_buf, spacer_buf, right_arm_buf);
        else
            snprintf(full_seq_buf, full_sz, "%s%s",
                     left_arm_buf, right_arm_buf);

        int total_len = h->arm_len * 2 + h->spacer_len;

        /* build hit dict */
        PyObject *d = PyDict_New();
        if (!d) {
            Py_DECREF(result);
            goto cleanup_err;
        }

#define _SET(key, val) \
    do { \
        PyObject *_v = (val); \
        if (!_v || PyDict_SetItemString(d, (key), _v) < 0) { \
            Py_XDECREF(_v); Py_DECREF(d); Py_DECREF(result); \
            goto cleanup_err; \
        } \
        Py_DECREF(_v); \
    } while (0)

        _SET("start",           PyLong_FromLong(h->start));
        _SET("end",             PyLong_FromLong(h->end));
        _SET("arm_length",      PyLong_FromLong((long)h->arm_len));
        _SET("spacer_length",   PyLong_FromLong((long)h->spacer_len));
        _SET("total_length",    PyLong_FromLong((long)total_len));
        _SET("ga_pct",          PyFloat_FromDouble((double)h->ga_pct));
        _SET("ct_pct",          PyFloat_FromDouble((double)h->ct_pct));
        _SET("mirror_identity", PyFloat_FromDouble((double)h->mirror_id));
        _SET("is_perfect",      PyBool_FromLong((long)h->is_perfect));
        _SET("left_arm",        PyUnicode_FromString(left_arm_buf));
        _SET("spacer",          PyUnicode_FromString(spacer_buf));
        _SET("right_arm",       PyUnicode_FromString(right_arm_buf));
        _SET("full_sequence",   PyUnicode_FromString(full_seq_buf));

#undef _SET

        PyList_SET_ITEM(result, (Py_ssize_t)i, d); /* steals ref to d */
    }

    /* free working buffers before potentially triggering a Python warning */
    free(left_arm_buf); free(spacer_buf);
    free(right_arm_buf); free(full_seq_buf);
    free(hits); free(dna);

    return result;

cleanup_err:
    free(left_arm_buf); free(spacer_buf);
    free(right_arm_buf); free(full_seq_buf);
    free(hits); free(dna);
    return NULL;
}

static PyObject *
py_scan_sequence_fast_gated(PyObject *self, PyObject *args, PyObject *kwargs)
{
    const char *raw_seq   = NULL;
    Py_ssize_t  raw_len   = 0;
    int         minrep    = 10;
    int         maxrep    = 1000;
    int         maxspacer = 10;
    double      purity    = 0.90;
    double      mismatch  = 0.10;
    int         gate_window = 0;
    int         gate_search_limit = 0;
    double      gate_mirror_frac = -1.0;
    double      gate_purity_frac = -1.0;
    int         fast_mode = 1;
    int         do_overlap = 1;
    int         use_purity_prefilter = 1;
    long        seq_offset = 1;

    static char *kwlist[] = {
        "seq",
        "minrep", "maxrep", "maxspacer",
        "purity", "mismatch",
        "gate_window", "gate_search_limit",
        "gate_mirror_frac", "gate_purity_frac",
        "fast_mode", "remove_overlaps",
        "use_purity_prefilter", "seq_offset",
        NULL
    };

    (void)self;
    if (!PyArg_ParseTupleAndKeywords(
            args, kwargs,
            "s#|iiiddiiddiiil",
            kwlist,
            &raw_seq, &raw_len,
            &minrep, &maxrep, &maxspacer,
            &purity, &mismatch,
            &gate_window, &gate_search_limit,
            &gate_mirror_frac, &gate_purity_frac,
            &fast_mode, &do_overlap,
            &use_purity_prefilter, &seq_offset))
        return NULL;

    if (minrep < 1 || maxrep < minrep || maxspacer < 0 ||
        purity < 0.0 || purity > 1.0 ||
        mismatch < 0.0 || mismatch > 1.0) {
        PyErr_SetString(PyExc_ValueError,
            "Invalid parameters: minrep>=1, maxrep>=minrep, maxspacer>=0, "
            "0<=purity<=1, 0<=mismatch<=1");
        return NULL;
    }

    if (raw_len == 0)
        return PyList_New(0);

    if (gate_window <= 0) gate_window = minrep;
    if (gate_search_limit <= 0)
        gate_search_limit = fast_mode ? 4 * minrep : maxrep;
    if (gate_search_limit > maxrep) gate_search_limit = maxrep;
    if (gate_mirror_frac < 0.0) gate_mirror_frac = 1.0 - mismatch;
    if (gate_purity_frac < 0.0) gate_purity_frac = purity;
    if (gate_window < 1 || gate_window > maxrep ||
        gate_search_limit < gate_window ||
        gate_mirror_frac < 0.0 || gate_mirror_frac > 1.0 ||
        gate_purity_frac < 0.0 || gate_purity_frac > 1.0) {
        PyErr_SetString(PyExc_ValueError,
            "Invalid gate parameters: 1<=gate_window<=maxrep, "
            "gate_search_limit>=gate_window, 0<=gate fractions<=1");
        return NULL;
    }

    char *dna = (char *)malloc((size_t)raw_len + 1);
    if (!dna) { PyErr_NoMemory(); return NULL; }
    for (Py_ssize_t i = 0; i < raw_len; i++)
        dna[i] = (char)tolower((unsigned char)raw_seq[i]);
    dna[raw_len] = '\0';

    int hit_capacity = INITIAL_HIT_CAPACITY;
    HDNA_HIT *hits = (HDNA_HIT *)malloc(sizeof(HDNA_HIT) * (size_t)hit_capacity);
    if (!hits) { free(dna); PyErr_NoMemory(); return NULL; }

    size_t arm_sz = (size_t)(maxrep + 2);
    size_t spacer_sz = (size_t)(maxspacer + 2);
    size_t full_sz = arm_sz * 2 + spacer_sz + 2;
    char *left_arm_buf = (char *)malloc(arm_sz);
    char *spacer_buf = (char *)malloc(spacer_sz);
    char *right_arm_buf = (char *)malloc(arm_sz);
    char *full_seq_buf = (char *)malloc(full_sz);
    if (!left_arm_buf || !spacer_buf || !right_arm_buf || !full_seq_buf) {
        free(left_arm_buf); free(spacer_buf);
        free(right_arm_buf); free(full_seq_buf);
        free(hits); free(dna);
        PyErr_NoMemory();
        return NULL;
    }

    int nhits;
    int raw_nhits;
    Py_BEGIN_ALLOW_THREADS
        nhits = findHDNA_fast_gated_core(
            dna, (int)raw_len,
            &hits, &hit_capacity,
            minrep, maxrep, maxspacer,
            (float)purity, (float)mismatch,
            gate_window, gate_search_limit,
            (float)gate_mirror_frac, (float)gate_purity_frac,
            use_purity_prefilter,
            seq_offset);
        raw_nhits = nhits;
        if (nhits >= 0 && do_overlap && nhits > 1)
            nhits = remove_overlaps_core(hits, nhits);
    Py_END_ALLOW_THREADS

    prof_raw_hits = raw_nhits >= 0 ? raw_nhits : 0;
    prof_final_hits = nhits >= 0 ? nhits : 0;

    if (nhits < 0) {
        free(hits); free(dna);
        free(left_arm_buf); free(spacer_buf);
        free(right_arm_buf); free(full_seq_buf);
        PyErr_NoMemory();
        return NULL;
    }

    PyObject *result = PyList_New((Py_ssize_t)nhits);
    if (!result) goto fg_cleanup_err;

    for (int i = 0; i < nhits; i++) {
        HDNA_HIT *h = &hits[i];
        strncpy(left_arm_buf, dna + h->left_start_idx, (size_t)h->arm_len);
        left_arm_buf[h->arm_len] = '\0';
        if (h->spacer_len > 0) {
            strncpy(spacer_buf, dna + h->spacer_start_idx, (size_t)h->spacer_len);
            spacer_buf[h->spacer_len] = '\0';
        } else {
            spacer_buf[0] = '.';
            spacer_buf[1] = '\0';
        }
        strncpy(right_arm_buf, dna + h->right_start_idx, (size_t)h->arm_len);
        right_arm_buf[h->arm_len] = '\0';
        if (h->spacer_len > 0)
            snprintf(full_seq_buf, full_sz, "%s%s%s",
                     left_arm_buf, spacer_buf, right_arm_buf);
        else
            snprintf(full_seq_buf, full_sz, "%s%s",
                     left_arm_buf, right_arm_buf);

        int total_len = h->arm_len * 2 + h->spacer_len;
        PyObject *d = PyDict_New();
        if (!d) {
            Py_DECREF(result);
            goto fg_cleanup_err;
        }

#define _SET_FG(key, val) \
    do { \
        PyObject *_v = (val); \
        if (!_v || PyDict_SetItemString(d, (key), _v) < 0) { \
            Py_XDECREF(_v); Py_DECREF(d); Py_DECREF(result); \
            goto fg_cleanup_err; \
        } \
        Py_DECREF(_v); \
    } while (0)

        _SET_FG("start",           PyLong_FromLong(h->start));
        _SET_FG("end",             PyLong_FromLong(h->end));
        _SET_FG("arm_length",      PyLong_FromLong((long)h->arm_len));
        _SET_FG("spacer_length",   PyLong_FromLong((long)h->spacer_len));
        _SET_FG("total_length",    PyLong_FromLong((long)total_len));
        _SET_FG("ga_pct",          PyFloat_FromDouble((double)h->ga_pct));
        _SET_FG("ct_pct",          PyFloat_FromDouble((double)h->ct_pct));
        _SET_FG("mirror_identity", PyFloat_FromDouble((double)h->mirror_id));
        _SET_FG("is_perfect",      PyBool_FromLong((long)h->is_perfect));
        _SET_FG("left_arm",        PyUnicode_FromString(left_arm_buf));
        _SET_FG("spacer",          PyUnicode_FromString(spacer_buf));
        _SET_FG("right_arm",       PyUnicode_FromString(right_arm_buf));
        _SET_FG("full_sequence",   PyUnicode_FromString(full_seq_buf));

#undef _SET_FG

        PyList_SET_ITEM(result, (Py_ssize_t)i, d);
    }

    free(left_arm_buf); free(spacer_buf);
    free(right_arm_buf); free(full_seq_buf);
    free(hits); free(dna);
    return result;

fg_cleanup_err:
    free(left_arm_buf); free(spacer_buf);
    free(right_arm_buf); free(full_seq_buf);
    free(hits); free(dna);
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Module definition                                                  */
/* ------------------------------------------------------------------ */

PyDoc_STRVAR(scan_sequence_doc,
"scan_sequence(seq, *, minrep=8, maxrep=3000, maxspacer=20,\n"
"              purity=0.90, mismatch=0.10, remove_overlaps=True,\n"
"              seq_offset=1) -> list[dict]\n"
"\n"
"Scan a raw DNA string for H-DNA / triplex mirror repeat motifs.\n"
"\n"
"Parameters\n"
"----------\n"
"seq : str\n"
"    DNA sequence (ACGTN; case-insensitive). 'N' bases break arm extension.\n"
"minrep : int\n"
"    Minimum arm length in bases (default 6).\n"
"maxrep : int\n"
"    Maximum arm length in bases (default 50).\n"
"maxspacer : int\n"
"    Maximum spacer between arms in bases (default 7).\n"
"purity : float\n"
"    Minimum fraction of GA or CT bases required in each arm (default 0.80).\n"
"    Set to 1.0 to reproduce strict non-B_gfa mirror-repeat results.\n"
"mismatch : float\n"
"    Maximum fraction of mirror-position mismatches allowed (default 0.20).\n"
"    Set to 0.0 for exact mirror only.\n"
"remove_overlaps : bool\n"
"    Remove overlapping hits, keeping the longest arm (default True).\n"
"seq_offset : int\n"
"    1-based genomic start coordinate of seq[0] (default 1).\n"
"    Pass the chromosomal start position when seq is a genomic slice.\n"
"\n"
"Returns\n"
"-------\n"
"list[dict]\n"
"    Each dict has keys: start, end, arm_length, spacer_length,\n"
"    total_length, ga_pct, ct_pct, mirror_identity, is_perfect,\n"
"    left_arm, spacer, right_arm, full_sequence.\n"
"    Coordinates are 1-based and inclusive.\n"
);

/* expose profiling counters to Python for diagnosis */
static PyObject *
py_profiling_info(PyObject *self, PyObject *noargs)
{
    (void)self;
    (void)noargs;
    return Py_BuildValue("{s:L, s:L, s:L, s:L, s:L, s:L, s:L, s:d}",
        "inner_iters",  prof_inner_iters,
        "ctr_sp_pairs", prof_ctr_sp_pairs,
        "skipped_invalid_kmax", prof_skipped_invalid_kmax,
        "skipped_purity_prefilter", prof_skipped_purity_prefilter,
        "gate_passed", prof_gate_passed,
        "raw_hits", prof_raw_hits,
        "final_hits", prof_final_hits,
        "runtime_seconds", prof_runtime_seconds);
}

static PyMethodDef HdnaMethods[] = {
    {
        "scan_sequence",
        (PyCFunction)py_scan_sequence,
        METH_VARARGS | METH_KEYWORDS,
        scan_sequence_doc,
    },
    {
        "scan_sequence_fast_gated",
        (PyCFunction)py_scan_sequence_fast_gated,
        METH_VARARGS | METH_KEYWORDS,
        "Scan a raw DNA string with FastGatedHSeeker."
    },
    {
        "profiling_info",
        (PyCFunction)py_profiling_info,
        METH_NOARGS,
        "Return profiling counters from the most recent detector call."
    },
    { NULL, NULL, 0, NULL }
};

PyDoc_STRVAR(module_doc,
"_hdna — low-level C extension for H-DNA / Triplex Mirror Repeat detection.\n\n"
"Use hseeker.scan_sequence() or hseeker.scan_fasta() instead of\n"
"calling this module directly.\n"
);

static struct PyModuleDef hdna_module = {
    PyModuleDef_HEAD_INIT,
    "_hdna",
    module_doc,
    -1,
    HdnaMethods
};

PyMODINIT_FUNC
PyInit__hdna(void)
{
    PyObject *m = PyModule_Create(&hdna_module);
    if (!m) return NULL;
    PyModule_AddIntConstant(m, "INITIAL_HIT_CAPACITY", INITIAL_HIT_CAPACITY);
    PyModule_AddStringConstant(m, "__version__", "0.1.0");
    return m;
}

/* ================================================================== */
/*  STANDALONE CLI mode  (gcc -DSTANDALONE _hdna.c -lm -o findHDNA)  */
/* ================================================================== */
#else  /* STANDALONE */

static char      sa_dna[MAX_DNA_STANDALONE + 1];
static HDNA_HIT *sa_hits     = NULL;
static int       sa_hits_cap = 0;
static char      sa_seq_id[MAX_SEQ_ID];
static long      sa_seq_offset;

static void sa_parse_header(const char *line) {
    char buf[MAX_SEQ_ID];
    long s = 0;
    const char *p = line;
    while (*p == '>' || *p == ' ') p++;
    int n = 0;
    while (*p && !isspace((unsigned char)*p) && n < MAX_SEQ_ID - 1)
        buf[n++] = *p++;
    buf[n] = '\0';
    char *colon = strchr(buf, ':');
    if (colon) {
        *colon = '\0';
        strncpy(sa_seq_id, buf, MAX_SEQ_ID - 1);
        sa_seq_id[MAX_SEQ_ID - 1] = '\0';
        s = atol(colon + 1);
    } else {
        strncpy(sa_seq_id, buf, MAX_SEQ_ID - 1);
        sa_seq_id[MAX_SEQ_ID - 1] = '\0';
        s = 1;
    }
    sa_seq_offset = s;
}

static int sa_read_fasta(FILE *fp) {
    int base, n = 0, truncated = 0;
    while ((base = getc(fp)) != EOF && base != '>') {}
    if (base == EOF) return 0;
    {
        char hdr[MAX_SEQ_ID + 64];
        int hi = 0;
        hdr[hi++] = '>';
        while ((base = getc(fp)) != EOF && base != '\n')
            if (hi < (int)sizeof(hdr) - 1) hdr[hi++] = (char)base;
        hdr[hi] = '\0';
        sa_parse_header(hdr);
    }
    while ((base = getc(fp)) != EOF) {
        if (base == '>') { ungetc(base, fp); break; }
        if (isalpha(base)) {
            if (n < MAX_DNA_STANDALONE)
                sa_dna[n++] = (char)tolower(base);
            else if (!truncated) {
                fprintf(stderr, "WARNING: '%s' truncated at %d bases.\n",
                        sa_seq_id, MAX_DNA_STANDALONE);
                truncated = 1;
            }
        }
    }
    sa_dna[n] = '\0';
    return n;
}

#define SA_ARM_BUF 100001

static void sa_print_header(FILE *fp) {
    fprintf(fp,
        "seq_id\tsource\tstart\tend\t"
        "arm_length\tspacer_length\ttotal_length\t"
        "ga_pct\tct_pct\tmirror_identity\tis_perfect\t"
        "left_arm\tspacer\tright_arm\tfull_sequence\n");
}

static void sa_print_hits(FILE *fp, int nhits) {
    char la[SA_ARM_BUF], sp_buf[SA_ARM_BUF], ra[SA_ARM_BUF];
    char fs[3 * SA_ARM_BUF + 2];
    for (int i = 0; i < nhits; i++) {
        HDNA_HIT *h = &sa_hits[i];
        strncpy(la, sa_dna + h->left_start_idx,  h->arm_len); la[h->arm_len] = '\0';
        if (h->spacer_len > 0) {
            strncpy(sp_buf, sa_dna + h->spacer_start_idx, h->spacer_len);
            sp_buf[h->spacer_len] = '\0';
        } else {
            strcpy(sp_buf, ".");
        }
        strncpy(ra, sa_dna + h->right_start_idx, h->arm_len); ra[h->arm_len] = '\0';
        if (h->spacer_len > 0)
            snprintf(fs, sizeof(fs), "%s%s%s", la, sp_buf, ra);
        else
            snprintf(fs, sizeof(fs), "%s%s", la, ra);
        fprintf(fp,
            "%s\tfindHDNA\t%ld\t%ld\t%d\t%d\t%d\t"
            "%.1f\t%.1f\t%.1f\t%d\t%s\t%s\t%s\t%s\n",
            sa_seq_id, h->start, h->end,
            h->arm_len, h->spacer_len, h->arm_len * 2 + h->spacer_len,
            h->ga_pct, h->ct_pct, h->mirror_id, h->is_perfect,
            la, sp_buf, ra, fs);
    }
}

static void sa_usage(const char *prog) {
    fprintf(stderr,
        "\nUsage: %s -seq <fasta> -out <prefix> [options]\n\n"
        "  -minrep    <int>    Minimum arm length       (default: 8)\n"
        "  -maxrep    <int>    Maximum arm length       (default: 3000)\n"
        "  -maxspacer <int>    Maximum spacer length    (default: 20)\n"
        "  -purity    <float>  Min GA or CT fraction    (default: 0.90)\n"
        "  -mismatch  <float>  Max mismatch fraction    (default: 0.10)\n"
        "  -skipoverlap        Skip overlap removal\n"
        "  -v                  Verbose\n\n", prog);
}

int main(int argc, char *argv[]) {
    char seq_fn[512] = {0}, out_pre[512] = {0}, out_fn[520] = {0};
    int  minrep = 8, maxrep = 3000, maxspacer = 20, do_overlap = 1, verbose = 0;
    float purity = 0.90f, mismatch = 0.10f;

    if (argc == 1) { sa_usage(argv[0]); return 1; }

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "-seq")        && argv[i+1]) strncpy(seq_fn,  argv[++i], 511);
        else if (!strcmp(argv[i], "-out")        && argv[i+1]) strncpy(out_pre, argv[++i], 511);
        else if (!strcmp(argv[i], "-minrep")     && argv[i+1]) minrep    = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-maxrep")     && argv[i+1]) maxrep    = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-maxspacer")  && argv[i+1]) maxspacer = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-purity")     && argv[i+1]) purity    = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "-mismatch")   && argv[i+1]) mismatch  = (float)atof(argv[++i]);
        else if (!strcmp(argv[i], "-skipoverlap"))              do_overlap = 0;
        else if (!strcmp(argv[i], "-v"))                        verbose    = 1;
        else if (!strcmp(argv[i], "-h") || !strcmp(argv[i], "--help"))
            { sa_usage(argv[0]); return 0; }
    }

    if (!seq_fn[0])  { fprintf(stderr, "ERROR: -seq required\n"); return 1; }
    if (!out_pre[0]) { fprintf(stderr, "ERROR: -out required\n"); return 1; }

    sa_hits = (HDNA_HIT *)malloc(sizeof(HDNA_HIT) * MAX_HITS_STANDALONE);
    if (!sa_hits) { fprintf(stderr, "ERROR: out of memory\n"); return 4; }
    sa_hits_cap = MAX_HITS_STANDALONE;

    FILE *fq = fopen(seq_fn, "r");
    if (!fq) { fprintf(stderr, "ERROR: Cannot open %s\n", seq_fn); return 2; }

    snprintf(out_fn, sizeof(out_fn), "%s_HDNA.tsv", out_pre);
    FILE *fo = fopen(out_fn, "w");
    if (!fo) { fclose(fq); fprintf(stderr, "ERROR: Cannot open %s\n", out_fn); return 3; }

    fprintf(stderr,
        "findHDNA — H-DNA / Triplex Mirror Repeat Detector\n"
        "  Input : %s\n  Output: %s\n"
        "  minrep=%d  maxrep=%d  maxspacer=%d\n"
        "  purity=%.2f  mismatch=%.2f\n\n",
        seq_fn, out_fn, minrep, maxrep, maxspacer, purity, mismatch);

    sa_print_header(fo);

    int total_records = 0, total_hits = 0, bases;
    while ((bases = sa_read_fasta(fq)) > 0) {
        total_records++;
        if (verbose)
            fprintf(stderr, "Processing %s (%d bases, offset %ld)...\n",
                    sa_seq_id, bases, sa_seq_offset);
        int nhits = findHDNA_core(
            sa_dna, bases, &sa_hits, &sa_hits_cap,
            minrep, maxrep, maxspacer, purity, mismatch, sa_seq_offset);
        if (nhits < 0) {
            fprintf(stderr, "ERROR: out of memory during scan of '%s'\n", sa_seq_id);
            fclose(fq); fclose(fo); free(sa_hits); return 5;
        }
        if (do_overlap && nhits > 1)
            nhits = remove_overlaps_core(sa_hits, nhits);
        if (verbose) fprintf(stderr, "  %d hits\n", nhits);
        sa_print_hits(fo, nhits);
        total_hits += nhits;
    }

    fclose(fq); fclose(fo); free(sa_hits);
    fprintf(stderr,
        "\nDone.\n  Records: %d\n  Hits: %d\n  Output: %s\n",
        total_records, total_hits, out_fn);
    return 0;
}

#endif /* STANDALONE */
