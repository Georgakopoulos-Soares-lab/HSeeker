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
 *   scan_sequence(seq, *, minrep=6, maxrep=50, maxspacer=7,
 *                 purity=0.80, mismatch=0.20, remove_overlaps=True,
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

/* ------------------------------------------------------------------ */
/*  Constants                                                          */
/* ------------------------------------------------------------------ */

/*
 * Hard cap on hits returned per scan_sequence() call.
 * If this is exceeded, a RuntimeWarning is issued and partial results
 * are returned.  For whole-chromosome scans, call scan_sequence once per
 * contig (as scan_fasta() does) rather than concatenating chromosomes.
 */
#define DEFAULT_MAX_HITS  1000000

/* ------------------------------------------------------------------ */
/*  Opt 2.2 — tolower lookup table                                     */
/*                                                                     */
/* Initialized once in PyInit__hdna (CPython mode).  In STANDALONE     */
/* mode it is initialized at the top of main().  Each entry maps an    */
/* ASCII byte to its lowercase equivalent — single array index, no     */
/* locale check, auto-vectorized by the compiler at -O2.               */
/* ------------------------------------------------------------------ */
static unsigned char lc_table[256];

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
 * hits         : pre-allocated output array (at least max_hits elements)
 * max_hits     : capacity of hits[]
 * minrep       : minimum arm length
 * maxrep       : maximum arm length
 * maxspacer    : maximum spacer between arms
 * purity_thresh: min fraction of GA or CT in each arm
 * mismatch_tol : max fraction of mismatched mirror positions
 * seq_offset   : 1-based genomic position of dna[0]
 *
 * Returns
 * -------
 * Number of hits written to hits[].  If equal to max_hits, the capacity
 * was reached and the caller should issue a warning.
 */
static int findHDNA_core(
    const char *dna, int total_bases,
    HDNA_HIT   *hits, int max_hits,
    int   minrep, int maxrep, int maxspacer,
    float purity_thresh, float mismatch_tol,
    long  seq_offset)
{
    int   ndx           = 0;
    float min_mirror_id = 1.0f - mismatch_tol;

    /* Opt 2.3: mismatch_budget is constant for the entire scan
     * (depends only on mismatch_tol and maxrep, not on ctr or sp).
     * Computing it once avoids a float multiply + truncation on every
     * inner-loop iteration, which executes billions of times at genome
     * scale. */
    int mismatch_budget = (int)(mismatch_tol * maxrep);

    for (int ctr = minrep - 1; ctr <= total_bases - minrep - 1; ctr++) {

        int max_sp = (ctr + maxspacer < total_bases - minrep)
                      ? maxspacer
                      : (total_bases - minrep - ctr - 1);

        for (int sp = 0; sp <= max_sp; sp++) {

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

                char lb = dna[left_i];
                char rb = dna[right_j];

                /* N bases cannot participate in H-DNA — stop extending */
                if (lb == 'n' || rb == 'n') break;

                if (lb != rb) mismatches++;
                k++;

                /* track composition of the right arm (read forward) */
                if (rb == 'g' || rb == 'a') ga_count++;
                if (rb == 'c' || rb == 't') ct_count++;

                if (k >= minrep) {
                    float mir_id = 1.0f - (float)mismatches / k;
                    float ga_f   = (float)ga_count / k;
                    float ct_f   = (float)ct_count / k;

                    if (mir_id >= min_mirror_id &&
                        (ga_f >= purity_thresh || ct_f >= purity_thresh)) {
                        best_k   = k;
                        best_mir = mir_id;
                        best_ga  = ga_f;
                        best_ct  = ct_f;
                    }

                    /* Opt 2.3: use pre-computed budget (hoisted above) */
                    if (mismatches > mismatch_budget) break;
                }

                left_i--;
                right_j++;
            }

            if (best_k >= minrep) {
                if (ndx >= max_hits) return ndx; /* cap reached */

                int left_start  = ctr - best_k + 1;
                int right_start = ctr + sp + 1;

                HDNA_HIT *h     = &hits[ndx];
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

/* ------------------------------------------------------------------ */
/*  Overlap removal                                                    */
/* ------------------------------------------------------------------ */

/* ------------------------------------------------------------------ */
/*  Opt 2.1 — O(n log n) overlap removal                              */
/*                                                                     */
/* Replacement for the original O(n²) all-pairs scan.                 */
/*                                                                     */
/* Algorithm:                                                           */
/*   1. Sort hits by (start ASC, arm_len DESC, spacer_len ASC).        */
/*      This places the "preferred" candidate first within any cluster  */
/*      of hits that share the same start coordinate, and ensures that  */
/*      a longer-arm hit encountered later in start order still wins    */
/*      over a shorter-arm hit that started earlier.                    */
/*   2. Greedy forward sweep with a single "last-kept" pointer.         */
/*      For each hit:                                                   */
/*        - No overlap with last kept → keep it (append).              */
/*        - Overlaps last kept → keep whichever has the longer arm      */
/*          (tiebreak: shorter spacer).  If the incoming hit wins,      */
/*          replace the last kept entry in-place.                       */
/*                                                                       */
/* Equivalence with the original O(n²) pass:                            */
/*   The original pass marks-and-compacts using "arm_len = -1".  The    */
/*   greedy sweep produces the same winning set because:                 */
/*   (a) Hits are sorted so the best candidate in any cluster is seen   */
/*       first.                                                          */
/*   (b) The replace logic handles the case where a later hit (higher   */
/*       start, longer arm) beats the current last-kept hit.             */
/*   (c) The tiebreak policy (longer arm, then shorter spacer) is        */
/*       preserved exactly.                                              */
/*                                                                       */
/* Complexity: O(n log n) sort + O(n) sweep = O(n log n) overall.       */
/* ------------------------------------------------------------------ */

/* qsort comparator: (start ASC, arm_len DESC, spacer_len ASC) */
static int cmp_hits(const void *a, const void *b)
{
    const HDNA_HIT *ha = (const HDNA_HIT *)a;
    const HDNA_HIT *hb = (const HDNA_HIT *)b;
    if (ha->start != hb->start)
        return (ha->start < hb->start) ? -1 : 1;
    if (ha->arm_len != hb->arm_len)
        return (ha->arm_len > hb->arm_len) ? -1 : 1; /* longer arm first */
    if (ha->spacer_len != hb->spacer_len)
        return (ha->spacer_len < hb->spacer_len) ? -1 : 1; /* shorter spacer first */
    return 0;
}

static int remove_overlaps_core(HDNA_HIT *hits, int nhits)
{
    if (nhits <= 1) return nhits;

    /* Step 1: sort so the preferred candidate leads each overlap cluster */
    qsort(hits, (size_t)nhits, sizeof(HDNA_HIT), cmp_hits);

    /* Step 2: greedy forward sweep */
    int n = 1; /* hits[0] is always kept */
    for (int i = 1; i < nhits; i++) {
        HDNA_HIT *prev = &hits[n - 1];
        HDNA_HIT *cur  = &hits[i];

        if (prev->start <= cur->end && cur->start <= prev->end) {
            /* intervals overlap — keep the better one */
            if (cur->arm_len > prev->arm_len ||
               (cur->arm_len == prev->arm_len &&
                cur->spacer_len < prev->spacer_len)) {
                /* cur wins: replace last kept entry in-place */
                *prev = *cur;
            }
            /* else prev wins: discard cur (do nothing) */
        } else {
            /* no overlap: append cur */
            hits[n++] = *cur;
        }
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
 *                 minrep: int = 6,
 *                 maxrep: int = 50,
 *                 maxspacer: int = 7,
 *                 purity: float = 0.80,
 *                 mismatch: float = 0.20,
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
    int         minrep    = 6;
    int         maxrep    = 50;
    int         maxspacer = 7;
    double      purity    = 0.80;
    double      mismatch  = 0.20;
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

    /* parameter validation
     * Upper bounds on maxrep / maxspacer prevent integer-overflow in
     * allocation-size arithmetic (opt 1.4).  100 000 bp arms are
     * biologically unrealistic and would produce absurdly large buffers.
     */
    if (minrep < 1 || maxrep < minrep || maxrep > 100000 ||
        maxspacer < 0 || maxspacer > 100000 ||
        purity  < 0.0 || purity  > 1.0 ||
        mismatch < 0.0 || mismatch > 1.0) {
        PyErr_SetString(PyExc_ValueError,
            "Invalid parameters: minrep>=1, maxrep>=minrep, "
            "maxrep<=100000, maxspacer>=0, maxspacer<=100000, "
            "0<=purity<=1, 0<=mismatch<=1");
        return NULL;
    }

    if (raw_len == 0)
        return PyList_New(0);

    /* ----------------------------------------------------------------
     * Allocate working buffers
     *
     * Opt 1.3: all extension-internal buffers use PyMem_Malloc /
     *          PyMem_Free so they are visible to Python memory tools
     *          (tracemalloc, valgrind wrappers, custom allocators).
     *
     * Opt 1.2: intermediate arm/spacer C string buffers are eliminated;
     *          arm strings are built directly from dna[] slices when
     *          constructing Python objects, so no extra heap allocation
     *          is needed here.
     * ---------------------------------------------------------------- */

    /* Opt 2.2 + 1.3: lowercase copy via lookup table (PyMem_Malloc).
     * lc_table[c] == tolower(c) for all c in 0..255, initialized once
     * at module load.  Single array index per byte — no locale check,
     * compiler-auto-vectorizable at -O2. */
    char *dna = (char *)PyMem_Malloc((size_t)raw_len + 1);
    if (!dna) { PyErr_NoMemory(); return NULL; }
    for (Py_ssize_t i = 0; i < raw_len; i++)
        dna[i] = (char)lc_table[(unsigned char)raw_seq[i]];
    dna[raw_len] = '\0';

    /* ----------------------------------------------------------------
     * Opt 1.1: dynamic hit buffer — start with a length-proportional
     * capacity, grow by 2× as needed up to DEFAULT_MAX_HITS.
     *
     * Previous code: malloc(sizeof(HDNA_HIT) * 1 000 000) = ~52 MB
     * unconditionally on every call, even for tiny inputs.
     *
     * Opt 1.1 (revised): initial capacity = raw_len / 32, clamped to
     * [HITS_MIN_CAP, DEFAULT_MAX_HITS].  For a 6 MB sequence this gives
     * ~187 000 entries (9.7 MB) which covers typical hit densities
     * (5–20 hits/KB) in a single findHDNA_core pass, eliminating the
     * re-scan overhead that previously caused 3–4× slowdown on genome-
     * scale inputs.  For short sequences the floor (HITS_MIN_CAP = 4096)
     * keeps the allocation negligible.
     *
     * Opt 1.3: PyMem_Raw* instead of PyMem_*; see note below.
     * ---------------------------------------------------------------- */
#define HITS_MIN_CAP   4096    /*  4 k × ~52 bytes ≈  208 KB — short seqs  */
#define HITS_GROW_FAC  2

    /* length / 32 ≈ 1 entry per 32 bp — covers 5–20 hits/KB in one pass */
    int hits_cap = (int)((size_t)raw_len >> 5);
    if (hits_cap < HITS_MIN_CAP)    hits_cap = HITS_MIN_CAP;
    if (hits_cap > DEFAULT_MAX_HITS) hits_cap = DEFAULT_MAX_HITS;
    /* Opt 1.3 note: hits uses PyMem_Raw* (malloc/realloc/free directly)
     * because it is reallocated inside Py_BEGIN_ALLOW_THREADS where the
     * GIL is released.  PyMem_Realloc internally calls
     * _PyInterpreterState_GET() which dereferences the thread state;
     * that pointer is NULL without the GIL, causing a SEGV on CPython
     * 3.12+.  PyMem_Raw* bypasses the interpreter state entirely. */
    HDNA_HIT *hits = (HDNA_HIT *)PyMem_RawMalloc(
                         sizeof(HDNA_HIT) * (size_t)hits_cap);
    if (!hits) { PyMem_Free(dna); PyErr_NoMemory(); return NULL; }

    /* ---- run C computation (GIL released) ---- */
    int nhits;
    int hit_cap_reached = 0;

    Py_BEGIN_ALLOW_THREADS
        /* Grow the hits buffer if findHDNA_core fills it, then re-scan.
         * Re-scanning is free because findHDNA_core is deterministic and
         * the sequence is already lowercased.  In practice re-scanning
         * only occurs for pathologically repetitive sequences. */
        for (;;) {
            nhits = findHDNA_core(
                dna, (int)raw_len,
                hits, hits_cap,
                minrep, maxrep, maxspacer,
                (float)purity, (float)mismatch,
                seq_offset);

            if (nhits < hits_cap || hits_cap >= DEFAULT_MAX_HITS)
                break; /* buffer was sufficient, or hard cap reached */

            /* buffer was too small — grow and retry */
            int new_cap = hits_cap * HITS_GROW_FAC;
            if (new_cap > DEFAULT_MAX_HITS) new_cap = DEFAULT_MAX_HITS;
            HDNA_HIT *tmp = (HDNA_HIT *)PyMem_RawRealloc(
                                hits, sizeof(HDNA_HIT) * (size_t)new_cap);
            if (!tmp) break; /* allocation failed: keep what we have */
            hits     = tmp;
            hits_cap = new_cap;
        }

        if (nhits >= DEFAULT_MAX_HITS) hit_cap_reached = 1;

        if (do_overlap && nhits > 1)
            nhits = remove_overlaps_core(hits, nhits);
    Py_END_ALLOW_THREADS

#undef HITS_MIN_CAP
#undef HITS_GROW_FAC

    /* Opt 1.2: full_sequence needs a contiguous buffer (concatenation of
     * up to 3 slices).  The other arm/spacer strings are built directly
     * from dna[] in the loop below — no heap allocation required.
     *
     * Maximum full_sequence length = 2*maxrep + maxspacer.
     * With the upper-bound validation (opt 1.4) maxrep <= 100 000 and
     * maxspacer <= 100 000, so the worst-case is 200 001 bytes.  We
     * allocate this once here and reuse it across all hits.
     *
     * Overflow-safe size computation (opt 1.4): validated parameters
     * guarantee maxrep <= 100 000, so 2*maxrep+maxspacer+2 <= 300 002,
     * well within size_t range on any supported platform.
     */
    size_t full_sz = (size_t)maxrep * 2 + (size_t)maxspacer + 2;
    char  *full_seq_buf = (char *)PyMem_Malloc(full_sz);
    if (!full_seq_buf) {
        PyMem_RawFree(hits); PyMem_Free(dna);
        PyErr_NoMemory();
        return NULL;
    }

    /* ---- build Python result list ---- */
    PyObject *result = PyList_New((Py_ssize_t)nhits);
    if (!result) goto cleanup_err;

    for (int i = 0; i < nhits; i++) {
        HDNA_HIT *h = &hits[i];

        /* Opt 1.2: build full_sequence directly in full_seq_buf from dna[]
         * slices — no intermediate heap buffers, no strncpy per arm.
         * left_arm and right_arm Python strings are created directly from
         * dna[] using PyUnicode_FromStringAndSize (one copy, not two). */
        int arm  = h->arm_len;
        int sp   = h->spacer_len;
        int total_len = arm * 2 + sp;

        /* Assemble full_sequence in the reusable stack-like buffer.
         * Buffer is sized to 2*maxrep + maxspacer + 2 (allocated above). */
        char *p = full_seq_buf;
        memcpy(p, dna + h->left_start_idx,  (size_t)arm); p += arm;
        if (sp > 0) {
            memcpy(p, dna + h->spacer_start_idx, (size_t)sp); p += sp;
        }
        memcpy(p, dna + h->right_start_idx, (size_t)arm); p += arm;
        *p = '\0';

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
        _SET("arm_length",      PyLong_FromLong((long)arm));
        _SET("spacer_length",   PyLong_FromLong((long)sp));
        _SET("total_length",    PyLong_FromLong((long)total_len));
        _SET("ga_pct",          PyFloat_FromDouble((double)h->ga_pct));
        _SET("ct_pct",          PyFloat_FromDouble((double)h->ct_pct));
        _SET("mirror_identity", PyFloat_FromDouble((double)h->mirror_id));
        _SET("is_perfect",      PyBool_FromLong((long)h->is_perfect));
        /* Opt 1.2: direct slice into dna[] — one copy instead of two */
        _SET("left_arm",  PyUnicode_FromStringAndSize(
                              dna + h->left_start_idx,  (Py_ssize_t)arm));
        _SET("spacer",    (sp > 0)
                              ? PyUnicode_FromStringAndSize(
                                    dna + h->spacer_start_idx, (Py_ssize_t)sp)
                              : PyUnicode_FromStringAndSize(".", 1));
        _SET("right_arm", PyUnicode_FromStringAndSize(
                              dna + h->right_start_idx, (Py_ssize_t)arm));
        _SET("full_sequence", PyUnicode_FromStringAndSize(
                              full_seq_buf, (Py_ssize_t)total_len));

#undef _SET

        PyList_SET_ITEM(result, (Py_ssize_t)i, d); /* steals ref to d */
    }

    /* free working buffers before potentially triggering a Python warning
     * Opt 1.3: hits uses PyMem_RawFree (matches PyMem_RawMalloc above);
     * dna and full_seq_buf use PyMem_Free (matches PyMem_Malloc). */
    PyMem_Free(full_seq_buf);
    PyMem_RawFree(hits);
    PyMem_Free(dna);

    if (hit_cap_reached) {
        if (PyErr_WarnFormat(PyExc_RuntimeWarning, 1,
                "scan_sequence: hit capacity reached (DEFAULT_MAX_HITS=%d); "
                "results may be incomplete. Consider splitting the input "
                "or recompiling with a larger DEFAULT_MAX_HITS.",
                DEFAULT_MAX_HITS) < 0) {
            Py_DECREF(result);
            return NULL;
        }
    }

    return result;

cleanup_err:
    /* Opt 1.3: hits uses PyMem_RawFree; others use PyMem_Free.
     * PyMem_RawFree(NULL) and PyMem_Free(NULL) are both no-ops. */
    PyMem_Free(full_seq_buf);
    PyMem_RawFree(hits);
    PyMem_Free(dna);
    return NULL;
}

/* ------------------------------------------------------------------ */
/*  Module definition                                                  */
/* ------------------------------------------------------------------ */

PyDoc_STRVAR(scan_sequence_doc,
"scan_sequence(seq, *, minrep=6, maxrep=50, maxspacer=7,\n"
"              purity=0.80, mismatch=0.20, remove_overlaps=True,\n"
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

static PyMethodDef HdnaMethods[] = {
    {
        "scan_sequence",
        (PyCFunction)py_scan_sequence,
        METH_VARARGS | METH_KEYWORDS,
        scan_sequence_doc,
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
    /* Opt 2.2: populate tolower lookup table once at module load */
    for (int i = 0; i < 256; i++)
        lc_table[i] = (unsigned char)tolower(i);

    PyObject *m = PyModule_Create(&hdna_module);
    if (!m) return NULL;
    PyModule_AddIntConstant(m, "DEFAULT_MAX_HITS", DEFAULT_MAX_HITS);
    PyModule_AddStringConstant(m, "__version__", "0.1.0");
    return m;
}

/* ================================================================== */
/*  STANDALONE CLI mode  (gcc -DSTANDALONE _hdna.c -lm -o findHDNA)  */
/* ================================================================== */
#else  /* STANDALONE */

static char      sa_dna[MAX_DNA_STANDALONE + 1];
static HDNA_HIT  sa_hits[MAX_HITS_STANDALONE];
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
                sa_dna[n++] = (char)lc_table[(unsigned char)base];
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
        "  -minrep    <int>    Minimum arm length       (default: 6)\n"
        "  -maxrep    <int>    Maximum arm length       (default: 50)\n"
        "  -maxspacer <int>    Maximum spacer length    (default: 7)\n"
        "  -purity    <float>  Min GA or CT fraction    (default: 0.80)\n"
        "  -mismatch  <float>  Max mismatch fraction    (default: 0.20)\n"
        "  -skipoverlap        Skip overlap removal\n"
        "  -v                  Verbose\n\n", prog);
}

int main(int argc, char *argv[]) {
    /* Opt 2.2: populate tolower lookup table once before first use */
    for (int i = 0; i < 256; i++)
        lc_table[i] = (unsigned char)tolower(i);

    char seq_fn[512] = {0}, out_pre[512] = {0}, out_fn[520] = {0};
    int  minrep = 6, maxrep = 50, maxspacer = 7, do_overlap = 1, verbose = 0;
    float purity = 0.80f, mismatch = 0.20f;

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
            sa_dna, bases, sa_hits, MAX_HITS_STANDALONE,
            minrep, maxrep, maxspacer, purity, mismatch, sa_seq_offset);
        if (do_overlap && nhits > 1)
            nhits = remove_overlaps_core(sa_hits, nhits);
        if (verbose) fprintf(stderr, "  %d hits\n", nhits);
        sa_print_hits(fo, nhits);
        total_hits += nhits;
    }

    fclose(fq); fclose(fo);
    fprintf(stderr,
        "\nDone.\n  Records: %d\n  Hits: %d\n  Output: %s\n",
        total_records, total_hits, out_fn);
    return 0;
}

#endif /* STANDALONE */
