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

                    /* early exit: mismatch budget exhausted even at maxrep */
                    if (mismatches > (int)(mismatch_tol * maxrep)) break;
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

/*
 * remove_overlaps_core
 *
 * For overlapping hits, keep the one with the longer arm; break ties by
 * keeping the shorter spacer.  Identical to the original remove_overlaps().
 *
 * Simple O(n²) pass — adequate for the typical hit counts per sequence.
 * Marks removed hits with arm_len = -1, then compacts the array.
 *
 * Returns the compacted hit count.
 */
static int remove_overlaps_core(HDNA_HIT *hits, int nhits) {
    for (int i = 0; i < nhits; i++) {
        if (hits[i].arm_len < 0) continue;
        for (int j = i + 1; j < nhits; j++) {
            if (hits[j].arm_len < 0) continue;

            long si = hits[i].start, ei = hits[i].end;
            long sj = hits[j].start, ej = hits[j].end;

            if (si <= ej && sj <= ei) { /* intervals overlap */
                if (hits[i].arm_len > hits[j].arm_len) {
                    hits[j].arm_len = -1;
                } else if (hits[j].arm_len > hits[i].arm_len) {
                    hits[i].arm_len = -1;
                    break;
                } else {
                    /* equal arm length: prefer shorter spacer */
                    if (hits[i].spacer_len <= hits[j].spacer_len)
                        hits[j].arm_len = -1;
                    else {
                        hits[i].arm_len = -1;
                        break;
                    }
                }
            }
        }
    }

    int n = 0;
    for (int i = 0; i < nhits; i++) {
        if (hits[i].arm_len >= 0)
            hits[n++] = hits[i];
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

    /* hits array */
    HDNA_HIT *hits = (HDNA_HIT *)malloc(sizeof(HDNA_HIT) * DEFAULT_MAX_HITS);
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
    int hit_cap_reached = 0;

    Py_BEGIN_ALLOW_THREADS
        nhits = findHDNA_core(
            dna, (int)raw_len,
            hits, DEFAULT_MAX_HITS,
            minrep, maxrep, maxspacer,
            (float)purity, (float)mismatch,
            seq_offset);

        if (nhits >= DEFAULT_MAX_HITS) hit_cap_reached = 1;

        if (do_overlap && nhits > 1)
            nhits = remove_overlaps_core(hits, nhits);
    Py_END_ALLOW_THREADS

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
    free(left_arm_buf); free(spacer_buf);
    free(right_arm_buf); free(full_seq_buf);
    free(hits); free(dna);
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
        "  -minrep    <int>    Minimum arm length       (default: 6)\n"
        "  -maxrep    <int>    Maximum arm length       (default: 50)\n"
        "  -maxspacer <int>    Maximum spacer length    (default: 7)\n"
        "  -purity    <float>  Min GA or CT fraction    (default: 0.80)\n"
        "  -mismatch  <float>  Max mismatch fraction    (default: 0.20)\n"
        "  -skipoverlap        Skip overlap removal\n"
        "  -v                  Verbose\n\n", prog);
}

int main(int argc, char *argv[]) {
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
