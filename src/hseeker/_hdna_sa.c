/*
 * _hdna_sa.c — Suffix array, Kasai LCP, sparse-table RMQ, and O(1) LCE query.
 *
 * This module is GENERIC — no HSeeker-specific concept (arms, spacers,
 * purity) appears here.  It is independently testable against a brute-force
 * O(n) LCE reference.
 *
 * The combined string is:
 *   S = dna[0..C) + SEP_BYTE + reverse(dna[0..C))
 *   len(S) = M = 2*C + 1
 *
 * SEP_BYTE is 0x01 — a byte that cannot appear in the lowercased DNA input
 * (which is restricted to a-z after tolower).  We assert this at build time.
 *
 * Suffix array construction: O(M log M) prefix-doubling (Manber & Myers).
 * LCP array: Kasai's O(M) algorithm.
 * Sparse table: O(M log M) build, O(1) RMQ.
 * LCE(a, b): O(1) via RMQ over LCP[rank[a]+1 .. rank[b]].
 *
 * Memory:
 *   SA:     M * sizeof(int32_t)
 *   rank:   M * sizeof(int32_t)
 *   LCP:    M * sizeof(int32_t)
 *   sparse: M * log2(M) * sizeof(int32_t)
 *   Total:  see _hdna_dp.h formula; ~175 MB for C=1e6.
 */

#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "_hdna_dp.h"

#define SEP_BYTE ((unsigned char)0x01)

/* ── Internal structs ────────────────────────────────────────────────────── */

struct SaRmqHandle {
    int32_t  *sa;       /* suffix array, length M */
    int32_t  *rank;     /* inverse of sa, length M */
    int32_t  *lcp;      /* LCP array, length M */
    int32_t **sparse;   /* sparse[j][i] = min LCP in [i, i+2^j), length log2(M)+1 rows */
    int       M;        /* 2*C+1 */
    int       C;        /* chunk size (length of dna half) */
    int       log2M;    /* floor(log2(M)) */
};

/* ── Projected memory calculator ────────────────────────────────────────── */

static size_t projected_bytes(int C)
{
    int M = 2 * C + 1;
    int log2M = 0;
    int tmp = M;
    while (tmp > 1) { log2M++; tmp >>= 1; }
    log2M++;  /* +1 for safety */

    size_t sa_bytes     = (size_t)M * 4;
    size_t rank_bytes   = (size_t)M * 4;
    size_t lcp_bytes    = (size_t)M * 4;
    size_t sparse_bytes = (size_t)M * (size_t)log2M * 4;
    return sa_bytes + rank_bytes + lcp_bytes + sparse_bytes;
}

/* ── Prefix-doubling suffix array construction ───────────────────────────── */

typedef struct { int32_t key1, key2, idx; } SaEntry;

static int sa_entry_cmp(const void *a, const void *b)
{
    const SaEntry *ea = (const SaEntry *)a;
    const SaEntry *eb = (const SaEntry *)b;
    if (ea->key1 != eb->key1) return (ea->key1 < eb->key1) ? -1 : 1;
    if (ea->key2 != eb->key2) return (ea->key2 < eb->key2) ? -1 : 1;
    return 0;
}

/* Build suffix array for string s of length M.
 * Returns 0 on success, -1 on OOM.
 * rank[] is filled with the inverse (rank of each suffix).
 */
static int build_sa(const unsigned char *s, int M, int32_t *sa, int32_t *rank)
{
    SaEntry *entries = (SaEntry *)malloc(sizeof(SaEntry) * (size_t)M);
    if (!entries) return -1;

    /* Initial rank = character value */
    for (int i = 0; i < M; i++) {
        sa[i]       = i;
        rank[i]     = (int32_t)s[i];
        entries[i].idx  = i;
        entries[i].key1 = (int32_t)s[i];
        entries[i].key2 = (i + 1 < M) ? (int32_t)s[i + 1] : -1;
    }

    for (int gap = 2; gap < M * 2; gap <<= 1) {
        qsort(entries, (size_t)M, sizeof(SaEntry), sa_entry_cmp);

        /* Assign new dense ranks */
        rank[entries[0].idx] = 0;
        for (int i = 1; i < M; i++) {
            rank[entries[i].idx] = rank[entries[i-1].idx];
            if (entries[i].key1 != entries[i-1].key1 ||
                entries[i].key2 != entries[i-1].key2)
                rank[entries[i].idx]++;
        }
        /* All ranks unique → done */
        if (rank[entries[M-1].idx] == M - 1) break;

        /* Update key2 using the new ranks */
        for (int i = 0; i < M; i++) {
            int j = entries[i].idx;
            entries[i].key1 = rank[j];
            int j2 = j + gap;
            entries[i].key2 = (j2 < M) ? rank[j2] : -1;
        }
    }

    /* Extract sorted SA from entries */
    for (int i = 0; i < M; i++)
        sa[i] = entries[i].idx;

    /* Recompute rank as true inverse of SA */
    for (int i = 0; i < M; i++)
        rank[sa[i]] = i;

    free(entries);
    return 0;
}

/* ── Kasai's LCP construction ─────────────────────────────────────────────── */

static int build_lcp(const unsigned char *s, int M,
                     const int32_t *sa, const int32_t *rank,
                     int32_t *lcp)
{
    lcp[0] = 0;
    int h = 0;
    for (int i = 0; i < M; i++) {
        int r = rank[i];
        if (r == 0) { h = 0; continue; }
        int j = sa[r - 1];
        while (i + h < M && j + h < M && s[i + h] == s[j + h])
            h++;
        lcp[r] = h;
        if (h > 0) h--;
    }
    return 0;
}

/* ── Sparse table for O(1) range-minimum query ───────────────────────────── */

static int build_sparse(const int32_t *lcp, int M, int32_t **sparse, int log2M)
{
    /* sparse[0][i] = lcp[i] */
    for (int i = 0; i < M; i++)
        sparse[0][i] = lcp[i];

    for (int j = 1; j <= log2M; j++) {
        int half = 1 << (j - 1);
        for (int i = 0; i + (1 << j) <= M; i++) {
            int32_t a = sparse[j-1][i];
            int32_t b = sparse[j-1][i + half];
            sparse[j][i] = (a < b) ? a : b;
        }
    }
    return 0;
}

/* O(1) range minimum query on sparse table */
static int32_t rmq(const SaRmqHandle *h, int l, int r)
{
    /* l and r are inclusive indices into the LCP array */
    if (l > r) return 0;
    int len = r - l + 1;
    int k   = 0;
    while ((1 << (k + 1)) <= len) k++;
    int32_t a = h->sparse[k][l];
    int32_t b = h->sparse[k][r - (1 << k) + 1];
    return (a < b) ? a : b;
}

/* ── Public API ──────────────────────────────────────────────────────────── */

SaRmqHandle *sa_rmq_build(const char *dna, int C, int *out_fallback)
{
    *out_fallback = 0;

    /* Memory-cap check */
    size_t projected = projected_bytes(C);
    if (projected > DP_MAX_CHUNK_BYTES) {
        prof_chunk_fallback_count++;
        *out_fallback = 1;
        return NULL;
    }

    int M = 2 * C + 1;
    assert(M < (int)2e9);  /* guard against INT32 overflow in SA indices */

    /* Build the combined string S = dna + SEP + rev(dna) */
    unsigned char *S = (unsigned char *)malloc((size_t)M);
    if (!S) return NULL;
    for (int i = 0; i < C; i++) {
        unsigned char b = (unsigned char)dna[i];
        /* Assert SEP_BYTE cannot appear in the input */
        assert(b != SEP_BYTE);
        S[i] = b;
    }
    S[C] = SEP_BYTE;
    for (int i = 0; i < C; i++)
        S[C + 1 + i] = (unsigned char)dna[C - 1 - i];

    SaRmqHandle *h = (SaRmqHandle *)calloc(1, sizeof(SaRmqHandle));
    if (!h) { free(S); return NULL; }
    h->M = M;
    h->C = C;

    /* Compute log2(M) */
    int log2M = 0;
    { int tmp = M; while (tmp > 1) { log2M++; tmp >>= 1; } }
    h->log2M = log2M;

    /* Allocate SA, rank, LCP */
    h->sa   = (int32_t *)malloc(sizeof(int32_t) * (size_t)M);
    h->rank = (int32_t *)malloc(sizeof(int32_t) * (size_t)M);
    h->lcp  = (int32_t *)malloc(sizeof(int32_t) * (size_t)M);
    if (!h->sa || !h->rank || !h->lcp) goto oom;

    /* Allocate sparse table rows */
    h->sparse = (int32_t **)calloc((size_t)(log2M + 1), sizeof(int32_t *));
    if (!h->sparse) goto oom;
    for (int j = 0; j <= log2M; j++) {
        h->sparse[j] = (int32_t *)malloc(sizeof(int32_t) * (size_t)M);
        if (!h->sparse[j]) goto oom;
    }

    /* Build */
    if (build_sa(S, M, h->sa, h->rank) != 0) goto oom;
    build_lcp(S, M, h->sa, h->rank, h->lcp);
    build_sparse(h->lcp, M, h->sparse, log2M);

    free(S);
    return h;

oom:
    free(S);
    sa_rmq_free(h);
    return NULL;
}

void sa_rmq_free(SaRmqHandle *h)
{
    if (!h) return;
    if (h->sparse) {
        for (int j = 0; j <= h->log2M; j++)
            free(h->sparse[j]);
        free(h->sparse);
    }
    free(h->sa);
    free(h->rank);
    free(h->lcp);
    free(h);
}

int sa_lce(const SaRmqHandle *h, int a, int b)
{
    if (a == b) return h->M - a;  /* same suffix — degenerate */
    int ra = h->rank[a];
    int rb = h->rank[b];
    if (ra > rb) { int tmp = ra; ra = rb; rb = tmp; }
    /* LCE = min(LCP[ra+1 .. rb]) */
    return (int)rmq(h, ra + 1, rb);
}

int sa_chunk_size(const SaRmqHandle *h)
{
    return h->C;
}
