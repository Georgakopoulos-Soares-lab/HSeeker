/*
 * _hdna_purity.c — Optimization A: prefix-sum GA/CT arrays + per-right0
 * purity memoization with LRU cache.
 *
 * This module has NO dependency on Python.h and NO dependency on the SA/LCE
 * code — it can be compiled and tested independently.
 *
 * Key data structures:
 *   PurityArrays     : prefix sums over the chunk buffer (built once per chunk)
 *   PerRight0Table   : best_purity_k[b] for a fixed right0 (built on demand)
 *   PurityCache      : small LRU (ring-buffer eviction) of PerRight0Tables
 *
 * All of these are declared in _hdna_dp.h.
 */

#include <stdlib.h>
#include <string.h>
#include "_hdna_dp.h"

/* The BT table is defined in _hdna.c and must be declared here as extern.
 * When compiling as a standalone test, define it locally. */
#ifndef HSEEKER_HDNA_C_INCLUDED
extern const unsigned char BT[256];
#endif

/* ── Profiling counters ──────────────────────────────────────────────────── */

long long prof_lce_queries         = 0;
long long prof_purity_cache_hits   = 0;
long long prof_purity_cache_misses = 0;
long long prof_chunk_fallback_count = 0;

void prof_dp_reset(void) {
    prof_lce_queries          = 0;
    prof_purity_cache_hits    = 0;
    prof_purity_cache_misses  = 0;
    prof_chunk_fallback_count = 0;
}

/* ── Feature flag ────────────────────────────────────────────────────────── */

int hdna_dp_enabled(void) {
    const char *val = getenv("HSEEKER_FAST_HDNA");
    return (val != NULL && val[0] == '1' && val[1] == '\0');
}

/* ── next_bad array ──────────────────────────────────────────────────────── */

int32_t *next_bad_build(const char *dna, int n, const unsigned char *BT_table)
{
    int32_t *nb = (int32_t *)malloc(sizeof(int32_t) * (size_t)(n + 1));
    if (!nb) return NULL;
    nb[n] = n;  /* sentinel */
    int last_bad = n;
    for (int i = n - 1; i >= 0; i--) {
        if (BT_table[(unsigned char)dna[i]] == 0)
            last_bad = i;
        nb[i] = last_bad;
    }
    return nb;
}

/* ── PurityArrays ────────────────────────────────────────────────────────── */

PurityArrays *purity_arrays_build(const char *dna, int n)
{
    PurityArrays *pa = (PurityArrays *)malloc(sizeof(PurityArrays));
    if (!pa) return NULL;

    pa->n         = n;
    pa->prefix_ga = (int32_t *)malloc(sizeof(int32_t) * (size_t)(n + 1));
    pa->prefix_ct = (int32_t *)malloc(sizeof(int32_t) * (size_t)(n + 1));
    if (!pa->prefix_ga || !pa->prefix_ct) {
        free(pa->prefix_ga);
        free(pa->prefix_ct);
        free(pa);
        return NULL;
    }

    pa->prefix_ga[0] = 0;
    pa->prefix_ct[0] = 0;
    for (int i = 0; i < n; i++) {
        unsigned char bt = BT[(unsigned char)dna[i]];
        pa->prefix_ga[i + 1] = pa->prefix_ga[i] + (bt == 1 ? 1 : 0);
        pa->prefix_ct[i + 1] = pa->prefix_ct[i] + (bt == 2 ? 1 : 0);
    }
    return pa;
}

void purity_arrays_free(PurityArrays *pa)
{
    if (!pa) return;
    free(pa->prefix_ga);
    free(pa->prefix_ct);
    free(pa);
}

/* ── PerRight0Table ──────────────────────────────────────────────────────── */

PerRight0Table *per_right0_build(
    const PurityArrays *pa,
    int right0, int maxrep,
    int purity_thresh_int)
{
    int K = pa->n - right0;
    if (K > maxrep) K = maxrep;
    if (K <= 0) {
        /* right0 at or past end of sequence: empty table */
        PerRight0Table *t = (PerRight0Table *)malloc(sizeof(PerRight0Table));
        if (!t) return NULL;
        t->best_k = NULL;
        t->K      = 0;
        t->right0 = right0;
        return t;
    }

    PerRight0Table *t = (PerRight0Table *)malloc(sizeof(PerRight0Table));
    if (!t) return NULL;
    t->right0 = right0;
    t->K      = K;
    t->best_k = (int32_t *)malloc(sizeof(int32_t) * (size_t)(K + 1));
    if (!t->best_k) {
        free(t);
        return NULL;
    }

    /* Backward pass: best_k[b] = largest k in [1,b] with purity true. */
    t->best_k[0] = 0;
    int running_best = 0;
    /* cumulative counts as we scan forward k=1..K */
    /* We need ga(right0,k) and ct(right0,k) for each k — use prefix sums */
    for (int k = 1; k <= K; k++) {
        int ga_k = purity_ga(pa, right0, k);
        int ct_k = purity_ct(pa, right0, k);
        int ga_ok = ga_k * 100 >= purity_thresh_int * k;
        int ct_ok = ct_k * 100 >= purity_thresh_int * k;
        if (ga_ok || ct_ok)
            running_best = k;
        t->best_k[k] = running_best;
    }
    return t;
}

void per_right0_free(PerRight0Table *t)
{
    if (!t) return;
    free(t->best_k);
    free(t);
}

/* ── PurityCache ─────────────────────────────────────────────────────────── */

PurityCache *purity_cache_create(int capacity)
{
    PurityCache *c = (PurityCache *)malloc(sizeof(PurityCache));
    if (!c) return NULL;
    c->entries = (PerRight0Table **)calloc((size_t)capacity, sizeof(PerRight0Table *));
    if (!c->entries) {
        free(c);
        return NULL;
    }
    c->capacity    = capacity;
    c->count       = 0;
    c->next_evict  = 0;
    return c;
}

void purity_cache_free(PurityCache *c)
{
    if (!c) return;
    for (int i = 0; i < c->capacity; i++)
        per_right0_free(c->entries[i]);
    free(c->entries);
    free(c);
}

const PerRight0Table *purity_cache_get(
    PurityCache *cache,
    const PurityArrays *pa,
    int right0, int maxrep,
    int purity_thresh_int)
{
    /* Linear scan: cache is tiny (≤2*maxspacer+4 ≈ 24 entries). */
    for (int i = 0; i < cache->count; i++) {
        if (cache->entries[i] && cache->entries[i]->right0 == right0) {
            prof_purity_cache_hits++;
            return cache->entries[i];
        }
    }
    prof_purity_cache_misses++;

    /* Build new table */
    PerRight0Table *t = per_right0_build(pa, right0, maxrep, purity_thresh_int);
    if (!t) return NULL;  /* OOM — caller must fall back to legacy */

    /* Evict oldest entry (ring-buffer) */
    int slot = cache->next_evict;
    per_right0_free(cache->entries[slot]);
    cache->entries[slot] = t;
    cache->next_evict = (slot + 1) % cache->capacity;
    if (cache->count < cache->capacity)
        cache->count++;

    return t;
}
