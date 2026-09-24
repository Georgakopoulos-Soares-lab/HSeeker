# RTR revision — session checkpoint

**Date:** 2026-09-24
**Branch:** `RTR`
**Benchmark:** `hdna_benchmark_balanced_v3.csv` (130 records, 65:65), sha256 `f4602d6453d0926c1eae10dd54311e9fad08afe2beb9e59a3a9a0153e3b2ad8c`
**Tests:** 150 passed, 3 skipped (chr1 regression tests skip without `benchmarks/data/chr1.fa`)

This file records where the revision stands so a fresh session can resume without
re-deriving anything. The machine-readable state is
`RTR/hseeker_revision_ledger.yaml`; this is the human summary.

---

## 1. What this branch contains

The Biomni planning session's deliverables (plan, ledger, dependency graph,
synthetic-negative generator), the curated benchmark, and the work executed since:
five repository fixes, a parameter-robustness study, a family-aware validation study,
and two comparator builds.

**Benchmark v3 is now the repository's benchmark.** Both analysis scripts default to
it, filter to `curation_decision == 'kept'`, accept
`--subset all|experimental|synthetic`, and record `input_csv_sha256` per run.

---

## 2. Status — 11 complete, 5 part-done, 18 pending

### Complete and verified

| action | item | outcome |
|---|---|---|
| WS1-A1 | R2.8 | `-march=native` now opt-in via `HSEEKER_NATIVE=1`; CI wheels built portably |
| WS1-A2 | R2.9 | Profiling counters made thread-local. Measured pre-fix corruption: 8 concurrent scans of identical input gave 8 different totals, worst **3.84× inflated** |
| WS1-A3 | R2.10 | Sequence fields normalised to uppercase at the Python boundary; C core untouched |
| WS1-A4 | R2.11 | **Streamlit removed entirely** (author decision superseded the planned rename); one front-end remains |
| WS1-A5 | R2.2 | Stale "54 sequences" docstring fixed; notebook's missing-input failure removed |
| AUDIT-1, CUR-1, SYN-1/2/3, WS2-A4 | R3.5, R2.1 | Provenance audit and benchmark v3 curation, re-verified in this repo |

### Part-done (evidence banked, not finished)

| action | item | what exists / what's missing |
|---|---|---|
| WS2-A2 | R2.7, R3.1 | Wilson + record-bootstrap + **family-cluster** CIs implemented. Missing: tool-vs-tool difference CIs (no R/Triplex here) |
| WS2-A3 | R3.2, R3.3 | LOFO + LOSO + both collapsed variants done for the **detection call**. Missing: the stability-score threshold |
| WS4-A3 | R3.9 | Full-factorial 576-config sweep done. Missing: pairing-only vs pairing+stacking ablation (WS4-A2) |
| WS3-A1 | R3.6 | NeSSie builds (commit `dbe6cb2`); interface characterised. Missing: the benchmark run |
| WS3-A2 | R3.6 | non-B_gfa builds (commit `a891b59`). Missing: the benchmark run |

---

## 3. Headline result

Under **every** validation scheme the winner is the configuration the manuscript
already states — `mismatch 0.15, minrep 8, purity 0.90, maxspacer 10, both filters`:

| scheme | sens | spec | MCC | record CI | family-cluster CI |
|---|---|---|---|---|---|
| leave-one-family-out | 0.954 | 0.985 | **0.939** | [0.877, 0.985] | [0.845, 0.988] |
| leave-one-study-out | 0.954 | 0.985 | 0.939 | [0.877, 0.985] | [0.836, 0.988] |
| sequence-level 5-fold CV | 0.954 | 0.985 | 0.939 | [0.877, 0.985] | — |
| collapsed per (family, class) | 0.933 | 0.966 | 0.889 | [0.772, 0.973] | [0.772, 0.974] |

It beats every tuned alternative (per-fold tuned 0.895, script default 0.908,
production CLI 0.856) — **and requires no selection at all**, so R3.3's circularity
objection does not apply to it.

### Why we can say HSeeker is not overfitted (R2.6)

1. Tuning on the benchmark makes results **worse** (0.939 pre-specified vs 0.895 tuned).
2. Over 200 repeated splits, tuning's median gain is exactly **+0.0000**.
3. **108 of 576** configurations lie within 0.05 MCC of the maximum — a plateau, not a peak.
4. `minrep 8` and `purity 0.85` selected in **54/54** family folds.
5. **Nothing is fitted.** LOFO, LOSO and full-data evaluation are numerically identical
   (TP62 FN3 TN64 FP1), because with fixed parameters the call depends only on the
   record's own sequence.

**Honest residual:** the AT/homopolymer thresholds *were* chosen with knowledge of
this benchmark's composition. That is where R2.6 retains genuine force; it must be
disclosed, not claimed as derived.

### On R3.2 specifically

Family-aware validation does not move the point estimate — it moves the *precision*.
The family-cluster CI is **1.26×** wider than the record-level CI. Both are reported
everywhere: the record-level interval answers R2.7, the cluster interval answers
R3.2's "inflate apparent statistical power".

---

## 4. Blocked on three author decisions

| decision | blocks | evidence |
|---|---|---|
| **DEC-2** FXN citation | WS0-A1 → 21 actions | `RTR/evidence/DEC-2_fxn_citation.md` |
| **DEC-3** add Hanvey 1988 `(TTC)8` | WS0-A1 → 21 actions | `RTR/evidence/DEC-3_pRW1406.md` |
| **DEC-5** canonical parameters | WS2-A1 → 16 actions | now answered empirically by §3 |

- **DEC-2:** the cited PMID (Romney 2008, a *C. elegans* ferritin paper — likely an
  FTN-1/FXN symbol confusion) and the cited DOI (Samuel 2007 innate-immunity
  minireview) are two *different* wrong papers. Both candidate replacements are
  reviews. No primary source was found for n = 10/20/33/66; Sakamoto 1999 used 75–270
  repeats and requires >59, so the `forming` label on GAA10/20/33 is itself untraced.
- **DEC-3:** Hanvey 1988's abstract establishes both the sequence and the forming
  label for `(TTC)8`. Recommend adding it as `Hanvey1988_(TTC)8` rather than
  `pRW1406`, because the plasmid-number mapping comes from the document already shown
  corrupted for this study.
- **DEC-5:** §3 shows the manuscript's values win out-of-fold under every scheme.

**WS0-A1 is not finished** even once decided: the 22 PMID and 2 DOI corrections sit
in `corrected_*` columns without being applied to `cited_*`, and the pCON terminal-C
fix (register D-9) is not applied — HDNA0029 is still 31 nt.

---

## 5. Next steps (planned, not started)

Four items need no decision. Two verified facts collapse the work:
`pairing_score`/`stacking_score` are already emitted per hit (`total == pairing +
stacking` exactly), and `--overlap-fraction` is consumed *after* the scan.

| # | item | ask | new compute |
|---|---|---|---|
| 1 | Score threshold out-of-fold | R3.3 | ~0 — one 27 ms scan, reuse fold machinery |
| 2 | Scoring ablation | R3.9 | **0** — re-threshold on components already returned |
| 3 | Overlap-selection rule | R3.10 | one chr1 scan, **shared with WS5** |
| 4 | Overlap-threshold sweep | R3.7 | **one** E. coli scan (genome already cached); sweep post-hoc |

**Item 1 is the priority** — it is the last live circularity, and the one R3.3 named.

### WS5 data is fetched and validated (not committed — 208 MB)

In the session scratchpad: S1-END-seq KM12 rep1 (pos+neg) plus the **matched END-seq
no-S1 control**, from GSE204808/GSE203632 (which are the same file set).

Three decisions the data forced:

1. **Everything is hg19**, verified on three files. Our chr1 predictions are hg38, so
   the plan's "human, so no liftover" is wrong about build. → **Re-run chr1 on hg19**
   rather than liftOver, which would distort exactly the repetitive regions H-DNA
   occupies.
2. **No called peaks exist** — only coverage tracks. → Use **continuous signal at
   loci vs composition-matched controls**, avoiding an invented peak threshold. R3.8
   explicitly accepts an enrichment analysis.
3. **A paired END-seq control exists for KM12** → use it, isolating S1-specific
   structural signal from DSB/accessibility background.

Risks: the chr1 `remove_overlaps=False` run may exhaust memory (cf. the OOM fix in
`cabf2f0`) — stream with `scan_fasta_iter`. And replicating the longest-arm rule in
Python must exactly reproduce production output before any comparison is reported.

---

## 6. Blocked on inputs we do not have

- **The manuscript and RTR `.docx` are not in this repository.** WS6-A2 (Limitations,
  R2.4/R3.4/R2.6), WS6-A3 (R2.5), WS6-A4 (R3.11), WS6-A5 (R3.12) cannot start.
  `algorithms.tex` is pseudocode panels only.
- **R / Bioconductor `triplex` is not installed** in this environment, so the direct
  arm has only ever been half-run and **no tool-vs-tool claim is currently supported
  by any run on benchmark v3.** This must be available on TACC.

---

## 7. Provenance note

The source archive `Biomni_lab_downloads_20260923_153955.zip` and
`biomni_history.txt` were present at the start of the session and are **no longer on
disk**; they were not deleted by any action taken here. All ten substantive files
they contained survive — eight in `RTR/`, two moved to the repository root — so
nothing is lost, but the originals are not recoverable from this branch.

Results in `analysis/results/` other than `robustness_v3/` and
`grouped_validation_v3/` are from **2026-06-24 against the old 80-row CSV** with
pre-filter code. They do not reflect benchmark v3 and should be treated as superseded.
