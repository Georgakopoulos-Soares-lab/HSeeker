# RTR revision — session checkpoint

**Date:** 2026-10-09 (previous checkpoints 2026-10-08, 2026-10-07, 2026-10-02, 2026-09-24)
**Branch:** `RTR-merge-main`, committed (the merge of `origin/main` into `origin/RTR` is `cbb4713`). Code as of `8e781ca`: Nikol's `fix/hdna-maximal-representation-v3` merged in `94742b3` (DEC-11), webapp filters `013be97` (DEC-12), FXN citations `8e781ca` (DEC-2). Upstream PRs #19 (memory fix) and #20 (score-mode speed-up) are open; Nikol's branch is merged only here, not in `main`
**Benchmark:** `hdna_benchmark_balanced_v3.csv` (128 records, 64 forming / 64 non-forming — 1:1 again after DEC-8; 69 experimental + 59 synthetic), sha256 `c5177b2479a5da4414aefa3c6ac6251f5645868260d2b5cefe9189de182f1ae4` since DEC-2 (2026-10-09; citation metadata of 4 FXN rows only). Superseded: `7bb49544…5907cf9` (still recorded in `analysis/results/*/metadata.json`, refreshed at the final freeze), `7f78555d…5bde89` (124 records at 64:60). Experimental CSV sha256 `dfa31654559e366fd666e98c0476349d9e6cb51c4e0fe36813f8caddbae0e273`
**Overlap default:** greedy (longest arm) in the library, CLI and webapp, identical to Nikol's main (DEC-10, 2026-10-08; supersedes DEC-9's score-informed default). Score mode is an option (`overlap_strategy="score"` / `-overlap-strategy score`). The detector includes Nikol's inward extension and partition deduplication (DEC-11); hg38 chr1 greedy default 96,731 loci (published, pre-merge: 96,729)
**Tests:** 349 passed, 3 skipped (2026-10-09) on the merged tree (chr1 regression tests skip without `benchmarks/data/chr1.fa`)
**Configuration labels:** the MCC 0.969 result with no filters is the **benchmark/script configuration** (minrep 8, mismatch 0.10, no filters), not the software default. The shipped defaults (minrep 10) give MCC 0.882 (see §0, 2026-10-08)

This file records where the revision stands so a fresh session can resume without
re-deriving anything. The machine-readable state is
`RTR/hseeker_revision_ledger.yaml`; this is the human summary.

---

## 0. Update log

### 2026-10-09 — Nikol's detector branch merged (DEC-11); webapp filters (DEC-12); FXN citations (DEC-2); pCON kept (DEC-13); Kouzine replaces E. coli (DEC-14)

Author decisions, all by Kimon and checked by him:

1. **DEC-11: Nikol's `fix/hdna-maximal-representation-v3` (`75c3494`) is merged into
   `RTR-merge-main`** (merge commit `94742b3`). `_hdna.c` now moves exact-matching spacer-end
   pairs into both arms (purity and maximal representation rechecked) and emits each
   arm/spacer partition once, before overlap removal and also with `-skipoverlap`. The
   benchmark is unchanged (0 of 73,728 detection calls across the 576-configuration grid
   differ). The hg38 chr1 greedy default (minrep 10, library, 16 workers) goes from 96,729 to
   **96,731** loci; 96,729 stays as the published, pre-merge number. `-skipoverlap`
   candidate sets roughly halve (golden file: e.g. 861 → 423 hits). The golden file
   `tests/data/batching_golden.json.gz` was regenerated with the pre-batching code (`4784c4b`)
   built against the merged detector: the 10 overlap-removal configurations are unchanged and
   the 10 `-skipoverlap` configurations changed. A stray file in the branch,
   `test_no.tsv_HDNA.tsv`, was not merged. Tests: 349 passed, 3 skipped. Not merged upstream
   (no PR exists).
2. **chr1 on the merged code** (16 workers, AMD EPYC 7763, single run each; WS6-A4 evidence
   §5):

   | mode | loci | wall |
   |---|---|---|
   | library greedy (default) | 96,731 | 25.9 s |
   | library score | 98,222 | 42.6 s (about 60 s before the merge: deduplication removes redundant candidates) |
   | CLI greedy | 41,939 | 22.9 s |
   | CLI score | 43,137 | 30.9 s |

   Greedy vs score, library: 100 % of greedy loci are overlapped by a score hit; 80,932
   (83.7 %) have identical coordinates; 1,075 loci are score-only; the reported motif differs
   at 15,799 loci. CLI: 91.5 % identical, 923 score-only. Benchmark grid, greedy vs score: 0 of
   73,728 calls differ. One merged greedy run peaked at 1,353 MB RSS. That is a single
   measurement: re-confirm it in X6 and do not quote it as final.
3. **DEC-12: the webapp applies the CLI composition filters** (commit `013be97`). The AT
   threshold is a form field (default 0.80, range 0–1.5), the homopolymer filter is always on,
   and the job page shows both. Checked with the FastAPI TestClient: webapp hits equal CLI hits
   for the defaults and for `-at-threshold 1.5`. The About page's "matches CLI defaults"
   statement is now true for the filters (plan register D-19).
4. **DEC-2 resolved: FXN records** (commit `8e781ca`). HDNA0060–0062 ((GAA)10/20/33) now cite
   Potaman VN et al. 2004, Nucleic Acids Res 32(3):1224-31 (PMID 14978261, DOI
   10.1093/nar/gkh274): in supercoiled plasmids at pH 7.4, (GAA)9 forms a stable
   intramolecular H-DNA, (GAA)23 a family of H-DNAs, and (GAA)42 a bi-triplex at higher
   supercoiling. HDNA0063 ((GAA)66) cites Sakamoto N et al. 1999, Mol Cell 3(4):465-75 (PMID
   10230399, DOI 10.1016/s1097-2765(00)80474-8): R·R·Y triplex / sticky DNA for more than 59
   repeats. The records stay secondary (class-level evidence; the exact constructs were not
   tested). The original PMID 18024960 / DOI 10.1074/jbc.R700013200 stay in the `cited_*`
   fields as an audit trail; the recorded title matches only the review Rajeswari 2012 (PMID
   22750988). Only `corrected_pmid`, `corrected_doi`, `resolved_origin`, `evidence` and
   `curation_justification` changed, in those 4 rows of both CSVs, so no detection call
   changes. New sha256: experimental `dfa31654…caddbae0e273`, balanced `c5177b24…de182f1ae4`.
   `analysis/results/*/metadata.json` still record the previous balanced sha; they are
   refreshed at the final benchmark freeze. This also settles the audit's open question that
   the forming label of (GAA)10/20/33 had no primary support.
5. **DEC-13: pCON (HDNA0029) is kept exactly as published** in the submitted manuscript's
   Supplementary Table: 31 nt `TTGCTGAACTTGATTCGACTCAGATGCAGTG`, identical in the original
   `hdna_experimental_sequences_final.csv` and in v3. The audit's terminal-C correction
   (register D-9, plan A.3 item 4) is **not** applied. HSeeker calls pCON non-forming with or
   without the terminal C (0 hits under the script, manuscript and shipped configurations), so
   no result changes. Every claim that the pCON sequence was "fixed" or "corrected" was removed
   (plan A.3, WS0-A1, R2.2, R3.5, D-9; this file §4).
6. **DEC-14: the E. coli insertion experiment is replaced by a genome-wide comparison with
   Kouzine et al. 2017** (Cell Systems 4:344-356, PMID 28237796; the dataset Reviewer 2 names
   in R2.3). Data: SRA072844 (KMnO4/S1 ssDNA-seq), raw reads only, because the paper's
   supplementary peak tables are not reachable from here (not open access in Europe PMC;
   cell.com returns 403). We reprocess mouse LPS+IL4-activated B-cell ssDNA (13 runs), resting
   B-cell ssDNA (4 runs) and activated-B input/sonicated DNA (2 runs) on mm10, and human Raji
   ssDNA (7 runs) on hg38 (no liftover; covers R3.8 on a human genome from the same study).
   bowtie2 2.5.5, samtools 1.24, iGenomes prebuilt indexes; MAPQ ≥ 10. The pipeline is running
   at `/workspaces/.hseeker-pr/kouzine` (`pipeline.sh`, `logs/pipeline.log`). The analysis
   design (ssDNA signal at HSeeker loci vs composition-matched controls, normalised to input;
   fraction of HSeeker loci in ssDNA-enriched regions) is fixed in
   `RTR/evidence/X8_kouzine_analysis_design_2026-10-09.md` (commit `a8d52e4`, before any
   result). **Dropped:** X5 (E. coli rerun) and WS4-A1 (overlap-threshold sweep); the
   experiment leaves the manuscript and Figure 3 will be replaced. R3.7 now explains the
   removal and that the overlap criterion no longer applies. X8 / WS5-A1 is now Kouzine (mouse)
   + Raji (human); WS5-A2 is superseded; the S1-END-seq KM12 plan (§5 below) is optional and
   not pursued.
7. **R3.10:** the final draft response is `RTR/evidence/R3.10_overlap_response_2026-10-09.md`.
   Plan Part C now points to it with a short summary, and its two `TODO(authors)` items are
   resolved by that text: the longest-arm default is kept (justified by geometry, the
   heuristic score's calibration and runtime), score mode is offered, a mode is recommended for
   each use, and near-homopolymer edge hits in score mode are reported (they are negligible on
   chr1).
8. **Phase 1 of `RTR/REMAINING_WORK_PLAN_2026-10-07.md`:** 1.1 code freeze **done** (DEC-11,
   DEC-12); 1.2 benchmark freeze partly done (DEC-2 and pCON/DEC-13 resolved; DEC-3 and the
   homopolymer panel pending); 1.3 configuration freeze pending (DEC-5, reporting
   configuration / minrep, comparator arguments); 1.4 genome-validation dataset **done**
   (DEC-14, Kouzine).
9. **Bookkeeping.** Ledger: DEC-2 resolved; DEC-11–DEC-14 added; WS4-A1 `cancelled`, WS5-A2
   `superseded`, WS5-A1 `in_progress`; WS0-A1, WS1-A4, WS4-A4, WS5-A3, WS6-A1, WS6-A4, CUR-1,
   SYN-3, WS4-A3 notes; artifact sha256s updated. Graph: DEC-11–DEC-14 and 6 `informs` edges
   (DEC-11 → WS4-A4, WS6-A4; DEC-12 → WS1-A4; DEC-13 → WS0-A1; DEC-14 → WS5-A1, WS4-A1) →
   80 nodes / 139 edges; DEC-2 resolved. Plan: 2026-10-09 header block; A.3, A.4, WS0-A1,
   WS0-A2, WS1-A4, WS4-A1, WS4-A4, WS5-A1–A3, WS6-A4, R2.2, R2.3, R3.5, R3.7, R3.8, R3.10,
   R3.11, D-5, D-9, D-10 revised; D-19 and D-20 added. Statements giving 96,729 as the current
   default chr1 count now give 96,731.

### 2026-10-08 — greedy is the overlap default again (DEC-10); configuration labels corrected

1. **DEC-10: the longest-arm (greedy) rule is the default again.** Author instruction
   2026-10-08: "make the greedy the default". `overlap_strategy` defaults to `"greedy"` in
   `scan_sequence`, `scan_fasta`, `scan_fasta_iter` and `scan_fasta_parallel`, whether or
   not scoring is on. The CLI's `-overlap-strategy` defaults to `greedy`, and the startup
   banner prints the overlap mode. The webapp form defaults to greedy, and choosing score
   with scoring off returns HTTP 400, as in Nikol's original code. These defaults are
   identical to Nikol's main (`d237bca`). **Why:** score mode's cost on a genome, measured on
   hg38 chr1 on an AMD EPYC 7763 with 16 vCPUs (the manuscript node's CPU): ~60 s vs 26.3 s
   for greedy at 16 workers (205.2 vs 172.1 s at 1 worker); published triplex 57.7 s. A
   further pure-Python optimisation gave identical output but no speedup and was reverted; a
   C-level selection was not done. Evidence:
   `RTR/evidence/WS6-A4_chr1_runtime_memory_2026-10-07.md`.
2. **Kept from DEC-9:** score mode as an option, with its verified result-preserving
   optimisation (DEC-9 evidence §3–4), and the PR #18 memory fix (chr1 peak RSS 2,704 →
   1,979 MB; output-neutral; not yet committed).
3. **Effects.** Benchmark detection metrics are unchanged (0 of 73,728 grid calls differ
   between the rules). The published speed numbers (25.2 s on chr1 with 16 workers; "2.3×
   faster than triplex") and the chr1 count 96,729 describe the **default** mode again:
   greedy reproduces 96,729 exactly and runs in 26.3 s on the same CPU model. They still need
   re-measuring on the final code (X6), and the memory fix still needs committing. The
   near-homopolymer edge-hit `TODO(authors)` now concerns only the optional score mode. The
   R3.10 draft was rewritten: it justifies **keeping** the longest-arm default (speed) and
   offering score mode, using the chr1 comparison (score mode retains every greedy locus;
   83.7 % identical coordinates; +1.5 % loci, 98,223 vs 96,729; ~2.3× runtime), and it
   recommends a mode (greedy for genome-wide scans, score where the most stable overlapping
   variant matters, `-skipoverlap` for exhaustive analysis). The R3.10 `TODO(authors)` on how
   to justify the choice was reworded for DEC-10. R3.11: the published runtime and hit
   numbers hold for the default again; the whole-genome figure still has to be labelled as
   an extrapolation or measured.
4. **Configuration labels corrected.** Measured 2026-10-07 on the 128-record benchmark:

   | configuration | where it is used | balanced (n=128) | experimental only (n=69) |
   |---|---|---|---|
   | minrep 8, mismatch 0.10, no filters | analysis scripts (`run_direct_sensitivity_analysis.py`, `REFERENCES["script_default"]`) | 62/2/64/0, MCC 0.969 | 62/2/5/0, MCC 0.832 |
   | minrep 10, no filters / both filters | **software defaults** (library / CLI) | 56/8/64/0, MCC 0.882 | 56/8/5/0, MCC 0.580 |
   | minrep 8, mismatch 0.15, both filters | manuscript | 63/1/63/1, MCC 0.969 (FP HDNA0028) | 63/1/4/1, MCC 0.784 |

   Earlier text here, in the plan, the ledger and the WS2-PRE/WS4-A3/WS2-A3 evidence called
   the minrep-8 result "library defaults" or "the shipped default"; that is relabelled
   everywhere as the benchmark/script configuration. New `TODO(authors)` (plan, under WS2;
   ledger DEC-5): decide which configuration the manuscript reports, and whether the shipped
   default minrep changes from 10 to 8 (nested CV picks minrep 8 in 100/100 folds; minrep is
   the dominant parameter, spread 0.203; minrep 10 reproduces the published chr1 count).
5. **Bookkeeping.** Ledger: new decision DEC-10; DEC-9 marked "default superseded by DEC-10"
   (its optimisation content kept); WS4-A4 set to `in_progress` (chr1 comparison measured,
   CSV not yet deposited). Graph: DEC-10 node and 3 `informs` edges (to WS4-A4, WS6-A4,
   WS4-A3) → 76 nodes / 133 edges. Plan: header bullet, WS4-A4, WS6-A4, R3.7, R3.10, R3.11,
   D-17 revised; D-18 added. Tests: 344 passed, 3 skipped.

### 2026-10-07 — 1:1 balance restored; score-informed overlap default; all executed results rerun

Two author decisions (recorded as **DEC-8** and **DEC-9** in the ledger), plus an
optimisation:

1. **DEC-8: four new synthetic non-forming records restore 1:1.** SYN0060
   `pAA32_mirrordisrupted` and SYN0061 `GA32_mirrordisrupted` (mirror-disrupted mutants),
   SYN0062 `pGC32_dinucshuffle` (dinucleotide shuffle) and SYN0063
   `random_matched_to_HDNA0004` (random GC/length-matched) follow the generator's
   perturbation-core split. They are drawn from an independent RNG
   (`REPLACEMENT_SEED = 20261002`) after the v3 pool, so the first 124 rows are
   byte-identical. Balanced set 124 (64:60) → **128 (64:64)**, 55 → 59 synthetic;
   experimental CSV unchanged. Validation passes (0 duplicates, 0 H-DNA-competent rows,
   max residual mirror identity 0.778). **None of the four is a homopolymer**, so the
   benchmark still can't test the filters, and the R2.1 poly-A `TODO(authors)` stands.
   Evidence: `RTR/evidence/DEC-8_rebalance_replacement_negatives.md`.
2. **DEC-9: overlap removal defaults to Nikol's score-informed strategy** (default
   superseded 2026-10-08 by DEC-10; see above). `overlap_strategy=None` resolves to
   `"score"` when scoring is on and to `"greedy"` when `score=False` / `-no-score`. An explicit `"score"` without scoring still raises
   `ValueError`. The webapp falls back to greedy when scoring is unchecked, and the CLI
   banner prints the mode. **No detection call changed on the benchmark** (0/128 differ
   greedy vs score at the reported configuration; 0 of 73,728 calls across the full WS4-A3
   grid, all four filter states). In general, with the homopolymer filter on, score can
   report a near-homopolymer hit where greedy reports none (e.g. ACGT×10 + G×100 + C×100 +
   ACGT×10 with the CLI filters: 0 hits greedy, 3 score); without that filter a call cannot
   change. Also **which hit is reported** in an overlapping cluster can change (e.g. GAMIR_SEQ:
   arm 6 under score, arm 7 under greedy). Evidence:
   `RTR/evidence/DEC-9_score_overlap_default_and_optimisation.md`.
3. **Score-based overlap removal optimised, with the logic kept the same:** incremental
   window search (O(arm³) → O(arm²) per hit), an O(arm) score upper bound, and lazy
   score-ordered selection. Poly-G 400 bp: 96 s → 0.06 s. A 2 Mb test genome: 39 s → 8.8 s
   (greedy 7.8 s). Score mode is still slower than greedy on long homopolymer-rich input
   (G3000+C2000+A5000: 5.4 s vs 1.6 s) and, as measured afterwards on hg38 chr1, about
   2.3× slower on a genome (~60 s vs 26.3 s at 16 workers; the reason for DEC-10). With
   the score default but before the
   optimisation, the suite had 35 failures (11 parallel-scan timeouts); after it, all
   passed (344 passed, 3 skipped as of 2026-10-08).
   New tests: `tests/test_score_overlap.py` (default semantics; lazy selection equals
   scoring every candidate; incremental window search equals naive recomputation). 24
   detector-geometry tests in `tests/test_hdna.py` pinned `overlap_strategy="greedy"`
   (no longer needed after DEC-10; the pins are gone).
   **The equivalence check against the pre-optimisation code has been run and passed**
   (DEC-9 evidence §4: 70,984 candidates × 4 flag combinations in the scorer, 19,710
   `scan_sequence` calls, all FASTA and parallel paths, long runs — all identical), so the
   optimisation is verified to give identical results to the pre-optimisation code.
   Remaining worst case: long G/C homopolymer runs with the homopolymer filter on
   (G400+C400: 40 s score vs 0.02 s greedy, growing roughly with the cube of the run
   length); `-overlap-strategy greedy` avoids it.

Effect on results (details in §3): whole-set metrics changed only through the four new
true negatives (TN 60 → 64); the random-split results moved more, because the split itself
changed (train 88 / eval 40). Fixed configurations still tie at the top (MCC 0.969) and still
beat per-fold tuning (nested CV 0.938, LOFO 0.923). The held-out split is balanced again
(eval 20:20). Tuning beats the fixed script configuration (minrep 8, no filters) in
**0/200** splits (tied 112, lost 88). Draft
answers R2.1, R2.2, R2.6, R3.1, R3.2, R3.5 and R3.9 were updated to 128; R3.7, R3.10 and
R3.11 were revised for DEC-9, with new `TODO(authors)` items under R3.10 and R3.11 (see
the plan).

Re-verified this round: WS1-A1 (default build 0× `-march=native`, opt-in 1×), WS1-A2/A3
tests, WS1-A4 (no streamlit), WS1-A5 (notebook load cells: 128 rows, 64/64), CUR-1
(experimental CSV unchanged), SYN-1/2/3, WS2-A4. WS3-A1/A2: only the HSeeker side was
re-run, because the comparator binaries aren't installed here; its values are unchanged.
Not rerun: the notebook beyond its load cells, `sensitivity_direct/`,
`sensitivity_injected/` (needs R), `ablation_eval.py`, chr1. The graph PNG is still stale.

### 2026-10-02 — main merged; homopolymers removed; all v3 results regenerated

Three author decisions (recorded as **DEC-6** and **DEC-7** in the ledger):

1. **Merge Nikol's `main` (PRs #16, #17) into RTR** (DEC-6).
2. **Where code overlaps, Nikol's methods win** (DEC-6). WS1-A1: opt-in is now
   `HSEEKER_NATIVE_BUILD=1`, and a native build targeting a different macOS
   architecture raises `RuntimeError`. WS1-A2: profiling counters are per-scan
   `ScanProfile` structs published into module `ProfileState` under the GIL.
   WS1-A3: **all sequence output, including `putative_triplex`, is lowercase**
   (RTR's uppercase normalisation is gone). Streamlit removal (WS1-A4) is kept.
3. **Keep RTR's data, minus the 6 pure single-base homopolymers** (DEC-7): HDNA0053
   pRW1405 (A×20, non-forming), HDNA0054 poly-(dC)30 (forming), SYN0053–SYN0056
   (poly-A15/A30, poly-T25/T35). Balanced set 130 (65:65) → **124 (64:60)**;
   experimental kept 71 (65/6) → **69 (64/5)**; 59 → 55 synthetic.

Effect on results (details in §3): all four false positives of the 130-record run were
homopolymers, so the AT and homopolymer filters now change **no** benchmark call, and
"no filter" is selected in every fold. The script default and the manuscript-stated
configuration now **tie** at the top (MCC 0.968), so the earlier "manuscript parameters
win under every scheme" is no longer true, and DEC-5 is no longer settled by data.
Tuning beats the fixed script configuration in **0/200** random splits. Draft letter text for R2.1, R2.6,
R3.5 and R3.9 now carries `TODO(authors)` flags in the plan.

Not rerun: `analysis/hdna_detection_benchmark.ipynb` (needs the E. coli genome + R
Triplex), `sensitivity_direct/`, `sensitivity_injected/`, and Nikol's
`analysis/scripts/ablation_eval.py`. `RTR/hseeker_revision_graph_preview.png` is
**stale** (pre-merge graph; it can't be regenerated here). The JSON/HTML graph was
hand-edited instead.

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

## 2. Status — 11 complete, 7 part-done, 14 pending, 2 dropped (2026-10-09, counted from the ledger; was 11 / 5 / 18 before WS4-A4 started on 2026-10-08 and DEC-14 on 2026-10-09; all complete items re-verified 2026-10-07)

### Complete and verified

| action | item | outcome |
|---|---|---|
| WS1-A1 | R2.8 | `-march=native` now opt-in via `HSEEKER_NATIVE_BUILD=1` (merged, Nikol's implementation); errors if `ARCHFLAGS` targets a different macOS arch; CI wheels built portably; `tests/test_build_flags.py`. Re-measured 2026-10-07 on the merged build: default 0×, opt-in 1× |
| WS1-A2 | R2.9 | Profiling counters now per-scan `ScanProfile` structs published under the GIL (merged, Nikol's implementation; replaced RTR's thread-local fix); exact totals across threads; `tests/test_profiling_threadsafe.py`. Measured pre-fix corruption: 8 concurrent scans of identical input gave 8 different totals, worst **3.84× inflated** |
| WS1-A3 | R2.10 | All sequence fields, including `putative_triplex`, are **lowercase** on every path and in the CLI TSV (merged, Nikol's implementation; replaced RTR's uppercase normalisation); `tests/test_output_case.py` |
| WS1-A4 | R2.11 | **Streamlit removed entirely** (author decision superseded the planned rename); one front-end remains |
| WS1-A5 | R2.2 | Stale "54 sequences" docstring fixed; notebook's missing-input failure removed; load cells give 128 rows (64/64) as of 2026-10-07 |
| AUDIT-1, CUR-1, SYN-1/2/3, WS2-A4 | R3.5, R2.1 | Provenance audit and benchmark v3 curation, re-verified in this repo; re-verified 2026-10-02 after the homopolymer removal (69 kept, 124 balanced) and 2026-10-07 after DEC-8 (69 kept, 128 balanced at 64:64) |

### Part-done (evidence banked, not finished)

| action | item | what exists / what's missing |
|---|---|---|
| WS2-A2 | R2.7, R3.1 | Wilson + record-bootstrap + **family-cluster** CIs implemented. Missing: tool-vs-tool difference CIs (no R/Triplex here) |
| WS2-A3 | R3.2, R3.3 | LOFO + LOSO + both collapsed variants done for the **detection call**. Missing: the stability-score threshold |
| WS4-A3 | R3.9 | Full-factorial 576-config sweep done. Missing: pairing-only vs pairing+stacking ablation (WS4-A2) |
| WS3-A1 | R3.6 | NeSSie builds (commit `dbe6cb2`); interface characterised. Missing: the benchmark run. Binary not installed in this environment, so only the HSeeker side was re-checked 2026-10-07 (unchanged) |
| WS3-A2 | R3.6 | non-B_gfa builds (commit `a891b59`). Missing: the benchmark run. Same 2026-10-07 caveat as WS3-A1 |
| WS4-A4 | R3.10 | chr1 greedy-vs-score comparison measured (2026-10-07 pre-merge; 2026-10-09 on the merged code) and the final R3.10 text written (`RTR/evidence/R3.10_overlap_response_2026-10-09.md`). Missing: `overlap_selection_comparison.csv` and the script that regenerates it (X7) |
| WS5-A1 | R2.3, R3.8 | Kouzine 2017 comparison (DEC-14): design pre-specified (`RTR/evidence/X8_kouzine_analysis_design_2026-10-09.md`), download and alignment running. Missing: HSeeker scans of mm10/hg38 in the frozen configuration, the analysis, every result |

### Dropped by an author decision (2026-10-09, DEC-14)

| action | item | why |
|---|---|---|
| WS4-A1 | R3.7 | `cancelled`: the E. coli insertion experiment leaves the manuscript, so the overlap-threshold sweep is not run (nor X5) |
| WS5-A2 | R3.8 | `superseded`: the mouse comparison is now part of WS5-A1 (Kouzine 2017 on mm10) |

---

## 3. Headline result (regenerated 2026-10-07 on the final code and the 128-record set)

Under **every** validation scheme the two **fixed** configurations beat per-fold tuning,
and they tie with each other: the script default (`purity 0.90, mismatch 0.10, minrep 8,
maxspacer 10, no filters` — the analysis scripts' benchmark configuration, **not** the
software default, which is minrep 10) and the manuscript-stated configuration (`mismatch 0.15,
minrep 8, purity 0.90, maxspacer 10, both filters`):

| scheme | config | sens | spec | MCC | record CI | family-cluster CI |
|---|---|---|---|---|---|---|
| leave-one-family-out | manuscript | 0.984 | 0.984 | **0.969** | [0.922, 1.000] | [0.901, 1.000] |
| leave-one-family-out | script default | 0.969 | 1.000 | **0.969** | [0.925, 1.000] | [0.909, 1.000] |
| leave-one-study-out | manuscript | 0.984 | 0.984 | 0.969 | [0.922, 1.000] | [0.892, 1.000] (study-cluster) |
| leave-one-study-out | script default | 0.969 | 1.000 | 0.969 | [0.925, 1.000] | [0.899, 1.000] (study-cluster) |
| sequence-level 5-fold CV | manuscript | 0.984 | 0.984 | 0.969 | [0.922, 1.000] | — |
| sequence-level 5-fold CV | script default | 0.969 | 1.000 | 0.969 | [0.924, 1.000] | — |
| collapsed per (family, class) | manuscript | 0.977 | 0.964 | 0.942 | [0.851, 1.000] | [0.8545, 1.000] |
| collapsed per (family, class) | script default | 0.955 | 1.000 | 0.944 | [0.865, 1.000] | [0.869, 1.000] |

(Precise MCCs: script default 0.9692, manuscript 0.9688. Collapsed per (family, class) is
now n = 72, 44/28; it was 71.)

Both beat per-fold tuning (LOFO/LOSO 0.923, nested CV 0.938) and the shipped software
defaults (production CLI and library alike, minrep 10: 0.882, whose 8 FNs come from
`minrep 10`, not from the filters). Neither requires any
selection, so R3.3's circularity objection does not apply to either. **History:** on 130
records the manuscript configuration won outright (0.939 vs tuned 0.895/0.909, script
default 0.908, production 0.856). On 124 records (2026-10-02) it tied the script default
at 0.968 (tuned 0.920/0.936, production 0.879). On 128 the tie holds, at 0.969. The
filters change no call, because the four homopolymer false positives they removed are out
of the benchmark and the DEC-8 replacements are not homopolymers. Neither DEC-9 nor
DEC-10 changes a detection call on the benchmark (0 of 73,728 grid calls differ between
the rules), so they move none of these numbers; in general, with the homopolymer filter on, score can report a near-homopolymer
hit where greedy reports none (see DEC-9 evidence). DEC-11 (Nikol's detector branch, merged
2026-10-09) changes none of the 73,728 grid calls either, and DEC-2 changed only citation
metadata, so every number in this section stands; only the input sha256 recorded in the
results' `metadata.json` is now out of date (refreshed at the final freeze).

### Why we can say HSeeker is not overfitted (R2.6)

1. Tuning on the benchmark makes results **worse** (0.969 for either fixed config vs
   0.938 nested-CV tuned and 0.923 LOFO tuned).
2. Over 200 repeated splits, tuning **never** beats the fixed script configuration
   (minrep 8, no filters): it wins **0/200**,
   ties 112 and loses 88 (median gain +0.000, 95% [−0.0955, +0.000]).
3. **104 of 576** configurations lie within 0.05 MCC of the held-out maximum (32 tie at
   it). The surface is a plateau, not a peak (124 records: 76 / 16; 130 records: 108 / 75).
   On this split the script configuration (minrep 8, no filters) reaches eval MCC 1.000
   and the train-tuned pick (maxspacer 5) 0.905; the shipped software defaults (minrep 10,
   `production_cli`) reach 0.951.
4. `minrep 8` and `purity 0.85` selected in **51/51** family folds.
5. **Nothing is fitted.** LOFO, LOSO and full-data evaluation are numerically identical
   for a fixed config (manuscript TP63 FN1 TN63 FP1; script default TP62 FN2 TN64 FP0),
   because with fixed parameters the call depends only on the record's own sequence.

**Honest residual (unchanged by DEC-8):** the AT/homopolymer thresholds *were* chosen
with knowledge of this benchmark's composition. After the 2026-10-02 removal the
benchmark holds no record that either filter acts on, and the four DEC-8 replacements
don't change that. The reported benchmark numbers therefore don't depend on the filters,
but the benchmark also can't test them: "no filter" wins 100/100 nested-CV folds and
51/51 family folds, through the fewer-filters tie-break. On 130 records the filters were
what removed all four FPs. R2.6 keeps its force here, and the removal must be disclosed,
not presented as neutral. See the `TODO(authors)` under R2.6 in the plan.

### On R3.2 specifically

Family-aware validation does not move the point estimate — it moves the *precision*.
For the manuscript configuration the family-cluster CI is **1.26×** wider than the
record-level CI (range 0.91–1.75×, median 1.23×, across all 14 reported rows; was
1.23× (median 1.20×) on 124 records and 1.26× on 130 records). Both are reported everywhere: the record-level interval
answers R2.7, the cluster interval answers R3.2's "inflate apparent statistical power".

---

## 4. Blocked on two author decisions (DEC-2 resolved 2026-10-09)

| decision | blocks | evidence |
|---|---|---|
| ~~**DEC-2** FXN citation~~ | resolved 2026-10-09 | `RTR/evidence/DEC-2_fxn_citation.md` |
| **DEC-3** add Hanvey 1988 `(TTC)8` | WS0-A1 → 21 actions | `RTR/evidence/DEC-3_pRW1406.md` |
| **DEC-5** canonical parameters | WS2-A1 → 16 actions | **no longer settled by data** (§3: tie) |

- **DEC-2 (resolved 2026-10-09):** the cited PMID (Romney 2008, a *C. elegans* ferritin
  paper — likely an FTN-1/FXN symbol confusion) and the cited DOI (Samuel 2007
  innate-immunity minireview) were two *different* wrong papers. The records now cite primary
  experimental studies: Potaman 2004 (PMID 14978261) for (GAA)10/20/33 and Sakamoto 1999
  (PMID 10230399) for (GAA)66. They stay secondary, and the original identifiers stay in
  `cited_*` as an audit trail. Potaman 2004 also supplies the primary support for the
  `forming` label of the short tracts that the audit had found missing.
- **DEC-3:** Hanvey 1988's abstract establishes both the sequence and the forming
  label for `(TTC)8`. Recommend adding it as `Hanvey1988_(TTC)8` rather than
  `pRW1406`, because the plasmid-number mapping comes from the document already shown
  corrupted for this study. If added, the balanced set becomes 129 records at 65:64
  (70 kept, 65/5), which breaks the 1:1 that DEC-8 restored; whether to add a matching
  synthetic negative is an author choice.
- **DEC-5:** on 130 records §3 showed the manuscript's values winning out-of-fold under
  every scheme. On 128 records (as on 124) they **tie** the script default (mismatch
  0.10, no filters) at MCC 0.969 (124: 0.968), and mismatch 0.10 and 0.15 have the same
  median held-out MCC across the grid (0.838; 124: 0.788). `TODO(authors):` decide DEC-5 on consistency with the
  reviewed manuscript. Benchmark performance no longer discriminates between the two.
  Related `TODO(authors)` (2026-10-08, plan under WS2): which configuration the manuscript
  reports, and whether the shipped default minrep changes from 10 to 8 (software defaults
  MCC 0.882 vs 0.969 for both minrep-8 configurations).

**WS0-A1 is not finished** even once decided: the 22 PMID and 2 DOI corrections sit
in `corrected_*` columns without being applied to `cited_*` (for the 4 FXN rows DEC-2
filled `corrected_*` and deliberately kept the original `cited_*` as an audit trail), and
the `metadata.json` sha256s have to be refreshed at the freeze. pCON is **not** an open item:
DEC-13 (2026-10-09) keeps HDNA0029 exactly as published (31 nt), so the audit's terminal-C
correction (register D-9) is not applied and no text may claim it was.

---

## 5. Next steps (planned, not started)

Four items needed no decision (item 4 was dropped 2026-10-09 by DEC-14, item 3 now has its final text). Two verified facts collapse the work:
`pairing_score`/`stacking_score` are already emitted per hit (`total == pairing +
stacking` exactly), and `--overlap-fraction` is consumed *after* the scan.

| # | item | ask | new compute |
|---|---|---|---|
| 1 | Score threshold out-of-fold | R3.3 | ~0 — one 27 ms scan, reuse fold machinery |
| 2 | Scoring ablation | R3.9 | **0** — re-threshold on components already returned |
| 3 | Overlap-selection rule | R3.10 | chr1 comparison **re-measured 2026-10-09 on the merged code** (WS6-A4 evidence §5: score keeps all 96,731 greedy loci, 83.7 % identical, 1,075 score-only loci, 42.6 s vs 25.9 s; pre-merge 2026-10-07: 96,729 / +1.5 % / ~60 s vs 26.3 s); the final R3.10 text exists (`RTR/evidence/R3.10_overlap_response_2026-10-09.md`), and its two `TODO(authors)` are resolved; still to deposit as `overlap_selection_comparison.csv` (X7). 0/128 benchmark calls differ (0 of 73,728 across the WS4-A3 grid), though with the homopolymer filter on score can in general report near-homopolymer hits greedy does not |
| 4 | ~~Overlap-threshold sweep~~ | R3.7 | **dropped 2026-10-09 (DEC-14):** the E. coli experiment leaves the manuscript, so there is nothing to sweep |

**Item 1 is the priority** — it is the last live circularity, and the one R3.3 named.

### WS5 data is fetched and validated (not committed — 208 MB)

**Superseded 2026-10-09 (DEC-14):** WS5 now uses Kouzine et al. 2017 ssDNA-seq (SRA072844;
mouse B cells on mm10 + human Raji on hg38), reprocessed from raw reads; the pipeline is
running at `/workspaces/.hseeker-pr/kouzine` and the design is fixed in
`RTR/evidence/X8_kouzine_analysis_design_2026-10-09.md`. The S1-END-seq KM12 plan below is
kept for the record; it is optional and not pursued.

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
`cabf2f0`) — stream with `scan_fasta_iter`. Both rules are now public library options
(`overlap_strategy="greedy"` / `"score"`), so the comparison no longer needs a Python
re-implementation of the longest-arm rule. The existing chr1 output
(`analysis/results/chr1/`, 96,729 loci) was produced under greedy, before DEC-9; greedy is
the default again (DEC-10), and a fresh greedy run on 2026-10-07 reproduced 96,729 exactly.
Since DEC-11 (2026-10-09) the merged code gives 96,731 loci with the same parameters; 96,729
is the published, pre-merge count.

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
`grouped_validation_v3/` (both regenerated 2026-10-07 on the final code and the
128-record set; previously 2026-10-02 on 124 records) are from **2026-06-24 against the old 80-row CSV** with
pre-filter code. They do not reflect benchmark v3 and should be treated as superseded.
`robustness_v3/` and `grouped_validation_v3/` record the balanced sha256 `7bb49544…`; since
DEC-2 (2026-10-09) the CSV is `c5177b24…` (citation metadata only, no call changes), so their
`metadata.json` is refreshed at the final benchmark freeze.
