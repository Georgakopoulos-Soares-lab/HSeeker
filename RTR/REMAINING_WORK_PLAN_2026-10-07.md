# HSeeker revision: remaining work plan (2026-10-07)

**Scope:** everything still needed to answer the 24 reviewer comments (R1, R2.1–R2.11,
R3.1–R3.12) and to remove stale claims from the code, docs, plots, RTR files and the manuscript.
**Basis:** verified on branch `RTR-merge-main` (uncommitted), which contains all of Nikol's pushed
work (origin/main `d237bca`, PR #16–#18). Three read-only audits (code/docs/plots; RTR
consistency/coverage; manuscript-claim map) were run and their material findings re-checked by
hand; the commands and numbers are cited inline. One audit claim was rejected as wrong ("library
default minrep is 8" — `inspect.signature(hseeker.scan_sequence)` gives `minrep=10`).
**Manuscript text is not in the repo**; manuscript claims are reconstructed from the reviewers'
verbatim comments (graph node titles), the RTR plan, `algorithms.tex`, `benchmarks/` and the
notebook outputs.

---

## Update 2026-10-08: DEC-10, greedy is the default again

Author instruction: "make the greedy the default". Library (all four `scan_*`), CLI and webapp defaults are `overlap_strategy="greedy"`, identical to Nikol's main; score mode stays optional with the verified optimisation. Tests: 344 passed, 3 skipped. Consequences for this plan: **F1b is resolved** (published speed claims 25.2 s / 2.3x vs triplex and the 96,729 chr1 count hold for the default mode again, pending re-measurement X6 and the memory fix F1c); the R3.10 edge-hit TODO no longer concerns the default; R3.10 now has to justify **keeping** the longest-arm default (speed; the chr1 comparison shows score mode retains every greedy locus) and recommend a mode; benchmark scripts must still pin `-overlap-strategy` for provenance.

## 0. Ground truth today

| Item | Verified state | Evidence |
|---|---|---|
| Code | origin/RTR + origin/main up to `d237bca` + DEC-9 + score-selection optimisation + PR #18 memory fix; pending merge parent = `d237bca` | `git rev-parse MERGE_HEAD`; line-by-line check that every line Nikol added is present or deliberately replaced |
| Tests | 343 passed, 3 skipped (hg38 data) | `pytest tests` |
| Equivalence of optimised score selection | identical outputs to pre-optimisation code (70,984 candidates × 4 flag combos; 19,710 `scan_sequence` calls; FASTA/parallel paths) **and** Nikol's golden outputs (`tests/test_parallel_batching.py`) | `RTR/evidence/DEC-9_…md` §4 |
| Not in main | Nikol's `fix/hdna-maximal-representation-v3` (75c3494, 2026-09-25: C spacer normalisation) | `git merge-base --is-ancestor` → not merged |
| Benchmark | balanced v3 = 128 (64:64; 69 experimental + 59 synthetic), sha `7bb49544…`; experimental 80 rows → 69 kept (64/5), 9 removed, 2 excluded, sha `7197c6df…` | `sha256sum`; CSV counts |
| Current results | `analysis/results/robustness_v3/`, `analysis/results/grouped_validation_v3/` (128 records, current code) | metadata.json |
| **Not current** | `analysis/results/sensitivity_direct/`, `sensitivity_injected/`, `sensitivity_summary.md`, `analysis/benchmark_outputs/*`, notebook outputs, `analysis/results/chr1/*`, `benchmarks/benchmark_comparison.md` + sweep plot — all June 2026, old 80-row CSV and/or pre-merge greedy code | audit A; e.g. chr1 TSV has uppercase `putative_triplex` (pre-merge) |
| chr1 runtime (same CPU model as manuscript node, EPYC 7763 16 vCPU) | greedy 16w 26.3 s / 96,729 hits (= published); **score default 16w ~60 s / 98,223**; 1w 172.1 vs 205.2 s; RSS 1,979 MB after the PR #18 fix (2,704 MB before) | `RTR/evidence/WS6-A4_chr1_runtime_memory_2026-10-07.md` |

### The configuration problem (affects every accuracy number)

Measured today on the 128-record benchmark (Wilson 95 % CIs):

| Configuration | Where it is the default | Balanced (n=128) | Experimental only (n=69) |
|---|---|---|---|
| minrep 8, mismatch 0.10, no filters | analysis scripts only (`run_direct_sensitivity_analysis.py:151`, `REFERENCES["script_default"]`) | 62/2/64/0, MCC **0.969**, sens 0.969 [0.893, 0.991] | 62/2/5/0, MCC 0.832 |
| minrep 8, mismatch 0.15, both filters | manuscript text | 63/1/63/1, MCC **0.969** (FP HDNA0028) | 63/1/4/1, MCC 0.784 |
| minrep 10, mismatch 0.10, no filters | **library** (`scan_*`) | 56/8/64/0, MCC **0.882**, sens 0.875 [0.772, 0.935] | 56/8/5/0, MCC 0.580 |
| minrep 10, mismatch 0.10, both filters | **CLI** | 56/8/64/0, MCC **0.882** | 56/8/5/0, MCC 0.580 |
| same as library (no filters) | **webapp** (`webapp/main.py:574-585` passes no filter args) | = library | = library |

The drafts quote 0.969 as the "default"/"library defaults" result — **wrong label**: the shipped
software defaults give 0.882. Either the manuscript reports the minrep-8 configuration explicitly
as "benchmark configuration", or the shipped default changes (decision F3 below). minrep is the
dominant parameter (largest spread of median held-out MCC in the grid, 0.203; nested CV picks
minrep 8 in 100/100 folds — `robustness_v3`).

---

## 1. Freeze decisions (blocking — every later run depends on them)

Order matters: code → benchmark → configuration. Re-running anything before these are settled
produces numbers that will be invalidated.

### F1. Code freeze
| # | Decision | Evidence | Recommendation |
|---|---|---|---|
| F1a | Merge Nikol's `fix/hdna-maximal-representation-v3`? | Changes `_hdna.c` hit representation (+90 lines) → can change calls, counts, coordinates | Ask Nikol now; if it is meant to ship, merge **before** any rerun |
| F1b | Score-default speed on genomes | chr1 16w: score ~60 s vs greedy 26.3 s; published triplex 57.7 s. Overhead ~33 s is serial Python on 1.39 M candidates (C scan ≈19–21 s both modes; exact scoring 18.4 s; rest per-candidate). A further pure-Python attempt (staged bounds) gave identical output and **no** speedup → reverted | Options: (a) move bound computation + selection loop into the C extension (bounds computed in the parallel scan threads; est. ~60 → ~35 s; needs cross-platform wheel rebuild); (b) keep score default for library/webapp, make the **CLI** genome-scan default greedy; (c) keep score everywhere and drop the speed claims. Recommend (a) if the speed claim matters to the paper, else (b) |
| F1c | Commit the PR #18 memory fix and tell Nikol | RSS 1,960 → 2,704 MB introduced by `d237bca`; fix restores 1,979 MB; output-neutral (batch contents unchanged; golden tests pass) | Commit; send to Nikol upstream |
| F1d | Webapp filter mismatch | Webapp applies no AT/homopolymer filter, but `about.html:169-171, :192` says it "matches CLI defaults exactly" | Either pass `at_threshold=0.8, filter_homopolymers=True` in `_run_job`, or fix the text. Decide which behaviour the paper describes |

### F2. Benchmark freeze (WS0-A1)
| # | Decision | Evidence |
|---|---|---|
| F2a | DEC-2 (FXN citation; 4 records) | All 4 FXN rows: `verification_status=citation_corrected` but `corrected_pmid`/`corrected_doi` empty (CSV check); DEC-2 evidence §3–4 |
| F2b | DEC-3 (add Hanvey 1988 (TTC)8; if yes, also +1 negative to keep 1:1 → 130) | `RTR/evidence/DEC-3_pRW1406.md` |
| F2c | Apply the WS0-A1 corrections the drafts already claim as done | pCON HDNA0029 is still 31 nt `…GCAGTG` in both CSVs (register D-9 lists the fix as the resolution); drafts R2.2/R3.5 say "found and fixed" → currently false |
| F2d | Homopolymer handling for R2.1/R3.5 (TODOs) | Separate panel (current code, minrep 8): A20/A30/T25/T35 hit without filters, rejected by either filter; poly-dC30 (forming, Kohwi 1988) kept by AT filter, rejected by homopolymer filter. Literature: Hanvey 1988 (A20 B-form), Fox 1990 (A69 forms, A23/A33 not), Kohwi 1988 |

### F3. Configuration freeze (WS2-A1)
| # | Decision | Evidence | Recommendation |
|---|---|---|---|
| F3a | DEC-5 mismatch 0.10 vs 0.15 | Tie out of fold (nested CV 0.9692 vs 0.9688; LOFO equal). 0.15 gains HDNA0058 (a tandem *direct* repeat, outside the model) and adds FP HDNA0028 (B-33), both at 87.5 % identity. Belotserkovskii 1990: single substitution in a 15-nt arm (~7 %) still forms H-DNA | 0.10; correct the manuscript |
| F3b | **Which configuration the manuscript reports, and whether the shipped default minrep changes 10 → 8** | Table in §0: 0.969 (minrep 8) vs 0.882 (shipped minrep 10). minrep 10 reproduces the published chr1 count (96,729) | Decide explicitly. If minrep stays 10, every benchmark claim must say "benchmark configuration minrep 8" and the 0.882 default result must be reported too; if it changes to 8, chr1/runtime must be re-measured |
| F3c | Comparator arguments | Notebook figures used triplex min_score 15 / min_len 8 / max_len 50 / prokaryotic and Triplexator `-l 10 -e 15 -fr off`; the deposited scripts and chr1 timing used package defaults (audit C4) | Pin one set per experiment in `analysis/benchmark_config.json` |

### F4. Author TODOs (text decisions; can run in parallel with F1–F3)
R2.1 poly-A; R2.6 (1) disclosure of the removal + pre-removal numbers, (2) justification for keeping
filters on; R3.5 one-sentence rationale (drafted in chat 2026-10-07); R3.10 (1) justify DEC-9 —
the reviewer explicitly asked to "compare … with selecting the highest-scoring candidate", and the
chr1 comparison now exists (100 % of greedy loci retained, 83.7 % identical coordinates, +1.5 %
loci), (2) edge hits (negligible on chr1: 5 of 936 CLI score-only loci ≥ 90 % one base; longest
G/C runs 27/21 bp); R3.11 which runtime numbers to publish (depends on F1b).

---

## 2. Runs required after the freeze

Feasibility checked today: GitHub (NeSSie, non-B_gfa), UCSC (hg38), NCBI GEO reachable;
`bioconductor-triplex` 1.50.0 (the published version) available on bioconda via `/opt/conda`;
Triplexator 1.3.2 is **not** on conda (published runs used a TACC build at `/scratch/10899/kimopro/...`).

| # | Run | Answers | Replaces / supports | Script & inputs | Est. effort |
|---|---|---|---|---|---|
| X1 | Direct benchmark, HSeeker + triplex + Triplexator, on the frozen benchmark (balanced and experimental subsets), with Wilson CIs and family-cluster bootstrap CIs on tool differences → `summary_metrics_v2.csv` | R2.1, R2.7, R3.1, R3.6, D-15 | every tool-comparison claim (currently **no v3 evidence**; the repo's own June direct run is a **tie**: HSeeker = triplex 68/3/7/2, MCC 0.703) | `run_direct_sensitivity_analysis.py` (+ Triplexator runner); needs R env + Triplexator build | 0.5–1 day |
| X2 | NeSSie and nBMST on the frozen benchmark, with a stated hit→call rule | R3.6 | removes the unsupported "nBMST lower than ours" sentence | builds exist only on another machine (WS3 evidence) | 0.5 day |
| X3 | Out-of-fold stability-score threshold (grouped/nested CV) | R3.3 | Fig 3B–D (AUC 0.980/0.981, Youden ≥ 76, MCC 0.729 → 0.875 — all in-sample, old data) | new code in `run_grouped_validation.py`; note: on v3 the AUC is degenerate (0.984 = (1+sens)/2: no negative has a hit), so expect "no threshold needed"; report honestly | 0.5 day |
| X4 | Component ablation → `ablation_metrics.csv` | R3.9 | — (Nikol's `analysis/scripts/ablation_eval.py` exists, never run on v3) | HSeeker only | 2 h |
| X5 | E. coli injected arm on v3, overlap mode pinned, + overlap-threshold sweep 0.5–0.95 → `overlap_sweep.csv` | R2.3, R3.7 | Figure 3A (old 79-seq run, mismatch 0.15, greedy, pre-filter code) | `run_injected_sensitivity_analysis.py` (already points at v3); needs K-12 genome + R | 0.5–1 day |
| X6 | Runtime/memory: chr1 1–16 worker sweep in the final mode(s) + full hg38 (or keep "extrapolation" label); optionally re-time triplex/Triplexator on the same CPU | R3.11 | `benchmarks/benchmark_comparison.md`, sweep plot, "25.2 s", "2.3× faster", "~5 min genome" | `benchmarks/parallel_sweep.py`, `benchmark.py` — **must pass `-overlap-strategy` explicitly** (they now run score + filters by default) | 0.5 day |
| X7 | Overlap-selection comparison deposit → `overlap_selection_comparison.csv` (benchmark + chr1) | R3.10 | — | numbers exist in WS6-A4 evidence; write a script that regenerates them | 2 h |
| X8 | S1-END-seq enrichment (re-fetch GSE204808/GSE203632; data are **hg19 coverage, no peaks**; previously fetched copy is gone) → re-run chr1 on hg19 in the final configuration | R2.3, R3.8 | WS5-A1 | new script | 1–2 days |
| X9 | Homopolymer panel table (if F2d = report separately) | R2.1 | — | trivial | 1 h |
| X10 | Regenerate Supplementary Table 1 from the frozen CSVs; regenerate chr1 dot plot (and commit its generator — the SVG has none) | R3.5, R2.2 | Supp. Table 1; `analysis/results/plots/*.svg` | — | 2 h |

---

## 3. Stale claims to fix (by location)

Severity A = wrong statement a user/reviewer would rely on. Line numbers as of today.

### 3.1 Code behaviour
- **Webapp ignores the composition filters** (`webapp/main.py:574-585`) while claiming CLI parity → F1d.
- `scan_fasta_parallel` returns hits **sorted by seq_id** when overlaps are removed (`__init__.py`, final sort), not in input order as README:307-308 says — fix the doc (or the code, if input order is intended).

### 3.2 README.md (A)
:86 overlap "longest arm first" (default is score) · :307-308 output order · :343 "-at-threshold Left-arm" (full motif) · :397-399/:472 "None if scoring disabled"/"removes the four columns" (API omits keys; TSV keeps empty columns) · :485 "costs little more than greedy (8.8 vs 7.8 s)" (chr1: ~60 vs 26 s) · :493-497 "~N× speedup", "chr1 ~78 s on 16 cores" (measured 26.3 s greedy / ~60 s score; 7.45× published) · :498 "capped at 1,000,000 hits" (no cap, `_hdna.c:78-83`) · :502 "N breaks arm extension" (right arm only; left-arm N → unscorable) · :400/:470 putative_triplex on the purine strand for CT-rich hits. (B/C) missing parameters (`parser`, `purity_rmq`, `at_threshold`, `filter_homopolymers`), "134 tests" (346 collected), quick-start example returns `[]` at minrep 10, example row not reproducible, HDNAhunter URLs.

### 3.3 Docstrings / help / developer docs
`__init__.py:382-384` and `ablation.py:55` "left-arm AT" (A) · `__init__.py:348` N handling (A) · `_hdna.c:806-818` stale C defaults (C) · `__main__.py:10-11` "drop-in replacement for findHDNA" (C) · CLI help does not say the homopolymer filter is always on · `.github/copilot-instructions.md`: scoring "after overlap removal" (:162), "longest arm" (:70, :81, :253, :257), non-existent `DEFAULT_MAX_HITS`, test counts, missing files.

### 3.4 Webapp text (`templates/about.html`, `job.html`)
CLI-parity claim (:169-171, :192); "longest arm" (:108, :240-241); the whole Validation / "Why minrep=8" section (39 sequences / 28 papers, 97.2 %, "100 % sensitivity across 35", minrep 6, "only 3 non-forming controls"; :12-14, :255-330, :393-395); "should not be raised above 0.80" vs shipped purity 0.90 (:228-230); spacer "" vs actual "." (:137); "Stable sites (score ≥ 50)" (`job.html:149, :236`, `main.py:518`) — **unvalidated threshold** (no out-of-fold evaluation; see X3).

### 3.5 `algorithms.tex` (manuscript algorithm panels) — all A
Phase 2 describes only greedy longest-arm removal **before** scoring (:99-119) → add score-ordered selection and reorder phases · :162 old window bound `R ← min(L+L_min, al)` (now never shorter than min(8, al)) · :180-181 "Return None" when no window fits (now scores the original arms) · no AT filter (full motif) and no homopolymer filter (PRE/POST) · :186-187 total not clamped (code: `max(total, 0)`) · :193-194 uppercase putative_triplex (now lowercase) · :10 RMQ prefilter presented as part of Phase 1 (optional, default off).

### 3.6 `benchmarks/`
`benchmark_comparison.md:33, :94, :115, :149-151` (25.2 s, 96,729, "2.3× faster than triplex") are **greedy, v0.1.0** — label or re-measure (X6) · `parallel_sweep.py:71-78`, `benchmark.py:164-170` must pin `-overlap-strategy` · sweep plot comparator numbers are hard-coded (`plot_parallel_sweep.py:32, 40`).

### 3.7 `analysis/`
`results/sensitivity_summary.md` "Current headline metrics" on the old 80-row CSV (A) · `sensitivity_direct/` + `sensitivity_injected/` summaries/metadata/plots (old CSV; in-sample Youden) · `chr1/` (greedy, other machine 10.37 s, pre-merge uppercase, negative scores no longer possible) · notebook saved outputs (79 sequences, AUC 0.977/0.961, MCC 0.729/0.702/0.649; `-mismatch 0.15`; "inject 5×" vs `N_REPLICATES=1`) · `run_parameter_robustness.py:10-11` docstring "124 … 64:60" · `plot_sensitivity_confusion_matrices.R:153` "80 curated". Action: rerun (X1, X5) or move to `analysis/results/archive_2026-06/` with a README stating data/code version.

### 3.8 RTR documents
- **Mislabel throughout:** "library defaults … MCC 0.969" → it is the minrep-8 script configuration (§0). Affects Plan R2.6 draft, WS2-PRE evidence tables, CHECKPOINT, ledger.
- **DEC-9 performance:** Plan:17, :139, :223; ledger DEC-9 (~:246-258); CHECKPOINT:47-60; DEC-9 evidence §3 — replace "8.8 vs 7.8 s / costs little" with the chr1 result (~60 vs 26 s).
- **Superseded by today's measurements:** Plan:230 "chr1 count under score not computed" (98,223 / 43,150); Plan:122 and ledger:898 "hit concordance not measured" (now measured); Plan:229 deposited timings not labelled greedy; Plan:128 WS5 "human, no liftover"/"triplex loci" (hg19 coverage, no peaks; data gone).
- **Claimed but not done:** pCON fix and "24 citation identifier corrections" (R2.2/R3.5 drafts; FXN rows uncorrected); R2.5 "done … verified"; R3.3 "threshold now selected within training folds" (not run).
- **Status drift:** WS2-A4 `done` without `stress_panel_metrics.csv`; WS4-A4/WS6-A4 `pending` though evidence exists; provenance_audit `verified` vs FXN; graph statuses vs ledger (WS1-A1…A5 pending in graph, verified in ledger; GATE_BM done while WS0-A1 pending); PNG stale.
- **Not recorded anywhere:** PR #18 merge, memory regression/fix, the minrep configuration problem, the unmerged Nikol branch, the speed decision F1b.
- **Test counts:** "225 passed" → 343 passed, 3 skipped (Plan:17, CHECKPOINT:7/50, ledger, DEC-9 evidence).
- **Stale code line references:** `_hdna.c:499-558` → 512-570; `__init__.py:191-220` → homopolymer filter 269-296; `__main__.py:104-121` → 107-122; `run_injected…:471` → 512; `run_direct…:99-141` → 145+; WS1-A3 evidence `_scoring.py`/`__init__.py` lines (list in audit B).
- **Cited output files that do not exist** (draft answers): `benchmark_config.json`, `summary_metrics_v2.csv`, `lofo/loso_metrics.csv` (content exists as `grouped_validation_v3/`), `stress_panel_metrics.csv`, `nessie/nbmst_benchmark.csv`, `overlap_sweep.csv`, `ablation_metrics.csv`, `parameter_sweep.csv` (exists as `robustness_v3/`), `overlap_selection_comparison.csv`, `chr1_vs_s1endseq_enrichment.csv`.
- Graph R-node titles contain RTR-author fragments appended to the reviewers' words (R2.3, R2.4, R2.7, R3.6 incl. the nBMST sentence, R3.7, R3.11, R3.12) — separate before any quote is reused.

### 3.9 Manuscript (reconstructed; highest risk first)
1. **"Outperformed existing tools" and every triplex/Triplexator accuracy number** — no v3 evidence; old direct run is a tie → X1/X2, then "at least comparable" unless CIs support more.
2. **"25.2 s on chr1; 2.3× faster than triplex; ~5 min genome"** — greedy only; default score mode ~60 s ≈ triplex; genome never measured → F1b + X6.
3. **AUC 0.981, Youden ≥ 76, MCC 0.729 → 0.875 (Fig. 3B–D)** — in-sample, old data; v3 AUC degenerate → X3.
4. **"79 sequences (70/9)", "7.8:1", "70 positives with strong evidence", specificity 66.7 %** — superseded (69 = 64/5; 57 primary; 128 balanced) → after F2.
5. **"96,729 chr1 loci under default parameters"; "minimum arm 8"** — greedy, minrep 10, no filters; current defaults give 98,223 (library) / 43,150 (CLI); 57,282 (59.2 %) of the 96,729 are ≥ 80 % A/T → F1b/F3b.
6. **Figure 3 (E. coli)** — old data/code → X5, reframe as recovery.
7. **Methods:** longest-arm overlap rule, mismatch 0.15, algorithm Phase 2/3 panels, filter definitions, Data availability "upon reasonable request", "thermodynamic" wording, Streamlit mentions, "previously computationally prohibitive" (triplex does chr1 in 57.7 s).
8. **PIT elements (if in the paper):** HSeeker 24/25 needs mismatch 0.3/purity 0.8; at the benchmark settings it finds 14/25 (= triplex) — state the settings or drop (audit C, re-run by the auditor; not re-verified here).

---

## 4. Corrections to work done in this session (to fix in RTR before reuse)
1. "Library defaults" label for the minrep-8 configuration (my R2.6 draft and the WS2-PRE rewrite) — wrong; see §0.
2. DEC-9 "fast" — true for pathological inputs and small files, not for genome scans (chr1 2.3× greedy).
3. README performance sentence (written by me) is misleading for genomes.
4. The R2.1 homopolymer panel in chat used minrep 8 — state the configuration.

---

## 5. Reviewer coverage (24 comments)

| Comment | Status | Remaining (blocking items) |
|---|---|---|
| R1 | complete | — |
| R2.1 | partial | F2d, X1, X9; PDB-structure suggestion unanswered |
| R2.2 | partial | F2c (claimed fixes not applied), F3 config file, X1 |
| R2.3 | placeholder | X5, X8; say why S1-END-seq rather than Kouzine data |
| R2.4 | partial | write Limitations (WS6-A2) |
| R2.5 | placeholder | manuscript edit (WS6-A3); remove "done … verified" |
| R2.6 | partial | F4 TODOs, F3, relabel configuration, fix line refs |
| R2.7 | placeholder | X1 |
| R2.8–R2.11 | complete | — |
| R3.1 | partial | F2, X1 (PPV/NPV CIs; primary-only n=57 metrics) |
| R3.2 | partial (results provisional) | rerun after F2/F3; cite real file names |
| R3.3 | placeholder (draft claims a run that was not done) | X3 |
| R3.4 | placeholder | wording/Limitations |
| R3.5 | partial | F2a–c, F4 rationale; account for 9 removals incl. HDNA0039 |
| R3.6 | placeholder | X1, X2; delete the nBMST sentence |
| R3.7 | placeholder | X5 |
| R3.8 | placeholder | X8 |
| R3.9 | partial | X4; "shipped default" wording (CLI minrep 10 = 0.882) |
| R3.10 | partial | F4, X7; recommend a mode explicitly |
| R3.11 | partial | F1b, X6, Data availability; memory (F1c) |
| R3.12 | placeholder | manuscript wording/figures |

## 6. Ordered task list (updated 2026-10-08)

Rule: nothing that produces a reported number runs before Phase 1 is closed; anything that does
not produce a reported number starts now, in parallel.

### Phase 0: protect and coordinate: DONE 2026-10-09
- 0.1 committed `cbb4713` (merge parents `007288b` RTR, `d237bca` main) and pushed `origin/RTR-merge-main`.
- 0.2 PRs to main for Nikol: **#19** memory fix (`fix/parallel-batch-memory`, chr1 peak RSS 2,708 → 1,978 MB greedy, 2,779 → 2,043 MB score; output sha256 identical) and **#20** optional score-mode speed-up (`perf/lazy-score-selection`, chr1 score 193.0 → 60.1 s, greedy 29.4 → 26.5 s; output identical on chr1 and on 5,152 `scan_sequence` calls). The question about `fix/hdna-maximal-representation-v3` still has to be asked (draft message given to Kimon).

### Phase 0: protect and coordinate (original task list)
| # | Task | Why first | Owner |
|---|---|---|---|
| 0.1 | Commit `RTR-merge-main` (merge of main `d237bca` + RTR + DEC-6…10 + memory fix) and push the branch | 65 uncommitted files; the session scratchpad was wiped twice, losing intermediate data | Kimon |
| 0.2 | Send Nikol: (a) the PR #18 memory-regression fix (2,704 → 1,979 MB, output-neutral) and, optionally, the score-mode optimisation (verified identical, incl. his golden tests); (b) the question whether `fix/hdna-maximal-representation-v3` ships | (b) changes the detector, so it gates every rerun | Kimon → Nikol |

### Phase 1: freeze decisions (authors; in this order)
| # | Decision | Blocks |
|---|---|---|
| 1.1 | **Code freeze**: merge or drop Nikol's maximal-representation branch; webapp filters (apply CLI filters or fix the About text) (F1a, F1d) | every number; `algorithms.tex` Phase 1 |
| 1.2 | **Benchmark freeze** (WS0-A1): DEC-2 (FXN), DEC-3 (Hanvey (TTC)8; if added, +1 negative), apply the corrections the drafts claim (pCON terminal C; citation identifiers), homopolymer panel yes/no (F2) | all accuracy numbers; Supp. Table 1 |
| 1.3 | **Configuration freeze** (WS2-A1): DEC-5 mismatch (recommend 0.10); reporting configuration vs shipped minrep (10 → 8?); comparator arguments → write `analysis/benchmark_config.json` (F3) | X1–X8 |
| 1.4 | **Genome-validation dataset** for R2.3/R3.8: Kouzine 2017 (mouse B cells, KMnO4/S1 footprinting; the dataset the reviewer names; PMID 28237796) as primary, S1-END-seq (human KM12, hg19) optional | X8 |

### Phase 2: start now, in parallel (no reported numbers)
| # | Task | Notes |
|---|---|---|
| 2.1 | Environment: conda R + `bioconductor-triplex` 1.50.0 (verified available); build Triplexator 1.3.2 from source (not on conda); build NeSSie and nBMST (GitHub reachable); download E. coli K-12, hg38, the mouse build used by Kouzine, Kouzine data | longest lead time |
| 2.2 | Write (not yet run for the record) the analysis scripts: out-of-fold score threshold (X3), ablation runner (X4), overlap-comparison deposit (X7), homopolymer panel (X9), Supp. Table 1 generator (X10), Kouzine enrichment pipeline with composition-matched controls (X8) | test on current data; rerun after Phase 1 |
| 2.3 | Provenance: pin `-overlap-strategy` (and every parameter) in `benchmarks/*.py`; record git sha + overlap mode + config in every results `metadata.json` | needed for R3.11 |
| 2.4 | Stale-claim sweep that does not depend on Phase 1 (§3.2–3.4, §3.6–3.7): README remainder, "left-arm AT" docstrings, CLI help, copilot-instructions, webapp About page (except the filter question), archive June results/notebook outputs to `analysis/results/archive_2026-06/` with a README | `algorithms.tex` waits for 1.1 |

### Phase 3: runs after the freeze (in this order)
| # | Run | Answers | Est. |
|---|---|---|---|
| 3.1 | Regenerate WS2-PRE, WS4-A3 (robustness), WS2-A3 (grouped) on the frozen benchmark/config | R2.6, R3.2, R3.9 | 1 h |
| 3.2 | **X1** direct benchmark HSeeker + triplex + Triplexator with Wilson and family-cluster bootstrap CIs → `summary_metrics_v2.csv` | R2.1, R2.7, R3.1, R3.6 (highest-risk claim) | 0.5–1 d |
| 3.3 | **X2** NeSSie, nBMST with stated hit→call rule | R3.6 | 0.5 d |
| 3.4 | **X3** out-of-fold threshold, **X4** ablation | R3.3, R3.9 | 0.5 d |
| 3.5 | **X5** E. coli injected arm (greedy pinned) + overlap-threshold sweep → regenerate Fig. 3 as a recovery test | R2.3, R3.7 | 0.5–1 d |
| 3.6 | **X8** Kouzine enrichment (mouse genome scan in the frozen config) | R2.3, R3.8 | 1–2 d |
| 3.7 | **X6** runtime/memory: chr1 1–16 workers + full hg38 (greedy default), optionally comparators on the same CPU | R3.11 | 0.5 d |
| 3.8 | **X7, X9, X10** deposits; chr1 dot plot regenerated with its generator committed | R3.10, R2.1, R3.5 | 0.5 d |

### Phase 4: text and release
| # | Task |
|---|---|
| 4.1 | One regeneration pass of RTR evidence/ledger/CHECKPOINT/graph with the final numbers; resync graph statuses; regenerate the graph PNG |
| 4.2 | Manuscript: Methods (configuration, filters, overlap rule, `algorithms.tex` panels), Results (benchmark with CIs, comparators, runtime labelled by mode, genome validation), Fig. 3 reframed, Limitations (R2.4/R3.4), Data availability (R3.11), wording (R3.12), Supp. Table 1 |
| 4.3 | Response letter from Part C; every number traced to a deposited file; delete the nBMST sentence unless X2 supports it |
| 4.4 | Final Part D consistency check; tests; release (version bump, wheels with the memory fix) and Zenodo archive of code + data + benchmark outputs |
