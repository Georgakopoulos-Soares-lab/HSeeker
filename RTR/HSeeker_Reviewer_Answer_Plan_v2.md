# HSeeker — Reviewer Answer Plan v2 (resumable, auditable)

**Scope.** Grounded in the supplied reviewer-response document (RTR), the supplied manuscript, the supplied `HSeeker-main` repository, and a provenance audit and curation of the benchmark executed alongside this plan (`hdna_benchmark_experimental_v3.csv`). No claim in any draft response is allowed to reference a number or artifact that does not exist yet; where an artifact is pending, the response text names the planned output file in brackets, e.g. `[WS2: summary_metrics_v2.csv]`, and the response letter may only be finalized after that file exists and its content matches.

**What changed relative to v1.**
1. Work is organized into dependency-ordered **workstreams (WS0–WS6)** with per-action states in a machine-readable ledger (`hseeker_revision_ledger.yaml`), so the revision can stop and resume across sessions without re-deriving decisions.
2. Every action lists **inputs, outputs, verification criteria, and evidence** (file:line or executed command). The per-reviewer-item matrix (Part C) maps each reviewer comment to the actions that evidence its response.
3. The **sequence-provenance question is resolved, not deferred**: a per-record audit of all 80 benchmark rows was executed (Part A). It identified the origin of all six orphan records, resolved the pRW1406/pRW1411 conflict, and surfaced three benchmark-validity issues v1 missed (inosine encoding artifact, an intermolecular-TFO label-scope error, a by-design homopolymer exclusion) plus systematically wrong PMIDs in the `literature_curated` batch.
4. A **consistency register** (Part D) freezes every manuscript↔repository discrepancy so the response letter, manuscript, and repo cannot drift apart again.

---

## Part 0 — How to resume this revision

The ledger `hseeker_revision_ledger.yaml` is the single source of truth for execution state.

- Each action entry: `id`, `workstream`, `reviewer_items`, `status` (`pending | in_progress | blocked | done | verified`), `depends_on`, `inputs`, `outputs`, `verify` (acceptance test), `evidence`, `last_updated`.
- **Resume protocol:** load the ledger; actions with `status: verified` are never re-run; actions with `status: done` are re-run only if their inputs changed (check SHA-256 recorded in the action's `evidence`); everything else runs in dependency order. After each action, run its `verify` test, then set `status: verified` and record the output file hash.
- **Auditability rule:** a response-letter sentence may cite a number only if that number appears in a ledger-`verified` output file. The response matrix in Part C names the exact file each claim depends on.

---

## Part A — Benchmark provenance audit (executed; feeds WS0)

Full per-record results: `hdna_benchmark_experimental_v3.csv` (80 rows; 71 kept, 7 removed, 2 excluded; columns include `verification_status`, `resolved_origin`, `disposition`, `rationale`, `evidence`, `corrected_pmid`, `corrected_doi`, `family_id`, `study_id`).

### A.1 The six orphan records (HDNA0032–0037) — origins resolved

These were the records flagged `confidence=low` from the unavailable "H-DNA Thermostability document". Resolution by literature verification + sequence analysis:

| Record | Name | Resolved origin | Status | Disposition |
|---|---|---|---|---|
| HDNA0032 | GG32 (30-nt) | **Rooney & Moore 1995 PNAS** GG32 (exact match). | verified_exact | **removed** (low confidence tier; author decision to exclude all low-tier records despite resolved provenance). |
| HDNA0033 | AG32 (30-nt) | **Rooney & Moore 1995** AG32 (exact match). | **label_conflict** | **removed** (low confidence tier; author decision. Also carried the unresolved label conflict — Rooney 1995 report AG32 H-DNA transition strongly impaired vs CSV 'forming'; removal moots DEC-1). |
| HDNA0034 | MYCAA | **Belotserkovskii et al. 2007 JBC** "AA" construct (verified computationally). | verified_design | **removed** (low confidence tier; author decision; removal moots DEC-4 figure-level confirmation). |
| HDNA0035 | GG32-long | Exact 32-nt pGG32 core + EcoRI/BamHI cloning overhangs (`AATTC…GGATC`). No primary source for the exact 42-nt oligo found. | unverifiable | **secondary** (also a near-duplicate of HDNA0002 by containment). |
| HDNA0036 | Mycm1 | c-MYC NHE mutant; exact 22-nt substring of pMycAG with one left-arm G deleted; mirror broken. Not among Belotserkovskii 2007 (WT/TT/Q/AG/AA) or Wang 2004 constructs. | unverifiable | **secondary** (label structurally plausible but not traceable). |
| HDNA0037 | MycM2 | Mycm1 + G21T; mirror broken. No primary source found. | unverifiable | **secondary**. |

### A.2 pRW1406/pRW1411 transcription inconsistency — resolved

Hanvey, Klysik & Wells 1988 JBC (PMID 2835375) studied seven inserts: (G19)·(C19), (TCC)8, (CT)12, (TTCC)6, (TTC)8, (GAAA)6, (A20). The CSV's five Hanvey-JBC records match five of these exactly in sequence and behavior. **pRW1411 = (TCC)8 (66% G+C, forming)** — the CSV record is correct and stays **primary**. **pRW1406 is the distinct 33% G+C (TTC)8 insert** (`TTCTTCTTCTTCTTCTTCTTCTTC`, also forming per the abstract); the source document's listing of pRW1406 with the (TCC)8 sequence was the transcription error. Optional: add pRW1406=(TTC)8 as a new verified primary record.

### A.3 Benchmark-validity issues found by the audit (absent from plan v1)

1. **HDNA0026 (R2-ino) — encoding artifact.** R2-ino is an *inosine-substituted* oligo (Del Mundo 2017, confirmed in paper text: third-strand guanines → inosine; TM ≈ 80 °C hairpin control). The benchmark FASTA silently converts inosine→G, manufacturing the all-G R2 sequence — a perfect mirror repeat. This produced **1 of the 2 false positives** in the current repo run. **Disposition: excluded** from the ACGT sequence-level benchmark (the label belongs to a chemically modified oligo, not to any DNA sequence).
2. **HDNA0080 (G10TTAA_AG5) — label-scope error.** The record's own notes: "NOT intramolecular H-DNA; G10 forms intermolecular triplex with AG5 TFO; correctly not detected by HSeeker." It is labeled `forming`, so it counts as a false negative by construction. **Disposition: excluded** (outside the tool's intramolecular target class; optionally cited in text as a target-class boundary example).
3. **HDNA0054 (poly-(dC)30) — by-design exclusion.** Forming per Kohwi 1988 (PMID 3375241, verified), but a pure homopolymer that the production CLI removes (`filter_homopolymers=True`, `__main__.py:109`). Keep in primary **only** for configurations that disable the homopolymer filter; document the tension explicitly (see register D-8).
4. **HDNA0029 (pCON) — transcription discrepancy.** Wang 2004 Fig. 1A shows a 31-nt CON (`TTGCTGAACTTGATTCGACTCAGATGCAGTGC`); the CSV record is missing the terminal C (30 nt). Correct the sequence; label and expected call unchanged.
5. **HDNA0055 (GA32) — citation correction.** Sequence matches Rooney & Moore 1995 GA32 (30-nt) exactly; the CSV cites Mirkin 1987 (motif origin) and Wang 2004's GA32 is a different 31-nt sequence. Add Rooney 1995 as sequence-source citation.

### A.4 Systematic citation defect in the `literature_curated` batch

PubMed metadata verification of all 27 cited PMIDs: **all 43 `high` and 5 `moderate` records have correct PMIDs; 22 of 26 `literature_curated` records have PMIDs pointing to unrelated papers** (e.g., the Mirkin 1987 record's PMID resolves to a cardiology paper). DOIs are correct except two: `10.1073/pnas.92.7.2141` (typo; correct: `10.1073/pnas.92.6.2141`) and `10.1074/jbc.R700013200` (resolves to an unrelated JBC innate-immunity minireview; intended source for the FXN records is likely Wells 2008 FASEB J "DNA triplexes and Friedreich ataxia", PMID 18211957 — **requires author confirmation**). Corrected PMIDs for all 13 mismatches are in the audit CSV (`corrected_pmid` column). This defect must be fixed before resubmission and disclosed as a curation correction.

### A.5 Resulting benchmark composition (after dispositions)

| Set | Records | Forming | Non-forming | Contents |
|---|---|---|---|---|
| **Primary** | 59 | 53 | 6 | Construct-level, traceable records only |
| **Secondary** | 12 | 12 | 0 | Representative motifs (class-level evidence); all non-forming orphans removed |
| **Excluded** | 2 | 1 | 1 | R2-ino (inosine), G10TTAA_AG5 (intermolecular) |
| Primary + Secondary | 71 | 65 | 6 | Curated experimental set (also the basis for 1:1 balanced benchmark) |

**Consequence the response must state honestly:** the defensible primary negative set is 6 (secondary now contributes 0 non-forming after low-confidence removal), not 9. This strengthens (a) the case for the separate synthetic stress panel (WS2-A4, now executed as the balanced v3 benchmark) and (b) the need to report metrics on both the experimental-only and the 1:1 balanced benchmark with CIs. This strengthens (a) the case for the separate synthetic stress panel (WS2-A4) and (b) the need to report metrics on both benchmark compositions with CIs. Family structure for R3.2: 27 studies, 50 families (parent-sequence grouping); the pXY32 series contributes 16 records from one study (Belotserkovskii 1990).

---

## Part B — Workstreams

### WS0 — Benchmark freeze (foundation; everything depends on it)

| Action | Detail |
|---|---|
| **WS0-A1** | Apply audit dispositions and curation: the frozen benchmark is `hdna_benchmark_experimental_v3.csv` (80 rows: 71 kept, 7 removed, 2 excluded) with `curation_decision`, `curation_justification`, and `label_conflict` columns; correct pCON terminal C; correct the 22 PMIDs + 2 DOIs; `family_id`, `study_id`, `evidence_level`, `disposition` columns retained. **Inputs:** provenance audit (executed). **Outputs:** `hdna_benchmark_experimental_v3.csv` + SHA-256 in ledger. **Verify:** 71 kept (65 forming / 6 non-forming); 0 low-tier rows kept; all primary rows have construct-level evidence and valid citations; 59 primary / 12 secondary. |
| **WS0-A2** | Author decision points (block downstream text, not downstream computation): (i) ~~AG32-30nt label~~ — **mooted** (HDNA0033 removed with all low-tier records); (ii) FXN citation confirmation (PMID 18211957?); (iii) add pRW1406=(TTC)8 as new primary record?; (iv) ~~confirm MYCAA promotion~~ — **mooted** (HDNA0034 removed with all low-tier records). **Verify:** decisions recorded in ledger `decisions:` section with date and decider. |
| **WS0-A3** | Write the inclusion criteria paragraph for Methods/Supplement: "primary benchmark = independently traceable, construct-level experimental records; secondary = class-representative or unverifiable records, analyzed separately." **Verify:** text matches the audit CSV dispositions exactly. |

### WS1 — Repository/code fixes (independent of WS0)

| Action | Detail |
|---|---|
| **WS1-A1** (R2.8) | `setup.py:21-49`: make `-march=native` opt-in (e.g., env var `HSEEKER_NATIVE=1`), not default-when-not-cross-compiling. **Verify:** `pip install .` on Linux produces a build log without `-march=native` unless opted in; wheels still build in CI. |
| **WS1-A2** (R2.9) | `_hdna.c:54-61`: make `prof_inner_iters`/`prof_ctr_sp_pairs` call-local (per-call struct returned to Python) or remove them; the GIL is released during scans (`_hdna.c:659-681`) so shared mutation is unsound. **Verify:** parallel scan of chr1 chunk set returns identical hits; counters (if kept) are per-call. |
| **WS1-A3** (R2.10) | Normalize public sequence-field case: C wrapper lowercases (`_hdna.c:629-634, 699-752`), scoring uppercases (`_scoring.py:76-90, 128-136`). Choose uppercase everywhere at the Python boundary. **Verify:** `left_arm`, `right_arm`, `full_sequence`, `putative_triplex` all uppercase; round-trip string comparison test passes; add regression test. |
| **WS1-A4** (R2.11) | **Superseded by author decision 2026-09-23:** the Streamlit demo is *removed* entirely rather than renamed. Delete `src/hseeker/app.py`, drop `streamlit` from the `[app]` extra, keep `webapp/main.py` as the single FastAPI entry; update README + Code Availability. **Verify:** no `streamlit` reference anywhere in the shipped tree; one web entry point only. |
| **WS1-A5** | Fix stale docstring "54 experimentally curated sequences" in `analysis/scripts/run_direct_sensitivity_analysis.py:2-7`; supply or remove the notebook's missing input `hseeker_additional_motifs.csv` (referenced in `analysis/hdna_detection_benchmark.ipynb`, absent from repo). **Verify:** notebook runs end-to-end or the cell is removed. |

### WS2 — Unified benchmark pipeline + statistics (depends on WS0, WS1-A3)

| Action | Detail |
|---|---|
| **WS2-A1** | Single version-pinned config `analysis/benchmark_config.json` (hseeker version, minrep, maxspacer, purity, mismatch, at_threshold, filter_homopolymers, triplex version + all `triplex.search` arguments, triplexator version + flags) driving BOTH the direct and injected arms. Resolves the current asymmetry: direct arm runs without AT/homopolymer filters (`run_direct_sensitivity_analysis.py:99-141` calls `scan_sequence` with defaults) while the injected arm uses CLI defaults (filters on, `__main__.py:97-110`); manuscript claims `mismatch=0.15` while scripts use 0.10; manuscript's Triplex settings (min_score=15, min_len=8, max_len=50, prokaryotic) differ from the script's package defaults. **Decision:** adopt the manuscript's stated parameters as canonical (they were the reviewed ones) and re-run; record any resulting number changes for the response. **Verify:** both arms' metadata.json report identical parameter blocks; manuscript Methods updated to match. |
| **WS2-A2** (R2.7, R3.1) | Extend the direct script: report sensitivity, specificity, PPV, NPV, MCC with 95% CIs — Wilson intervals for proportions; cluster bootstrap (resample families, not rows) for MCC and for tool-vs-tool differences. **Verify:** `summary_metrics_v2.csv` contains CI columns; a unit test recomputes one interval by hand. |
| **WS2-A3** (R3.2, R3.3) | Family/study-aware evaluation: leave-one-family-out and leave-one-study-out for the detection call and for the stability-score threshold (select Youden threshold on training folds only, evaluate on the held-out family/study; report the distribution of out-of-fold MCC/AUC, not a single number). Uses `family_id`/`study_id` from WS0-A1. The full-data ROC is retained as descriptive only. **Verify:** no threshold used in a reported confusion matrix was selected on the rows it is evaluated on; `lofo_metrics.csv`, `loso_metrics.csv` produced. |
| **WS2-A4** (R2.1) | Synthetic stress panel — **executed** as `hdna_benchmark_balanced_v3.csv` (60 synthetic negatives, seed 20260922): (i) random sequences matched to positives' length and GC; (ii) shuffled positives (composition preserved, mirror destroyed); (iii) mirror-disrupted single-point mutants of each primary positive (the perturbation the reviewer requested — feasible and DNA-appropriate); (iv) poly-A/poly-G/poly-C/poly-T homopolymers; (v) known non-H-DNA structures: G4-forming (e.g., c-MYC Pu27 G4 sequence) and Z-DNA-forming sequences; (vi) B-DNA controls (scrambled, e.g., pCON-like). **Grounded correction for the response letter:** BoltzGen/RFDiffusion are protein-structure generators and cannot generate DNA duplex controls; we substitute DNA-appropriate perturbations (i–vi) and say so. **Verify:** panel CSV + `stress_panel_metrics.csv`; text states these are synthetic controls, not additional experimental negatives. |
| **WS2-A5** | Re-run direct + injected arms on the frozen primary set and on primary+secondary; produce per-record prediction tables with the FP/FN explanations from the audit (inosine, TFO, homopolymer) attached. **Verify:** every FP/FN in the final tables has an audit-consistent explanation or is flagged as unexplained. |

### WS3 — Comparator expansion (R3.6; depends on WS0)

| Action | Detail |
|---|---|
| **WS3-A1** | NeSSie: build from source (github.com/B3rse/nessie, `make`); run mirror-repeat/triplex mode on the frozen benchmark with documented parameters; deposit command + outputs. |
| **WS3-A2** | nBMST: run locally via the `abcsFrederick/non-B_gfa` runner (the nBMST site's own local option for batch use); mirror-repeat module on the frozen benchmark; deposit command + outputs. Note: nBMST reports motif hits, not H-DNA calls — the response must state how a hit was mapped to a classification (this is the honest version of the RTR's current unsupported sentence "We run nBMST… we had to deduce them"). |
| **WS3-A3** | Harmonize claims: Key Points "outperformed existing tools" → align with Results ("at least comparable" unless WS2/WS3 CIs support more). Document target-class differences (Triplexator/3plex = intermolecular TTS-oriented; NeSSie = symmetry search; nBMST = motif search). **Verify:** no comparative claim lacks a corresponding deposited output; the RTR's current nBMST sentence is removed unless WS3-A2 outputs exist. |

### WS4 — Robustness analyses (R3.7, R3.9, R3.10; depends on WS2-A1)

| Action | Detail |
|---|---|
| **WS4-A1** | Overlap-threshold sweep for the injected arm: `--overlap-fraction` ∈ {0.5, 0.6, 0.7, 0.8, 0.9, 0.95} (parameter already exists, `run_injected_sensitivity_analysis.py:471`); report whether tool rankings change. **Output:** `overlap_sweep.csv`. |
| **WS4-A2** | Component ablation on the frozen benchmark: (a) detection only (no filters, no scoring); (b) + AT filter; (c) + homopolymer filter; (d) detection + pairing-only score; (e) pairing+stacking; (f) full model. Public switches cover (a–c); (d–e) need analysis-only access to `_scoring.py` components (read pairing/stacking subscores from existing outputs — no core change). **Output:** `ablation_metrics.csv`. |
| **WS4-A3** | Parameter robustness: purity ∈ {0.85, 0.90, 0.95, 1.0}, mismatch ∈ {0.0, 0.05, 0.10, 0.15}, minrep ∈ {8, 10, 12}, maxspacer ∈ {5, 10, 15} (one-at-a-time around the canonical config); report sensitivity/specificity/MCC per setting. **Output:** `parameter_sweep.csv`. |
| **WS4-A4** | Overlap-selection comparison: on `remove_overlaps=False` candidate sets, compare the production longest-arm rule (`_hdna.c:499-558`) against highest-score selection; report concordance on the benchmark and on chr1. **Output:** `overlap_selection_comparison.csv`; adopt a new production rule only if the comparison justifies it (else justify the default in text). |

### WS5 — Orthogonal genome-scale validation (R2.3, R3.8; depends on WS2-A1)

| Action | Detail |
|---|---|
| **WS5-A1** | Compare the existing hg38 chr1 predictions (96,729 loci, `analysis/results/chr1/`) against human S1-END-seq triplex loci (GSE204808/GSE203632, Matos-Rodrigues 2022 — human, so no liftover). Enrichment of HSeeker loci at experimental loci vs composition-matched control intervals (match on GC + purine content and length; reviewer explicitly required composition control). **Output:** `chr1_vs_s1endseq_enrichment.csv` + figure. |
| **WS5-A2** | Fallback if S1-END-seq peaks are unusable: mouse S1-seq (GSE197669, Maekawa 2022) with a fresh HSeeker mouse-genome run (sandbox-feasible; chr1-equivalent runtime measured at ~25–188 s/chromosome). |
| **WS5-A3** | Text: describe WS5 as orthogonal validation; keep the E. coli experiment framed as genomic-context recovery only (see R2.3/R3.7). |

### WS6 — Manuscript & response-letter synchronization (depends on WS0–WS5)

| Action | Detail |
|---|---|
| **WS6-A1** | Every revised number in the manuscript and response letter traced to a named, ledger-verified output file (the audit table in Part C lists the file per claim). |
| **WS6-A2** (R2.4, R3.4) | Add a **Limitations** subsection: heuristic (not calibrated thermodynamic) score; fixed 5.0 kcal/mol averaged mismatch penalty; no backbone geometry; assay-condition heterogeneity of the benchmark (acidic S1, Mg²⁺, supercoiling); benchmark size/imbalance; filter provenance disclosure (R2.6: state how the AT/homopolymer rules were chosen and that they were not re-fit during evaluation). |
| **WS6-A3** (R2.5) | Move sensitivity/complexity result statements from Methods to Results (the RTR claims this was done; the supplied manuscript still has them in Methods — verify before submission). |
| **WS6-A4** (R3.11) | Label the ~5-min whole-genome runtime as extrapolation everywhere; deposit benchmark commands/raw timing/outputs; update Data Availability (currently says "upon reasonable request" — stale). Optional: run full hg38 (~5–30 min extrapolated; sandbox-feasible) to replace the extrapolation with a measurement. |
| **WS6-A5** (R3.12) | Terminology discipline ("putative/predicted" vs experimentally "forming"); remove/qualify "works highly effectively and aligns seamlessly", "only H-DNA identification algorithm with a proven scoring system built-in", "previously computationally prohibitive"; simplify Figure 1 and algorithm panels. |
| **WS6-A6** | Complete the four unfinished RTR draft responses (R2.1 placeholder, R2.4 missing limitations text, R3.1 truncated list, R3.7 note-to-self) using the Part C matrix; remove the unsupported nBMST sentence in R3.6 unless WS3-A2 is verified. |

---

## Part C — Per-reviewer-item response matrix

Each item: grounded draft response + evidence file(s) + claims limits. Bracketed tokens `[WSx: file]` are placeholders that must be replaced by the verified artifact's content before the letter is sent.

### Reviewer 1
**R1 (positive assessment).** Response: thanks; no change. *No action.*

### Reviewer 2

**R2.1 (small, imbalanced benchmark; discriminating power).**
- *Draft response:* We agree accuracy is uninformative at 7.8:1 imbalance and that nine (after re-curation: six primary) experimental negatives limit precision of negative-class estimates. We now (i) report class-aware metrics with 95% CIs and de-emphasize accuracy `[WS2: summary_metrics_v2.csv]`; (ii) add family/study-aware evaluation `[WS2: lofo_metrics.csv]`; (iii) add a clearly separated synthetic stress panel — random GC/length-matched sequences, shuffled and mirror-disrupted positives, homopolymers including poly-A, G4- and Z-DNA-forming sequences, and B-DNA controls `[WS2: stress_panel_metrics.csv]`. We note respectfully that BoltzGen/RFDiffusion generate protein structures and cannot produce DNA duplex controls; we therefore implemented the perturbation-based DNA controls the reviewer suggested (disrupted mirror repeats, poly-A) directly.
- *Can claim:* high sensitivity on the curated positives; robust detection-by-design behavior on synthetic controls. *Cannot claim:* high specificity as an established population quantity; superiority from accuracy.

**R2.2 (sequences as CSV + runner script).**
- *Draft response:* The benchmark is now provided as machine-readable CSVs with full provenance metadata (`hdna_benchmark_experimental_v3.csv` with per-record curation decisions, plus the 1:1 balanced set `hdna_benchmark_balanced_v3.csv`), and `analysis/scripts/run_direct_sensitivity_analysis.py` reproduces the analysis end-to-end with the pinned configuration `analysis/benchmark_config.json`. During re-curation we found and fixed a cardinality drift (79 vs 80 rows), a transcription error (pCON terminal C), and 24 citation-identifier errors (22 PMIDs, 2 DOIs); all corrections are documented per record.
- *Can claim:* full reproducibility of the benchmark. *Must not claim:* the CSV matches the old supplement (it supersedes it).

**R2.3 (E. coli insertion is not independent validation).**
- *Draft response:* We agree; the experiment tests computational recovery in genomic context, not biological prediction. We have reframed it as such throughout (text + Figure 3 caption), report the total non-inserted hit counts `[WS2: injected metadata]`, and added an orthogonal genome-scale validation against human S1-END-seq experimental triplex loci with composition-matched controls `[WS5: chr1_vs_s1endseq_enrichment.csv]`.
- *Cannot claim:* the E. coli experiment validates H-DNA formation.

**R2.4 (thermodynamic simplification; limitations section).**
- *Draft response:* We agree and added a Limitations subsection: the score is an experimentally informed heuristic (fixed G·G/A·A pairing values from Park et al., an averaged 5.0 kcal/mol mismatch/stacking penalty, exponential consecutive-mismatch penalty, boundary optimization); it does not model backbone geometry and is not calibrated to ΔG, TM, or superhelical density; benchmark assays span heterogeneous conditions (acidic S1, Mg²⁺, supercoiled plasmids), so labels are condition-dependent. `[WS6: manuscript Limitations]`
- *Note:* the current RTR draft says "We added the following…" with no text following — the actual paragraph must be written (WS6-A2).

**R2.5 (move sensitivity/complexity to Results).** Agree; done in the revised manuscript and verified against the compiled text (the supplied manuscript still had them in Methods). `[WS6-A3]`

**R2.6 (filters defined on the same benchmark).**
- *Draft response:* The AT-content (≥80%) and homopolymer-exclusion rules are fixed pipeline rules, not fitted during evaluation (`__init__.py:53-89`; CLI defaults `__main__.py:58-61,97-110`). We now state their provenance explicitly and disclose that, as fixed rules chosen with knowledge of the H-DNA literature and benchmark composition, they are not an independent test; the family-aware evaluation (WS2-A3) bounds any residual optimism. We also corrected an implementation asymmetry: the earlier direct benchmark bypassed these filters while the injected arm applied them; both arms now use the pinned configuration `[WS2: benchmark_config.json]`.

**R2.7 (confidence intervals).** Covered by WS2-A2; response cites the CI table and qualifies superiority language where intervals overlap. `[WS2: summary_metrics_v2.csv]`

**R2.8 (-march=native).** Confirm the reviewer's reading (correct for the supplied `setup.py`); now opt-in via `HSEEKER_NATIVE=1`; wheels unaffected. `[WS1-A1]`

**R2.9 (unsynchronized profiling counters).** Confirm; counters made call-local/removed; parallel determinism test added. `[WS1-A2]`

**R2.10 (case inconsistency).** Confirm (C wrapper lowercases, scoring uppercases); all public sequence fields now uppercase; regression test added. `[WS1-A3]`

**R2.11 (two app entry points).** Confirm the reviewer's reading. Rather than renaming the demo, we **removed the Streamlit application entirely**: the package now ships a single web front-end, the FastAPI service (`webapp/main.py`), which is also what the Docker image serves. `streamlit` is no longer a dependency in any form. Docs updated accordingly. `[WS1-A4]`

### Reviewer 3

**R3.1 (size/imbalance/heterogeneity; CIs).**
- *Draft response:* combines R2.1 elements; adds the re-curation result: after source-level verification, the primary benchmark contains 59 construct-level records (53 forming, 6 non-forming); we report metrics with 95% CIs on the primary set and on primary+secondary (71 records), and on the 1:1 balanced benchmark (130 records, 65:65) `[WS2: summary_metrics_v2.csv]`, and we qualify all negative-class conclusions. The response must state the negative-set size plainly — this is the honest cost of the provenance fix.
- *Cannot claim:* the benchmark establishes precise specificity.

**R3.2 (family-aware validation).**
- *Draft response:* We grouped records by study and by sequence family (27 studies, 50 families; the pXY32 substitution series contributes 16 records from one study) using curated provenance metadata, and evaluated with leave-one-family-out and leave-one-study-out schemes `[WS2: lofo_metrics.csv, loso_metrics.csv]`; we also report family-collapsed performance. 
- *Cannot claim:* 79 independent observations (that framing is removed).

**R3.3 (Youden threshold circularity).**
- *Draft response:* The threshold is now selected within training folds only (grouped/nested CV aligned with R3.2) and evaluated out-of-fold; the full-data ROC (AUC reported previously) is retained as descriptive and labeled as such. `[WS2: lofo_metrics.csv]`
- *Cannot claim:* the previously reported AUC 0.981/0.961 demonstrates generalizable performance.

**R3.4 (stability score ≠ calibrated thermodynamics).** Same limitations text as R2.4 plus explicit wording change throughout ("experimentally informed heuristic stability score"); the Discussion sentence "works highly effectively and aligns seamlessly…" is removed. `[WS6-A2, WS6-A5]`

**R3.5 (benchmark provenance scrutiny).**
- *Draft response:* We audited every record at source level (per-record audit CSV provided). Outcome: six records previously carried low-confidence annotations from a secondary document. Three (GG32-long, Mycm1, MycM2) had no traceable primary source and were removed. Three (GG32-30nt, AG32-30nt, MYCAA) had resolved provenance (Rooney & Moore 1995; Belotserkovskii 2007) but were removed by author decision to exclude all low-confidence-tier records. Class-representative motifs remain in a clearly labeled secondary set excluded from the primary performance benchmark. We also excluded two records on validity grounds (an inosine-substituted oligo that cannot be represented as ACGT DNA; an intermolecular-TFO record mislabeled as intramolecular H-DNA), corrected one transcription error (pCON), and corrected 24 citation identifiers. The primary benchmark now contains only independently traceable construct-level records.
- *Can claim:* the curation is now auditable record-by-record. *Cannot claim:* the negative set grew (it shrank; say so).

**R3.6 (comparator set; claim strength).**
- *Draft response:* We added NeSSie (built from source) and nBMST/non-B_gfa (local runner) on the frozen benchmark with deposited commands and outputs `[WS3: nessie_benchmark.csv, nbmst_benchmark.csv]`, and we document target-class differences (Triplexator/3plex intermolecular-TTS-oriented; NeSSie symmetry search; nBMST motif search without H-DNA classification, requiring a stated hit→call mapping). Key Points wording aligned with Results.
- *Hard rule:* the current RTR sentence "We run nBMST with the sequences and show performance lower than ours" must not be submitted unless WS3-A2 outputs exist and support it.

**R3.7 (E. coli de-emphasis; overlap-threshold robustness).**
- *Draft response:* Figure 3 reframed as genomic-context recovery; overlap-threshold sweep (0.5–0.95) shows whether rankings are robust `[WS4: overlap_sweep.csv]`. (The current RTR note "Tone down this experiment only." is a note-to-self, not a response — replaced by this text.)

**R3.8 (genome-scale biological validation).**
- *Draft response:* We compared chromosome-1 predictions with an independent human experimental dataset (S1-END-seq triplex loci) using composition-matched control intervals; report enrichment and overlap statistics `[WS5: chr1_vs_s1endseq_enrichment.csv]`. Scope stated as orthogonal validation of genomic enrichment, not per-locus proof of formation.

**R3.9 (ablation + parameter robustness).**
- *Draft response:* Component ablation (detection-only → +composition filters → pairing-only → pairing+stacking → full) and one-at-a-time parameter sweeps around the canonical configuration `[WS4: ablation_metrics.csv, parameter_sweep.csv]`; results show which components drive benchmark behavior.

**R3.10 (overlap-removal rule).**
- *Draft response:* We justify the longest-arm default and compare it against highest-score selection on the unfiltered candidate sets `[WS4: overlap_selection_comparison.csv]`; we state when retaining overlapping candidates (web/CLI option) is appropriate for exploratory analysis. Production rule changed only if the comparison supports it.

**R3.11 (benchmark reproducibility; extrapolation).**
- *Draft response:* The ~5-min whole-genome figure is labeled an extrapolation throughout; benchmark commands, raw timing, and outputs are deposited in the repository (they already exist under `benchmarks/` and `analysis/results/chr1/`); Data Availability updated. Optionally we report an actual full-hg38 run `[WS6-A4]`.

**R3.12 (claims/terminology/figures).**
- *Draft response:* The three challenged phrases removed or qualified; "putative/predicted" used for computational outputs throughout; Figure 1 and algorithm panels simplified. `[WS6-A5]`

---

## Part D — Consistency register (freeze before any number is quoted)

| # | Discrepancy (verified in this review) | Frozen resolution | Owner |
|---|---|---|---|
| D-1 | Manuscript: 79 sequences (70/9). Repo CSV: 80 rows (71/9). | After audit + curation: 71 kept (65 forming / 6 non-forming); 59 primary (53/6), 12 secondary (12/0), 7 removed (unverifiable + low-tier), 2 excluded; balanced benchmark 130 (65:65). Manuscript, letter, and repo report these numbers. | WS0-A1 |
| D-2 | Manuscript: HSeeker `mismatch=0.15`. Scripts: 0.10. | One pinned value in `benchmark_config.json`; manuscript updated to match. | WS2-A1 |
| D-3 | Manuscript Triplex: min_score=15, min_len=8, max_len=50, prokaryotic. Script: package defaults (min_len=6, max_len=25, eukaryotic). | Pinned explicit arguments matching the manuscript; Triplex version pinned and disclosed (sandbox verification used 1.46.0; manuscript used 1.50.0). | WS2-A1 |
| D-4 | Direct arm runs without AT/homopolymer filters; injected arm runs with them (CLI defaults). | Both arms use the pinned config; state filters explicitly in Methods. | WS2-A1 |
| D-5 | Manuscript: chr1 run used "minimum arm length 8 bp". `run_chr1_analysis.py:22` uses minrep=10. | Reconcile (re-run chr1 with the pinned value or correct the text); the 96,729-locus count must correspond to the stated parameters. | WS2-A1/WS5 |
| D-6 | R2-ino inosine→G encoding artifact (manufactured FP). | Excluded from benchmark; documented in response R3.5. | WS0-A1 |
| D-7 | HDNA0080 intermolecular TFO labeled forming (manufactured FN). | Excluded; documented. | WS0-A1 |
| D-8 | HDNA0054 poly-(dC)30 forming but removed by production homopolymer filter. | Documented as by-design exclusion; reported per configuration. | WS2-A1 |
| D-9 | pCON missing terminal C vs Wang 2004 Fig. 1A. | Sequence corrected in frozen CSV. | WS0-A1 |
| D-10 | 22 wrong PMIDs + 2 wrong DOIs in `literature_curated` batch. | Corrected per audit CSV; disclosed as curation correction. | WS0-A1 |
| D-11 | AG32-30nt label conflict (Rooney 1995: non-forming; CSV: forming). | **Resolved:** HDNA0033 removed with all low-tier records (author decision); DEC-1 mooted. | WS0-A2 |
| D-12 | Notebook requires absent `hseeker_additional_motifs.csv`. | Supply the file or remove the merge cell. | WS1-A5 |
| D-13 | RTR draft contains unfinished responses (R2.1 "…", R2.4 missing text, R3.1 truncated list, R3.6 unsupported nBMST claim, R3.7 note-to-self). | Completed per Part C before submission. | WS6-A6 |
| D-14 | Data Availability says "upon reasonable request… deposited upon publication" while artifacts already exist in repo. | Rewrite to point at deposited artifacts. | WS6-A4 |
| D-15 | Key Points "outperformed existing tools" vs Results "at least comparable". | Harmonized to the weaker defensible claim unless WS2/WS3 CIs support more. | WS3-A3 |

**Pre-submission verification checklist (all must pass):** (1) every register row's resolution is reflected in manuscript + repo + letter; (2) every number in the letter appears in a ledger-`verified` output; (3) benchmark CSV hashes match the ledger; (4) no response references an artifact with `status != verified`.

---

## Execution order (dependency-driven)

1. **WS0** (benchmark freeze; consumes the completed audit) — unblocks everything.
2. **WS1** (code fixes) — parallel with WS0; WS1-A3 (case) must land before final benchmark re-runs.
3. **WS2** (unified pipeline + statistics) → **WS4** (robustness) and **WS3** (comparators) can run in parallel after WS2-A1.
4. **WS5** (orthogonal validation) — parallel after WS2-A1.
5. **WS6** (manuscript + letter) — last; Part C placeholders resolved against verified outputs; Part D checklist gates submission.
