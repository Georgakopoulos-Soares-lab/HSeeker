# Plan: bringing the submitted manuscript to the revised state

**Source files in this folder**

| File | Role |
|---|---|
| `Revised HSeeker_ an algorithm for systematic H-DNA sequence identification .docx` | The submitted manuscript. 177 paragraphs, 5 figures, 1 supplementary figure, 95 embedded images (inline equations and the three algorithm panels). |
| `HSEEKER_RTR (1).docx` | The response letter. |

**Output:** one edited copy of the manuscript with every change marked as a suggestion, plus
replacement figure files and supplementary tables. The letter is handled separately and is not
changed by this plan.

---

## 1. Guiding principle: minimal edit

The manuscript's structure, argument and voice stay as they are. Most paragraphs are untouched. Of
the 177 paragraphs, this plan edits **41**, deletes **6**, and adds **19**. Everything else —
the whole Introduction except two sentences, the entire candidate-discovery and scoring Methods, the
command-line section, both web-interface sections, the funding and conflicts statements, and all 95
inline equation images — is left exactly as submitted.

Where a reviewer's point can be met by changing a clause rather than a paragraph, the plan changes
the clause. Where a number has changed, only the number is replaced. New text is added only where a
reviewer explicitly asked for something that does not yet exist: a limitations section, a
family-aware validation section, an ablation section, a tool-comparison section, and a genome-wide
comparison section.

---

## 2. How the changes will be marked

Every edit goes in as a **Word tracked change** (`w:ins` / `w:del`), so the authors can walk the
document with Accept / Reject and nothing is applied silently. All changes are attributed to a
single author, `HSeeker revision`, dated, so Word renders them in one consistent colour and the
whole set can be accepted or rejected in one action if desired.

On the light blue: Word assigns tracked-change colours per author, and that assignment lives in the
reader's Word installation, not in the file. To guarantee the changes read as light blue regardless
of who opens the file, inserted runs also carry an explicit character colour of **`#4A90D9`**.

There is one consequence worth stating before we start, because it is a decision for the authors:

| Option | Behaviour | Trade-off |
|---|---|---|
| **A — tracked changes + explicit light blue** *(recommended)* | Light blue in every viewer; Accept/Reject works | After accepting, the text stays light blue until the author selects all and sets font colour to Automatic — one action for the whole document |
| **B — tracked changes only** | Accept/Reject works; colour is whatever Word assigns the author | Not guaranteed to be light blue |
| **C — light blue text only, no tracking** | Always light blue | No Accept/Reject; changes must be applied by hand |

The plan assumes **Option A**. Deleted text appears struck through in the same colour, as Word
renders tracked deletions.

Comment balloons are **not** used: a reviewer comment is not a change, and the reasoning for each
edit lives in this plan and in the response letter.

---

## 3. Figures: the renumbering

This is the largest structural change and must be settled before any text is edited, because six
in-text citations depend on it.

| Submitted | Content | Becomes | Action |
|---|---|---|---|
| Figure 1 | Workflow, panels A and B | **Figure 1** | Edit in place (§6.1) |
| Figure 2 | Parallel scaling, 3 tools, 4 panels | **Figure 2** | Replace image and caption |
| Figure 3 | Sensitivity/specificity: A E. coli confusion matrices, B ROC, C raw detection, D Youden threshold | **Figure 3** | Replace image and caption; panels B and D are withdrawn |
| — | *(new)* | **Figure 4** | Insert: comparison with single-stranded DNA maps |
| Figure 4 | Web interface | **Figure 5** | Renumber only |
| Figure 5 | Interactive visualisation | **Figure 6** | Renumber only |
| Supplementary Figure 1 | Direct-sequence confusion matrices, 3 tools | **Supplementary Figure 1** | Replace image and caption, now 5 tools |

Replacement files, already generated:

| Target | File |
|---|---|
| Figure 2 | `analysis/results/plots/figure2_performance.pdf` |
| Figure 3 | `analysis/results/plots/figure3_ecoli_recovery.pdf` |
| Figure 4 | `analysis/results/plots/figure3_kouzine.pdf` |
| Supplementary Figure 1 | `analysis/results/plots/suppfig1_confusion_matrices.pdf` |

Figure 1 has no source file in the repository and is the one figure that needs hand editing in
whatever tool drew it (§6.1).

**In-text citations to update:** `Figure 4` → `Figure 5` in paragraphs 106 and 108; `Figure 4` →
`Figure 5` in the Figure 4 caption itself; `Figure 5` → `Figure 6` in its caption. New citations to
Figure 3 and Figure 4 are introduced by the new Results text.

---

## 4. Section-by-section edit list

Paragraph numbers are positions in the submitted `.docx`, counting all paragraphs including empty
ones, so they can be found deterministically. Each row names the reviewer comment that requires it.

### 4.1 Abstract — paragraph 10

Four clause-level edits inside one paragraph. The paragraph is otherwise untouched.

| Edit | Comment |
|---|---|
| "a state-of-the-art computational tool" → "a computational tool" | R3.12 |
| "The greedy overlap removal phase resolves overlapping candidates by retaining the longest and most compact motif within each overlapping region." → describes all three modes | R3.10 |
| "the thermodynamic stability scoring algorithm" → "an experimentally informed stability scoring model"; add that the score is a relative ranking, not a calibrated free energy | R2.4, R3.4 |
| Add one validation sentence before "HSeeker is also user-friendly": the balanced benchmark result and the single-stranded DNA comparison | R2.1, R3.8 |

### 4.2 Key Points — paragraphs 14, 15, 16, 17

| Para | Edit | Comment |
|---|---|---|
| 14 | "thermodynamic scoring" → "stability scoring" | R3.4 |
| 15 | "greedy overlap removal" → "three selectable overlap-resolution rules" | R3.10 |
| 16 | Replace in full. "On 79 experimentally validated sequences, HSeeker achieved 98.6% sensitivity and outperformed existing tools in benchmarking" is wrong on all three counts: the benchmark is 130, the sensitivity is 96.9%, and the superiority claim must carry the triplex qualification | R2.1, R2.7, R3.1, R3.6 |
| 17 | "~25 seconds" → "16 s", and add the measured whole-genome figure | R3.11 |

### 4.3 Introduction — paragraph 28 only

Two sentences inside one paragraph. The rest of the Introduction, including all of paragraphs
19–27 and 30, is untouched.

| Edit | Comment |
|---|---|
| nBMST: "only includes mirror repeats of at least 90% purines or pyrimidines, without tolerating any mismatch in mirror symmetry" — the 90% figure is not what the tool applies. Replace with a description of what it does: no mismatch tolerance, and a purine/pyrimidine count on the reported strand rather than a purity threshold | R3.6 |
| NeSSie: "optimized for detecting the mirror repeats characteristic of H-DNA" — it returns symmetric or purine-rich k-mers, which is a different object from an arm–spacer–arm motif. Qualify the sentence | R3.6 |

### 4.4 Methods

| Para | Action | Comment |
|---|---|---|
| 42 | Add two sentences on inward extension and deduplication, to match the corrected Algorithm 1 panel | consistency with §6.2 |
| 45 | Replace the composition-filters paragraph. Keep both filter descriptions as written; add the disclosure that the rules were chosen with knowledge of this benchmark, that neither filter changes any benchmark call, and that evidence for them comes from outside the benchmark | R2.6 |
| 47 | Keep the Park et al. attribution. Add that the combined score is a ranking without validated physical units and is not reported in kcal/mol | R2.4, R3.4 |
| 49 | Correct the spacer constraint. The text says the spacer is capped at 20 bp; the implementation requires the trimmed spacer not to exceed the selected window | correctness |
| 53–62 | Untouched except one added sentence after 62: negative totals are clamped to zero | correctness |
| **new, after 63** | **Insert: *Benchmark dataset*.** The 130-record benchmark, the inclusion criteria, the designed negatives, the separate homopolymer panel, and the assay conditions | R2.1, R3.1, R3.4, R3.5 |
| **new** | **Insert: *Detection settings*.** The default parameters in one place, and the statement that the package, CLI and web application share them | R2.11, R3.9 |
| 64 | Retitle *Sensitivity analysis* → *Validation design* and cut it back to design only: the call rule, the interval methods, and the grouping scheme. The results move to Results | R2.5, R2.7, R3.1, R3.2 |
| 65 | Correct the command line. It currently shows `-mismatch 0.15` and an injected `.fna` filename in the direct-sequence section | correctness, R3.9 |
| 66–67 | Keep. Add the two new comparators and the hit rule for each | R3.6 |
| 68 | **Delete.** The Youden threshold and the ROC are withdrawn | R3.3 |
| 69 | Retitle *Genomic context validation* → *Genomic context recovery* and replace. Keep the design; remove the biological framing; state that the overlap criterion is varied as a sensitivity analysis | R2.3, R3.7 |
| 70 | Replace. AUC-ROC is withdrawn; MCC is kept and the interval methods are named | R3.3, R3.1 |
| **new** | **Insert: *Comparison with existing tools*.** Versions, command lines, hit rules, and that the design was fixed before running | R3.6 |
| **new** | **Insert: *Genome-wide comparison with single-stranded DNA maps*.** Data, matching procedure, the two analyses and their controls | R2.3, R3.8 |
| **new** | **Insert: *Computational benchmarking*.** Hardware, three runs, external peak-memory measurement | R3.11 |

### 4.5 Results

| Para | Action | Comment |
|---|---|---|
| 85 | Replace the numbers only. Sentence structure is kept | R3.11 |
| 87 | Replace the comparator paragraph. Four tools at 16 cores, same node and day | R3.6, R3.11 |
| 89 | Replace. The linear extrapolation to "approximately 5 min" is withdrawn in favour of the measured whole-genome run | R3.11 |
| **new, after 89** | **Insert: *Algorithmic complexity*.** The complexity analysis moved verbatim out of Methods | R2.5 |
| 94 | Retitle *Comparison of HSeeker with other H-DNA detection tools using experimental data* → *Detection performance on the revised benchmark* | — |
| 95 | Replace in full. This paragraph carries 98.6% / 66.7%, the 79-sequence set, the E. coli recovery claim and MCC 0.73, all of which change | R2.1, R3.1, R3.3, R3.7 |
| 96 | Replace. "at least comparable to Triplex and Triplexator" is superseded by the four-tool comparison | R3.6 |
| **new** | **Insert: *Family- and study-aware validation*** | R3.2 |
| **new** | **Insert: *Component ablation and parameter robustness*** | R3.9 |
| **new** | **Insert: *Comparison with existing tools*** + **Table 1** | R2.7, R3.1, R3.6 |
| **new** | **Insert: *Genomic context recovery*** (Figure 3) | R2.3, R3.7 |
| **new** | **Insert: *Comparison with single-stranded DNA maps*** (Figure 4) | R2.3, R3.8 |
| 104 | Add one sentence: a single web interface sharing the command-line defaults | R2.11 |
| 106 | Update the preset values, which changed with the defaults; renumber Figure 4 → 5 | R3.9 |
| 108 | Renumber Figure 4 → 5. Remove the "Stable Sites (score ≥ 50)" tile from the results-page description if present in the final text | R3.3 |

### 4.6 Discussion

| Para | Action | Comment |
|---|---|---|
| 125 | Untouched | — |
| 127 | Two edits. "two existing H-DNA and mirror-repeat detection tools" → four; and the closing sentence "works highly effectively and aligns seamlessly with previously published thermodynamic parameters" → the moderated wording agreed in the letter | R3.6, R3.12 |
| 129 | Two edits. "HSeeker is the only H-DNA identification algorithm with a proven scoring system built-in" → moderated; "~25 seconds" → 16 s; and the trailing "previously computationally prohibitive" clause → the agreed rephrasing | R3.11, R3.12 |
| **new, after 129** | **Insert: *Limitations*.** The section Reviewer 2 asked for | R2.4, R3.4 |
| 130 | Untouched | — |

### 4.7 Availability statements

| Para | Action | Comment |
|---|---|---|
| 131 | Add the release tag | R2.8, R2.9 |
| 132 | Replace. "79 experimentally validated" → 130; add the Kouzine data source; remove "available from the corresponding authors upon reasonable request" | R2.2, R3.11 |

### 4.8 Global terminology pass

Reviewer 3's comment 12 asks for four terms to be kept distinct. This is a find-and-check pass, not
a rewrite: predictions become *putative H-DNA-forming sequences* or *predicted H-DNA motifs*;
*H-DNA-forming sequence* unqualified is reserved for sequences with experimental support; *mirror
repeat* is used for the sequence feature, never as a synonym for the structure. Each substitution
is tracked individually so the authors can reject any that read awkwardly.

---

## 5. Supplementary material

The submitted manuscript cites "Supplementary Table 1" for two different things: the benchmark
sequences (paragraphs 64, 95, 127, 132) and the stability-score pairing values (paragraph 51). This
has to be resolved before anything else is numbered.

| Number | Content | Status |
|---|---|---|
| Supplementary Table 1 | Benchmark records: sequence, label, provenance, curation decision and justification, family and study identifiers, assay conditions with supporting quotation | Replace — now 130 records with the new provenance and conditions columns |
| **Supplementary Table 2** | Stability-score pairing and stacking values | **New number** for existing content; update the citation in paragraph 51 |
| **Supplementary Table 3** | Comparator versions, command lines and hit rules | New |
| **Supplementary Table 4** | Per-record calls for all five tools | New |
| **Supplementary Table 5** | Homopolymer panel across tools and filter settings | New |
| Supplementary Figure 1 | Confusion matrices, now five tools on 130 records | Replace image and caption |

---

## 6. Figures and algorithm panels in detail

### 6.1 Figure 1 — the only figure needing hand editing

No source file exists in the repository, so these are instructions for whoever redraws it.

- **Panel B** contains a "Sensitivity Analysis / AUC" inset. That analysis is withdrawn
  (R3.3) — **remove the inset.** This also relieves the density Reviewer 3 noted in comment 12.
- **Panel B** labels the stability balance "H-DNA Stability" → "Stability score (heuristic
  ranking)", so the figure does not reintroduce the free-energy reading (R2.4, R3.4).
- **Panel B**'s "H-DNA-forming ✓ / non-forming ✗" decision tree implies a score-based decision →
  relabel as detection only (R3.3).
- **Panel A** omits two steps now in the algorithm: inward extension and deduplication. Add both
  before "Return candidate list", and thin the decision-diamond chain (R3.12).

### 6.2 Algorithm panels — paragraphs 71–79

The three panels are embedded images rendered from LaTeX. The submitted panels no longer match the
code: the purity prefilter is shown as mandatory, inward extension and deduplication are missing,
overlap removal shows only the greedy rule, and the scorer is shown returning nothing when no window
fits.

Corrected LaTeX exists in `main.tex` in the repository root and supersedes `algorithms.tex`, which
is the submitted version. `scoring.tex` is a byte-identical duplicate of the Phase 3 panel already
inside `main.tex` and should be deleted rather than maintained.

The revised set is four panels, which also answers the density point in R3.12:

| Panel | Content | Change |
|---|---|---|
| **Algorithm 1** | Pipeline: order of stages | **New** — fixes the stage order, which the submitted Methods described incorrectly |
| **Algorithm 2** | Candidate enumeration | Replaces submitted Algorithm 1; adds the `rmq` flag, inward extension, deduplication |
| **Algorithm 3** | Overlap removal | Replaces submitted Algorithm 2; adds the score-based branch |
| **Algorithm 4** | Stability scoring | Replaces submitted Algorithm 3; symbolic constants, fallback window, clamp, lower-case output |

Caption wording: "Thermodynamic Stability Scoring" → "Heuristic stability scoring" (R2.4, R3.4).

Each panel is re-rendered from `main.tex` to PNG at the same width as the current images and
swapped in place, so the surrounding layout does not move.

---

## 7. Coverage: every comment to a location in the paper

| Comment | Where it lands | Paper change? |
|---|---|---|
| Reviewer 1 | — | None requested |
| R2.1 benchmark size | Methods *Benchmark dataset*; Results 95; Abstract; Key Point 16 | Yes |
| R2.2 reproducibility | Data Availability 132 | Yes |
| R2.3 E. coli / Kouzine | Methods *Genomic context recovery* 69 + *Genome-wide comparison*; Results new sections; Figures 3, 4 | Yes |
| R2.4 limitations | Discussion *Limitations*; Abstract; Methods 47 | Yes |
| R2.5 move analyses | Section moves; Methods 64 cut to design | Yes |
| R2.6 filter circularity | Methods 45 | Yes |
| R2.7 confidence intervals | Results *Comparison with existing tools*; Table 1 | Yes |
| R2.8 `-march=native` | Code; release tag in Code Availability 131 | Tag only |
| R2.9 profiling counters | Code; release tag in Code Availability 131 | Tag only |
| R2.10 output casing | Code | No |
| R2.11 two `app.py` | Results 104 | Yes |
| R3.1 size and intervals | Table 1; Results 95 | Yes |
| R3.2 family-aware validation | Methods *Validation design*; Results new section | Yes |
| R3.3 threshold circularity | Methods 68 deleted, 70 replaced; Results 95; Figure 3 panels B, D removed; Figure 1 inset removed | Yes |
| R3.4 thermodynamic claim | Methods 47; Discussion *Limitations*; assay conditions in *Benchmark dataset* | Yes |
| R3.5 curation | Methods *Benchmark dataset*; Supplementary Table 1 | Yes |
| R3.6 comparators | Introduction 28; Methods *Comparison with existing tools*; Results; Key Point 16; Discussion 127 | Yes |
| R3.7 E. coli, overlap threshold | Methods 69; Results *Genomic context recovery*; Figure 3 | Yes |
| R3.8 genome-scale validation | Methods *Genome-wide comparison*; Results; Figure 4 | Yes |
| R3.9 ablation | Methods *Detection settings*; Results *Component ablation* | Yes |
| R3.10 overlap strategy | Abstract; Key Point 15; Methods 69 and overlap description; Results; Algorithm 3 | Yes |
| R3.11 benchmarking | Results 85, 87, 89; Figure 2; Data Availability 132 | Yes |
| R3.12 claims, terms, figures | Discussion 127, 129; global terminology pass; Figures 1 and algorithm panels | Yes |

Three comments (R2.8, R2.9, R2.10) are answered by code changes and need nothing in the paper beyond
the release tag. Every other comment has at least one tracked change attached to it.

---

## 8. Order of work

1. **Settle the figure numbering** (§3). Six in-text citations depend on it, and doing it first
   means they are edited once.
2. **Resolve the Supplementary Table 1 collision** (§5). Two Methods citations depend on it.
3. **Re-render the algorithm panels** from `main.tex` (§6.2) and swap the images.
4. **Apply the text edits** section by section in document order (§4), each as a tracked change.
5. **Swap Figures 2, 3, Supplementary Figure 1; insert Figure 4** with captions.
6. **Terminology pass** (§4.8) last, so it also covers the newly inserted text.
7. **Figure 1** is handed off separately, since it needs a drawing tool.

---

## 9. Decisions needed before step 1

1. **Marking mechanism** — Option A, B or C from §2. The plan assumes A.
2. **Figure 1** — who redraws it, and whether the four changes in §6.1 are accepted as described.
3. **Figure 2 panel count** — the submitted figure has four panels sweeping 1–16 workers for three
   tools. That cannot be reproduced: two of the four comparators have no parallel mode, so a sweep
   is meaningless for them. The replacement has three panels. If the four-panel design is preferred,
   sweeps for triplex and Triplexator can be run, and the other two shown as single points.
4. **Biographical Note, paragraph 18** — names only two authors. Not a reviewer point; flagged
   because it will be noticed at proof stage.
