#!/usr/bin/env python3
"""Apply the HSeeker revision to the submitted manuscript as tracked changes.

Every edit here is listed in REVISION_PLAN.md with the reviewer comment that
requires it. Paragraph indices are positions in the submitted document and do
not shift as the script inserts or deletes.

Re-runnable: it always starts from the submitted file and writes a new copy.
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from docx_track import Doc  # noqa: E402

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
SRC = HERE / "Revised HSeeker_ an algorithm for systematic H-DNA sequence identification .docx"
OUT = HERE / "HSeeker_revised_tracked.docx"
PLOTS = REPO / "analysis" / "results" / "plots"
SCRATCH = Path("/tmp/claude-903286/-scratch-10899-kimopro-HSeeker/"
              "2ace8cff-952f-44f6-83eb-34d39660be40/scratchpad")


def main() -> None:
    d = Doc(SRC)

    # =====================================================================
    # Abstract — paragraph 10
    # =====================================================================
    d.replace_text(10, "a state-of-the-art computational tool",
                   "a computational tool", "R3.12")
    d.replace_text(
        10,
        "The greedy overlap removal phase resolves overlapping candidates by "
        "retaining the longest and most compact motif within each overlapping region.",
        "Overlapping candidates are resolved by one of three user-selectable rules: "
        "retaining the longest and most compact motif within each overlapping region, "
        "retaining the highest-scoring motif, or reporting every candidate.", "R3.10")
    d.replace_text(
        10,
        "Finally, the thermodynamic stability scoring algorithm evaluates the "
        "candidate motifs using an experimentally informed scoring model that "
        "incorporates Hoogsteen G-G and A-A bonds, consecutive-pair stacking, and "
        "imposes penalties for mismatches and disrupted stacks.",
        "Finally, an experimentally informed stability scoring model evaluates the "
        "candidate motifs using Hoogsteen G-G and A-A bonds, consecutive-pair "
        "stacking, and penalties for mismatches and disrupted stacks; the resulting "
        "score ranks candidates relative to one another within this scheme and is not "
        "a calibrated free energy.", "R2.4/R3.4")
    d.replace_text(
        10,
        "HSeeker is also user-friendly,",
        "On a balanced benchmark of 130 experimentally characterised and designed "
        "sequences, HSeeker detected 63 of 65 H-DNA-forming sequences and rejected all "
        "65 non-forming sequences, and its genome-wide predictions were enriched for "
        "mirror repeats found single-stranded in cells by permanganate/S1 footprinting. "
        "HSeeker is also user-friendly,", "R2.1/R3.8")

    # =====================================================================
    # Key Points — paragraphs 14-17
    # =====================================================================
    d.replace_text(14, "thermodynamic scoring", "stability scoring", "R3.4")
    d.replace_text(
        15,
        "a center-outward mirror-repeat search, greedy overlap removal, and an "
        "experimentally informed stability-scoring model",
        "a center-outward mirror-repeat search, three selectable overlap-resolution "
        "rules, and an experimentally informed stability-scoring model", "R3.10")
    d.replace_text(
        16,
        "On 79 experimentally validated sequences, HSeeker achieved 98.6% sensitivity "
        "and outperformed existing tools in benchmarking.",
        "On a balanced benchmark of 130 experimentally characterised and designed "
        "sequences, HSeeker achieved higher accuracy than four existing tools, matching "
        "triplex on experimentally characterised sequences while rejecting the "
        "mirror-disrupted controls that triplex accepted.", "R2.1/R2.7/R3.1/R3.6")
    d.replace_text(
        17,
        "analyzing human chromosome 1 in ~25 seconds",
        "analysing human chromosome 1 in 16 s and the complete human genome in 3.4 min "
        "on 16 cores", "R3.11")

    # =====================================================================
    # Introduction — paragraph 28
    # =====================================================================
    d.replace_text(
        28,
        "but its implementation only includes mirror repeats of at least 90% purines "
        "or pyrimidines, without tolerating any mismatch in mirror symmetry",
        "but its mirror-repeat module tolerates no mismatch in mirror symmetry, and the "
        "subset it labels as triplex-forming is selected by a purine/pyrimidine count on "
        "the reported strand rather than by a purity threshold, so mirror repeats written "
        "on the pyrimidine-rich strand are not labelled as triplex-forming", "R3.6")
    d.replace_text(
        28,
        "NeSSie is a C/C++ library and tool optimized for detecting the mirror repeats "
        "characteristic of H-DNA, utilizing an optimal global alignment based on the "
        "Needleman-Wunsch algorithm",
        "NeSSie is a C/C++ library and tool for the detection of approximate sequence "
        "symmetries, including mirror symmetry and triplex-forming potential, utilizing "
        "an optimal global alignment based on the Needleman-Wunsch algorithm", "R3.6")
    d.replace_text(
        28,
        "and allowing gaps, mismatches, or asymmetric loops to be present.",
        "and allowing gaps, mismatches, or asymmetric loops to be present. Its search "
        "returns symmetric or purine-rich k-mers rather than arm-spacer-arm motifs, so a "
        "mirror repeat whose two arms flank a non-symmetric spacer is not returned as a "
        "single motif.", "R3.6")

    # =====================================================================
    # Methods
    # =====================================================================
    # candidate discovery: the two steps added since submission
    # Appended after the prefilter sentence, not spliced into it.
    d.replace_text(
        42,
        "eliminates center-spacer pairs whose right arm cannot satisfy the composition "
        "constraint.",
        "eliminates center-spacer pairs whose right arm cannot satisfy the composition "
        "constraint. After enumeration, matching base pairs at the inner ends of the "
        "spacer are absorbed into the arms where the purity criterion permits, so that a "
        "motif is reported with its maximal arm representation, and candidates arising "
        "from different centre and spacer combinations that normalise to the same arm "
        "and spacer partition are reported once.", "panel consistency")

    # composition filters: the disclosure Reviewer 2 asked for
    d.replace_text(
        45,
        "as these do not form stable inter-strand triplex interactions despite "
        "satisfying homopurine/homopyrimidine requirements.",
        "as these do not form stable inter-strand triplex interactions despite "
        "satisfying homopurine/homopyrimidine requirements. Both rules are fixed and are "
        "not fitted during evaluation. They were, however, chosen with knowledge of the "
        "H-DNA literature and of the composition of our benchmark, and the benchmark is "
        "therefore not an independent test of them: because the single-base homopolymers "
        "are reported separately, neither filter changes any benchmark call, so the "
        "reported metrics neither depend on the filters nor support them. Evidence "
        "bearing on the filters is drawn instead from published single-stranded DNA maps "
        "and is presented in the Results.", "R2.6")

    # stability score: no longer a free energy
    # Placed after the kcal/mol inputs are described, so it reads as a statement
    # about the combined score rather than contradicting the sentence that follows.
    d.replace_text(
        47,
        "Taken together, the HSeeker stability scoring was calculated as follows:",
        "The combination of these terms is a ranking score, not a free energy: the inputs "
        "are drawn from reported interaction energies, but the total has no validated "
        "physical units and we do not report it in kcal/mol. Taken together, the HSeeker "
        "stability scoring was calculated as follows:", "R2.4/R3.4")

    # spacer constraint: the text did not match the implementation
    d.replace_text(
        49,
        "However, H-DNA is thermodynamically unfavourable with very long spacer lengths; "
        "thus, we imposed a constraint on the total spacer length that we consider at a "
        "maximum of 20 bp.",
        "However, H-DNA is thermodynamically unfavourable with very long spacer lengths; "
        "thus, we require that the spacer formed by folding the trimmed arm tails inward "
        "does not exceed the length of the selected arm window.", "correctness")

    # The displayed constraint in this paragraph is an inline image still showing the
    # superseded m <= min(20, L); replace it with the constraint the code applies.
    d.swap_image(49, SCRATCH / "spacer_constraint.png", "spacer constraint", n=5)

    # The scoring-weights table and the benchmark both claimed "Supplementary
    # Table 1" in the submitted manuscript. The weights become Table 2.
    d.replace_text(51, "we consider the values in Supplementary Table 1",
                   "we consider the values in Supplementary Table 2", "numbering")

    # total score is clamped
    d.replace_text(
        62,
        "The procedure yields a stacking score, a pairing score, their sum as the total "
        "score,",
        "The procedure yields a stacking score, a pairing score, their sum as the total "
        "score, clamped at zero so that the reported score is non-negative,", "correctness")

    # ---- Methods: new subsections ---------------------------------------
    d.insert_paragraph_after(63, BENCHMARK_DATASET, style_from=45, label="R2.1/R3.1/R3.4/R3.5")
    d.insert_paragraph_after(63, DETECTION_SETTINGS, style_from=45, label="R2.11/R3.9")

    # Sensitivity analysis -> Validation design
    d.replace_text(
        64,
        "Sensitivity analysisTo assess the ability of HSeeker to detect H-DNA-forming "
        "sequences, we conducted a sensitivity analysis using a benchmark set of 79 "
        "experimentally validated sequences: 70 known to form H-DNA and 9 negative "
        "controls. Positive controls were selected based on strong experimental evidence "
        "of H-DNA formation, including nucleotide-level structural fine-mapping via "
        "chemical modification assays in a plasmid context (Supplementary Table 1).",
        "Validation designA sequence is called H-DNA-forming when at least one candidate "
        "survives detection and filtering; no score threshold is applied. Sensitivity, "
        "specificity, positive and negative predictive values are reported with Wilson "
        "intervals and the Matthews correlation coefficient with percentile bootstrap "
        "intervals. Because many records are closely related variants from the same "
        "experimental series, intervals are computed both by resampling records and by "
        "resampling whole sequence families, and both are reported. Records were grouped "
        "by sequence family and by study from the curated provenance, and designed "
        "negatives inherit the group of the record from which they were derived. "
        "Parameter selection, where performed, was confined to the training portion of "
        "each split, under repeated stratified five-fold cross-validation, "
        "leave-one-family-out and leave-one-study-out schemes.", "R2.5/R2.7/R3.1/R3.2")

    # The three panels become four, so every in-text reference moves up one.
    d.replace_text(42, "as described in Algorithm 1",
                   "as described in Algorithm 2", "panel renumber")
    d.replace_text(42, "The candidate enumeration step is described in Algorithm 2",
                   "The overlap-removal step is described in Algorithm 3", "panel renumber")
    d.replace_text(62, "The scoring procedure is described in Algorithm 3.",
                   "The scoring procedure is described in Algorithm 4.", "panel renumber")

    # the quoted command line was wrong in two ways
    d.replace_text(
        65,
        "We first benchmarked HSeeker against two existing tools, Triplex and "
        "Triplexator, running all three independently on the 79 benchmark sequences.",
        "We benchmarked HSeeker against four existing tools, triplex, Triplexator, "
        "NeSSie and nBMST, running all five independently on the 130 benchmark "
        "sequences.", "R3.6")
    d.replace_text(
        65,
        "hseeker -seq hdna_benchmark_injected.fna -minrep 8 -maxspacer 10 -mismatch 0.15 "
        "-purity 0.9 -maxrep 1000 -out hdna_benchmark",
        "hseeker -seq benchmark.fna -out hdna_benchmark", "correctness/R3.9")

    # the Youden threshold is withdrawn
    d.delete_paragraph(68, "R3.3")

    # genomic context validation -> recovery
    d.replace_text(
        69,
        "Genomic context validation. We next tested whether each tool could recover the "
        "benchmark sequences when embedded in a realistic genomic background.",
        "Genomic context recovery. We next tested whether embedding a sequence in a "
        "genome changes the call made on it.", "R2.3/R3.7")
    d.replace_text(
        69,
        "randomly inserted each of the 79 sequences into it",
        "randomly inserted each of the 130 sequences into it", "R2.1")
    # Appended after the existing justification, so that clause keeps its referent.
    d.replace_text(
        69,
        "a threshold chosen to accommodate cases where tools extend hits into flanking "
        "sequence or detect only the triplex core within a longer insertion.",
        "a threshold chosen to accommodate cases where tools extend hits into flanking "
        "sequence or detect only the triplex core within a longer insertion. The "
        "criterion is varied from 0.50 to 0.95 as a sensitivity analysis, and the "
        "unmodified genome was scanned under the same settings to establish the "
        "background locus rate. This experiment measures recovery and the stability of a "
        "call across sequence context; it does not test whether H-DNA forms, and we draw "
        "no biological conclusion from it.", "R3.7")

    # AUC-ROC is withdrawn; the interval methods are named
    d.replace_text(
        70,
        "For both the direct and genomic-context analyses, performance was evaluated "
        "using AUC-ROC and the Matthews Correlation Coefficient (MCC), the latter chosen "
        "to account for the substantial class imbalance in the dataset (~7.8:1 ratio of "
        "H-DNA-forming to non-forming sequences).",
        "For both the direct and genomic-context analyses, performance was evaluated "
        "using the Matthews correlation coefficient, reported with percentile bootstrap "
        "intervals computed over records and over whole sequence families, together with "
        "sensitivity, specificity and predictive values reported with Wilson intervals.",
        "R3.1/R3.3")

    d.insert_paragraph_after(70, TOOL_COMPARISON_METHODS, style_from=45, label="R3.6")
    d.insert_paragraph_after(70, GENOMEWIDE_METHODS, style_from=45, label="R2.3/R3.8")
    d.insert_paragraph_after(70, BENCHMARKING_METHODS, style_from=45, label="R3.11")

    # =====================================================================
    # Results
    # =====================================================================
    d.replace_text(
        85,
        "All runs produced an identical set of 96,729 non-overlapping H-DNA loci, "
        "confirming deterministic output across all parallelisation levels. At a single "
        "worker, the full pipeline was completed in 187.7 s (approximately 3.1 min). "
        "Increasing the thread count to 2, 4, 8, and 16 yielded wall-clock times of "
        "105.3 s, 59.5 s, 35.9 s, and 25.2 s, corresponding to speedups of 1.78-fold, "
        "3.15-fold, 5.22-fold, and 7.44-fold (Figure 2). Parallel efficiency remained "
        "above 79% through 4 workers and was 65% at 8 workers and 47% at 16 workers,",
        "All runs produced an identical set of 45,927 non-overlapping putative H-DNA "
        "loci, "
        "confirming deterministic output across all parallelisation levels. At a single "
        "worker, the full pipeline was completed in 173.8 s (approximately 2.9 min). "
        "Increasing the thread count to 2, 4, 8, and 16 yielded wall-clock times of "
        "89.6 s, 47.4 s, 26.7 s, and 16.3 s, corresponding to speedups of 1.94-fold, "
        "3.67-fold, 6.51-fold, and 10.65-fold (Figure 2). Parallel efficiency remained "
        "above 90% through 4 workers and was 81% at 8 workers and 67% at 16 workers,",
        "R3.11")
    d.replace_text(
        85,
        "with all other parameters held at their defaults (minimum arm length 8 bp, "
        "maximum spacer 10 bp, purity threshold 0.90, RMQ purity prefilter enabled, "
        "scoring enabled).",
        "with all other parameters held at their defaults (minimum arm length 8 bp, "
        "maximum spacer 10 bp, purity threshold 0.90, mismatch tolerance 0.10, "
        "AT-content threshold 0.80 with homopolymer exclusion, RMQ purity prefilter "
        "enabled, scoring enabled).", "R3.9")

    d.replace_paragraph(87, COMPARATOR_TIMINGS, "R3.6/R3.11")
    d.replace_paragraph(89, WHOLE_GENOME, "R3.11")

    # Reviewer 2, comment 5: the complexity analysis belongs in the Results.
    # Moved whole, so its 34 inline equations and citation links travel with it.
    d.move_paragraphs([76, 77], dest=89, label="R2.5 complexity -> Results")

    d.replace_text(
        94,
        "Comparison of HSeeker with other H-DNA detection tools using experimental data",
        "Detection performance on the revised benchmark", "section title")
    d.replace_paragraph(95, DETECTION_RESULTS, "R2.1/R3.1/R3.3/R3.7")
    d.replace_paragraph(96, DETECTION_ERRORS, "R3.6")

    # Reversed, because each call inserts immediately after paragraph 96.
    for text, label in ((KOUZINE_RESULTS, "R2.3/R3.8"), (ECOLI_RESULTS, "R2.3/R3.7"),
                        (TOOL_RESULTS, "R2.7/R3.1/R3.6"), (ABLATION, "R3.9"),
                        (FAMILY_AWARE, "R3.2")):
        d.insert_paragraph_after(96, text, style_from=95, label=label)

    # web interface: one interface, updated presets
    d.replace_text(
        104,
        "we developed a web interface that enables users to identify putative "
        "H-DNA-forming sequences without requiring command-line experience.",
        "we developed a single web interface that enables users to identify putative "
        "H-DNA-forming sequences without requiring command-line experience. It uses the "
        "same scanning library and the same default settings as the command-line tool, "
        "so that the two return identical results for identical input.", "R2.11")

    # =====================================================================
    # Discussion
    # =====================================================================
    d.replace_text(
        127,
        "We also performed thorough benchmarking against two existing H-DNA and "
        "mirror-repeat detection tools, Triplex and Triplexator",
        "We also performed thorough benchmarking against four existing H-DNA and "
        "mirror-repeat detection tools, triplex, Triplexator, NeSSie and nBMST", "R3.6")
    d.replace_text(
        127,
        "The successful discrimination of experimentally validated sequences confirms "
        "that this biophysically informed scoring system works highly effectively and "
        "aligns seamlessly with previously published thermodynamic parameters.",
        "The stability score ranks candidates in a manner consistent with the "
        "biochemical and genetic evidence available for the sequences in our benchmark; "
        "it is an experimentally informed heuristic and has not been calibrated against "
        "measured thermodynamic quantities.", "R3.12")
    d.replace_text(
        129,
        "HSeeker is the only H-DNA identification algorithm with a proven scoring system "
        "built-in, and can achieve high-throughput processing speeds.",
        "Among currently available tools for intramolecular H-DNA, HSeeker is the only "
        "one that reports a per-motif stability ranking alongside the detected "
        "coordinates, and it achieves high-throughput processing speeds.", "R3.12")
    d.replace_text(
        129,
        "it was able to screen the entire human chromosome 1 in approximately 25 seconds "
        "using a standard multi-threaded environment.",
        "it was able to screen the entire human chromosome 1 in approximately 16 seconds, "
        "and the complete human genome in 3.4 minutes, using a standard multi-threaded "
        "environment.", "R3.11")
    d.insert_paragraph_after(129, LIMITATIONS, style_from=127, label="R2.4/R3.4")

    # =====================================================================
    # Availability
    # =====================================================================
    d.replace_paragraph(132, DATA_AVAILABILITY, "R2.2/R3.11")

    # =====================================================================
    # Figures: images, captions, renumbering
    # =====================================================================
    # Figure 1 is left exactly as submitted, on the authors' instruction. The
    # redrawn version is in analysis/results/plots/figure1_workflow.png if it is
    # wanted later; see the note in REVISION_PLAN.md about the AUC inset.

    # Algorithm panels: the submitted three are superseded by four
    d.swap_image(71, PLOTS / "algorithm2.png", "Algorithm 2 (was 1)")
    d.swap_image(73, PLOTS / "algorithm3.png", "Algorithm 3 (was 2)")
    d.swap_image(78, PLOTS / "algorithm4.png", "Algorithm 4 (was 3)")
    d.insert_paragraph_after(70, ALG1_CAPTION, style_from=72,
                             label="new Algorithm 1 caption")
    d.insert_image_after(70, PLOTS / "algorithm1.png", template=71,
                         label="new Algorithm 1 (pipeline)")
    d.replace_paragraph(72, ALG2_CAPTION, "Algorithm 2 caption")
    d.replace_paragraph(74, ALG3_CAPTION, "Algorithm 3 caption")
    d.replace_paragraph(79, ALG4_CAPTION, "Algorithm 4 caption")

    d.swap_image(91, PLOTS / "figure2_performance.png", "Figure 2")
    d.replace_paragraph(92, FIG2_CAPTION, "Figure 2 caption")

    d.swap_image(97, PLOTS / "figure3_ecoli_recovery.png", "Figure 3")
    d.replace_paragraph(98, FIG3_CAPTION, "Figure 3 caption")
    d.insert_paragraph_after(98, FIG4_CAPTION, style_from=98, label="new Figure 4")
    d.insert_image_after(98, PLOTS / "figure3_kouzine.png", template=97,
                         label="new Figure 4")

    # web-interface figures move to 5 and 6
    d.replace_text(106, "(Figure 4)", "(Figure 5)", "renumber")
    d.replace_text(108, "(Figure 4)", "(Figure 5)", "renumber")
    d.replace_text(119, "Figure 4: HSeeker web interface.",
                   "Figure 5: HSeeker web interface.", "renumber")
    d.replace_text(121, "Figure 5: Interactive visualization",
                   "Figure 6: Interactive visualization", "renumber")

    # New reference for the Mg2+-dependent poly-(dC) triplex cited in Limitations.
    # Every other work cited by the new text is already in the list: the triplex
    # algorithm (21), nBMST (23), NeSSie (24), Triplexator (25) and Kouzine (4).
    d.insert_paragraph_after(175, REFERENCE_28, style_from=175, label="new reference")

    d.swap_image(144, PLOTS / "suppfig1_confusion_matrices.png", "Supp Figure 1")
    d.replace_paragraph(145, SUPPFIG1_CAPTION, "Supp Figure 1 caption")

    # Table 1
    d.insert_table_after(96, TABLE1_HEADER, TABLE1_ROWS, "Table 1")
    d.insert_paragraph_after(96, TABLE1_CAPTION, style_from=92, label="Table 1 caption")

    # Terminology: predictions are putative unless experimentally demonstrated
    for para, old, new in TERMINOLOGY:
        try:
            d.replace_text(para, old, new, "R3.12 terminology")
        except KeyError:
            print(f"   terminology: p{para} {old[:44]!r} not found, skipped")

    d.save(OUT)
    print()
    for line in d.log:
        print("  ", line)


# =========================================================================
# Inserted blocks
# =========================================================================
BENCHMARK_DATASET = (
    "Benchmark dataset. The benchmark contains 130 sequences, 65 reported to form H-DNA "
    "and 65 reported or designed not to form it. Seventy are experimental records, 65 "
    "forming and five non-forming. Each was traced to its primary publication; records "
    "whose annotation derived only from a secondary compilation, and records with no "
    "traceable primary source, were removed. Fifty-seven of the retained records are "
    "construct-level, meaning the published experiment tested the sequence in the "
    "benchmark, and are used for the primary performance figures. The remaining thirteen "
    "are class-level, meaning the published experiment tested a closely related construct "
    "of the same repeat family, and are reported separately. Two records were excluded on "
    "validity grounds: an inosine-substituted oligonucleotide that cannot be represented "
    "in A/C/G/T, and a record describing an intermolecular triplex-forming "
    "oligonucleotide rather than an intramolecular structure. Because only five "
    "experimental non-forming constructs survived curation, we added 60 designed "
    "negatives, generated by a documented procedure under a fixed random seed: 22 "
    "positives whose mirror symmetry was destroyed by purine-preserving substitutions, 11 "
    "dinucleotide-preserving shuffles of positives, 12 random sequences matched to "
    "positives in length and GC content, four G-quadruplex-forming sequences, three "
    "Z-DNA-forming (CG)n repeats, five random mixed-sequence controls, and three perfect "
    "mirror repeats of balanced base composition. Apart from the last group, every "
    "designed negative was verified to contain no homopurine mirror repeat satisfying the "
    "detection criteria; the three balanced perfect mirrors are included precisely because "
    "they do contain such a repeat but are not purine-rich, and they therefore test "
    "whether mirror symmetry alone triggers a call. Single-base homopolymers are reported "
    "as a separate panel rather than within the benchmark, because the composition filters "
    "determine their calls by construction. For each experimental record we recorded the "
    "assay used, the DNA context, and the pH, temperature, ionic conditions and "
    "supercoiling state under which the structure was assessed, together with the passage "
    "in the source supporting those values; values that the source does not report are "
    "marked as such. Sixty-four of the 70 experimental records have a reported pH. Because "
    "the labels are outcomes under heterogeneous conditions, the pooled performance "
    "figures should be read as discrimination of reported sequence labels across assay "
    "conditions rather than as propensities under a single common condition "
    "(Supplementary Table 1)."
)

DETECTION_SETTINGS = (
    "Detection settings. Unless stated otherwise, all analyses used the default settings: "
    "minimum arm length 8 bp, maximum arm length 1,000 bp, maximum spacer 10 bp, purity "
    "threshold 0.90, mismatch tolerance 0.10, AT-content threshold 0.80 with homopolymer "
    "exclusion enabled, and overlapping candidates resolved by the longest-arm rule. These "
    "are the defaults of the Python package, the command-line tool and the web application "
    "alike, so that the three interfaces return identical results for identical input, and "
    "parameters are given explicitly wherever an analysis departs from these values."
)

TOOL_COMPARISON_METHODS = (
    "Comparison with existing tools. HSeeker was compared with triplex (v1.50.0), "
    "Triplexator (v1.3.2), NeSSie and nBMST on the same 130 records. Each tool was run at "
    "its documented defaults; with its geometry matched to HSeeker's where it exposes the "
    "corresponding parameters; and with a small parameter grid selected inside training "
    "folds only. A record was scored as positive when a tool reported at least one "
    "qualifying hit within it, the qualifying hit type being a triplex of any type for "
    "triplex, a duplex target site for Triplexator, a triplex-potential motif for NeSSie, "
    "and a mirror repeat classified as triplex-forming for nBMST; for NeSSie and nBMST the "
    "alternative hit types are reported alongside. Tool versions, command lines and the hit "
    "rule for each tool are given in Supplementary Table 3, and the per-record calls in "
    "Supplementary Table 4. The settings and the statistical comparison were fixed before "
    "the tools were run."
)

GENOMEWIDE_METHODS = (
    "Genome-wide comparison with single-stranded DNA maps. We compared HSeeker's "
    "predictions with the permanganate/S1 footprinting data of Kouzine and colleagues, "
    "using their published non-B DNA motif calls and single-stranded DNA signal tracks for "
    "activated and resting mouse B cells and for human Raji cells, on the genome assemblies "
    "used in that work. Both genomes were scanned with HSeeker at its default settings. "
    "Recall was computed as the fraction of the motifs they called formed that are "
    "overlapped by an HSeeker locus. Selectivity was computed among their own predicted "
    "H-DNA motifs, as the odds that a motif is called formed when HSeeker reports it "
    "against when HSeeker does not, with the mappability of each motif to 36-nt reads used "
    "as a sensitivity analysis. Enrichment of the measured signal was assessed against "
    "composition-matched controls. For each HSeeker locus, control windows were drawn on "
    "the same chromosome with the same length, the same A/T count and purine/pyrimidine "
    "asymmetry within one count, containing no ambiguous base and overlapping no HSeeker "
    "locus or predicted H-DNA motif. Two analyses were performed. The first used all loci "
    "with controls drawn genome-wide. The second restricted both loci and controls to "
    "sequence that 36-nt reads can map to, excluded intervals within 100 bp of a predicted "
    "G-quadruplex, Z-DNA, cruciform or stress-induced melting motif and within 2 kb of an "
    "annotated transcription start site, and drew ten controls per locus within 50 kb. "
    "Enrichment is reported as the ratio of mean signal with a bootstrap interval over "
    "loci, together with the probability that a locus exceeds one of its own controls. Two "
    "checks accompany each comparison: one control per locus was compared against that "
    "locus's remaining controls, and the motifs called formed by the original authors were "
    "passed through the identical procedure."
)

BENCHMARKING_METHODS = (
    "Computational benchmarking. Runtime and peak resident set size were measured on a "
    "single node of the TACC Lonestar6 cluster with 16 virtual cores of an AMD EPYC 7763 "
    "processor, 32 GB of RAM and Rocky Linux 8.10, with no other workload on the node. Each "
    "configuration was run three times and wall-clock time and peak memory were recorded "
    "externally. Inputs were human chromosome 1 and the hg38 primary assembly, both "
    "obtained from the UCSC Genome Browser. Comparators were timed on the same node under "
    "the same protocol."
)

COMPARATOR_TIMINGS = (
    "For comparison we timed four existing tools on the same chromosome, the same node and "
    "the same day. Triplexator, which provides native multi-threading, completed in 45.0 s "
    "at 16 threads and reported 275,476 duplex target sites. The triplex package is "
    "single-threaded; split into 1 Mbp windows and run as 16 parallel processes it "
    "completed in 64.7 s and reported 81,906 hits. nBMST and NeSSie provide no parallel "
    "mode and completed in 96.5 s and 295.7 s on one core, reporting 147,222 mirror "
    "repeats, of which 33,519 were classified as triplex-forming, and 720,811 motif "
    "occurrences respectively. Peak memory was 1,714 MB for Triplexator, 660 MB for each "
    "triplex worker process, 866 MB for nBMST and 1,056 MB for NeSSie, against 1,334 MB for "
    "HSeeker. The hit counts reported by triplex and Triplexator reproduce our earlier "
    "measurements exactly. The comparison should be read with the differences in target "
    "class in mind: Triplexator enumerates duplex target sites for a separately supplied "
    "third strand and NeSSie enumerates symmetric k-mers, neither of which is the same "
    "object as an intramolecular mirror repeat, and the loci counts are correspondingly not "
    "comparable."
)

WHOLE_GENOME = (
    "Rather than extrapolating from chromosome 1, we ran the complete hg38 primary "
    "assembly. At 16 workers it completed in 204.5 s, that is 3.4 min, with a peak memory "
    "of 2,296 MB, reporting 559,601 non-overlapping H-DNA loci; three independent runs "
    "produced identical output. Score-based overlap selection, which defers overlap "
    "resolution until after scoring, took 23.9 s on chromosome 1 and reported 47,217 loci."
)

DETECTION_RESULTS = (
    "HSeeker detected 63 of the 65 H-DNA-forming sequences and rejected all 65 non-forming "
    "sequences, corresponding to a sensitivity of 0.969 (95% CI 0.895-0.992), a specificity "
    "of 1.000 (0.944-1.000) and a Matthews correlation coefficient of 0.970 (0.925-1.000) "
    "(Table 1). The two sequences that were not detected, the c-MYC nuclease-sensitive "
    "element construct (ACCCTCCCC)4 and the gamma-globin -228/-189 element, are tandem "
    "direct repeats rather than mirror repeats and therefore fall outside the structural "
    "model that HSeeker implements. Restricting the evaluation to the 70 experimental "
    "records gives a sensitivity of 0.969 (0.895-0.992) and a specificity of 1.000 "
    "(0.566-1.000); the width of the latter interval reflects the five experimental "
    "non-forming constructs that survived curation, and conclusions about specificity rest "
    "largely on the designed negatives. We also asked whether the stability score could "
    "improve classification beyond detection. Within each training fold we selected the "
    "score threshold maximising Youden's J statistic and applied it only to the held-out "
    "records. Out of fold this produced a coefficient of 0.949 (0.937-0.961) under repeated "
    "five-fold cross-validation and 0.955 (0.898-1.000) under leave-one-family-out and "
    "leave-one-study-out validation, below the 0.970 obtained with no threshold. The "
    "threshold removed no false positive, because no non-forming record retains a candidate "
    "once the purity requirement is applied, and it discarded forming records whose scores "
    "fell below a threshold learned without them. For the same reason the area under the "
    "ROC curve on this benchmark, 0.985, restates the detection result rather than adding "
    "to it. We therefore report no score threshold: a sequence is called H-DNA-forming when "
    "it contains a candidate satisfying the detection criteria, and the stability score is "
    "used to rank candidates."
)

DETECTION_ERRORS = (
    "The errors are informative. Mirror-repeat detection alone, before the purity "
    "requirement is applied, calls seven of the 65 non-forming sequences: the "
    "experimentally non-forming c-MYC control B-33, one G-quadruplex-forming sequence, two "
    "Z-DNA-forming (CG)n repeats and all three perfect mirror repeats of balanced "
    "composition. The purity requirement removes all seven without losing a forming "
    "sequence, which is the clearest demonstration in the benchmark that mirror symmetry "
    "alone is not sufficient and that the composition criterion is doing the discriminating "
    "work. No designed negative of any class is called at the default settings."
)

FAMILY_AWARE = (
    "Family- and study-aware validation. Because many benchmark records are closely related "
    "variants from the same experimental series, we grouped records by sequence family and "
    "by study and repeated the evaluation holding out one whole group at a time. The 130 "
    "records span 52 sequence families and 27 studies, the largest family, the pXY32 "
    "substitution series and its derivatives, accounting for 37 of them. Under "
    "leave-one-family-out and leave-one-study-out validation the fixed default "
    "configuration reached a coefficient of 0.970, against 0.924 when the parameters were "
    "tuned within each training fold. Collapsing the benchmark to one representative per "
    "family gave 0.862 over 52 records, and to one representative per family and class "
    "0.945 over 74 records. Because the parameters are fixed, a record's call depends only "
    "on its own sequence, and group-aware validation does not shift the point estimate; it "
    "widens the confidence interval, and both the record-level and the family-level "
    "intervals are reported throughout."
)

ABLATION = (
    "Component ablation and parameter robustness. To establish which components of the "
    "pipeline account for its behaviour, we evaluated the stages cumulatively, calling a "
    "record forming when at least one candidate survives the stage in question. "
    "Mirror-repeat detection alone detected 63 of 65 forming records and seven of 65 "
    "non-forming records (sensitivity 0.969, specificity 0.892, coefficient 0.864). Adding "
    "the purity requirement removed all seven without losing a forming record (specificity "
    "1.000, coefficient 0.970). Adding the AT-content filter, pairing-only scoring, pairing "
    "plus stacking, and the complete model with boundary optimisation and homopolymer "
    "filtering changed no record-level call. Retaining all candidates, longest-arm selection "
    "and score-based selection gave identical record-level calls. Varying one parameter at a "
    "time around the defaults: purity 0.80-0.90 gave 0.970, while 0.95 and 1.00 lost nine "
    "forming records (0.843); mismatch tolerances of 0 and 0.05 detected 49 of 65 forming "
    "records (0.778), 0.10 detected 63 (0.970), and 0.15 and 0.20 detected 64 but added a "
    "false positive (0.969); a minimum arm length of 6 added ten false positives (0.839), 7 "
    "and 8 gave 0.985 and 0.970, 9 and 10 lost eight forming records (0.884) and 12 lost "
    "sixteen (0.778); a maximum spacer of 3 detected 41 forming records (0.679), 5 detected "
    "55 (0.856), and 8 to 15 detected 63 (0.970). Across all 576 combinations with selection "
    "confined to training folds, the fixed default configuration performed as well as or "
    "better than per-fold tuning under nested cross-validation, leave-one-family-out and "
    "leave-one-study-out validation (0.970 against 0.939, 0.924 and 0.924), and tuning did "
    "not beat the fixed configuration in any of 200 random splits. A minimum arm length of 8 "
    "was selected in 99% of nested folds and in every family and study fold, and a mismatch "
    "tolerance of 0.10 in 48 of 52 family folds. Performance therefore does depend on these "
    "choices, but the reported results do not rest on a single hand-selected "
    "parameterisation. The scorer does not expose an independent mismatch-penalty switch, so "
    "we do not claim to have ablated that penalty separately."
)

TOOL_RESULTS = (
    "Comparison with existing tools. On the balanced benchmark HSeeker reached a coefficient "
    "of 0.970 (family-level 95% CI 0.909-1.000), followed by triplex at 0.819 (0.744-0.936), "
    "NeSSie at 0.583 (0.421-0.780), nBMST at 0.449 (0.346-0.560) and Triplexator at 0.024 "
    "(-0.137-0.249) (Table 1, Supplementary Figure 1). Using paired differences with "
    "confidence intervals, the difference was +0.15 (0.03-0.22) against triplex, with 10 "
    "records called correctly by HSeeker and none by triplex alone (exact McNemar p = 0.002), "
    "+0.39 (0.18-0.55) against NeSSie, +0.52 (0.41-0.62) against nBMST and +0.95 (0.71-1.11) "
    "against Triplexator. The advantage over triplex comes from the designed negatives: on "
    "the 70 experimental records alone the difference was +0.06 with an interval reaching "
    "zero (0.00-0.23), and on that set the two tools are comparable. These differences follow "
    "from what each tool searches for. Triplexator enumerates duplex target sites for a "
    "separately supplied third strand and does not require a mirror repeat; at the settings "
    "used in our original manuscript it detects 64 of 65 forming records but also calls 42 of "
    "the 65 negatives, including all 22 mirror-disrupted mutants. NeSSie reports k-mers that "
    "are themselves mirror-symmetric rather than arm-spacer-arm motifs; it made no false call "
    "but detected none of the 16 pXY32 constructs. nBMST's mirror-repeat module matches "
    "HSeeker's geometry, but its triplex classification accepts only mirror repeats with more "
    "purines than pyrimidines on the reported strand, so it detected none of the 19 forming "
    "records written on the pyrimidine-rich strand, whereas its unclassified mirror-repeat "
    "output detected 15 of them. triplex, the only other tool designed for intramolecular "
    "triplexes, detected nearly the same forming records as HSeeker but called nine designed "
    "negatives forming, seven of them mirror-disrupted mutants. In the homopolymer panel, "
    "triplex and NeSSie call four of the five A/T homopolymers forming, and both call "
    "poly-(dC)30 forming; Triplexator calls none; nBMST calls poly-A30 only; and HSeeker "
    "calls none, which is correct for the five A/T homopolymers and incorrect for "
    "poly-(dC)30 (Supplementary Table 5)."
)

ECOLI_RESULTS = (
    "Genomic context recovery. Sequences are rarely scanned in isolation, so we asked "
    "whether embedding a sequence in a genome changes the call made on it. Inserting the 130 "
    "benchmark sequences into the E. coli K-12 genome left HSeeker's calls unchanged for "
    "every record: the same 63 forming sequences were detected, no non-forming sequence was "
    "called, and the resulting coefficient of 0.970 is identical to the value obtained on the "
    "sequences alone (Figure 3A, B). The matching hit covered the whole inserted sequence for "
    "the majority of records and at least half of it in every case (Figure 3D). triplex, run "
    "on the same modified genome, returned a different call for 44 of the 130 records, in "
    "each case a detection present when the sequence was scanned alone and absent once it was "
    "embedded, and its coefficient fell from 0.819 to 0.512. Varying the overlap required to "
    "match a hit to an insertion from 0.50 to 0.95 changed neither result (Figure 3C). The "
    "unmodified genome yielded 83 loci across its 4.64 Mb, or 17.9 per Mb, and the hits in "
    "the modified genome are accounted for by those background loci together with the "
    "recovered insertions. This experiment measures the stability of a call across sequence "
    "context, not whether H-DNA forms. We draw no biological conclusion from it; the "
    "comparison with measured single-stranded DNA follows."
)

KOUZINE_RESULTS = (
    "Comparison with single-stranded DNA maps. We compared HSeeker's predictions with the "
    "permanganate/S1 footprinting data of Kouzine and colleagues. HSeeker loci overlap 79.0% "
    "of the H-DNA structures found to be formed in mouse (13,517 of 17,109) and 82.3% in "
    "human Raji cells (6,879 of 8,356). Among their predicted H-DNA motifs, those HSeeker "
    "calls are formed more often than those it does not: 4.5% against 2.8% in mouse (odds "
    "ratio 1.65, 95% CI 1.59-1.71) and 5.7% against 1.7% in human (3.39, 3.20-3.59); "
    "restricting to motifs that the 36-nt reads can map to raises these to 1.68-2.07 and "
    "4.16-5.70 (Figure 4A, B). Across all HSeeker loci with genome-wide composition-matched "
    "controls, loci in activated mouse B cells carried 14% more single-stranded signal than "
    "their controls (ratio of means 1.137, 1.127-1.147; 825,749 loci); in resting B cells the "
    "ratio was 0.822 (0.818-0.827) and in human Raji cells 0.977 (0.969-0.984). Restricting "
    "loci and controls to mappable sequence more than 100 bp from other predicted non-B "
    "motifs and more than 2 kb from annotated transcription start sites, activated mouse B "
    "cells carried 25% more signal than their controls (1.254, 1.179-1.332; a locus exceeded "
    "one of its own controls with probability 0.514, p = 2 x 10-17), while human Raji cells "
    "showed no detectable excess (1.023, 0.988-1.061; p = 0.08). Two built-in checks behaved "
    "as required: comparing one control against the remaining controls of the same locus gave "
    "no difference (0.500-0.502), and the motifs the authors called formed were enriched "
    "7.4-fold in mouse and 2.1-fold in human under the identical procedure (Figure 4C). The "
    "excess increased with mirror-repeat arm length, from no detectable effect at 8-9 nt "
    "(1.12, 0.93-1.35) to 1.58 (1.34-1.83) at 12-14 nt and 1.6-1.8 at longer arms (Spearman "
    "p = 3 x 10-21), and was absent in resting B cells (0.858, 0.831-0.885), consistent with "
    "the activation- and supercoiling-dependent formation the authors report (Figure 4D). "
    "Most formed H-DNA that HSeeker does not report are short A/T-rich or single-base motifs "
    "excluded by the composition filters; in these data that class is formed at about half "
    "the rate of other predicted H-DNA in mouse and a quarter in human, which supports "
    "excluding it by default, though it is 17-19% of the formed H-DNA. These results "
    "associate HSeeker loci in activated mouse B cells with single-stranded DNA beyond what "
    "composition, mappability, local context, the other predicted non-B structures and "
    "promoter transcription account for, and the association strengthens with arm length, as "
    "H-DNA formation should. We did not observe the same enrichment in human Raji cells. The "
    "per-locus effect is modest even in mouse: in a given cell state most HSeeker loci, like "
    "most control windows, carry no single-stranded signal, so the loci are candidates rather "
    "than sites demonstrated to be structured."
)

LIMITATIONS = (
    "Limitations. HSeeker has several limitations that bear on how its output should be "
    "read. The stability score is an experimentally informed heuristic, not a calculated free "
    "energy. It assigns fixed weights to matched G and A positions, a uniform penalty to "
    "every other arm pairing, a constant reward to adjacent matched positions and an "
    "escalating penalty to runs of disrupted ones, and it then selects the highest-scoring "
    "arm window. The underlying interaction and stacking energy parameters are drawn from the "
    "ab initio calculations of Park et al., which isolated nucleobase interactions by "
    "replacing backbones with hydrogen atoms and left out other factors, such as torsional "
    "stress, backbone structural distortions, superhelical stress, pH-dependent protonation, "
    "ionic composition, temperature, loop or spacer energetics, or competition with duplex "
    "and other non-canonical conformations. It is therefore a relative ranking of candidates "
    "within the scheme defined here, and it carries no physical units. We report no score "
    "threshold, and a candidate is called on detection alone. The benchmark labels come from "
    "assays performed under different conditions, including acidic buffers, Mg2+-containing "
    "buffers and negatively supercoiled plasmids. We record those conditions for every "
    "record, but a sequence characterised as H-DNA-forming under one set of conditions does "
    "not necessarily form the structure under physiological conditions, and the pooled "
    "performance figures should be read as discrimination of reported sequence labels across "
    "conditions. The composition filters are fixed rules rather than fitted parameters, but "
    "they were chosen with knowledge of the H-DNA literature and of the composition of this "
    "benchmark, so the benchmark is not an independent test of them. They are not free of "
    "cost: the homopolymer filter rejects poly-(dC)-(dG) tracts, which form an "
    "Mg2+-dependent triplex28, and the AT-content filter removes a class of motifs that "
    "constitutes 17-19% of the H-DNA found to be formed in published single-stranded DNA "
    "maps. Finally, HSeeker detects antiparallel, purine-motif H-DNA at mirror repeats. "
    "Tandem direct repeats, parallel triplexes requiring cytosine protonation, and "
    "intermolecular triplexes are outside its model, and the two forming sequences it misses "
    "in our benchmark are direct repeats. Its predictions identify sequences with the "
    "potential to adopt H-DNA; they are not evidence that the structure forms at a given "
    "locus in vivo."
)

DATA_AVAILABILITY = (
    "Data AvailabilityThe benchmark sequences, their provenance and the assay conditions "
    "under which each construct was characterised are provided as machine-readable files in "
    "the repository and as Supplementary Table 1. Scripts that regenerate every analysis and "
    "figure in this work, together with the parameter settings used for each, are provided in "
    "the same repository. The human chromosome 1 and hg38 primary assembly sequences used for "
    "performance benchmarking were obtained from the UCSC Genome Browser "
    "(http://hgdownload.soe.ucsc.edu/downloads.html#human). The E. coli K-12 reference genome "
    "used for genomic-context recovery was downloaded from NCBI. The single-stranded DNA maps "
    "and non-B DNA motif calls used for the genome-wide comparison are those published by "
    "Kouzine and colleagues. Raw benchmark commands, timing output and result files are "
    "included in the repository and will be archived in a public repository with a permanent "
    "identifier on publication."
)

REFERENCE_28 = (
    "28.Kohwi, Y. & Kohwi-Shigematsu, T. Magnesium ion-dependent triple-helix structure formed "
    "by homopurine-homopyrimidine sequences in supercoiled plasmid DNA. Proc Natl Acad Sci U S A "
    "85, 3781-3785 (1988)."
)

ALG1_CAPTION = (
    "Algorithm 1: HSeeker pipeline, order of stages. Candidate enumeration runs first. If the "
    "longest-arm overlap rule is selected, overlaps are resolved before any filter or score is "
    "applied; the AT-content filter is applied next, then scoring, then the homopolymer filter, "
    "which inspects the scored putative triplex when scoring is enabled. Score-based overlap "
    "resolution, if selected, runs last, after filtering and scoring."
)

ALG2_CAPTION = (
    "Algorithm 2: Center-outward mirror repeat candidate enumeration with optional purity RMQ "
    "prefilter. For each hinge centre and spacer length, HSeeker computes the maximum reachable "
    "arm length and, when the prefilter is enabled, discards centre-spacer pairs whose right arm "
    "cannot satisfy the composition constraint; the prefilter is exact and does not change the "
    "output. Candidates are then extended outward, matching base pairs at the inner ends of the "
    "spacer are absorbed into the arms where purity permits, and candidates that normalise to the "
    "same start, end and arm length are reported once."
)

ALG3_CAPTION = (
    "Algorithm 3: Overlap removal. Under the longest-arm rule, hits are sorted by genomic start "
    "coordinate and a candidate with a longer arm replaces an overlapping incumbent, with a "
    "shorter spacer breaking ties. Under the score-based rule, candidates are considered in "
    "descending order of total score and retained when they do not overlap a candidate already "
    "retained, with ties resolved by longer arm, then shorter spacer, then earlier coordinate."
)

ALG4_CAPTION = (
    "Algorithm 4: Heuristic stability scoring. Each candidate is oriented onto the purine-rich "
    "strand, the right arm is reversed for positional comparison, and a per-position pairing score "
    "and match indicator are computed. The highest-scoring arm window satisfying the spacer "
    "constraint is selected; if no window satisfies it, the detected arms are scored as they "
    "stand. The total is the sum of pairing and stacking contributions, clamped at zero."
)

TABLE1_CAPTION = (
    "Table 1: Detection performance on the revised benchmark. Values are point estimates with 95% "
    "confidence intervals: Wilson intervals for proportions and percentile bootstrap intervals for "
    "the Matthews correlation coefficient. The balanced benchmark combines 70 experimental records "
    "with 60 designed negatives. Construct-level records are those for which the published "
    "experiment tested the exact sequence in the benchmark."
)

TABLE1_HEADER = ["Evaluation set", "n", "TP", "FN", "TN", "FP",
                 "Sensitivity", "Specificity", "PPV", "NPV", "MCC"]

TABLE1_ROWS = [
    ["Balanced benchmark", "130", "63", "2", "65", "0",
     "0.969 (0.895-0.992)", "1.000 (0.944-1.000)", "1.000 (0.943-1.000)",
     "0.970 (0.898-0.992)", "0.970 (0.925-1.000)"],
    ["Experimental records", "70", "63", "2", "5", "0",
     "0.969 (0.895-0.992)", "1.000 (0.566-1.000)", "1.000 (0.943-1.000)",
     "0.714 (0.359-0.918)", "0.832 (0.560-1.000)"],
    ["Construct-level records", "57", "50", "2", "5", "0",
     "0.962 (0.870-0.989)", "1.000 (0.566-1.000)", "1.000 (0.929-1.000)",
     "0.714 (0.359-0.918)", "0.829 (0.556-1.000)"],
]

# (paragraph, old, new) - predictions described as putative unless demonstrated
TERMINOLOGY = [
    (10, "identify and score potential H-DNA-forming sequences",
         "identify and score putative H-DNA-forming sequences"),
    (30, "scalable annotation of highly probable candidate H-DNA loci",
         "scalable annotation of putative H-DNA-forming loci"),
    (101, "export predicted H-DNA-forming loci",
          "export putative H-DNA-forming loci"),
    (104, "retrieve annotated candidate H-DNA motifs",
          "retrieve annotated putative H-DNA motifs"),
    (108, "the predicted H-DNA motifs", "the putative H-DNA motifs"),
]

FIG2_CAPTION = (
    "Figure 2: Computational performance of HSeeker. A) Wall-clock time and speedup for the "
    "complete pipeline on human chromosome 1 (hg38, 248.9 Mbp) across 1-16 worker threads, "
    "mean of three runs; the dashed line is ideal linear speedup. B) Peak resident set size "
    "across the same sweep, with the complete hg38 primary assembly shown for reference. C) "
    "Parallel efficiency across the same sweep. The comparator tools are not plotted because "
    "two of the four provide no parallel mode, so a worker sweep is undefined for them; at 16 "
    "cores on the same node and the same day, chromosome 1 took 45.0 s with Triplexator and "
    "64.7 s with triplex run as 16 parallel processes, against 16.3 s for HSeeker, while "
    "nBMST and NeSSie took 96.5 s and 295.7 s on one core."
)

FIG3_CAPTION = (
    "Figure 3: Recovery of benchmark sequences embedded in a genomic background. The 130 "
    "benchmark sequences were inserted into the E. coli K-12 MG1655 genome at positions drawn "
    "under a fixed random seed, and the modified genome was scanned with HSeeker and with "
    "triplex at the settings used throughout. A) Percentage of records for which a tool "
    "returns the same call when the sequence is embedded as when it is scanned on its own. B) "
    "Detection performance in the genomic background. C) Matthews correlation coefficient as "
    "a function of the positional overlap required to match a reported hit to an inserted "
    "sequence; the dashed line marks the 0.80 criterion used elsewhere. D) Fraction of each "
    "inserted sequence covered by its matching hit, for the sequences HSeeker recovered; the "
    "unmodified genome yields 83 loci across 4.64 Mb. This analysis measures recovery and the "
    "stability of a call across sequence context; it does not test whether H-DNA forms."
)

FIG4_CAPTION = (
    "Figure 4: HSeeker predictions compared with single-stranded DNA mapped in cells. HSeeker "
    "was run at its default settings on the genome assemblies used by Kouzine and colleagues. "
    "A) Recall: the percentage of the H-DNA motifs they found single-stranded that are "
    "overlapped by an HSeeker locus, in activated mouse B cells and in human Raji cells. B) "
    "Selectivity: the odds ratio that one of their predicted H-DNA motifs is found "
    "single-stranded when HSeeker reports it against when it does not, for all analysed "
    "motifs and for motifs that 36-nt reads can map to. C) Single-stranded signal at HSeeker "
    "loci relative to control windows matched for length, A/T content and purine/pyrimidine "
    "asymmetry. All loci uses genome-wide controls; isolated loci restricts loci and controls "
    "to mappable sequence more than 100 bp from other predicted non-B motifs and more than 2 "
    "kb from annotated transcription start sites, with controls drawn within 50 kb; authors' "
    "formed motifs applies the identical procedure to the motifs they called formed; placebo "
    "compares one control against the remaining controls of the same locus. D) Isolated loci "
    "stratified by mirror-repeat arm length. Points are ratios of mean signal, bars are 95% "
    "bootstrap confidence intervals, and the dashed line marks no difference. Scales in B-D "
    "are logarithmic."
)

SUPPFIG1_CAPTION = (
    "Supplementary Figure 1: Classification of the balanced benchmark by each tool. Confusion "
    "matrices for HSeeker and the four comparators on the 130-record benchmark (65 forming, "
    "65 non-forming), each tool at its documented default settings. Cells give the number of "
    "records and the percentage of the corresponding reported class. NeSSie is shown for its "
    "triplex-potential output and nBMST for the mirror repeats it classifies as "
    "triplex-forming."
)


if __name__ == "__main__":
    main()
