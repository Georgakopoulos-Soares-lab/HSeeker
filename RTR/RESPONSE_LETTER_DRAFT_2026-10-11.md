# Response to reviewers — draft of 2026-10-11

Manuscript: HSeeker (Briefings in Bioinformatics). One answer per comment, in the order of the
decision letter. All analyses are complete (no placeholders remain); text in the final section,
"Notes for the authors", is not for the reviewers.

All numbers refer to the revised benchmark (`hdna_benchmark_balanced_v3.csv`, 130 records: 65 forming /
65 non-forming; 70 experimental records and 60 designed negatives) and the revised default settings
(minimum arm length 8 nt, maximum spacer 10 nt, purity 0.90, mismatch tolerance 0.10, AT filter 0.80,
homopolymer filter on, longest-arm overlap selection), which are now the same in the Python package, the
command-line tool and the web application.

---

## Reviewer 1

> HSeeker, an algorithm for systematic H-DNA sequence identification, is well presented, … nothing to
> change in this reviewer's opinion.

We thank the reviewer for the positive assessment. Although no changes were requested, the revisions
made for the other reviewers (an expanded and balanced benchmark, family-aware validation, a genome-wide
comparison with experimental data and revised runtime measurements) also apply to the points the
reviewer highlighted.

---

## Reviewer 2

### Comment 1 — discriminating power, imbalanced benchmark, negative controls

We agree that a benchmark with nine negatives cannot establish discriminating power, and that accuracy
is uninformative at that imbalance. We rebuilt the benchmark as follows.

- **Experimental records.** Every record was re-checked against its primary source (see Reviewer 3,
  comment 5). The experimental set now has 70 records (65 forming, 5 non-forming).
- **Designed negatives, as the reviewer suggested.** We added 60 negatives, generated with a documented
  recipe and a fixed random seed, so that the benchmark is balanced (65:65): 22 positives whose mirror
  repeat was disrupted by purine-preserving substitutions (the "perturbed positive cases" the reviewer
  proposed), 11 dinucleotide-preserving shuffles of positives, 12 random sequences matched to positives
  in length and GC content, 12 sequences known to form other structures (four G-quadruplex sequences,
  including the c-KIT and telomeric G4s, three Z-DNA (CG)n repeats and five random B-DNA sequences), and
  three perfect mirror repeats with balanced base composition, which test that mirror symmetry alone is
  not called H-DNA. Every designed negative was checked to contain no H-DNA-competent mirror repeat.
- **Structure-generation tools.** BoltzGen and RFdiffusion design or predict macromolecular structures;
  they cannot establish whether a DNA sequence adopts H-DNA under given conditions, so their output could
  not provide negative labels. The perturbation-based controls above implement the reviewer's
  alternative suggestion directly.
- **Poly-A and similar cases.** Single-base homopolymers are reported as a separate panel rather than in
  the benchmark metrics, because the composition filters decide their calls by construction and a run of
  one base is mirror-symmetric at every position. With the default filters HSeeker rejects A20 (shown not
  to form H-DNA by Hanvey et al. 1988) and poly-A15/A30 and poly-T25/T35; without filters it calls four
  of them. The panel also shows a limitation that we now state: the homopolymer filter rejects
  poly-(dC)30, which forms an Mg²⁺-dependent triplex (Kohwi & Kohwi-Shigematsu 1988).

**Result.** HSeeker detects 63 of 65 forming records and none of the 65 non-forming records
(sensitivity 0.969, 95 % CI 0.895–0.992; specificity 1.000, 0.944–1.000; MCC 0.970, 0.926–1.000).
Mirror-repeat detection alone, without the purity requirement, calls 7 negatives (B-33, one G4, two Z-DNA
repeats and the three balanced perfect mirrors); the purity requirement removes all seven. The two
forming records HSeeker misses are a tandem direct repeat, (ACCCTCCCC)4, which lies outside the
mirror-repeat model, and the γ-globin −228/−189 element.

**Comparison with other tools on the same records.** We ran triplex, Triplexator, and the two
mirror-repeat tools Reviewer 3 names (NeSSie and nBMST) on the same 130 records, each at its documented
defaults, with geometry matched to HSeeker, and tuned within training folds, with settings fixed before
any tool was run (details under Reviewer 3, comment 6). At default settings HSeeker had the highest MCC
(0.970, 95 % CI 0.909–1.000, family-cluster bootstrap), followed by triplex (0.819), NeSSie (0.583),
nBMST (0.449) and Triplexator (0.024). triplex detected almost the same forming records (62 of 65) but
called 9 designed negatives forming, 7 of them mirror-disrupted mutants, so its lower specificity comes
from exactly the perturbed controls the reviewer proposed. On the 70 experimental records alone, HSeeker
and triplex differ by a single record and perform comparably (MCC 0.832 vs 0.772). In the homopolymer
panel, triplex and NeSSie call four of the five A/T homopolymers forming; HSeeker calls none.

### Comment 2 — sequences as CSV and a script to reproduce the analysis

We agree. The benchmark is now provided as CSV files with full provenance:
`hdna_benchmark_experimental_v3.csv` (all 81 curated records, with the curation decision and its
justification, original and corrected citations, family and study identifiers, and the experimental
conditions of each assay) and `hdna_benchmark_balanced_v3.csv` (the 130-record benchmark). The analyses
are reproduced by scripts in `analysis/scripts/`, with all settings pinned in
`analysis/benchmark_config.json`: `run_direct_sensitivity_analysis.py` (benchmark metrics),
`run_parameter_robustness.py` and `run_grouped_validation.py` (cross-validation and family/study-aware
validation), `run_ablation_comment9.py`, `run_score_threshold_oof.py`, `run_overlap_comparison.py`,
`run_homopolymer_panel.py` and the Kouzine comparison scripts. The designed negatives are regenerated by
`RTR/generate_synthetic_negatives.py`.

### Comment 3 — E. coli insertion; validation on a genome tested for non-B DNA

We agree that the insertion experiment tested recovery of known sequences in a genomic background, not
H-DNA formation, and we have removed it from the manuscript (see also Reviewer 3, comment 7). Following
the reviewer's suggestion, we compared HSeeker with the genome-wide data of Kouzine et al. (2017), who
mapped single-stranded DNA in cells with permanganate and S1 nuclease footprinting in activated mouse B
cells and human Raji cells, predicted non-B DNA motifs genome-wide, and called a motif formed when it
carried single-stranded signal. We used their published calls and signal tracks on the genome builds
they used, and scanned both genomes with HSeeker at its default settings.

**Recall.** HSeeker loci overlap 79.0 % of the H-DNA structures the authors found formed in mouse
(13,517 of 17,109) and 82.3 % in human Raji cells (6,879 of 8,356).

**Selectivity.** Among the authors' own predicted H-DNA motifs, those HSeeker calls are formed more often
than those it does not: 4.5 % vs 2.8 % in mouse (odds ratio 1.65, 95 % CI 1.59–1.71) and 5.7 % vs 1.7 %
in human (3.39, 3.20–3.59); 1.68–2.07 and 4.16–5.70 on motifs that the 36-nt reads can map to.

**On the subset criterion.** A strict subset is not attainable by any sequence-based predictor on these
data, because formation depends on cell state: only about 3–5 % of the authors' own predicted H-DNA
motifs were single-stranded in a given cell type. The informative questions are whether HSeeker's calls
are enriched for structures that form, which they are, and how many formed structures it finds, about
four in five.

**What HSeeker misses.** Most formed H-DNA that HSeeker does not report are short A/T-rich or
single-base motifs, which the composition filters exclude. In the authors' data this class is formed at
about half the rate of other predicted H-DNA in mouse and a quarter in human, which supports excluding
it by default; it is nevertheless 17–19 % of the formed H-DNA, and users can disable the filters.

The enrichment of the measured single-stranded signal at HSeeker loci is described in our response to
Reviewer 3, comment 8 (new Figure 3).

### Comment 4 — thermodynamic simplification; limitations section

We agree and have added a Limitations section. It states that the stability score is an experimentally
informed heuristic, not a free-energy model: it combines fixed pairing weights for matched G and A
positions (7.6 and 3.83, from triplex interaction energies), a uniform penalty for other pairings
(−5.715), a stacking term of 5.0 for adjacent matches and an escalating penalty for consecutive
mismatches, with boundary optimisation. It does not model backbone geometry, supercoiling, pH, ionic
conditions, temperature or competition with other structures, and it is not calibrated against measured
ΔG, melting temperatures or superhelical transition points. Benchmark labels come from assays performed
under different conditions (acidic S1 digestion, Mg²⁺-containing buffers, supercoiled plasmids), which
we now list for every record. The section also covers the homopolymer filter's rejection of poly-(dC)
tracts, the dependence of genome-wide formation on cell state, and that HSeeker predicts sequence
potential, not formation at a given locus. We removed kcal/mol units and all language describing the
score as a free energy (see also Reviewer 3, comment 4).

### Comment 5 — complexity and sensitivity analyses belong to Results

We agree and have moved both analyses to the Results.

### Comment 6 — filters defined on the same dataset (overfitting)

We thank the reviewer for this point and address it in three ways.

1. **Disclosure.** The AT-content filter (motifs ≥ 80 % A/T) and the homopolymer filter are fixed rules
   that were chosen with knowledge of the H-DNA literature and of the original benchmark's composition.
   In the original benchmark, the only false positives at the analysis settings were four homopolymers,
   and the filters removed them. In the revised benchmark the homopolymers are reported separately (see
   comment 1), so the filters change no benchmark call: the reported metrics do not depend on them, and
   equally the benchmark provides no evidence for them.
2. **Independent evidence for the filters.** The genome-wide data of Kouzine et al. provide such
   evidence, independent of our benchmark: predicted H-DNA motifs that the filters exclude are formed at
   about half the rate of other predicted H-DNA in mouse (odds ratio 0.45–0.58) and a quarter in human
   (0.15–0.30). Without the filters HSeeker calls almost all predicted motifs and loses its selectivity
   (odds ratio 0.54).
3. **Selection out of sample.** Every parameter was evaluated with selection inside training folds only
   (576 combinations of purity, mismatch, arm length, spacer and filter state). Under nested
   cross-validation the per-fold tuned configuration reached MCC 0.939 (95 % CI 0.877–0.985), below the
   fixed default configuration (0.970), and tuning beat the fixed configuration in none of 200 random
   train/test splits. Under leave-one-family-out and leave-one-study-out validation the results were the
   same (tuned 0.924, fixed 0.970).

### Comment 7 — confidence intervals for the tool comparison

We agree. All metrics are now reported with 95 % confidence intervals: Wilson intervals for sensitivity,
specificity, PPV and NPV, and bootstrap intervals for MCC, at the record level and with families or
studies as the resampling unit (see Reviewer 3, comments 1 and 2). For the tool comparison we report
paired differences (HSeeker minus each tool) in sensitivity, specificity and MCC, with bootstrap
confidence intervals that resample records or whole sequence families, and exact McNemar tests on the
records where the tools disagree. We fixed in advance that we would call HSeeker better only where the
family-level interval of the MCC difference excludes zero and the difference has the same sign with
matched and with tuned settings. By that rule HSeeker outperforms all four tools on the balanced
benchmark: ΔMCC +0.15 (95 % CI 0.03–0.22; McNemar 10 vs 0 discordant records, p = 0.002) against triplex,
+0.39 (0.18–0.56) against NeSSie, +0.52 (0.40–0.63) against nBMST and +0.95 (0.71–1.11) against
Triplexator. Against triplex the difference comes from the designed negatives; on experimental records
alone the interval reaches zero (+0.06, 0.00–0.22) and we describe the two tools as comparable there. The
manuscript's wording now follows these results.

### Comment 8 — `-march=native`

*(Keep the answer in the current letter; it is correct for the code.)* We thank the reviewer for this
comment. The old cross-compilation check only inspected macOS `ARCHFLAGS`; on an ordinary Linux build it
allowed `-march=native` by default. That was not appropriate for redistributable binaries. Linux and
macOS builds now use portable compiler targets by default (`-O3 -Wall`); native tuning is available only
when a user sets `HSEEKER_NATIVE_BUILD=1` for a local build, and the option rejects a different macOS
target architecture. Our wheels do not enable it. A regression test checks the compiler flags.

### Comment 9 — unsynchronised profiling counters

*(Keep the answer in the current letter.)* We removed the shared counters from the scan loops. Each scan
counts in its own structure while the GIL is released and publishes the counts after reacquiring it, so
concurrent scans neither lose nor mix counts; a regression test runs 40 concurrent scans and checks the
totals. The counters are diagnostic only and no reported number depends on them.

### Comment 10 — upper-case `putative_triplex`

*(Keep the answer in the current letter.)* All sequence fields, including `putative_triplex`, are now
lower case in the Python API and the command-line output, with regression tests. `putative_triplex` can
still differ from `full_sequence` because the scorer may optimise the arm boundaries or represent a
CT-rich motif on the reverse-complement strand.

### Comment 11 — two `app.py` files

We agree that two front-ends were confusing. We removed the Streamlit demo (`src/hseeker/app.py`) and its
dependency from the package. HSeeker now has a single web interface, the FastAPI application in
`webapp/main.py`, which is the hosted service and uses the same Python scanning API as the command-line
tool; for the same settings the web application and the command-line tool now give identical results.
We thank the reviewer for the kind words about the software engineering.

---

## Reviewer 3

### Comment 1 — benchmark size, imbalance, confidence intervals

We agree. The benchmark is now balanced (130 records, 65:65; see Reviewer 2, comment 1), and we report
all metrics with 95 % confidence intervals on three sets, stating the size of the negative class
plainly:

| Set | n | TP / FN / TN / FP | Sensitivity | Specificity | PPV | NPV | MCC |
|---|---|---|---|---|---|---|---|
| Balanced benchmark | 130 | 63 / 2 / 65 / 0 | 0.969 (0.895–0.992) | 1.000 (0.944–1.000) | 1.000 (0.943–1.000) | 0.970 (0.898–0.992) | 0.970 (0.926–1.000) |
| Experimental records | 70 | 63 / 2 / 5 / 0 | 0.969 (0.895–0.992) | 1.000 (0.566–1.000) | 1.000 (0.943–1.000) | 0.714 (0.359–0.918) | 0.832 (0.569–1.000) |
| Primary experimental records | 58 | 51 / 2 / 5 / 0 | 0.962 (0.872–0.990) | 1.000 (0.566–1.000) | 1.000 (0.930–1.000) | 0.714 (0.359–0.918) | 0.829 (0.556–1.000) |

Only five experimental non-forming constructs remain after source-level verification. The specificity
estimate rests mainly on the designed negatives, and we say so in the text; conclusions about specificity
are qualified accordingly.

### Comment 2 — related sequence families; family-aware validation

We agree that the records are not independent. We grouped them by sequence family and by study using the
curated provenance (the 70 experimental records span 48 families and 26 studies; designed negatives
inherit the family and study of the record they were derived from, giving 52 family and 27 study groups;
the largest family, the pXY32 substitution series and its derivatives, holds 37 records). Results:

| Validation | Fixed default configuration: MCC (95 % CI) | Per-fold tuned |
|---|---|---|
| Leave-one-family-out | 0.970 (record CI 0.925–1.000; family-cluster CI 0.907–1.000) | 0.924 |
| Leave-one-study-out | 0.970 (study-cluster CI 0.901–1.000) | 0.924 |
| One record per family | 0.862 (52 records) | — |
| One record per family and class | 0.945 (74 records) | — |

With the parameters fixed, the call for each record depends only on its own sequence, so family-aware
validation leaves the point estimate unchanged; what it changes is the precision: the family-cluster
interval is wider than the record-level one, and we report both.

### Comment 3 — circularity of the score threshold

We agree, and we re-evaluated the threshold exactly as suggested. Within each training fold we selected
the threshold that maximised Youden's J and applied it only to the held-out records (repeated stratified
5-fold cross-validation, leave-one-family-out and leave-one-study-out). A threshold chosen this way never
improved performance: out of fold it gave MCC 0.949 (5-fold, 95 % CI 0.937–0.961) and 0.955
(leave-one-family-out and leave-one-study-out, 0.898–1.000), against 0.970 (0.925–1.000) without any
threshold. It removed no false positive, because after the purity requirement no non-forming record has
any candidate, and it lost forming records whose score fell below a threshold learned without them. For
the same reason the score's AUC on this benchmark (0.985) only restates the detection result. We
therefore no longer use or report a score threshold: a sequence is called H-DNA-forming when it contains
a candidate that passes the detection criteria, and the stability score ranks candidates. The
Youden-selected threshold, the in-sample AUC and the threshold-based confusion matrix were removed (former
Figure 3, now replaced), and the web interface no longer labels sites with score ≥ 50 as "stable".

### Comment 4 — thermodynamic interpretation of the score; assay conditions

*(Keep the answer in the current letter, which is correct, with one addition.)* We agree and describe
the output as an experimentally informed, sequence-based heuristic stability score, not a calculated free
energy or a quantitative prediction of H-DNA formation. [Current letter text on the weights 7.6 / 3.83 /
−5.715, the stacking term 5.0, the removal of kcal/mol labels, and the revised main-text sentence.]
**Addition:** As promised, the benchmark now records the experimental conditions of every experimental
record (assay, DNA context, pH, temperature, ionic conditions, supercoiling), taken from the cited papers
with the supporting passage, and marks values that the source does not report. The pooled results are
described as discrimination across assay conditions, and the Limitations section states that a label
obtained under one condition does not imply formation under physiological conditions.

### Comment 5 — benchmark provenance and quality control

We audited every record against its primary source and deposit the per-record audit. The outcome:

- Six records carried low-confidence annotations from a secondary document. Three (GG32-long, Mycm1,
  MycM2) had no traceable primary source and were removed. Three (GG32 and AG32 30-nt versions, MYCAA)
  had traceable sources but were removed with all other low-confidence-tier records.
- Two records were excluded on validity grounds: an inosine-substituted oligonucleotide that cannot be
  represented as A/C/G/T DNA, and an intermolecular triplex-forming-oligonucleotide record mislabelled as
  intramolecular H-DNA.
- The two pure single-base homopolymers (A20, non-forming; poly-(dC)30, forming) were moved to a separate
  panel, because the composition filters decide their calls by construction (see Reviewer 2, comment 1).
- Citation identifiers were corrected for 25 records (PMIDs) and 9 records (DOIs); the original citations
  are kept in the files. For the four Friedreich ataxia (GAA)n records, whose cited identifiers pointed to
  unrelated papers, we now cite the primary experimental studies (Potaman et al. 2004 for n = 10, 20, 33;
  Sakamoto et al. 1999 for n = 66); because those studies tested related constructs rather than these
  exact sequences, the records are kept in a clearly labelled secondary set.
- **pRW1406 / pRW1411.** The primary abstract (Hanvey et al. 1988) resolves the inconsistency: the
  (TCC)8 insert is pRW1411, and the (TTC)8 insert, which the abstract lists among the five inserts that
  form the intramolecular triplex, is a distinct construct. We kept pRW1411 and added (TTC)8 as a new
  primary record, named after its verified sequence rather than the plasmid number, whose mapping came
  from the unreliable secondary document.
- The benchmark now records, for every experimental record, the experimental conditions under which it
  was tested (see comment 4).

The primary benchmark set consists of 58 construct-level records (53 forming, 5 non-forming); 12
class-level records form the secondary set.

### Comment 6 — comparator set; claim strength

We agree. We added the two mirror-repeat tools the reviewer names, NeSSie and nBMST (non-B_gfa), alongside
triplex and Triplexator, and ran all of them on the revised 130-record benchmark. The settings, the rule
for turning each tool's output into a call, the statistics and the wording rule were fixed and deposited
before any tool was run (`RTR/TOOL_COMPARISON_PLAN_2026-10-11.md`; code `analysis/scripts/
run_tool_comparison.py`; builds `benchmarks/tools/`). Each tool was run at its documented defaults, with
its geometry matched to HSeeker's (minimum arm 8 nt, spacer up to 10 nt) where it has such parameters, and
with a small parameter grid selected inside training folds only, the same opportunity HSeeker's defaults
had. A record was called forming when the tool reported at least one hit in it.

| Tool, default settings (hit type) | Sensitivity | Specificity | MCC (95 % CI) | Tuned out of fold, MCC |
|---|---|---|---|---|
| HSeeker | 63/65 | 65/65 | 0.970 (0.909–1.000) | 0.939 |
| triplex (any hit) | 62/65 | 56/65 | 0.819 (0.743–0.935) | 0.819 |
| NeSSie (triplex-potential motif) | 33/65 | 65/65 | 0.583 (0.421–0.779) | 0.547 |
| nBMST (mirror repeats classified as triplex) | 24/65 | 64/65 | 0.449 (0.351–0.561) | 0.597 |
| Triplexator (triplex target site) | 8/65 | 58/65 | 0.024 (−0.132–0.250) | 0.520 |

The differences reflect what each tool searches for, which we now describe in the manuscript.
Triplexator searches for target sites for a third strand supplied by a separate oligonucleotide and does
not require a mirror repeat: at its defaults it finds few forming records, and at the settings used in our
original manuscript it finds 64 of 65 but also calls all 22 mirror-disrupted mutants. NeSSie reports
k-mers that are themselves mirror-symmetric and purine- or pyrimidine-rich, not arm–spacer–arm motifs; it
made no false call but found none of the 16 pXY32 constructs. nBMST's mirror-repeat module matches
HSeeker's geometry, but its triplex classification counts purines and pyrimidines with integer division,
so it accepts only mirror repeats with more purines than pyrimidines on the reported strand: it found none
of the 19 forming records written on the pyrimidine-rich strand, whereas its unclassified mirror-repeat
output found 15 of them. triplex, the only other tool designed for intramolecular triplexes, detected
nearly the same forming records as HSeeker but called 9 designed negatives forming, mostly
mirror-disrupted mutants; on the experimental records alone the two are comparable.

Using paired differences with confidence intervals (Reviewer 2, comment 7), HSeeker outperforms all four
tools on the balanced benchmark, and Triplexator, NeSSie and nBMST also on the experimental records alone;
against triplex the advantage is confined to the designed negatives. We have revised the Key Points
accordingly: "On a balanced benchmark of experimentally characterised and designed sequences, HSeeker
achieved higher accuracy than four existing tools, matching triplex on experimentally characterised
sequences while rejecting mirror-disrupted controls that triplex accepted."

### Comment 7 — E. coli insertion experiment; overlap threshold

We agree. The insertion experiment tests computational recovery in a genomic background, not H-DNA
formation, and we have removed it, together with the former Figure 3, from the manuscript rather than
toning it down. In its place the manuscript compares HSeeker's predictions with experimental
single-stranded DNA maps (Reviewer 2, comment 3, and comment 8 below). Because the insertion experiment
is no longer part of the paper, the 80 % overlap criterion it used no longer enters any reported result.

### Comment 8 — genome-scale biological validation

We thank the reviewer for this suggestion. We compared HSeeker predictions with the structure-sensitive
nuclease mapping of Kouzine et al. (2017): permanganate/S1 footprinting of single-stranded DNA in
activated and resting mouse B cells and in human Raji cells. Recall and selectivity against the authors'
own calls are described in our response to Reviewer 2, comment 3. Here we describe two enrichment
analyses of the measured single-stranded signal (new Figure 3). The first was specified, with its
controls and statistics, before any result was computed. The second was added after the first human
result and was likewise specified before it was computed. We report both and deposit the design
document with every change and its reason.

**Composition matching.** Controls were matched to each HSeeker locus on length, A/T count and
purine/pyrimidine asymmetry, and did not overlap any HSeeker locus or predicted H-DNA motif. Matching on
asymmetry is essential: controls matched on A/T content alone carry 1.4–3 times more signal than HSeeker
loci.

**First analysis: all HSeeker loci, genome-wide controls.** In activated mouse B cells HSeeker loci
carried 14 % more single-stranded signal than matched controls (ratio of means 1.14, 95 % CI 1.13–1.15;
825,749 loci); in resting B cells 0.82 (0.82–0.83); in human Raji cells 0.98 (0.97–0.98). This analysis
does not separate H-DNA from other sources of single-stranded DNA and does not account for the poor
mappability of repetitive purine tracts.

**Second analysis: loci isolated from other sources of single-stranded DNA.** We removed, from loci and
controls alike, intervals the reads cannot map to, intervals within 100 bp of motifs predicted to form
G-quadruplex, Z-DNA, cruciform or stress-induced melting, and intervals within 2 kb of transcription start
sites, and drew controls within 50 kb. In activated mouse B cells isolated HSeeker loci carried 25 % more
signal than their controls (1.25, 1.18–1.33; a locus exceeded one of its own controls with probability
0.514, p = 2 × 10⁻¹⁷). The excess increased with arm length, from none detectable at 8–9 nt to 1.6–1.8 at
12 nt and longer (Spearman p = 3 × 10⁻²¹), as expected for H-DNA but not for the other structures, which
do not depend on mirror symmetry; it was absent in resting B cells (0.86), in which even the authors'
formed motifs showed no excess. In human Raji cells there was no detectable excess (1.02, 0.99–1.06;
p = 0.08). A placebo comparison of controls against controls showed no difference in either species, and
the motifs the authors called formed were enriched 7.4-fold (mouse) and 2.1-fold (human) under the same
procedure.

**Interpretation and limitations.** In activated mouse B cells HSeeker loci are associated with
single-stranded DNA beyond what composition, mappability, local context, the other predicted structures
and promoter transcription explain, and the association grows with the mirror-repeat arm length. We did
not observe this in human Raji cells, where HSeeker's calls are nevertheless strongly enriched among the
motifs the authors found formed; the Raji data have no input control and a smaller dynamic range. The
per-locus effect is small even in mouse: most HSeeker loci, like most controls, carry no signal in a given
cell state, so HSeeker loci are candidates, not sites shown to be structured, and single-stranded signal
does not prove triplex formation at an individual locus. These limitations are stated in the Discussion.

The chromosome 1 count quoted in the comment (96,729) used the previous default minimum arm length of
10 nt. With the revised defaults HSeeker reports 45,927 loci on hg38 chromosome 1 (106,481 without the
composition filters).

### Comment 9 — ablation and parameter robustness

We added a reproducible ablation and parameter analysis on the balanced benchmark (notebook
`analysis/hseeker_ablation_comment9.ipynb`, script `analysis/scripts/run_ablation_comment9.py`, tables in
`analysis/results/ablation_comment9/`). A record is called forming when at least one candidate survives;
no score threshold is used.

**Component ablation.** Mirror-repeat detection alone detected 63 of 65 forming and 7 of 65 non-forming
records (sensitivity 0.969, specificity 0.892, MCC 0.864); the seven were B-33, one G-quadruplex, two
Z-DNA repeats and the three perfect mirrors with balanced composition. Adding the purity requirement
removed all seven (specificity 1.000, MCC 0.970) without losing a forming record. Adding the AT filter,
pairing-only scoring, pairing plus stacking and the complete model (boundary optimisation and homopolymer
filtering) changed no call. Pairing-only and pairing-plus-stacking scores rank the same candidates, and no
non-forming record has any candidate once purity is applied, so the exploratory score AUROC (0.985) only
restates the detection result. On the 70 experimental records alone, MCC is 0.708 for detection alone and
0.832 from the purity step onward. The scorer has no separate switch for the mismatch penalty, so we do
not claim to have ablated it separately. Retaining all candidates, greedy selection and score-based
selection gave identical calls.

**Parameters, one at a time.** Purity 0.80–0.90 gave MCC 0.970, while 0.95–1.00 lost 9 forming records
(0.843). Mismatch tolerance 0–0.05 detected 49 of 65 forming records (0.778), 0.10 detected 63 (0.970),
and 0.15–0.20 detected 64 but added a false positive (0.969). Minimum arm length 6 added 10 false
positives (0.839), 7–8 gave 0.985 and 0.970, 9–10 lost 8 forming records (0.884) and 12 lost 16 (0.778).
Maximum spacer 3 detected 41 forming records (0.679), 5 detected 55 (0.856) and 8–15 detected 63 (0.970).
The single record that minimum arm 7 and mismatch 0.15 gain is the tandem direct repeat (ACCCTCCCC)4,
which is outside the mirror-repeat model.

**Joint variation, selected out of sample.** Across 576 combinations with selection inside training folds
only, the fixed default configuration performed as well as or better than per-fold tuning under nested
cross-validation, leave-one-family-out and leave-one-study-out validation (0.970 vs 0.939, 0.924, 0.924).
Minimum arm 8 was selected in 99 % of nested folds and in every family and study fold, and mismatch 0.10
in 48 of 52 family folds. On this basis we changed the default minimum arm length from 10 to 8 nt.
Performance therefore depends on the parameters as described, but the reported results do not rest on
one hand-picked parameterisation.

### Comment 10 — overlap-removal strategy

We thank the reviewer for the comment. We agree that arm length is not a validated proxy for biological
activity, and we now offer the alternative the reviewer proposes and compare it with the longest-arm
rule. We keep the longest-arm rule as the default for three reasons. First, among overlapping candidates
that share a mirror-repeat region, the longest arm describes the full extent of the repeat, which is the
region that can take part in strand exchange. Second, it is a purely geometric criterion that does not
depend on the stability score, which is a heuristic rather than a calibrated free energy (comment 4), so
small score differences between overlapping candidates should not decide what a user sees by default.
Third, it is faster on genomes, because overlaps are resolved before scoring.

HSeeker now offers three modes in the Python API, the command-line tool and the web interface: greedy
selection (the default); score-based selection (`overlap_strategy="score"`, `-overlap-strategy score`,
web: *Highest scores first*), which scores all candidates and keeps them in descending score order when
they do not overlap a candidate already kept (ties favour longer arms, then shorter spacers, then earlier
coordinates); and retaining all candidates (`remove_overlaps=False`, `-skipoverlap`, web: *Keep all raw
hits*). The AT filter runs before scoring; greedy selection runs before the scored homopolymer check and
score-based selection after it. We recommend the default for genome-wide annotation, score-based selection
for locus-level work such as choosing a variant to test, and retaining all candidates for exhaustive or
statistical analyses. A higher score helps prioritise candidates but is not evidence of formation.

*Examples* (`minrep=7, maxrep=25, maxspacer=5, purity=0.9, mismatch=0.2`, filters on). For
GAGAAAAGAGGGAGAAAAAGGAA both rules keep 1–21 (arm 10; the overlapping 3–19 candidate, arm 8, has the same
score, 53.635). For GGAAGAAGAAGGAGGAGAGAGGGAAGGAA greedy keeps 9–29 (arm 10, score 20.36) and score-based
selection keeps 1–15 (arm 7, score 38.58). An overlap chain need not collapse to one locus: for (GAA)20,
both rules keep 3–32 and 35–60.

*Comparison.* On the benchmark, 43 of the 70 experimental records have overlapping candidates; the two
rules return identical intervals for 38 and different intervals for 5, and the detection call is the same
for all 130 records and for all 74,880 record-by-configuration calls across our parameter grid. On hg38
chromosome 1 with the default settings, greedy selection reports 45,927 loci and score-based selection
47,217; every greedy locus is overlapped by a score-selected locus, 91.2 % have identical coordinates, and
score-based selection adds 991 loci, at a modest cost in time (23.9 s vs 16.3 s on chromosome 1 with
16 workers). The comparison and its scripts are deposited
(`analysis/results/overlap_selection/`).

### Comment 11 — computational benchmarking; extrapolation; data deposition

We agree and re-measured runtime and memory on the final code and default settings, on the same hardware
as in the original benchmark (AMD EPYC 7763, 16 vCPUs, 32 GB, Rocky Linux 8.10), with three runs per point
and peak memory from `/usr/bin/time -v`. On hg38 chromosome 1, HSeeker took 173.8 s with one worker and
16.3 s with 16 workers (a 10.7-fold speedup), with a peak memory of 1,334 MB. Instead of extrapolating, we
ran the complete hg38 primary assembly: 204.5 s (3.4 min) with 16 workers and a peak memory of 2,296 MB,
reporting 559,601 loci; the three runs gave identical output. The optional score-based overlap selection
took 23.9 s on chromosome 1. All commands, raw timing output, hardware descriptions and input checksums
are deposited in the repository (`benchmarks/results_2026-10-11/`, `benchmarks/parallel_sweep.py`), and
the code, data and benchmark outputs will be archived on Zenodo at publication; the Data Availability
statement now says so instead of "upon reasonable request". We also moderated the wording of the
performance claims (see comment 12). To compare fairly, we re-timed the other tools on the same node on
the same day. With 16 cores, chromosome 1 took 16.3 s with HSeeker, 45.0 s with Triplexator and 64.7 s
with triplex (run as 16 parallel processes); nBMST and NeSSie have no parallel mode and took 96.5 s and
295.7 s on one core. Each tool's output was identical across the three runs, and triplex and Triplexator
reproduced the hit counts of our original measurement exactly.

### Comment 12 — claims, terminology and figures

We agree. We removed or rephrased the three statements the reviewer quotes: the scoring system "works
highly effectively and aligns seamlessly", HSeeker as "the only H-DNA identification algorithm with a
proven scoring system", and analyses "previously computationally prohibitive". We now use "putative
H-DNA-forming sequence" or "predicted H-DNA motif" for computational output, reserve "H-DNA-forming
sequence" for experimentally demonstrated cases, and use these terms consistently. Figure 1 and the
algorithm panels were simplified, and the manuscript was language-edited.

---

## Notes for the authors (not for the reviewers)

1. **No runs pending.** Every number in this letter comes from a deposited script and output.
2. **HDNA0065 (XDP (CCCTCT)10)** is kept in the benchmark (author decision 2026-10-11) but its cited
   review states that no experimental data exist for this repeat; its label is inferred. It is currently
   marked `primary`, which conflicts with the R3.5 sentence about the primary set being construct-level
   records with experimental evidence. Either move it to the secondary set or adjust that sentence.
3. **Weaker evidence recorded in the condition columns, labels unchanged:** HDNA0055/0056 (the cited
   paper is a recombination assay; structure inferred from Mirkin 1987), HDNA0059 (triplex shown on
   oligonucleotides at pH 8), HDNA0066–0068 (S1-seq at pH 4.5 with Zn²⁺; the authors cannot exclude
   formation after lysis), HDNA0058 (full text not accessible).
4. **Manuscript edits** that the answers promise: Methods (defaults, filters, overlap modes, algorithm
   panels), Results (new Figure 3, benchmark table with CIs, ablation, runtime), Limitations, Data
   Availability, terminology. The manuscript file is not in the repository.
5. **R2.8–R2.10:** the current letter cites commit `063e3b5`; that commit is in `main`. The final release
   commit should be cited instead once tagged.
