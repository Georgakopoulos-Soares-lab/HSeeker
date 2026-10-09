# X8: genome-wide comparison with Kouzine et al. 2017 (pre-specified design)

**Status:** design fixed 2026-10-09, **before** any result was computed. Data download and
alignment are running. Any later deviation from this design is to be recorded below with its
reason.
**Decision:** DEC-14, which replaces the E. coli insertion experiment (R2.3, R3.7, R3.8).

## Question

R2.3 asks whether HSeeker's predictions are a subset of the non-B DNA that has actually been
observed in a genome tested for it. R3.8 asks for an orthogonal genome-wide validation that
accounts for sequence composition. Translated into testable statements:

1. **Enrichment:** HSeeker loci carry more single-stranded-DNA (ssDNA) signal than
   composition-matched control intervals.
2. **Subset:** a larger fraction of HSeeker loci falls in ssDNA-enriched regions than of the
   matched controls.
3. **Score (secondary):** within HSeeker loci, ssDNA signal increases with the stability score.
   This is an out-of-sample check of the heuristic score (R3.3, R3.4).
4. **Activity (secondary):** enrichment is larger in activated than in resting B cells. Kouzine et
   al. report transcription- and supercoiling-dependent non-B formation, and H-DNA requires
   negative supercoiling.

## Data

Kouzine F. et al. 2017, *Cell Systems* 4:344-356 (PMID 28237796). Raw reads are from SRA072844
(ENA). The authors' processed peak tables are not accessible from our environment.

| Group | Runs | Genome |
|---|---|---|
| mouse LPS+IL4-activated B cells, KMnO4/S1 ssDNA-seq | 13 | mm10 |
| mouse resting B cells, ssDNA-seq | 4 | mm10 |
| mouse activated B cells, input / sonicated DNA (control) | 2 | mm10 |
| human Raji (Burkitt lymphoma), ssDNA-seq | 7 | hg38 |

**Processing.** Alignment uses bowtie2 2.5.5 with default end-to-end settings against the iGenomes
UCSC mm10/hg38 prebuilt indexes. Reads are kept at MAPQ ≥ 10 and merged per group with samtools
1.24. Analysis covers chromosomes 1–19/22, X and Y only; random, Un and chrM contigs are excluded.
Human Raji has no matched input in this study, so its controls are composition-matched intervals
only (stated as a limitation).

## HSeeker predictions

The genome is scanned with HSeeker in two configurations:
- the **shipped CLI defaults:** minrep 10, mismatch 0.10, purity 0.90, AT and homopolymer
  filters on, greedy overlap removal;
- the **benchmark configuration:** minrep 8, no filters.

The **primary** configuration is the one the authors freeze in Phase 1.3; the other is reported as a
sensitivity analysis. The exact configuration is recorded in the outputs.

## Signal

For each interval, ssDNA signal is the number of MAPQ ≥ 10 reads overlapping the interval
(`samtools bedcov`), expressed as reads per kb per million mapped reads (RPKM). For mouse, an
enrichment value log2((ssDNA RPKM + 0.5) / (input RPKM + 0.5)) is also computed. The input
normalisation corrects for mappability and copy number.

## Controls (composition matching)

For every HSeeker locus, **10** control intervals are sampled. Each matches the locus on:
- the same chromosome and the same length;
- GC fraction within ±0.02;
- purine-strand asymmetry, max(GA, CT) fraction, within ±0.02;
- no overlap with any HSeeker locus (±1 kb);
- no N bases;
- input coverage present (mouse only), so unmappable regions are excluded symmetrically.

Sampling uses a fixed seed (20261009) from a genome-wide pool of candidate windows. A locus with
fewer than 10 matches after 1,000 draws keeps the matches it has; the count is reported.

## Tests and reporting

1. **Enrichment:** median and mean signal, loci vs controls. Effect is reported as the ratio of
   means with a 95 % bootstrap CI (1,000 resamples of loci with their own controls), plus a
   Wilcoxon rank-sum p-value.
2. **Subset:** ssDNA-enriched regions are called from 200 bp windows, as a Poisson excess over the
   input-scaled local background (mouse) or over the genome-wide mean (Raji), Benjamini–Hochberg
   FDR 1 %. The fraction of loci overlapping an enriched region is compared with the controls'
   fraction (odds ratio with 95 % CI). The fraction of enriched regions containing an HSeeker
   locus is reported too. That second fraction is not expected to be high, because ssDNA also
   arises from promoters, G4, Z-DNA and other structures.
3. **Score:** Spearman correlation between stability score and enrichment within loci, plus signal
   by score quartile.
4. **Activity:** activated vs resting enrichment ratio, with bootstrap CI.

No threshold, window size or control parameter will be tuned after looking at the results. If a
step proves infeasible, the change and its reason are recorded here.

## Limitations (to be stated)

- ssDNA-seq detects single-strandedness from any cause (H-DNA's single-stranded loop and
  displaced strand, but also transcription bubbles, G4, Z-DNA, cruciforms). It shows association,
  not proof of triplex formation at a given locus.
- Mouse B cells and human Raji are specific cell states.
- Raji has no input control.
- The reprocessing (alignment, MAPQ, enrichment calling) is ours, not the authors'.

## Deviations log

**Deviation 1 (2026-10-09, before any ssDNA signal was computed): control sampling.**
- *What failed.* The pre-specified sampling (1,000 random draws per locus on the same chromosome,
  GC and purine asymmetry within ±0.02) produced only 0.18 controls per locus on mm10, and
  828,757 of 828,939 loci had fewer than 10.
- *Why.* HSeeker loci are strongly purine-asymmetric: on chr19 the asymmetry quartiles are
  0.87 / 0.95 / 1.00. Windows that asymmetric are rare among random genomic positions. On chr19
  only 7,213 of 3.38 M unblocked 30-bp windows (stride 10) reach asymmetry ≥ 0.88, and the ±1 kb
  exclusion blocks 39 % of chr19. Random draws therefore almost never hit a match.
- *Change.* Controls are now drawn from a genome-wide pool:
  - every N-free window outside the ±1 kb exclusion is enumerated at a 10 bp stride, for each
    locus-length class (window length = median locus length in the class);
  - windows are binned by GC (0.02) × asymmetry (0.02);
  - each locus draws 10 controls with replacement from its own cell;
  - only if a cell is empty is the asymmetry bin widened by one, then two, bins, and every widening
    is counted;
  - unmatched loci are excluded and counted.
- *What this gives up.* The same-chromosome requirement and sampling without replacement. Both
  are reported in `controls_meta.json`.
- *What is unchanged.* Matching on GC, purine asymmetry and length, and every test and threshold.
- *Why the controls are meaningful.* They are purine-rich windows that HSeeker does not call,
  i.e. purine richness without mirror-repeat structure. The comparison therefore asks whether H-DNA
  potential adds single-stranded signal beyond composition.

**Deviation 2 (2026-10-09, before any ssDNA signal was computed): second control set.**
- *Why.* With pool matching, 525,092 of 828,939 mm10 loci (63 %) found asymmetry-matched controls
  (62,133 after widening by 1 bin, 9,564 by 2 bins). The 303,847 unmatched loci are mostly the
  purest purine/pyrimidine runs. Excluding them would drop HSeeker's most H-DNA-like calls and
  bias the comparison.
- *Change.* A second control set (names `G{i}_{k}`, 10 per locus) is matched on length class and
  GC bin only, drawn from the same pool aggregated over asymmetry, and covers all loci.
- *Which is primary.* The asymmetry-matched set remains **primary**: it answers "beyond
  composition?". The GC-matched set answers "across all loci?". Both are reported, with every test
  run identically on each.

**Deviation 3 (2026-10-09, before any ssDNA signal was computed): implementation bug in deviation 1.**
- *The bug.* Pool windows of a class length cannot take every GC or asymmetry bin. Their fractions
  come in steps of 1/length, e.g. 1/22 ≈ 0.045, which is coarser than the 0.02 bins. A 21-bp locus
  with GC 7/21 (bin 16) could never match a 22-bp window (bins 15 or 18). That is why 280,021 loci
  found no GC-only match.
- *The fix.* Loci up to 60 bp (≈ 90 % of loci; p90 length 73) are matched against windows of
  **their exact length**, with the same GC count and the same purine-asymmetry count. Widening is
  by ±1, then ±2 asymmetry counts. This is stricter than ±0.02 and has no discretisation mismatch.
  Loci longer than 60 bp keep the length-class pools with 0.02 bins.
- *Unchanged.* All tests and thresholds.

**Deviation 4 (2026-10-09, before any ssDNA signal was computed): the authors' own data become primary.**
- *What we found.* The authors' processed data are public on the Przytycka lab page cited in the
  paper's Key Resources table:
  https://www.ncbi.nlm.nih.gov/CBBresearch/Przytycka/software/nonbdna.html
  It holds:
  - `nonB_DNA_predicted.tar` — predicted non-B motifs (SMnB), mm9 / hg19;
  - `nonB_DNA_ssDNA_enriched.tar` — the motifs the authors called as formed non-B structures
    (ssDNA+ SMnB);
  - `ssDNA_wiggle.tar` — their ssDNA-seq signal for activated B, resting B and Raji.
- *Verified against the paper.* Mouse activated-B ssDNA+ H-DNA = 17,109 = 16,876 + 233 (Table 1).
  Raji ssDNA+ H-DNA = 8,356 = 8,008 + 348 (Table S1B). Predicted mouse H-DNA SMnB = 728,355
  (main text).
- *Change.* Because these are the authors' published calls, they are independent of any processing
  choice of ours. They therefore become the **primary** comparison. HSeeker scans mm9 and hg19
  directly (no liftover), with the same configurations.
- *Primary tests (fixed now):*
  - **P1, recall:** fraction of the authors' ssDNA+ H-DNA structures overlapped by an HSeeker locus
    (mouse activated B; Raji).
  - **P2, selectivity:** among the authors' predicted H-DNA motifs, the odds that a motif is ssDNA+
    when HSeeker calls it vs when it does not (odds ratio, 95 % CI). Within HSeeker-called motifs,
    stability score of ssDNA+ vs ssDNA− motifs (Mann–Whitney).
  - **P3, enrichment:** the authors' activated-B and resting-B ssDNA signal (wiggle) at HSeeker loci
    vs the composition-matched controls, with exactly the control construction and statistics
    specified above.
- *Demoted.* Our reprocessing of the raw reads (mm10 / hg38) becomes a secondary reproducibility
  check, with the analysis as originally specified.
- *Note on comparability.* The authors' H-DNA motif definition (Inverted Repeats Finder with mirror
  option; ≥ 90 % purine or pyrimidine; loop ≤ 4 bp; motif ≤ 80 bp; ≥ 90 % matching; SINE-proximal
  motifs removed) differs from HSeeker's. P2 compares the two directly.

