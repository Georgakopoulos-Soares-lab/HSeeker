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

(none yet)
