# H-DNA Sensitivity Analysis Overview

## Primary Validation

The primary validation runs HSeeker and Triplex directly on the 80 curated
experimental sequences:

- Summary: `analysis/results/sensitivity_direct/summary.md`
- Per-sequence comparison: `analysis/results/sensitivity_direct/method_comparison_by_sequence.csv`
- ROC plot: `analysis/results/sensitivity_direct/plots/hseeker_score_roc.pdf`
- Confusion matrices: `analysis/results/sensitivity_direct/plots/confusion_matrices_hseeker_triplex.pdf`

This is the main benchmark for comparing HSeeker and Triplex sequence-level
classification.

## Secondary Validation

The secondary/context validation injects the same 80 sequences into the
E. coli K-12 `NC_000913.3` genome using seed 42 and a minimum inter-insertion
distance of 500 bp, then matches hits back to insertion intervals using the
80% overlap rule:

- Summary: `analysis/results/sensitivity_injected/summary.md`
- Per-insertion comparison: `analysis/results/sensitivity_injected/method_comparison_by_insertion.csv`
- ROC plot: `analysis/results/sensitivity_injected/plots/hseeker_score_roc.pdf`
- Confusion matrices: `analysis/results/sensitivity_injected/plots/confusion_matrices_hseeker_triplex.pdf`

This benchmark evaluates whether hits can be recovered from realistic genomic
context and flanking sequence. It should be reported as secondary to the direct
sequence validation.

## Current Headline Metrics

| analysis | method | TP | FN | TN | FP | sensitivity | specificity | F1 | MCC |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| direct primary | HSeeker direct hit | 68 | 3 | 7 | 2 | 0.958 | 0.778 | 0.965 | 0.703 |
| direct primary | HSeeker score Youden | 66 | 5 | 8 | 1 | 0.930 | 0.889 | 0.957 | 0.701 |
| direct primary | Triplex direct hit | 68 | 3 | 7 | 2 | 0.958 | 0.778 | 0.965 | 0.703 |
| injected secondary | HSeeker overlap | 69 | 2 | 7 | 2 | 0.972 | 0.778 | 0.972 | 0.750 |
| injected secondary | HSeeker score Youden | 69 | 2 | 8 | 1 | 0.972 | 0.889 | 0.979 | 0.822 |
| injected secondary | Triplex overlap | 31 | 40 | 9 | 0 | 0.437 | 1.000 | 0.608 | 0.283 |
