# H-DNA Sensitivity Analysis Overview

## Primary Validation

The primary validation runs HSeeker and Triplex directly on the 54 curated
experimental sequences:

- Summary: `analysis/results/sensitivity_direct/summary.md`
- Per-sequence comparison: `analysis/results/sensitivity_direct/method_comparison_by_sequence.csv`
- ROC plot: `analysis/results/sensitivity_direct/plots/hseeker_score_roc.pdf`
- Confusion matrices: `analysis/results/sensitivity_direct/plots/confusion_matrices_hseeker_triplex.pdf`

This is the main benchmark for comparing HSeeker and Triplex sequence-level
classification.

## Secondary Validation

The secondary/context validation injects the same 54 sequences into the
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
| direct primary | HSeeker direct hit | 45 | 0 | 7 | 2 | 1.000 | 0.778 | 0.978 | 0.863 |
| direct primary | HSeeker score Youden | 44 | 1 | 8 | 1 | 0.978 | 0.889 | 0.978 | 0.867 |
| direct primary | Triplex direct hit | 44 | 1 | 7 | 2 | 0.978 | 0.778 | 0.967 | 0.793 |
| injected secondary | HSeeker overlap | 45 | 0 | 7 | 2 | 1.000 | 0.778 | 0.978 | 0.863 |
| injected secondary | HSeeker score Youden | 44 | 1 | 8 | 1 | 0.978 | 0.889 | 0.978 | 0.867 |
| injected secondary | Triplex overlap | 12 | 33 | 9 | 0 | 0.267 | 1.000 | 0.421 | 0.239 |
