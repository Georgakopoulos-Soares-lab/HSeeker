# Primary Direct H-DNA Sensitivity Analysis

## Dataset

- Experimental sequences: 80 (71 forming, 9 non-forming).
- HSeeker and Triplex were run directly on each curated sequence. No genomic insertion or flanking sequence was used in this primary analysis.

## Methods

- HSeeker parameters: `minrep=8`, `maxrep=1000`, `maxspacer=10`, `purity=0.90`, `mismatch=0.10`, scoring enabled, exact `purity_rmq=True`.
- Triplex was run with installed package defaults via `triplex.search(DNAString(seq))`.
- HSeeker score discrimination was evaluated by ROC/AUC, and the classification threshold was selected by Youden's J.

## Runtime

- HSeeker direct runtime: 0.002197 sec.
- Triplex direct runtime: 0.269000 sec.

## Classification Metrics

| method | TP | FN | TN | FP | sensitivity | specificity | precision | F1 | accuracy | MCC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HSeeker direct hit | 68 | 3 | 7 | 2 | 0.958 | 0.778 | 0.971 | 0.965 | 0.938 | 0.703 |
| HSeeker score Youden | 66 | 5 | 8 | 1 | 0.930 | 0.889 | 0.985 | 0.957 | 0.925 | 0.701 |
| Triplex direct hit | 68 | 3 | 7 | 2 | 0.958 | 0.778 | 0.971 | 0.965 | 0.938 | 0.703 |

## HSeeker Score ROC

- AUC-ROC: 0.961.
- Youden-optimal threshold: 86.495.
- Youden J: 0.818; MCC at this threshold: 0.701.

## Secondary Analysis

- The E. coli injected-genome benchmark is retained separately in `analysis/results/sensitivity_injected/` as a secondary/context validation.

## Output Files

- `hseeker_direct_predictions.csv`: best HSeeker hit/score per sequence.
- `triplex_direct_predictions.csv`: best Triplex hit/score per sequence.
- `method_comparison_by_sequence.csv`: per-sequence HSeeker/Triplex comparison.
- `summary_metrics.csv`: classification metrics.
- `hseeker_score_roc.csv`: ROC thresholds and metrics.
- `plots/hseeker_score_roc.pdf`: ROC plot.
- `plots/confusion_matrices_hseeker_triplex.pdf`: HSeeker vs Triplex confusion matrices.
