# HSeeker Validation / Sensitivity Summary

Dataset: `/Users/archit/Projects/repeat_finder/HSeeker/hdna_experimental_sequences_final.csv`
Dataset size: 54
Forming: 45
Non-forming: 9
Score threshold: `60.0`
Purity RMQ enabled: `True`

## Method Comparison

| method | parameters | TP | FN | TN | FP | sensitivity | specificity | precision | F1 | accuracy |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HSeeker default | minrep=10, maxrep=100, maxspacer=10, purity=0.90, mismatch=0.10 | 38 | 7 | 8 | 1 | 0.844 | 0.889 | 0.974 | 0.905 | 0.852 |
| HSeeker tuned | minrep=8, maxrep=100, maxspacer=10, purity=0.90, mismatch=0.10 | 45 | 0 | 7 | 2 | 1.000 | 0.778 | 0.957 | 0.978 | 0.963 |
| Triplex default | min_score=15, p_value=0.05, min_len=6, max_len=25, min_loop=3, max_loop=10 | 44 | 1 | 7 | 2 | 0.978 | 0.778 | 0.957 | 0.967 | 0.944 |

The tuned HSeeker row uses the same parameters for every sequence: only `minrep` is changed from 10 to 8.

Tuned HSeeker and Triplex agreement: 53/54 (0.981).

Triplex is treated as a binary presence/absence caller here; its score is not directly comparable to HSeeker's thermodynamic score.

## Runtime Comparison

HSeeker default total runtime: 0.000895 sec; mean per sequence: 0.000017 sec.
HSeeker tuned total runtime: 0.000866 sec; mean per sequence: 0.000016 sec.
Triplex total runtime: 0.134000 sec; mean per sequence: 0.002481 sec.
Triplex/tuned-HSeeker runtime ratio on this validation set: 154.79x.
These are per-sequence search timings on very short sequences; HSeeker timings include scoring, while Triplex timings exclude R package startup.

## Default HSeeker Failure Cases

- HDNA0022 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0023 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0024 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0025 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0027 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0030 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0034 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0053 (non-forming): score=83.30, arm=10, spacer=0, GA=100.0, CT=0.0, mirror=100.0

## Tuned HSeeker Failure Cases

- HDNA0026 (non-forming): score=92.03, arm=8, spacer=7, GA=100.0, CT=0.0, mirror=100.0
- HDNA0053 (non-forming): score=83.30, arm=10, spacer=0, GA=100.0, CT=0.0, mirror=100.0

## Tuned HSeeker vs Triplex Disagreements

- HDNA0042 (forming): tuned HSeeker=True score=86.50; Triplex=False count=0 best_score=NA
