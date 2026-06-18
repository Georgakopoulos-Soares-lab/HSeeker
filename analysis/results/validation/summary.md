# HSeeker Validation / Sensitivity Summary

Dataset: `/Users/archit/Projects/repeat_finder/HSeeker/hdna_experimental_sequences_final.csv`
Dataset size: 45
Forming: 45
Non-forming: 0

## Default HSeeker

Score threshold: `60.0`
Purity RMQ enabled: `True`
TP=38 FN=7 TN=0 FP=0
Sensitivity=0.844, specificity=0.000, precision=1.000, F1=0.916, accuracy=0.844

## Best Grid Row

{'minrep': 6, 'maxrep': 100, 'maxspacer': 10, 'purity': 0.85, 'mismatch': 0.1, 'TP': 45, 'FN': 0, 'TN': 0, 'FP': 0, 'sensitivity': 1.0, 'specificity': 0.0, 'precision': 1.0, 'F1': 1.0, 'accuracy': 1.0, 'n': 45}

## Triplex Presence/Absence

TP=44 FN=1 TN=0 FP=0
Sensitivity=0.978, specificity=0.000, precision=1.000, F1=0.989, accuracy=0.978

Triplex is treated as a binary presence/absence caller here; its score is not directly comparable to HSeeker's thermodynamic score.

## HSeeker Failure Cases

- HDNA0022 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0023 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0024 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0025 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0027 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0030 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
- HDNA0034 (forming): score=0.00, arm=0, spacer=0, GA=0.0, CT=0.0, mirror=0.0
