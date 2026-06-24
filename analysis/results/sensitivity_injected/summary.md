# Injected E. coli H-DNA Sensitivity Analysis

## Dataset

- Experimental sequences: 54 (45 forming, 9 non-forming).
- Reference genome: NC_000913.3; length 4,641,652 bp.
- Injected genome length: 4,643,455 bp.
- Insertion seed: 42; minimum inter-insertion distance: 500 bp.
- Insertions used sequences exactly as provided in `hdna_experimental_sequences_final.csv`.

## Methods

- HSeeker command: `hseeker -seq hdna_benchmark_injected.fna -minrep 8 -maxspacer 10 -mismatch 0.1 -purity 0.9 -maxrep 1000 -out hdna_benchmark`.
- This run used the exact purity RMQ prefilter (`-purity-rmq`) to preserve HSeeker output while reducing runtime.
- Triplex was run on the same injected FASTA using `triplex.search(DNAString(seq))` with the installed package defaults.
- A detected hit matched an insertion if overlap was at least 80% of either hit length or insertion length.

## Runtime

- HSeeker runtime: 0.404 sec; hits: 139.
- Triplex runtime: 8.852 sec; hits: 24.

## Classification Metrics

| method | TP | FN | TN | FP | sensitivity | specificity | precision | F1 | accuracy | MCC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HSeeker overlap | 45 | 0 | 7 | 2 | 1.000 | 0.778 | 0.957 | 0.978 | 0.963 | 0.863 |
| HSeeker score Youden | 44 | 1 | 8 | 1 | 0.978 | 0.889 | 0.978 | 0.978 | 0.963 | 0.867 |
| Triplex overlap | 12 | 33 | 9 | 0 | 0.267 | 1.000 | 1.000 | 0.421 | 0.389 | 0.239 |

## HSeeker Score ROC

- AUC-ROC: 0.980.
- Youden-optimal threshold: 86.495.
- Youden J: 0.867; MCC at this threshold: 0.867.

## Output Files

- `data/NC_000913.3.fna`: downloaded E. coli K-12 MG1655 reference.
- `data/hdna_benchmark_injected.fna`: injected benchmark genome.
- `insertion_manifest.csv`: inserted sequence coordinates.
- `hseeker_insert_matches.csv`: best HSeeker match per insertion.
- `triplex_insert_matches.csv`: best Triplex match per insertion.
- `method_comparison_by_insertion.csv`: per-insertion HSeeker/Triplex comparison.
- `hseeker_score_roc.csv`: ROC thresholds and metrics.
- `plots/hseeker_score_roc.svg`: ROC plot.
