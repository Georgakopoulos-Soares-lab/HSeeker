# Secondary Injected E. coli H-DNA Sensitivity Analysis

## Dataset

- This injected-genome benchmark is a secondary/context validation. The primary validation is the direct-on-sequences benchmark.
- Experimental sequences: 130 (65 forming, 65 non-forming).
- Reference genome: NC_000913.3; length 4,641,652 bp.
- Injected genome length: 4,646,348 bp.
- Insertion seed: 42; minimum inter-insertion distance: 500 bp.
- Insertions used sequences exactly as provided in `hdna_benchmark_balanced_v3.csv`.

## Methods

- HSeeker command: `hseeker -seq hdna_benchmark_injected.fna -minrep 8 -maxspacer 10 -mismatch 0.1 -purity 0.9 -maxrep 1000 -out hdna_benchmark`.
- This run used the exact purity RMQ prefilter (`-purity-rmq`) to preserve HSeeker output while reducing runtime.
- Triplex was run on the same injected FASTA using `triplex.search(DNAString(seq))` with the installed package defaults.
- A detected hit matched an insertion if overlap was at least 80% of either hit length or insertion length.

## Runtime

- HSeeker runtime: 1.285 sec; hits: 147.
- Triplex runtime: 15.311 sec; hits: 89.

## Classification Metrics

| method | TP | FN | TN | FP | sensitivity | specificity | precision | F1 | accuracy | MCC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HSeeker overlap | 63 | 2 | 65 | 0 | 0.969 | 1.000 | 1.000 | 0.984 | 0.985 | 0.970 |
| HSeeker score Youden | 63 | 2 | 65 | 0 | 0.969 | 1.000 | 1.000 | 0.984 | 0.985 | 0.970 |
| Triplex overlap | 27 | 38 | 65 | 0 | 0.415 | 1.000 | 1.000 | 0.587 | 0.708 | 0.512 |

## HSeeker Score ROC

- AUC-ROC: 0.985.
- Youden-optimal threshold: 76.375.
- Youden J: 0.969; MCC at this threshold: 0.970.

## Output Files

- `data/NC_000913.3.fna`: downloaded E. coli K-12 MG1655 reference.
- `data/hdna_benchmark_injected.fna`: injected benchmark genome.
- `insertion_manifest.csv`: inserted sequence coordinates.
- `hseeker_insert_matches.csv`: best HSeeker match per insertion.
- `triplex_insert_matches.csv`: best Triplex match per insertion.
- `method_comparison_by_insertion.csv`: per-insertion HSeeker/Triplex comparison.
- `hseeker_score_roc.csv`: ROC thresholds and metrics.
- `plots/hseeker_score_roc.pdf`: ROC plot.
- `plots/confusion_matrices_hseeker_triplex.pdf`: HSeeker vs Triplex confusion matrices.
