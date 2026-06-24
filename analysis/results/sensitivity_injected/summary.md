# Secondary Injected E. coli H-DNA Sensitivity Analysis

## Dataset

- This injected-genome benchmark is a secondary/context validation. The primary validation is the direct-on-sequences benchmark.
- Experimental sequences: 80 (71 forming, 9 non-forming).
- Reference genome: NC_000913.3; length 4,641,652 bp.
- Injected genome length: 4,644,767 bp.
- Insertion seed: 42; minimum inter-insertion distance: 500 bp.
- Insertions used sequences exactly as provided in `hdna_experimental_sequences_final.csv`.

## Methods

- HSeeker command: `hseeker -seq hdna_benchmark_injected.fna -minrep 8 -maxspacer 10 -mismatch 0.1 -purity 0.9 -maxrep 1000 -out hdna_benchmark`.
- This run used the exact purity RMQ prefilter (`-purity-rmq`) to preserve HSeeker output while reducing runtime.
- Triplex was run on the same injected FASTA using `triplex.search(DNAString(seq))` with the installed package defaults.
- A detected hit matched an insertion if overlap was at least 80% of either hit length or insertion length.

## Runtime

- HSeeker runtime: 0.417 sec; hits: 163.
- Triplex runtime: 8.807 sec; hits: 100.

## Classification Metrics

| method | TP | FN | TN | FP | sensitivity | specificity | precision | F1 | accuracy | MCC |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| HSeeker overlap | 69 | 2 | 7 | 2 | 0.972 | 0.778 | 0.972 | 0.972 | 0.950 | 0.750 |
| HSeeker score Youden | 69 | 2 | 8 | 1 | 0.972 | 0.889 | 0.986 | 0.979 | 0.963 | 0.822 |
| Triplex overlap | 31 | 40 | 9 | 0 | 0.437 | 1.000 | 1.000 | 0.608 | 0.500 | 0.283 |

## HSeeker Score ROC

- AUC-ROC: 0.971.
- Youden-optimal threshold: 76.375.
- Youden J: 0.861; MCC at this threshold: 0.822.

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
