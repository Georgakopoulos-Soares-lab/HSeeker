# Genomic context recovery — E. coli insertion analysis (2026-10-15)

The 130 balanced benchmark sequences inserted into *E. coli* K-12 MG1655
(NC_000913.3) and rescanned at the current default settings.

Produced by:

    analysis/scripts/run_injected_sensitivity_analysis.py \
        --outdir analysis/results/sensitivity_injected_2026-10-15 --workers 16
    analysis/scripts/run_ecoli_overlap_sweep.py

Figure: `analysis/scripts/plot_figure3_ecoli.py`.

`data/` holds the reference genome and the modified genome. Both are omitted
from version control because they are reconstructible: the reference is the
public accession NC_000913.3, and `insertion_manifest.csv` gives the exact
coordinate of every inserted sequence, drawn under seed 42 with at least 500 bp
between successive sites. Rerunning the command above downloads the reference
and rebuilds an identical modified genome.

## Results

| | HSeeker | triplex |
|---|---|---|
| In genomic background | 63 TP, 2 FN, 65 TN, 0 FP; MCC 0.970 | 27 TP, 38 FN, 65 TN, 0 FP; MCC 0.512 |
| Calls unchanged vs scanned alone | 130/130 | 86/130 (44 detections lost) |

Overlap criterion varied 0.50–0.95: HSeeker MCC 0.954–0.970, triplex 0.512
throughout. Unmodified genome: 83 loci in 4,641,652 bp (17.9 per Mb).
