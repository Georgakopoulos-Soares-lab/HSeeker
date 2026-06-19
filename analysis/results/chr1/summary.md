# chr1 hg38 HSeeker Summary

FASTA: `/Users/archit/Projects/repeat_finder/HSeeker/benchmarks/data/chr1.fa`

HSeeker output: `/Users/archit/Projects/repeat_finder/HSeeker/analysis/results/chr1/chr1_hseeker_16threads_scored.tsv`

HSeeker hits: 96,729

Run settings: 16 worker threads, scoring enabled, `purity_rmq=True`, minrep=10, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10.

HSeeker runtime: 10.37 seconds.

The hg38 chr1 FASTA was downloaded as a temporary cache for the run and removed afterward.

Score range: min=-9780.04, median=71.415, max=1785.99

Triplex status: not run for full chr1 in this analysis.

## Triplex Runtime Estimate

Full chr1 Triplex was not run. Runtime was extrapolated from Triplex runs on hg38 chr1 prefix samples using the same default Triplex settings as the validation comparison (`min_score=15`, `p_value=0.05`, `min_len=6`, `max_len=25`, `min_loop=3`, `max_loop=10`).

Sample timing table: `/Users/archit/Projects/repeat_finder/HSeeker/analysis/results/chr1/triplex_chr1_runtime_samples.csv`

Estimate JSON: `/Users/archit/Projects/repeat_finder/HSeeker/analysis/results/chr1/triplex_chr1_runtime_estimate.json`

Fit method: least_squares_through_origin_using_largest_3_samples.

Predicted full chr1 Triplex runtime: 455.2 seconds (7.6 minutes; 0.13 hours).

Relative to the measured HSeeker chr1 runtime (10.37 sec), this is about 43.9x slower.

This is an extrapolation, not a completed full-chromosome Triplex run.
