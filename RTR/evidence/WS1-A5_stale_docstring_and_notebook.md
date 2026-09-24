# WS1-A5 — stale docstring + notebook's missing input (reviewer R2.2)

**Status:** done, verified
**Date:** 2026-09-23
**Reviewer item R2.2 (verbatim):** "The results are non-reproducible, as the
sequences cannot be retrieved from the supplement. It would be much better if the
sequences would be given as a csv file… Also, it would be good if in the repo there
was a script to run the analysis on the 79 test cases."

## Part 1 — stale "54 sequences" docstring

`analysis/scripts/run_direct_sensitivity_analysis.py:2-7` described the benchmark as
"the 54 experimentally curated sequences" while the script actually ran on 80 rows.
Rewritten during the benchmark-v3 replacement to describe the real input (balanced
v3: 130 records at 65:65, with `--subset experimental` for the 71-record
composition).

Verified: no occurrence of "54 experimentally curated" or "54 sequences" remains
anywhere in the repository outside the revision archive.

## Part 2 — notebook could not run: missing `hseeker_additional_motifs.csv`

`analysis/hdna_detection_benchmark.ipynb` cell 4 did
`pd.read_csv("hseeker_additional_motifs.csv")`. That file is **absent from the
repository** and absent from every commit, so the notebook failed at cell 4 for
anyone who tried to reproduce it — exactly the reproducibility failure R2.2 reports.

**Resolved by removing the cell, not by inventing the file.** The merge is now
obsolete: benchmark v3 is the single curated source and already contains every
retained record with its provenance, curation decision and label, plus the
synthetic negatives. Cell 4 is replaced by a markdown note stating this, so a reader
sees why it disappeared rather than finding an unexplained gap.

## A defect this uncovered in my own earlier change

Repointing the notebook at v3 (done during the CSV replacement) **silently broke
cell 6**: it wrote FASTA headers from `row['first_author']` and `row['year']`, and
those two columns exist in the old CSV but **not** in v3. The notebook would have
raised `KeyError` at cell 6.

Caught by tracing column dependencies before declaring the action done, and fixed:
the header now uses `record_id, sequence_name, record_type, label` — all present in
v3 — with a comment recording that provenance moved to `cited_pmid` / `study_id`.
Cells 4 and 6 were the only two cells referencing the dropped columns.

## Verification — the data path was executed, not merely parsed

Ran the notebook's load cell and FASTA cell against the real file:

```
Wrote 130 sequences to .../hdna_experimental_sequences.fasta
records written: 130
label counts: {'forming': 65, 'non_forming': 65}

first FASTA lines:
  >HDNA0001, pAA32, experimental, forming
  AAGGGAGAAAGGGGTATAGGGGAAAGAGGGAA
  >HDNA0002, pGG32, experimental, forming
  AAGGGAGAAGGGGGTATAGGGGGAAGAGGGAA
```

Notebook validates under `nbformat` (49 cells).

**Scope limit, stated honestly:** only the data-loading and FASTA-writing cells were
executed. Later cells download the E. coli K-12 genome and run the full injection
benchmark, which is the TACC-bound work (WS2-A1) and was not run here. So this
action verifies that the notebook's *missing-input* failure is gone and its data
path is correct — not that every downstream cell now succeeds.

## Files changed

- `analysis/scripts/run_direct_sensitivity_analysis.py` (docstring, during CSV replacement)
- `analysis/hdna_detection_benchmark.ipynb` (cell 4 → markdown note; cell 6 header fields)
