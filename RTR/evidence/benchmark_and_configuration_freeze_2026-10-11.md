# Benchmark and configuration freeze (2026-10-11)

Author instructions, 2026-10-11:
- Hanvey 1988 (TTC)8: "do it", with one negative to keep 1:1.
- Homopolymers: "whatever is most correct".
- Citation corrections: "copy the 22 PMID and 2 DOI corrections into the citation columns. Do that."
- Configuration: "prefer the better one"; mismatch 0.10 vs 0.15 "based on data"; comparator
  settings "when the comparison time comes".

## 1. Benchmark

**New files.** `hdna_benchmark_experimental_v3.csv` sha256
`19418318b28b0e4c07ed635b5abe6b65ddf22689ab3751abd541bed35d1dfee0`;
`hdna_benchmark_balanced_v3.csv` sha256
`02cdd769c382a84f687d41dfa80d1f70c1d37f778c8f38c685facd19985ea20b`.
Balanced: **130 records, 65 forming / 65 non-forming** (70 experimental + 60 synthetic).
Experimental kept: 70 (65 forming / 5 non-forming).

### 1a. Hanvey 1988 (TTC)8 added (HDNA0081)

`Hanvey1988_(TTC)8`, `TTCTTCTTCTTCTTCTTCTTCTTC`, forming, primary, confidence high, PMID 2835375,
study `10.1016/S0021-9258(18)68654-5`, family `fam:Hanvey1988_(TTC)8`. The abstract lists
(TTC)8·(GAA)8 among the five inserts that adopt the intramolecular triplex (evidence:
`DEC-3_pRW1406.md`). Named from the verified insert, not as `pRW1406`. Its reverse complement,
(GAA)8, is contained in HDNA0041 and the four FXN (GAA)n records; this is recorded in
`relationships` (`revcomp_substring_of:…`). Family assignment follows the existing practice of one
family per repeat-length construct for these records.

### 1b. One matched negative (SYN0064)

`RTR/generate_synthetic_negatives.py` now leaves records added after the v3 pool
(`ADDED_AFTER_V3 = ('HDNA0081',)`) out of the pool computation, so the 59 earlier synthetic records are
unchanged, and draws one random GC/length-matched negative per added record from its own RNG stream
(`ADDED_SEED = 20261011`). Result: SYN0064 `random_matched_to_HDNA0081`,
`AAAGCCAAACGAATGAGAAAAGAA` (24 nt; residual arm ≥ 8 mirror: arm 9, identity 0.667; purine 0.833).

Checks: the 128 earlier records are unchanged except the citation fields of 21 records (1c);
0 exact or reverse-complement duplicates; H-DNA-competent synthetic rows `[]`.

### 1c. Citation corrections copied into the citation columns

The audit counted "22 PMID and 2 DOI" corrections. The current file holds 25 records with a corrected
PMID (21 kept, 3 removed, 1 excluded) and 9 with a corrected DOI (6 kept, 3 removed); the difference
from the audit's count is the 4 FXN records corrected on 2026-10-09 and the records removed or
excluded since. All of them were applied, to both CSVs:
- `cited_pmid` ← `corrected_pmid`;
- `cited_doi` ← `corrected_doi`, or kept when no corrected DOI exists (in every such case it already
  equals the DOI PubMed gives for the corrected PMID; checked via E-utilities `esummary`);
- `cited_paper` ← the PubMed title when the old title was a different paper (FXN rows, HDNA0055,
  HDNA0056, removed HDNA0032–0034);
- `study_id` ← the corrected DOI where `study_id` held the old DOI.

The original values are kept in four new columns, `original_cited_paper`, `original_cited_doi`,
`original_cited_pmid` and `original_study_id`, and each record's `curation_justification` notes the
change. Effects on grouping: HDNA0055 and HDNA0056 now share the study Rooney 1995
(`10.1073/pnas.92.6.2141`); the FXN records split into Potaman 2004 (HDNA0060–0062) and Sakamoto 1999
(HDNA0063). Leave-one-study-out validation has 27 study groups.

### 1d. Homopolymers: separate panel, not in the benchmark

Most correct: keep them out of the benchmark metrics and report them as a separate panel, because
(1) the composition filters decide their calls by construction and were designed with them in view;
(2) a single-base run is mirror-symmetric at every position, so it does not test mirror-repeat
detection; (3) their labels depend on conditions. Panel at the shipped configuration
(`analysis/scripts/run_homopolymer_panel.py` → `analysis/results/homopolymer_panel/`):

| Record | Sequence | Label (evidence) | no filters | AT only | homopolymer only | both (CLI) |
|---|---|---|---|---|---|---|
| HDNA0053 pRW1405 | A20 | non-forming (Hanvey 1988) | hit ✗ | no hit ✓ | no hit ✓ | no hit ✓ |
| HDNA0054 poly-(dC)30 | C30 | forming (Kohwi 1988, Mg²⁺-dependent) | hit ✓ | hit ✓ | no hit ✗ | no hit ✗ |
| SYN0053 poly-A15 | A15 | non-forming (synthetic) | no hit ✓ | no hit ✓ | no hit ✓ | no hit ✓ |
| SYN0054 poly-A30 | A30 | non-forming (synthetic) | hit ✗ | no hit ✓ | no hit ✓ | no hit ✓ |
| SYN0055 poly-T25 | T25 | non-forming (synthetic) | hit ✗ | no hit ✓ | no hit ✓ | no hit ✓ |
| SYN0056 poly-T35 | T35 | non-forming (synthetic) | hit ✗ | no hit ✓ | no hit ✓ | no hit ✓ |
| correct | | | 2/6 | 6/6 | 5/6 | 5/6 |

To disclose: the shipped homopolymer filter rejects poly-(dC)30, which forms an Mg²⁺-dependent
triplex. Poly-A15 gets no hit even without filters because 15 nt cannot hold two 8-nt arms.

## 2. Configuration

Re-run on the 130-record benchmark (`robustness_v3/`, `grouped_validation_v3/`, regenerated
2026-10-11):

| | nested CV (5×20) | leave-one-family-out | leave-one-study-out | collapsed per family |
|---|---|---|---|---|
| minrep 8, mismatch 0.10, no filters (benchmark configuration) | 0.970 | 0.970 | 0.970 | 0.862 |
| minrep 8, mismatch 0.15, both filters (manuscript) | 0.969 | 0.969 | 0.969 | 0.835 |
| minrep 10, mismatch 0.10, both filters (old CLI default) | 0.884 | 0.884 | 0.884 | 0.650 |
| per-fold tuned | 0.939 | 0.924 | 0.924 | — |

Selection frequencies inside training folds: minrep 8 in 99 % (nested), 52/52 (family), 27/27
(study); mismatch 0.10 in 52 % (nested), 48/52 (family), 24/27 (study); maxspacer 10 in 80 %, 50/52,
26/27; filters "none" everywhere (the filters change no benchmark call; "none" wins the
fewer-filters tie-break).

**Decisions.**
- **Minimum arm 8 (was 10).** Better under every scheme (0.970 vs 0.884) and selected in every
  grouped fold; on the authors' Kouzine calls it raises recall (mouse 74.1 → 79.0 %, human 68.0 →
  82.3 %) with an unchanged per-call ssDNA+ rate. Minimum arm 7 scores 0.985 on all 130 records
  but only by detecting HDNA0058, a tandem direct repeat outside the model, and minimum arm 6 adds 10
  false positives; 7 was not in the pre-specified grid. 8 is kept.
- **Mismatch 0.10.** Ties 0.15 (0.9697 vs 0.9692). The two differ on two records: 0.15 detects
  HDNA0058 (the direct repeat) and adds a false positive on HDNA0028 B-33 (experimentally
  non-forming). Grouped validation selects 0.10 in 48/52 and 24/27 folds; collapsed MCC 0.862 vs
  0.835.
- **Purity 0.90** (unchanged). 0.85 and 0.90 give identical calls on all 130 records; the folds'
  "selection" of 0.85 is the grid order breaking the tie.
- **Maximum spacer 10** (unchanged).
- **Filters on in the CLI and web app** (unchanged; AT 0.80 + homopolymer). The benchmark cannot
  test them; the Kouzine comparison supports them (the excluded A/T-rich class forms at half the
  rate in mouse and a quarter in human; without filters selectivity is lost). The Python library
  still defaults to no filters (`at_threshold=None`, `filter_homopolymers=False`), as in Nikol's
  `main`.
- **Comparator settings:** deferred to the tool-comparison run, as instructed.

**Code changed:** the `minrep` default 10 → 8 in `scan_sequence`, `scan_fasta`, `scan_fasta_iter`,
`scan_fasta_parallel`, `hseeker.ablation`, the CLI (`-minrep`), the web form and its presets, the
About page, the README and `.github/copilot-instructions.md`. The standalone C `main()` in `_hdna.c`
already defaulted to 8. Tests: 349 passed, 3 skipped (chr1 tests skip without
`benchmarks/data/chr1.fa`).

**Consequences.**
- hg38 chr1 at the new default (greedy): **106,481** loci with the library defaults (no filters) and
  **45,927** with the CLI defaults (was 96,731 / 41,939 at minimum arm 10; 96,729 published).
- The published speed figures (25.2 s, 2.3× triplex, ~5 min genome) were measured at minimum arm 10;
  they must be re-measured at 8 on the manuscript CPU.
- Benchmark results quoted anywhere at 128 records (MCC 0.969, 64:64, 73,728 grid calls) are now at
  130 records (MCC 0.970, 65:65, 74,880 grid calls).
