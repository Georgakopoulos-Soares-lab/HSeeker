# Plan: comparing HSeeker with Triplexator, NeSSie and nBMST (and triplex) — 2026-10-11

**Purpose.** Answer the reviewers' comparator comments with a fair, pre-specified, reproducible comparison
on the frozen benchmark, and fill the four places marked PENDING in
`RTR/RESPONSE_LETTER_DRAFT_2026-10-11.md`.

**Comments this answers.**
- **R2.1:** with the old 70:9 imbalance, "not worse than other models" was not convincing.
- **R2.7:** confidence intervals are needed for the tool-vs-tool comparison before claiming superiority.
- **R3.1:** CIs for sensitivity, specificity, MCC, PPV and NPV.
- **R3.6:** Triplexator/3plex target intermolecular triplexes. Add NeSSie and nBMST, or document why they
  cannot be included. Harmonise "outperformed existing tools" (Key Points) with "at least comparable"
  (Results).
- **R3.11:** speed claims, raw benchmark data deposited.
- Reviewer 1 praised the comparison, so it must stay in the paper.

**triplex (Bioconductor)** is included although the request names the other three. It is the only
comparator designed for intramolecular triplexes (H-DNA) and was in the original comparison, so dropping
it would look selective.

---

## 1. What exists today

| Tool | Accuracy on a benchmark | Speed | Build / environment on TACC |
|---|---|---|---|
| **triplex** 1.50.0 (R/Bioconductor) | June 2026, old 80-record set: **tied** with HSeeker (68/3/7/2, MCC 0.703 each). Settings: `min_score 15, min_len 8, max_len 50, seq_type "prokaryotic"` (chosen for the E. coli insertion arm) | chr1 57.7 s at 16 cores (mclapply), measured against the old HSeeker run | Conda env `r_triplex` (`/home1/10899/kimopro/WORK/miniconda3/envs/r_triplex`) loads triplex 1.50.0 — **ready** |
| **Triplexator** 1.3.2 | Run in the June notebook only for the E. coli arm (`-ds <fasta> -l 10 -L -1 -e 15 -fr off`, duplex target-site mode); no direct benchmark result | chr1 67.5 s (`-ds -rm 1 -p 16 -of 0`, defaults) | `/scratch/10899/kimopro/triplexator` is a git checkout of `Gurado/triplexator` at `4505bba` whose **files and binary were removed by the scratch purge**; the two local C++11 patches to the bundled SeqAn headers were not committed and are lost |
| **NeSSie** (`B3rse/nessie`, `dbe6cb2`) | Never run on the benchmark. Smoke test (Sept, other machine): `-M -k 8` gives **no** mirror motif in the canonical H-DNA construct pGG32, but one in poly-A20; `-T` reports purine-rich k-mers | never | Built on another machine only; GitHub reachable from TACC |
| **nBMST / non-B_gfa** (`abcsFrederick/non-B_gfa`, `a891b59`) | Never run on the benchmark. Smoke test: hits pGG32 (repeat 14, spacer 4), **misses pcMyc** (a verified former), hits poly-A20; `Score` column is `NA` | never | Built on another machine only; **crashes with a buffer overflow on long path arguments** (run with short relative paths) |
| **HSeeker** (final defaults) | 130 records: 63/2/65/0, MCC 0.970 (0.926–1.000) | chr1 16.3 s at 16 workers; whole hg38 204.5 s | ready |

---

## 2. Principles fixed before any comparator is run

1. **Same input, same labels.** All tools read the same FASTA of the frozen benchmark
   (`hdna_benchmark_balanced_v3.csv`, sha256 `c3f3b842…`; 130 records, 65:65; length 15–198 nt, median
   32). Results are reported on three sets: balanced (130), experimental (70: 65/5) and primary (58:
   53/5).
2. **One call rule for every tool.** A record is called forming if the tool reports at least one
   qualifying hit anywhere in it, on either strand. No score threshold is used for any tool (as for
   HSeeker, see R3.3). The qualifying hit type is fixed per tool in §4 before running.
3. **Three settings per tool, all reported, decided now:**
   - **A, documented defaults** (what a user gets) — *primary comparison*;
   - **B, geometry matched to HSeeker** where the tool has the parameters (arm/half-repeat ≥ 8 nt,
     spacer ≤ 10 nt, ≤ 10 % mismatch, ≥ 90 % purine/pyrimidine), to separate the algorithm from the
     parameter choice;
   - **C, tuned out of fold:** a small grid per tool, selected inside training folds exactly as for
     HSeeker (nested 5-fold × 20, leave-one-family-out, leave-one-study-out). This answers the concern
     that HSeeker's defaults were chosen on this benchmark: every tool gets the same chance.
   The settings used in the original manuscript figures are run as well, as a reference, so the old
   claims can be checked.
4. **Wording rule, decided now.** "HSeeker performs better than X" only if the 95 % CI of the paired
   difference in MCC (family-cluster bootstrap) excludes 0 in the primary comparison **and** the
   difference has the same sign under settings B and C. Otherwise "comparable". This sets the Key Points
   wording.
5. **Pinned and deposited.** Tool versions, commits, build flags, exact command lines and every raw output
   file are written to `analysis/results/tool_comparison/` and the settings into
   `analysis/benchmark_config.json` (the `comparators` block) before the first run. A design note with
   §2–§6 is committed **before** results exist, as was done for the Kouzine analysis.

---

## 3. Phase 1 — environment and builds (≈ 3 h)

All builds go to `$WORK/tools/` (not `/scratch`, which is purged) with a build script per tool in
`benchmarks/tools/`, so they can be rebuilt.

### 3.1 triplex
- Use the existing `r_triplex` env; record `sessionInfo()` (R, triplex 1.50.0, Biostrings versions).
- Smoke test on pGG32, pcMyc and a random negative.

### 3.2 Triplexator 1.3.2
1. Fresh clone of `https://github.com/Gurado/triplexator` at `4505bba` into `$WORK/tools/triplexator`
   (or `git checkout -- .` in the purged copy).
2. `module load gcc/9.4.0 cmake/4.1.1`; `cmake -DCMAKE_BUILD_TYPE=Release`, `make`.
3. Expect the two C++11 incompatibilities in the bundled SeqAn headers that the published build patched.
   Fix them minimally, save the change as `benchmarks/tools/triplexator_seqan_cxx11.patch` and apply it
   from the build script, so the build is reproducible.
4. Check: `triplexator --version` = 1.3.2; run the bundled `demos/` example and compare with
   `demos/reference/`.

### 3.3 NeSSie
1. Clone `https://github.com/B3rse/nessie` at `dbe6cb2` into `$WORK/tools/nessie`; `make`.
2. Its output parsers need Python 2.7, so write a small Python 3 parser for the raw output (`>name`,
   `$|length|motif[|alignment]`, `@counts`, `@indexes`) in `analysis/scripts/tool_parsers.py`, with unit
   tests against the README examples.
3. Check: the README examples reproduce; pGG32 / poly-A20 reproduce the Sept smoke test.

### 3.4 nBMST (non-B_gfa)
1. Clone `https://github.com/abcsFrederick/non-B_gfa` at `a891b59` into `$WORK/tools/non-B_gfa`; `make`.
2. Always run from inside the output directory with short relative paths (`./gfa -seq s.fa -out o`),
   with `-skipWGET`.
3. Check: the bundled `gfa_test.fasta` reproduces the 14 reference output files.

---

## 4. Phase 2 — settings and hit types (written down before running)

### 4.1 Triplexator
Triplexator models a third strand from a separate oligonucleotide (TFO) binding a duplex target site
(TTS); it does not model intramolecular H-DNA.
- **Hit type:** a triplex target site (TTS) in duplex mode (`-ds`), as in the manuscript.
- **A (defaults):** `-ds <fa> -of 0` (minimum length 16, error rate 5 %, guanine ≥ 10 %, low-complexity
  filter on, all motifs).
- **B (matched):** `-ds <fa> -l 16 -e 10 -of 0`. A TTS spanning both arms of an 8+8 mirror repeat is 16
  nt; the error rate is matched to HSeeker's 10 % mismatch.
- **Manuscript reference:** `-ds <fa> -l 10 -L -1 -e 15 -fr off`.
- **C (tuned):** grid `-l {10,12,16}` × `-e {5,10,15}` × `-fr {on,off}` (18 settings).
- **Exploratory, not part of the claims:** the same FASTA as both `-ss` and `-ds` (TFO–TTS pairs within
  one molecule) as the closest Triplexator gets to an intramolecular triplex. Report its calls, but
  state that the mode is not designed for H-DNA.
- Minimum length 10 is a hard limit; records shorter than the minimum length simply get no hit, which
  is reported.

### 4.2 NeSSie
NeSSie searches k-mers that are mirror-symmetric (`-M`), or purine-rich with triplex-forming potential
(`-T`); it does not model an arm–spacer–arm motif.
- **Primary hit type:** a `-T` motif (the authors' "DNA-triplex forming potential" mode).
- **A (defaults):** `-T -k 16` (the shortest motif covering two 8-nt arms; NeSSie requires `-k`), other
  parameters at their defaults.
- **B (matched):** `-T -k 16 -K 60 -MAX -p 10 -m 10`.
- **Secondary hit type, reported alongside:** a degenerate mirror motif, `-M -k 16 -K 60 -MAX -m 10 -g 10`
  (gaps allowed so that a non-symmetric spacer can be absorbed).
- **C (tuned):** `-T` with `-k {12,16,20}` × `-p {0,10,20}`, and `-M` with `-k {12,16,20}` × `-m {0,10}`
  × `-g {0,10}`.
- Report the Sept observation (no `-M` motif in pGG32, one in poly-A20) in the letter as the documented
  difference in target class, with the numbers from this run.

### 4.3 nBMST (non-B_gfa)
nBMST's mirror-repeat module matches HSeeker's geometry, and its "triplex" subset (mirror repeats with
enough purine/pyrimidine content and a short spacer) is the Non-B DB H-DNA prediction.
- **Primary hit type:** the triplex subset of mirror repeats (`*_MR.tsv` rows flagged as triplex /
  the triplex output file; confirm the file name in the Phase 1 test).
- **Secondary hit type:** any mirror repeat.
- **A (defaults):** `-skipWGET -skipAPR -skipSTR -skipDR -skipIR -skipGQ -skipZ` (only the mirror-repeat
  module; defaults `-minMRrep 10 -maxMRspacer 100 -minTriplexYRpercent 10 -maxTriplexSpacer 8`).
- **B (matched):** as A with `-minMRrep 8 -maxMRspacer 10 -maxTriplexSpacer 10`.
- **C (tuned):** `-minMRrep {8,10,12}` × `-maxMRspacer {5,10}` × `-minTriplexYRpercent {10,50,90}`.
- nBMST gives no score, so no ranking metric is computed for it.

### 4.4 triplex
- **Hit type:** any `triplex.search` hit (all triplex types 0–7).
- **A (defaults):** `triplex.search(dna)` (min_score 15, p_value 0.05, min_len 6, max_len 25,
  eukaryotic model).
- **B (matched):** `min_len = 8, max_len = 50`, other defaults.
- **Manuscript reference:** `min_score 15, min_len 8, max_len 50, seq_type "prokaryotic"`. The
  prokaryotic model was chosen for the E. coli arm, which is now removed, so it is not used as primary.
- **C (tuned):** `min_score {10,15,20,25}` × `p_value {0.01,0.05}` × `min_len {6,8}`.

### 4.5 HSeeker
- **A = B:** shipped defaults.
- **C:** the existing 576-configuration out-of-fold results (`analysis/results/robustness_v3/`,
  `grouped_validation_v3/`).

---

## 5. Phase 3 — benchmark runs and statistics (≈ 4 h)

**Script:** `analysis/scripts/run_tool_comparison.py`, extending `run_direct_sensitivity_analysis.py`.
It writes the FASTA, runs each tool × setting (subprocess, timeouts, return codes recorded), parses the
outputs into one table and computes everything below. Each tool run is cached by a hash of its command, so
reruns are cheap.

**Outputs** (`analysis/results/tool_comparison/`):
- `per_record_calls.csv`: record × tool × setting → call, number of hits, best score where one exists, hit
  coordinates;
- `summary_metrics_v2.csv`: per tool, setting and set (130 / 70 / 58) — TP, FN, TN, FP, sensitivity,
  specificity, PPV and NPV with Wilson CIs, MCC with record-bootstrap and family-cluster-bootstrap CIs;
- `paired_differences.csv`: HSeeker minus tool for sensitivity, specificity and MCC — paired bootstrap
  CIs (records and families as units) and exact McNemar tests on the discordant records;
- `by_negative_class.csv`: false-positive rate per negative class (mirror-disrupted mutants, shuffles,
  random GC/length-matched, G-quadruplex, Z-DNA, B-DNA, balanced perfect mirrors) and false-negative rate
  per positive family;
- `tuned_out_of_fold.csv`: setting C, nested / leave-one-family-out / leave-one-study-out;
- `homopolymer_panel_tools.csv`: the six homopolymers per tool (the panel reported for R2.1);
- `raw/`: every raw tool output and command line; `metadata.json` with versions, commits, sha256 of
  inputs, node and dates.

**Checks before reading any result:**
- every record appears for every tool and setting (no silently missing records);
- tools that fail on very short records (15–20 nt) are reported, not dropped;
- the June manuscript settings reproduce the June triplex calls on the 69 old records still in the
  benchmark (a consistency check of the pipeline).

**Exploratory, clearly labelled:** score-based ROC curves for tools that report a score (triplex score,
Triplexator error/score, NeSSie score), with any threshold chosen out of fold only.

---

## 6. Phase 4 — speed and memory (≈ 2 h, same node as HSeeker's re-measurement)

Same node type and protocol as `RTR/evidence/R3.11_runtime_memory_2026-10-11.md` (AMD EPYC 7763,
16 vCPUs; `/usr/bin/time -v`; three runs; hg38 chr1, sha256 `04ee6db2…`), with nothing else running:
- **triplex:** `benchmarks/triplex_benchmark.R`, 249 windows of 1 Mb with 1 kb padding,
  `mclapply(mc.cores = 16)`, setting A. Report R-internal time and total wall time.
- **Triplexator:** `benchmarks/triplexator_benchmark.py` (`-ds -rm 1 -p 16 -of 0`), setting A.
- **nBMST:** single-threaded, mirror-repeat module only; report as such (no parallel mode exists), and
  also with the chromosome split into 16 chunks run in parallel if a fair 16-core figure is wanted.
- **NeSSie:** single-threaded, `-T` setting A, with a 2-hour limit; if it does not finish, report the
  limit, not an extrapolated time.

Table: tool, version, settings, cores, wall time, peak RSS, loci on chr1. This replaces the June
comparison table in `benchmarks/benchmark_comparison.md`.

---

## 7. Phase 5 (optional, strongest independent evidence) — genome-wide comparison against Kouzine 2017

Run each tool (setting A) on mm9 and hg19 and compute the same recall and selectivity used for HSeeker
(R2.3): the fraction of the authors' formed H-DNA each tool overlaps, and the odds ratio of formation for
motifs each tool calls. This compares the tools against experimental data that none of them was tuned
on. Cost: triplex and Triplexator ≈ 10–20 min per genome at 16 cores, nBMST ≈ 30–60 min single-threaded
per genome, NeSSie possibly many hours (run per chromosome in parallel, or skip if infeasible and say
so). Reuse `run_kouzine_published.py` with a `--loci` option for any BED of tool hits.

---

## 8. Phase 6 — writing (≈ 2 h)

Fill the PENDING places in `RTR/RESPONSE_LETTER_DRAFT_2026-10-11.md`:
- **R2.1:** per-tool results on the balanced benchmark, and where each tool's false positives come from
  (`by_negative_class.csv`).
- **R2.7:** paired differences with CIs and McNemar p-values; the wording rule of §2.4.
- **R3.6:**
  - the target-class differences, with numbers: Triplexator duplex TTS; NeSSie k-mer symmetry and purine
    content; nBMST mirror-repeat geometry without a score;
  - the results under settings A, B and C;
  - Key Points wording set by §2.4.
- **R3.11:** the speed table from Phase 4.
- **Methods:** a paragraph per tool with version, commit, command lines and hit rule; **Supplementary
  table:** `summary_metrics_v2.csv`.

---

## 9. Risks and how they are handled

| Risk | Handling |
|---|---|
| Triplexator does not compile with the current toolchain | Fix the SeqAn headers minimally (as the published build did) and commit the patch; if it still fails, use the `gcc/9.4.0` module the published build used |
| Scratch purge removes tools or data again | Build in `$WORK`; keep build scripts and patches in the repository; record sha256 of inputs |
| nBMST buffer overflow on long paths | Run from inside the output directory with relative paths (documented) |
| NeSSie parsers need Python 2.7 | Own Python 3 parser with tests |
| A tool cannot be scored as an H-DNA classifier (e.g. NeSSie `-M`) | Report it with the stated hit rule *and* say in R3.6 why the mapping is approximate; that is the "documented technical limitation" R3.6 allows |
| Results favour or disfavour HSeeker | Settings, hit rules and the wording rule are fixed and committed before running; all settings reported |
| Short negatives (15–20 nt) below a tool's minimum length | Count as "no hit" (the correct call for negatives) and report how many records were below each tool's minimum |

---

## 10. Order and time

1. Commit this plan's §2–§5 as the design note, before any comparator result (15 min).
2. Phase 1 builds and smoke tests (≈ 3 h).
3. Write `run_tool_comparison.py` and the parsers, with unit tests (≈ 2 h).
4. Phase 3 runs and statistics (≈ 1 h of compute on this node).
5. Phase 4 speed runs, node otherwise idle (≈ 1 h).
6. Phase 6 writing (≈ 2 h).
7. Optional Phase 5 (≈ half a day of compute, mostly unattended).

About one working day without Phase 5. Everything runs on the interactive node; no `sbatch` is needed.

---

## 11. Changes recorded during Phase 1 (2026-10-11, before any benchmark run)

**Builds (all in `$WORK/tools`, script `benchmarks/tools/build_comparators.sh`):**
- Triplexator 1.3.2 (`4505bba`): builds with `gcc/9.4.0`, `cmake/4.1.1`
  (`-DCMAKE_POLICY_VERSION_MINIMUM=3.5`), `-std=gnu++98`, and a two-line patch
  (`benchmarks/tools/triplexator_seqan_guessFileFormat.patch`) that declares the never-instantiated
  `seqan::guessFileFormat` as returning a const reference. Its own smoke test: 14 passed, 0 failed.
- NeSSie (`dbe6cb2`) and non-B_gfa (`a891b59`): plain `make`; non-B_gfa reproduces its 14 reference
  output files exactly.

**NeSSie.** Searches over a length range (`-K`, with or without `-MAX`) stop silently after the first
input sequence, so ranges are not used. Instead each setting is the **union of fixed lengths k … k+4**
(both parities, so alternating repeats such as (GA)n, which are mirror-symmetric only at odd lengths, can
be found). `-T` reports motifs that are mirror-symmetric *and* purine- or pyrimidine-rich (e.g. (GA)n and
(CT)n at k = 17, poly-A; not pGG32 at k = 16). Settings:
- A: `-T`, k = 16–20; secondary `-M`, k = 16–20;
- B: `-T -p 10 -m 10`, k = 16–20; secondary `-M -m 10 -g 10`, k = 16–20;
- C: `-T` with k = {12, 16, 20} + 0…4 × `-p {0, 10, 20}`; `-M` with k = {12, 16, 20} + 0…4 × `-m {0, 10}` ×
  `-g {0, 10}`.

**nBMST.** The triplex subset (`Subset` = 1 in `*_MR.tsv`) is decided by `(nY / nR) * 100 <= limit`
with **integer division**, so it reduces to "more purines than pyrimidines on the + strand (or no
pyrimidines)" plus spacer `< maxTriplexSpacer`. Consequence, verified: (GA)12 is flagged as triplex but
the same duplex written as (CT)12 is not. This is the tool's behaviour and is scored as is; R3.6 reports
it. Settings:
- A: mirror-repeat module only (`-skipWGET -skipAPR -skipSTR -skipDR -skipIR -skipGQ -skipZ`), defaults,
  hit = triplex subset; secondary hit = any mirror repeat;
- B: as A with `-minMRrep 8 -maxMRspacer 10 -maxTriplexSpacer 11` (the test is "less than", so 11 allows
  spacers up to 10 as in HSeeker);
- C: `-minMRrep {8, 10, 12}` × `-maxMRspacer {5, 10}` × hit type {triplex subset, any mirror repeat}
  (`-minTriplexYRpercent` values between 0 and 99 are equivalent because of the integer division, so it is
  not varied).

**Triplexator.** At its defaults it reports no target site in (GA)12, (CT)12 or pGG32 (its
low-complexity filter and 16-nt minimum); with the manuscript settings it reports all three. Settings as
in §4.1. Output coordinates are 0-based.
