# HSeeker

**HSeeker** is a fast, cross-platform Python package for detecting **H-DNA (intramolecular triplex DNA)** and **mirror repeat** structures in genomic sequences. It ships a compiled C extension (`_hdna.c`) for core scanning performance.

[![PyPI version](https://img.shields.io/pypi/v/hseeker)](https://pypi.org/project/hseeker/)
[![Python](https://img.shields.io/pypi/pyversions/hseeker)](https://pypi.org/project/hseeker/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![CI — Build Wheels](https://github.com/Georgakopoulos-Soares-lab/HDNAhunter/actions/workflows/build_wheels.yml/badge.svg)](https://github.com/Georgakopoulos-Soares-lab/HDNAhunter/actions/workflows/build_wheels.yml)

**Package:** [`hseeker` on PyPI](https://pypi.org/project/hseeker/)

---

**HSeeker Web App:** https://hseeker-production.up.railway.app/
 
## Table of Contents

1. [Biological Background](#1-biological-background)
2. [Algorithm Overview](#2-algorithm-overview)
3. [Installation](#3-installation)
4. [Quick Start](#4-quick-start)
5. [Python API](#5-python-api)
6. [Command-Line Interface](#6-command-line-interface)
7. [Output Format](#7-output-format)
8. [Parameter Reference](#8-parameter-reference)
9. [Understanding the Results](#9-understanding-the-results)
10. [Performance Notes](#10-performance-notes)
11. [Development & Testing](#11-development--testing)
12. [Benchmarks](#12-benchmarks)
13. [Citation](#13-citation)
14. [License](#14-license)

---

## 1. Biological Background

**H-DNA** (also called intramolecular triplex DNA) forms when a single DNA strand folds back onto the Watson–Crick duplex and inserts into the major groove via **Hoogsteen** or **reverse-Hoogsteen** pairing. The result is a triple-helical structure over a local region of the duplex.

```
  5'─────[left arm]─────[spacer]─────[right arm]─────3'
                    ╲       (hinge)      ╱
                     ╲__triplex core___╱
```

### Structural requirements

| Requirement | Rationale |
|---|---|
| **Mirror repeat geometry** | Both arms must be approximate reverse complements of each other (reading the sequence in the same 5′→3′ direction, the two halves mirror each other). |
| **Compositional asymmetry (purity)** | Stable H-DNA requires a dominantly purine (GA-rich) or dominantly pyrimidine (CT-rich) arm to enable Hoogsteen pairing. |
| **Short spacer / hinge** | Loops >10 bp introduce thermodynamic strain that precludes stable folding. Typical genomic H-DNA has spacers of 0–7 bp. |
| **Arm length ≥ 6 bp** | Shorter tracts lack the binding energy for stable triplex formation. |

### Biological significance

H-DNA is found throughout eukaryotic and prokaryotic genomes and has been associated with:
- **Transcription regulation** — GA-rich mirror repeats in promoters modulate RNA Pol II stalling.
- **Genomic instability** — H-DNA structures are linked to DNA breakage and large-scale rearrangements.
- **Epigenetic features** — Triplex-forming regions are enriched at nucleosome-depleted regulatory elements.
- **Human disease** — Polypurine tracts implicated in cancer, neurodegenerative diseases, and trinucleotide expansion disorders.

---

## 2. Algorithm Overview

HSeeker implements a **center-outward biologically informed heuristic** that supports imperfect mirrors and compositional tolerance, enabling detection of biologically realistic H-DNA candidates.

### How it works

1. **Candidate generation** — For every position `ctr` (potential center of the hinge) and every spacer length `sp` from 0 to `maxspacer`:
   - A left pointer starts at `ctr` and moves left (toward 5′).
   - A right pointer starts at `ctr + sp + 1` and moves right (toward 3′).
2. **Outward extension** — At each step `k`, the base pair `(dna[left], dna[right])` is compared:
   - If equal: mirror match.
   - If different: mismatch counter is incremented.
   - Arm composition (fraction GA or CT) is tracked continuously.
3. **Threshold evaluation** — After every base pair beyond `minrep`, the algorithm evaluates:

$$\text{mirror identity} = 1 - \frac{\text{mismatches}}{k}$$

$$\text{purity} = \max\!\left(\frac{\text{GA count}}{k},\ \frac{\text{CT count}}{k}\right)$$

If both $\text{mirror identity} \geq 1 - \texttt{mismatch}$ and $\text{purity} \geq \texttt{purity}$, the current arm length is recorded as a valid candidate.

4. **Longest-arm retention** — For each `(ctr, sp)` pair, only the longest valid arm is kept.
5. **Inward extension** — Matching bases at both ends of a spacer move into the arms while purity and `maxrep` still permit it. Duplicate representations of the same arm and spacer are removed. A nonmatching pair remains in the spacer; the mirror mismatch rule is unchanged.
6. **Overlap removal** — After inward extension, overlapping hits are resolved by keeping the **longest arm first**; ties are broken by **shorter spacer**.

### `is_perfect` classification

A hit is classified as `is_perfect = True` only when:
- The arm is **100% pure** (either all GA or all CT), **and**
- The mirror identity is **100%** (zero mismatches).

---

## 3. Installation

### From PyPI (recommended)

```bash
pip install hseeker
```

Binary wheels are provided for **CPython 3.9, 3.10, 3.11, 3.12, and 3.13** on Linux (x86_64, aarch64), macOS (x86_64, arm64, universal2), and Windows (AMD64). No compiler is needed.

### From source

```bash
git clone https://github.com/Georgakopoulos-Soares-lab/HDNAhunter.git
cd HDNAhunter
pip install .
```

A C compiler (GCC, Clang, or MSVC) is required when building from source. The C extension compiles automatically during `pip install`.

#### Optional: CPU-native optimization

By default the extension is compiled **portably**, so the result can be copied to
or shared with any machine of the same platform. To tune the binary to the CPU of
the machine you are building on, opt in explicitly:

```bash
HSEEKER_NATIVE_BUILD=1 pip install .
```

This adds `-march=native`. The resulting binary is **not portable** — running it on
a host lacking the same instruction-set extensions can abort with an illegal
instruction — so use it only for a build that stays on the machine that produced
it, never for a wheel you intend to distribute. The build fails with an error if
`ARCHFLAGS` targets a different macOS architecture.

### Development install

```bash
pip install -e ".[dev]"
```

This installs pytest, build, twine, and cibuildwheel in addition to the package itself.

### Dependencies

| Package | Version | Purpose |
|---|---|---|
| Python | 3.9–3.13 | Core runtime |
| pandas | ≥ 1.5 | DataFrame output helper |

---

## 4. Quick Start

### In Python

```python
import hseeker

# Scan a raw sequence string
hits = hseeker.scan_sequence(
    "GGGAAAGGGTTTTCCCAAACCC",
    minrep=10,
    maxspacer=10,
    purity=0.90,
    mismatch=0.10,
)
for h in hits:
    print(h["start"], h["end"], h["arm_length"], h["is_perfect"], h["total_score"])

# Scan an entire FASTA file (all hits collected into a list)
hits = hseeker.scan_fasta("genome.fa", minrep=10, purity=0.90)

# Memory-efficient streaming alternative — yields one hit at a time
for hit in hseeker.scan_fasta_iter("genome.fa", minrep=10, purity=0.90):
    print(hit["seq_id"], hit["start"], hit["end"])

# Multi-core parallel scan — all chromosomes scanned concurrently
hits = hseeker.scan_fasta_parallel("genome.fa", minrep=10, purity=0.90)
hits = hseeker.scan_fasta_parallel("genome.fa", minrep=10, workers=4)
```

### From the command line

```bash
# Minimal — scan test.fa and write test_HDNA.tsv
hseeker -seq test.fa -out test

# Strict mode — exact mirror, 100 % pure
hseeker -seq genome.fa -out strict -purity 1.0 -mismatch 0.0

# Relaxed mode — allow up to 10 % mismatch, 90 % purity
hseeker -seq genome.fa -out relaxed -purity 0.90 -mismatch 0.10

# Verbose output, longer arms only
hseeker -seq genome.fa -out long -minrep 10 -maxrep 1000 -v
```

### Web interface

There is a single web front-end: the **FastAPI** application in `webapp/main.py`,
which serves the online tool and is what the Docker image runs.

```bash
uvicorn webapp.main:app --host 0.0.0.0 --port 8000
```

It is not required to use the library or the CLI.

---

## 5. Python API

### `hseeker.scan_sequence`

```python
hseeker.scan_sequence(
    seq: str,
    *,
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    remove_overlaps: bool = True,
    overlap_strategy: str = "greedy",
    seq_offset: int = 1,
    score: bool = True,
) -> list[dict]
```

Scan a single DNA string for H-DNA mirror repeat motifs. The C core runs with the Python GIL released, so other threads remain unblocked during the scan.

**Parameters** — see [Parameter Reference](#8-parameter-reference).

**Returns** a list of hit dictionaries. Each dictionary contains all fields described in [Output Format](#7-output-format).

---

### `hseeker.scan_fasta`

```python
hseeker.scan_fasta(
    path: str | Path,
    *,
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    remove_overlaps: bool = True,
    overlap_strategy: str = "greedy",
    score: bool = True,
) -> list[dict]
```

Scan every record in a FASTA file. Adds a `seq_id` key to each hit dict. Genomic offsets are parsed automatically from UCSC/Ensembl-style headers (`>seqid:start-end`); other headers default to offset 1.

---

### `hseeker.scan_fasta_iter`

```python
hseeker.scan_fasta_iter(
    path: str | Path,
    *,
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    remove_overlaps: bool = True,
    overlap_strategy: str = "greedy",
    score: bool = True,
) -> Generator[dict, None, None]
```

Streaming alternative to `scan_fasta`. Yields one hit dict at a time, keeping only a single record's worth of hits in memory at any moment. Use this for large FASTA files (e.g. whole-genome scans).

---

### `hseeker.scan_fasta_parallel`

```python
hseeker.scan_fasta_parallel(
    path: str | Path,
    *,
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    remove_overlaps: bool = True,
    overlap_strategy: str = "greedy",
    workers: int | None = None,
    chunk_size: int = 1_000_000,
    score: bool = True,
) -> list[dict]
```

Scans all FASTA records in parallel using a thread pool.  Because the C scan releases
the GIL, threads achieve true CPU parallelism on the heavy computation.

Records longer than *chunk_size* are split into overlapping chunks and scanned in
parallel.  With the default `chunk_size=1_000_000`, even single-chromosome genomes
like *E. coli* (4.6 Mb) are split into multiple chunks for thread-level parallelism.
On a 5-core machine this achieves **~4.6× speedup** over the sequential scan on
a single-chromosome genome.

`workers` defaults to ``os.cpu_count()``.  Pass an explicit integer to cap thread
count (e.g. ``workers=4``).  Results are returned in the same record order as the
input FASTA file.

---

## 6. Command-Line Interface

```
hseeker -seq <FASTA> -out <PREFIX> [options]
```

or equivalently:

```
python -m hseeker -seq <FASTA> -out <PREFIX> [options]
```

### Positional / required arguments

| Argument | Description |
|---|---|
| `-seq FASTA` | Path to the input FASTA file. May contain one or multiple records. Both single-line and wrapped FASTA are supported. |
| `-out PREFIX` | Output file prefix. The tool writes `<PREFIX>_HDNA.tsv`. |

### All optional arguments

| Argument | Type | Default | Description |
|---|---|---|---|
| `-minrep INT` | int | `10` | **Minimum arm length** in base pairs. Arms shorter than this are not reported, even if they satisfy all other criteria. Increasing this value reduces noise and focuses on longer, likely more stable motifs. |
| `-maxrep INT` | int | `1000` | **Maximum arm length** cap. The algorithm stops extending an arm beyond this length. Set higher for repetitive regions; setting it very high increases runtime without improving signal in most genomes. |
| `-maxspacer INT` | int | `10` | **Maximum spacer / hinge loop length** in base pairs. The spacer is the single-stranded loop between the two mirror arms. H-DNA with spacers > 10 bp is thermodynamically unfavoured. |
| `-purity FLOAT` | float | `0.90` | **Minimum compositional purity** of each arm (0.0–1.0). The arm must be ≥ `purity` fraction GA (purine) **or** ≥ `purity` fraction CT (pyrimidine). Set to 1.0 for perfect purine or pyrimidine tracts. |
| `-mismatch FLOAT` | float | `0.10` | **Maximum mirror mismatch fraction** (0.0–1.0). Fraction of base positions in the arm where `dna[left] ≠ dna[right]` (i.e. the mirror is broken). `0.0` requires perfect mirror symmetry; `0.10` allows up to 10 % divergence. |
| `-skipoverlap` | flag | *(off)* | **Skip overlap removal**. Keep distinct overlapping candidates after inward extension and removal of duplicate representations. |
| `-overlap-strategy greedy\|score` | choice | `greedy` | Choose the longest arm before scoring (`greedy`, default) or score all candidates and retain the highest-scoring nonoverlapping hits (`score`). Score mode requires scoring to be enabled and is slower on genome-scale input (see Overlap removal behaviour). |
| `-score` | flag | *(on)* | **Apply thermodynamic stability scoring**. Enables stacking and pairing score computation for each hit. Pass `-no-score` to disable and keep only the core detection columns. |
| `-at-threshold FLOAT` | float | `0.80` | **Left-arm AT-content filter**, applied before scoring. Hits whose left arm has AT content ≥ this value are dropped — AT-rich arms are dominated by weak A·A stacking and are unlikely to fold into stable H-DNA. Raise it (e.g. `1.5`) to effectively disable it. |
| `-workers INT` | int | *(all cores)* | **Parallel worker threads**. The CLI uses `scan_fasta_parallel` internally and defaults to all available CPU cores. Set this to a lower value to cap CPU usage. |
| `-v` | flag | *(off)* | **Verbose mode**. Prints progress information to `stderr`. Useful for monitoring large jobs. |

### CLI examples

```bash
# 1. Basic scan
hseeker -seq genome.fa -out results

# 2. Strict mode — exact mirror, 100 % pure
hseeker -seq genome.fa -out strict -purity 1.0 -mismatch 0.0

# 3. Whole-genome scan with minimum arm length 10 for high confidence
hseeker -seq hg38.fa -out hg38_hdna -minrep 10 -purity 0.85 -v

# 4. Exploratory scan — very permissive, keep overlapping candidates
hseeker -seq region.fa -out explore -purity 0.80 -mismatch 0.20 \
            -minrep 8 -maxspacer 10 -skipoverlap

# 5. Long perfect H-DNA only
hseeker -seq genome.fa -out perfect -purity 1.0 -mismatch 0.0 \
            -minrep 12 -maxspacer 3

# 6. Verbose output to monitor progress on a large genome
hseeker -seq hg38.fa -out hg38 -minrep 8 -v 2>progress.log

# 7. Use 8 worker threads instead of all cores
hseeker -seq hg38.fa -out hg38 -workers 8
```

---

## 7. Output Format

The CLI writes a **tab-separated file** (`<prefix>_HDNA.tsv`). `scan_fasta()` returns the same data as a list of dictionaries. Column descriptions:

| Column | Type | Description |
|---|---|---|
| `seq_id` | str | FASTA record identifier |
| `source` | str | Always `findHDNA` (for GFF3 / database compatibility) |
| `start` | int | **1-based, inclusive** genomic start of the left arm |
| `end` | int | **1-based, inclusive** genomic end of the right arm |
| `arm_length` | int | Length of each arm in bp (both arms are always equal length) |
| `spacer_length` | int | Length of the hinge loop between the two arms in bp |
| `total_length` | int | `arm_length × 2 + spacer_length` — full motif span |
| `ga_pct` | float | Percentage of G+A (purine) bases in the right arm (0–100) |
| `ct_pct` | float | Percentage of C+T (pyrimidine) bases in the right arm (0–100) |
| `mirror_identity` | float | Fraction of mirror positions that match × 100 (0–100) |
| `is_perfect` | bool | `True` if arm is 100 % pure **and** mirror identity is 100 % |
| `left_arm` | str | Sequence of the left arm (5′→3′) |
| `spacer` | str | Sequence of the hinge loop |
| `right_arm` | str | Sequence of the right arm (5′→3′) |
| `full_sequence` | str | Complete motif sequence: `left_arm + spacer + right_arm` |
| `stacking_score` | float | Stacking energy score from the triplex stability model (`None` if scoring was disabled) |
| `pairing_score` | float | Pairing energy score from the triplex stability model (`None` if scoring was disabled) |
| `total_score` | float | Combined stacking + pairing score (`None` if scoring was disabled) |
| `putative_triplex` | str | Adjusted triplex sequence in `left_arm[spacer]right_arm` notation after score-based boundary optimisation |

### Example output

```
seq_id  source    start  end  arm_length  spacer_length  total_length  ga_pct   ct_pct   mirror_identity  is_perfect  left_arm  spacer  right_arm  full_sequence  stacking_score  pairing_score  total_score putative_triplex
chr1    findHDNA  1001   1020  7           6              20            85.71    14.29    100.0            False       gggaaat   tttttt  taaaggg    gggaaattttttttaaaggg  35.0            42.3            77.3      gggaaa[tttttt]taaaggg
```

---

## 8. Parameter Reference

| Parameter | API name | CLI flag | Type | Default | Valid range | Notes |
|---|---|---|---|---|---|---|
| Minimum arm length | `minrep` | `-minrep` | int | 10 | ≥ 1 | Arms shorter than this are discarded entirely |
| Maximum arm length | `maxrep` | `-maxrep` | int | 1000 | ≥ `minrep` | Hard cap on extension; rarely needs changing |
| Maximum spacer | `maxspacer` | `-maxspacer` | int | 10 | ≥ 0 | Set to 0 for zero-loop (adjacent arms) only |
| Purity threshold | `purity` | `-purity` | float | 0.90 | 0.0–1.0 | Fraction GA or CT required in arm |
| Mismatch tolerance | `mismatch` | `-mismatch` | float | 0.10 | 0.0–1.0 | Fraction of arm positions allowed to mismatch the mirror |
| Overlap removal | `remove_overlaps` | `-skipoverlap` (inverts) | bool | True | — | When True, filters overlapping hits using `overlap_strategy` |
| Overlap strategy | `overlap_strategy` | `-overlap-strategy` | str | `greedy` | `greedy`, `score` | Greedy keeps the original arm-length rule; score mode filters after scoring all candidates |
| Sequence offset | `seq_offset` | *(automatic in CLI)* | int | 1 | ≥ 1 | 1-based start coordinate of `seq[0]`; used for correct genomic coordinates |
| Worker threads | *(CLI only)* | `-workers` | int | all cores | ≥ 1 | Number of parallel threads used by the CLI (via `scan_fasta_parallel`) |
| Stability scoring | `score` | `-score` / `-no-score` | bool | True | — | When True, computes stacking and pairing scores and an optimised triplex sequence for each hit. Scoring adds `stacking_score`, `pairing_score`, `total_score`, and `putative_triplex` columns. |

### Choosing parameters for your use case

**Genome-wide survey (high confidence)**
```
minrep=10, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10
```

**Exploratory / maximum sensitivity**
```
minrep=8, maxrep=1000, maxspacer=10, purity=0.80, mismatch=0.20
```

**Downstream ML or statistical analysis (overlapping candidates)**
```
remove_overlaps=False, skipoverlap (CLI)
```

---

## 9. Understanding the Results

### Interpreting `ga_pct` and `ct_pct`

These two values always sum to ≤ 100 %. The difference is due to ambiguous or masked bases:
- **GA-rich arm** (`ga_pct ≈ 80–100`): purine-rich tract → **H-r triplex** (Pu·Pu·Py)
- **CT-rich arm** (`ct_pct ≈ 80–100`): pyrimidine-rich tract → **H-y triplex** (Py·Pu·Py)
- **Mixed** (`neither > purity`): filtered out unless `purity` is set very low

### `mirror_identity` vs purity

`mirror_identity` measures how closely the two arms are mirror images of each other. A value of 100 means every base on the left arm is identical to its mirror partner on the right arm. A value of 90 means 10 % of positions have a mismatch. `purity`, by contrast, measures the fraction of the arm that is GA-rich or CT-rich, regardless of the mirror match.

### `is_perfect` flag

`is_perfect = True` is a stringent classification reserved for hits that are both:
1. **100 % pure** — every base in the arm is GA, or every base is CT.
2. **100 % mirror-symmetric** — no mismatches between the two arms.

These are the most likely candidates for stable H-DNA structures.

### Scoring results (`stacking_score`, `pairing_score`, `total_score`)

When scoring is enabled (default), each hit is passed through a thermodynamic stability model that evaluates the energetic favourability of triplex formation. The model assigns a **pairing score** (based on hydrogen bonding interactions) and a **stacking score** (base-stacking interactions), summed to give `total_score`. Higher scores indicate stronger predicted stability.

The scoring pass also produces a `putative_triplex` field — the motif sequence in `left_arm[spacer]right_arm` notation — whose arm/spacer boundaries may differ from the original detection because scoring re-optimises the boundaries within a local window to maximise predicted stability.

Scoring can be disabled at the API level (`score=False`) or via the CLI flag `-no-score`, which removes the four scoring columns from the TSV output.

### AT-content and homopolymer filtering

Two post-detection filters remove hits that are unlikely to form real H-DNA:

- **AT-content filter** (before scoring, `at_threshold` / `-at-threshold`, default `0.80`): drops hits whose full detected motif — both arms and the spacer — is ≥ 80% A/T.
- **Homopolymer filter** (`filter_homopolymers`, always on in the CLI): drops hits whose detected `full_sequence` is a single repeated base (poly-A, poly-G, ...). With scoring on, it also drops hits whose `putative_triplex` — arms plus spacer, brackets stripped — collapses to a single base. These carry no real mirror-repeat structure.

Both are opt-in at the Python API level (`hseeker.scan_fasta(...)`, off by default) so existing integrations are unaffected; the CLI enables the AT filter by default at `0.80` and always applies the homopolymer filter, with or without scoring.

### Overlap removal behaviour

Inward extension runs first, even when overlap removal is disabled. This keeps only one copy of an identical arm/spacer representation while allowing distinct overlapping candidates to remain.

When `remove_overlaps=True` (default), two hits overlap if their 1-based, inclusive genomic ranges share at least one base. The default `overlap_strategy="greedy"` keeps the longer arm, breaking ties by shorter spacer, before scoring. Set `overlap_strategy="score"` (CLI: `-overlap-strategy score`) to keep, instead, the highest-scoring candidate of each overlapping group: candidates are visited in descending `total_score` order and a hit is accepted when it does not overlap a previously accepted hit (ties favor longer arms, shorter spacers, then earlier coordinates); this requires scoring. `remove_overlaps=False` (CLI: `-skipoverlap`) retains every distinct candidate (after inward extension and removal of duplicate representations) regardless of strategy. Score mode is slower on genome-scale input because it has to consider every raw candidate rather than only the overlap-free set: on hg38 chromosome 1 with 16 workers it took about 60 s versus 26 s for greedy (AMD EPYC 7763), and it reports 1.5% more loci, all greedy loci being retained.

---

## 10. Performance Notes

- The C extension releases the Python GIL during scanning, enabling multi-threaded use.
- `scan_fasta_parallel(path, workers=N)` with default `chunk_size=1_000_000` achieves
  **~N× speedup on N-core machines** even for single-chromosome genomes like *E. coli*
  (4.6 Mb), which are automatically split into overlapping chunks.
- The **CLI uses `scan_fasta_parallel` by default**, automatically parallelising
  across all available CPU cores. A single human chromosome (≈ 250 Mb, chr1) processes
  in ~78 seconds on 16 cores. Use `-workers N` to limit parallelism.
- The internal hit buffer is capped at **1,000,000 hits per `scan_sequence` call**.
  For whole-genome scans, use `scan_fasta_parallel` (which splits each contig into
  overlapping chunks) rather than concatenating chromosomes into a single string.
- Memory usage is approximately **O(n_hits × 200 bytes)** plus the sequence string itself.
- N bases in the sequence (`N`, `n`) break arm extension, as per the original algorithm.

---

## 11. Development & Testing

### Running the test suite

```bash
# Install dev dependencies first
pip install -e ".[dev]"

# Run all 134 tests
pytest -v tests/
```

The test suite (`tests/test_hdna.py`) covers:
- Empty / trivial inputs
- Known GA and CT mirror motifs with exact expected values
- Strict mode (purity=1.0, mismatch=0.0) and relaxed mode
- Coordinate offset propagation (`seq_offset`)
- Multi-record FASTA parsing, `scan_fasta_iter` streaming generator, and `scan_fasta_parallel`
- Parallel chunking edge cases: boundary hits, tiny chunks, N-bases near boundaries, determinism
- Overlap removal correctness
- Parameter validation (out-of-range inputs)
- Exact field values for reference sequences
- Purity and maxspacer boundary conditions
- Case insensitivity (mixed-case input)
- Output field self-consistency (`start + total_length - 1 == end`)
- Flanking N base handling
- `is_perfect` flag logic
- Mirror identity with known mismatch counts
- Determinism (identical calls produce identical results)
- CLI end-to-end integration

### Project structure

```
hseeker/
├── src/
│   └── hseeker/
│       ├── _hdna.c          # C extension — core algorithm
│       ├── __init__.py      # Python API (scan_sequence, scan_fasta, scan_fasta_iter, scan_fasta_parallel)
│       └── __main__.py      # CLI entry point (hseeker / python -m hseeker)
├── tests/
│   └── test_hdna.py         # 134 comprehensive tests (pytest)
├── .github/
│   └── workflows/
│       └── build_wheels.yml # CI: build binary wheels + sdist, publish to PyPI
├── pyproject.toml           # Build config (setuptools, cibuildwheel, pytest)
├── setup.py                 # C extension definition
├── MANIFEST.in              # sdist file inclusion rules
└── README.md
```

### Building wheels locally

```bash
# Build the current platform wheel
pip install build
python -m build --wheel

# Build all platform wheels (requires Docker on Linux/Windows)
pip install cibuildwheel
cibuildwheel --platform linux
```

### Adding a new test

Tests live in `tests/test_hdna.py`. Each section focuses on one concern. Add new sequences as module-level constants and group related assertions under a descriptively named function. Run `pytest -v tests/` to validate.

---

## 12. Benchmarks

`benchmarks/benchmark.py` is a CLI-only end-to-end benchmark that downloads real
hg38/GRCh38 chromosomes from NCBI and measures wall time, peak RSS, and TSV
throughput for the full `hseeker` pipeline (FASTA read → scan → overlap removal
→ scoring → TSV write).

### Setup

```bash
pip install -e ".[dev]"
```

### Running

```bash
# All 24 chromosomes (~3 GB download on first run)
python benchmarks/benchmark.py

# Single chromosome quick test
python benchmarks/benchmark.py --chromosomes chr1

# Save machine-readable results to JSON
python benchmarks/benchmark.py --chromosomes chr1 --json results.json

# Write a Markdown report
python benchmarks/benchmark.py --chromosomes chr1 --report report.md

# Force re-download
python benchmarks/benchmark.py --no-cache
```

Chromosome FASTA files are cached in `benchmarks/data/` (gitignored).
When chr1, chr2, and chr3 are all present, the benchmark automatically
builds a concatenated `chr1_2_3.fa` for multi-record throughput testing.

---

## 13. Citation

If you use HSeeker in published research, please cite:

> Georgakopoulos-Soares Lab, UT Austin.  
> **HSeeker**: fast, cross-platform detection of H-DNA and mirror repeat structures in genomic sequences.  
> https://github.com/Georgakopoulos-Soares-lab/HDNAhunter

---

## 14. License

MIT License — see [LICENSE](LICENSE) for full text.

Copyright © 2024 Georgakopoulos-Soares Lab, UT Austin.
