# H-DNA Hunter

**H-DNA Hunter** is a fast, cross-platform Python package for detecting **H-DNA (intramolecular triplex DNA)** and **mirror repeat** structures in genomic sequences. It ships a compiled C extension (`_hdna`) as its computational core and exposes both a clean Python API and a drop-in CLI replacement for the original `findHDNA` binary.

[![PyPI version](https://img.shields.io/pypi/v/hdna-hunter)](https://pypi.org/project/hdna-hunter/)
[![Python](https://img.shields.io/pypi/pyversions/hdna-hunter)](https://pypi.org/project/hdna-hunter/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![CI — Build Wheels](https://github.com/Georgakopoulos-Soares-lab/HDNAhunter/actions/workflows/build_wheels.yml/badge.svg)](https://github.com/Georgakopoulos-Soares-lab/HDNAhunter/actions/workflows/build_wheels.yml)

---

## Table of Contents

1. [Biological Background](#1-biological-background)
2. [Algorithm Overview](#2-algorithm-overview)
3. [Installation](#3-installation)
4. [Quick Start](#4-quick-start)
5. [Python API](#5-python-api)
6. [Command-Line Interface](#6-command-line-interface)
7. [Output Format](#7-output-format)
8. [Parameter Reference](#8-parameter-reference)
9. [Streamlit Web App](#9-streamlit-web-app)
10. [Understanding the Results](#10-understanding-the-results)
11. [Performance Notes](#11-performance-notes)
12. [Development & Testing](#12-development--testing)
13. [Citation](#13-citation)
14. [License](#14-license)

---

## 1. Biological Background

**H-DNA** (also called intramolecular triplex DNA) forms when a single DNA strand folds back onto the Watson–Crick duplex and inserts into the major groove via **Hoogsteen** or **reverse-Hoogsteen** hydrogen bonds. The result is a three-stranded helical segment joined by a short single-stranded hinge loop.

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

H-DNA Hunter implements a **center-outward biologically informed heuristic** — a significant improvement over strict mirror-only approaches such as non-B_gfa, which require exact purity.

### How it works

1. **Candidate generation** — For every position `ctr` (potential center of the hinge) and every spacer length `sp` from 0 to `maxspacer`:
   - A left pointer starts at `ctr` and moves left (toward 5′).
   - A right pointer starts at `ctr + sp + 1` and moves right (toward 3′).
2. **Outward extension** — At each step `k`, the base pair `(dna[left], dna[right])` is compared:
   - If equal: mirror match.
   - If different: mismatch counter is incremented.
   - Arm composition (fraction GA or CT) is tracked continuously.
3. **Threshold evaluation** — After every base pair beyond `minrep`, the algorithm evaluates:

$$\text{mirror\_identity} = 1 - \frac{\text{mismatches}}{k}$$

$$\text{purity} = \max\left(\frac{\text{GA count}}{k},\, \frac{\text{CT count}}{k}\right)$$

If both `mirror_identity ≥ 1 − \texttt{mismatch}` and `purity ≥ \texttt{purity}`, the current arm length is recorded as a valid candidate.

4. **Longest-arm retention** — For each `(ctr, sp)` pair, only the longest valid arm is kept. This ensures no redundant shorter arms at the same position.
5. **Overlap removal** — After scanning the full sequence, overlapping hits are resolved by keeping the **longest arm first**; ties are broken by **shorter spacer**. This is the behaviour of the original non-B_gfa suite.

### `is_perfect` classification

A hit is classified as `is_perfect = True` only when:
- The arm is **100% pure** (either all GA or all CT), **and**
- The mirror identity is **100%** (zero mismatches).

---

## 3. Installation

### From PyPI (recommended)

```bash
pip install hdna-hunter
```

Binary wheels are provided for CPython 3.9–3.13 on Linux (x86_64, aarch64), macOS (x86_64, arm64, universal2), and Windows (AMD64). No compiler is needed.

### From source

```bash
git clone https://github.com/Georgakopoulos-Soares-lab/HDNAhunter.git
cd HDNAhunter
pip install .
```

A C compiler (GCC, Clang, or MSVC) is required when building from source. The C extension compiles automatically during `pip install`.

### Development install

```bash
pip install -e ".[dev]"
```

This installs pytest, build, twine, and cibuildwheel in addition to the package itself.

### Dependencies

| Package | Version | Purpose |
|---|---|---|
| Python | ≥ 3.9 | Core runtime |
| pandas | ≥ 1.5 | DataFrame output helper |
| streamlit | ≥ 1.30 | Web app (optional, `[app]` extra) |
| plotly | ≥ 5.0 | Charts in web app (optional, `[app]` extra) |

---

## 4. Quick Start

### In Python

```python
import hdna_hunter

# Scan a raw sequence string
hits = hdna_hunter.scan_sequence(
    "GGGAAAGGGTTTTCCCAAACCC",
    minrep=6,
    maxspacer=7,
    purity=0.80,
    mismatch=0.20,
)
for h in hits:
    print(h["start"], h["end"], h["arm_length"], h["is_perfect"])

# Scan an entire FASTA file
for hit_dict in hdna_hunter.scan_fasta("genome.fa", minrep=10, purity=0.85):
    print(hit_dict)
```

### From the command line

```bash
# Minimal — scan test.fa and write test_HDNA.tsv
hdna-hunter -seq test.fa -out test

# Strict mode — exact mirror, 100 % pure, like non-B_gfa
hdna-hunter -seq genome.fa -out strict -purity 1.0 -mismatch 0.0

# Relaxed mode — allow up to 20 % mismatch, 80 % purity
hdna-hunter -seq genome.fa -out relaxed -purity 0.80 -mismatch 0.20

# Verbose output, longer arms only
hdna-hunter -seq genome.fa -out long -minrep 10 -maxrep 50 -v
```

---

## 5. Python API

### `hdna_hunter.scan_sequence`

```python
hdna_hunter.scan_sequence(
    seq: str,
    *,
    minrep: int = 6,
    maxrep: int = 50,
    maxspacer: int = 7,
    purity: float = 0.80,
    mismatch: float = 0.20,
    remove_overlaps: bool = True,
    seq_offset: int = 1,
) -> list[dict]
```

Scan a single DNA string for H-DNA mirror repeat motifs. The C core runs with the Python GIL released, so other threads remain unblocked during the scan.

**Parameters** — see [Parameter Reference](#8-parameter-reference).

**Returns** a list of hit dictionaries. Each dictionary contains all fields described in [Output Format](#7-output-format).

---

### `hdna_hunter.scan_fasta`

```python
hdna_hunter.scan_fasta(
    path: str | Path,
    *,
    minrep: int = 6,
    maxrep: int = 50,
    maxspacer: int = 7,
    purity: float = 0.80,
    mismatch: float = 0.20,
    remove_overlaps: bool = True,
) -> list[dict]
```

Scan every record in a FASTA file. Adds a `seq_id` key to each hit dict. Genomic offsets are parsed automatically from UCSC/Ensembl-style headers (`>seqid:start-end`); other headers default to offset 1.

---

### `hdna_hunter.parse_fasta`

```python
hdna_hunter.parse_fasta(path: str | Path) -> Generator[tuple[str, str, int], None, None]
```

Low-level FASTA parser. Yields `(seq_id, sequence, offset)` tuples where:
- `seq_id` — the identifier after `>` up to the first whitespace.
- `sequence` — uppercase sequence string.
- `offset` — 1-based genomic start of the first base (parsed from `>id:start-end` headers).

---

## 6. Command-Line Interface

```
hdna-hunter -seq <FASTA> -out <PREFIX> [options]
```

or equivalently:

```
python -m hdna_hunter -seq <FASTA> -out <PREFIX> [options]
```

### Positional / required arguments

| Argument | Description |
|---|---|
| `-seq FASTA` | Path to the input FASTA file. May contain one or multiple records. Both single-line and wrapped FASTA are supported. |
| `-out PREFIX` | Output file prefix. The tool writes `<PREFIX>_HDNA.tsv`. |

### All optional arguments

| Argument | Type | Default | Description |
|---|---|---|---|
| `-minrep INT` | int | `6` | **Minimum arm length** in base pairs. Arms shorter than this are not reported, even if they satisfy all other criteria. Increasing this value reduces noise and focuses on structurally stable, longer H-DNA elements (recommended ≥ 10 for whole-genome scans). |
| `-maxrep INT` | int | `50` | **Maximum arm length** cap. The algorithm stops extending an arm beyond this length. Set higher for repetitive regions; setting it very high increases runtime without improving sensitivity for typical H-DNA. |
| `-maxspacer INT` | int | `7` | **Maximum spacer / hinge loop length** in base pairs. The spacer is the single-stranded loop between the two mirror arms. H-DNA with spacers > 10 bp is thermodynamically unfavourable in vivo; the default of 7 reflects experimentally validated structures. Setting this to 0 requires the two arms to be immediately adjacent. |
| `-purity FLOAT` | float | `0.80` | **Minimum compositional purity** of each arm (0.0–1.0). The arm must be ≥ `purity` fraction GA (purine) **or** ≥ `purity` fraction CT (pyrimidine). Set to `1.0` to require a perfectly pure homopurine/homopyrimidine arm (equivalent to strict non-B_gfa mode). Lower values (e.g. `0.70`) increase sensitivity at the cost of more false positives. |
| `-mismatch FLOAT` | float | `0.20` | **Maximum mirror mismatch fraction** (0.0–1.0). Fraction of base positions in the arm where `dna[left] ≠ dna[right]` (i.e. the mirror is broken). `0.0` requires a perfect mirror; `0.20` allows 1 mismatch per 5 bp. This parameter is independent of purity — both must be satisfied simultaneously. |
| `-skipoverlap` | flag | *(off)* | **Skip overlap removal**. By default, overlapping hits are collapsed to the longest arm (ties broken by shortest spacer). Pass this flag to disable overlap removal and receive every raw hit at every `(center, spacer)` combination that passes the thresholds. Useful for statistical analyses or when you want to inspect the full hit landscape. |
| `-v` | flag | *(off)* | **Verbose mode**. Prints per-sequence statistics (record name, length, offset, hit count) to `stderr` as each FASTA record is processed. Useful for monitoring progress on large genomes. |

### Parameter interaction summary

```
Sensitivity ──────────────────────────────────── Specificity
  ↑ lower purity            higher purity ↑
  ↑ higher mismatch         lower mismatch ↑
  ↑ smaller minrep          larger minrep ↑
  ↑ larger maxspacer        smaller maxspacer ↑
```

### CLI examples

```bash
# 1. Basic scan
hdna-hunter -seq genome.fa -out results

# 2. Strict mode (exact mirror, 100 % pure) — matches non-B_gfa output
hdna-hunter -seq genome.fa -out strict -purity 1.0 -mismatch 0.0

# 3. Whole-genome scan with minimum arm length 10 for high confidence
hdna-hunter -seq hg38.fa -out hg38_hdna -minrep 10 -purity 0.85 -v

# 4. Exploratory scan — very permissive, keep all raw hits
hdna-hunter -seq region.fa -out explore -purity 0.60 -mismatch 0.30 \
            -minrep 6 -maxspacer 10 -skipoverlap

# 5. Long perfect H-DNA only
hdna-hunter -seq genome.fa -out perfect -purity 1.0 -mismatch 0.0 \
            -minrep 12 -maxspacer 3

# 6. Verbose output to monitor progress on a large genome
hdna-hunter -seq hg38.fa -out hg38 -minrep 8 -v 2>progress.log
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

### Example output

```
seq_id  source    start  end  arm_length  spacer_length  total_length  ga_pct   ct_pct   mirror_identity  is_perfect  left_arm  spacer  right_arm  full_sequence
chr1    findHDNA  1001   1020  7           6              20            85.71    14.29    100.0            False       gggaaat   tttttt  taaaggg    gggaaattttttttaaaggg
```

---

## 8. Parameter Reference

| Parameter | API name | CLI flag | Type | Default | Valid range | Notes |
|---|---|---|---|---|---|---|
| Minimum arm length | `minrep` | `-minrep` | int | 6 | ≥ 1 | Arms shorter than this are discarded entirely |
| Maximum arm length | `maxrep` | `-maxrep` | int | 50 | ≥ `minrep` | Hard cap on extension; rarely needs changing |
| Maximum spacer | `maxspacer` | `-maxspacer` | int | 7 | ≥ 0 | Set to 0 for zero-loop (adjacent arms) only |
| Purity threshold | `purity` | `-purity` | float | 0.80 | 0.0–1.0 | Fraction GA or CT required in arm |
| Mismatch tolerance | `mismatch` | `-mismatch` | float | 0.20 | 0.0–1.0 | Fraction of arm positions allowed to mismatch the mirror |
| Overlap removal | `remove_overlaps` | `--skipoverlap` (inverts) | bool | True | — | When True, keeps only the longest non-overlapping hit |
| Sequence offset | `seq_offset` | *(automatic in CLI)* | int | 1 | ≥ 1 | 1-based start coordinate of `seq[0]`; used for correct genomic coordinates |

### Choosing parameters for your use case

**Genome-wide survey (high confidence)**
```
minrep=10, maxrep=50, maxspacer=7, purity=0.85, mismatch=0.10
```

**Reproduce non-B_gfa exact behaviour**
```
minrep=6, maxrep=50, maxspacer=7, purity=1.0, mismatch=0.0
```

**Exploratory / maximum sensitivity**
```
minrep=6, maxrep=50, maxspacer=10, purity=0.70, mismatch=0.30
```

**Downstream ML or statistical analysis (all raw candidates)**
```
remove_overlaps=False, skipoverlap (CLI)
```

---

## 9. Streamlit Web App

An interactive web application is included in the package. It provides:
- Paste-in or file-upload of a FASTA sequence
- Interactive parameter sliders
- Annotated sequence viewer highlighting arms and spacers
- Summary statistics table with Plotly charts

### Launch the app

```bash
pip install "hdna-hunter[app]"
streamlit run $(python -c "import hdna_hunter, pathlib; print(pathlib.Path(hdna_hunter.__file__).parent / 'app.py')")
```

Or from the cloned repository:

```bash
streamlit run src/hdna_hunter/app.py
```

---

## 10. Understanding the Results

### Interpreting `ga_pct` and `ct_pct`

These two values always sum to ≤ 100 %. The difference is due to ambiguous or masked bases:
- **GA-rich arm** (`ga_pct ≈ 80–100`): purine-rich tract → **H-r triplex** (Pu·Pu·Py)
- **CT-rich arm** (`ct_pct ≈ 80–100`): pyrimidine-rich tract → **H-y triplex** (Py·Pu·Py)
- **Mixed** (`neither > purity`): filtered out unless `purity` is set very low

### `mirror_identity` vs purity

`mirror_identity` measures how closely the two arms are mirror images of each other. A value of 100 means every base on the left arm is identical to its mirror partner on the right arm. A value of 80 means 80 % of positions match (1 mismatch per 5 bp in a 6-mer arm). This is independent of purity: a CT-rich arm with one GA interruption could still have 100 % mirror identity.

### `is_perfect` flag

`is_perfect = True` is a stringent classification reserved for hits that are both:
1. **100 % pure** — every base in the arm is GA, or every base is CT.
2. **100 % mirror-symmetric** — no mismatches between the two arms.

These are the most likely candidates for stable H-DNA structures.

### Overlap removal behaviour

When `remove_overlaps=True` (default), two hits overlap if their genomic ranges on the DNA share at least one base. Among all overlapping hits, the one with the **longest arm** is retained. Ties are broken by choosing the **shorter spacer**. This is identical to the behaviour of the original non-B_gfa suite.

---

## 11. Performance Notes

- The C extension releases the Python GIL during scanning, enabling multi-threaded use.
- A single human chromosome (≈ 250 Mb) completes in under 30 seconds on a single core with default parameters.
- The internal hit buffer is capped at **1,000,000 hits per `scan_sequence` call**. For whole-genome scans, call `scan_fasta()` (which iterates per-contig) rather than concatenating chromosomes into a single string.
- Memory usage is approximately **O(n_hits × 200 bytes)** plus the sequence string itself.
- N bases in the sequence (`N`, `n`) break arm extension, as per the original algorithm.

---

## 12. Development & Testing

### Running the test suite

```bash
# Install dev dependencies first
pip install -e ".[dev]"

# Run all 110 tests
pytest -v tests/
```

The test suite covers:
- Empty / trivial inputs
- Known GA and CT mirror motifs with exact expected values
- Strict mode (purity=1.0, mismatch=0.0) and relaxed mode
- Coordinate offset propagation (`seq_offset`)
- Multi-record FASTA parsing
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
hdna_hunter/
├── src/
│   └── hdna_hunter/
│       ├── _hdna.c          # C extension — core algorithm
│       ├── __init__.py      # Python API (scan_sequence, scan_fasta, parse_fasta)
│       ├── __main__.py      # CLI entry point (hdna-hunter / python -m hdna_hunter)
│       └── app.py           # Streamlit web application
├── tests/
│   └── test_hdna.py         # 110 comprehensive tests (pytest)
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

Tests live in `tests/test_hdna.py`. Each section focuses on one concern. Add new sequences as module-level constants and group related assertions under a descriptively named function. Run `pytest -v tests/test_hdna.py -k <your_test_name>` to run just your new test.

---

## 13. Citation

If you use H-DNA Hunter in published research, please cite:

> Georgakopoulos-Soares Lab, Penn State University.  
> **H-DNA Hunter**: fast, cross-platform detection of H-DNA and mirror repeat structures in genomic sequences.  
> https://github.com/Georgakopoulos-Soares-lab/HDNAhunter

The core algorithm is inspired by and validated against:

> Cer, R.Z. et al. (2013). **Non-B DB v2.0: a database of predicted non-B DNA-forming motifs and its associated tools.** *Nucleic Acids Research*, 41(D1), D94–D100. https://doi.org/10.1093/nar/gks955

---

## 14. License

MIT License — see [LICENSE](LICENSE) for full text.

Copyright © 2024 Georgakopoulos-Soares Lab, Penn State University.
