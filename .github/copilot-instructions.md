# GitHub Copilot Instructions — HSeeker

This file provides Copilot with deep context about the HSeeker codebase so that it can suggest accurate, idiomatic changes. Read this before proposing any modifications.

---

## Project at a glance

**HSeeker** (`hseeker` on PyPI) detects H-DNA / triplex mirror repeat structures in DNA sequences. It is a Python package whose computational core is a **compiled CPython C extension** (`_hdna.c`). The Python layer provides a clean API and a CLI that is a drop-in replacement for the original `findHDNA` binary from MirrorHunter / TriplexDetector.

- Repository: https://github.com/Georgakopoulos-Soares-lab/HDNAhunter
- PyPI: https://pypi.org/project/hseeker/
- Python 3.9–3.13, Linux / macOS / Windows
- MIT license

---

## Repository structure

```
hseeker/
├── src/
│   └── hseeker/
│       ├── _hdna.c          ← C extension (DO NOT refactor lightly)
│       ├── __init__.py      ← Public Python API
│       ├── __main__.py      ← CLI entry point
├── tests/
│   └── test_hdna.py         ← 120 pytest tests
├── benchmarks/
│   ├── benchmark.py         ← performance benchmark suite
│   └── data/                ← generated FASTA files (gitignored)
├── .github/
│   ├── copilot-instructions.md   ← this file
│   └── workflows/
│       └── build_wheels.yml ← cibuildwheel CI
├── pyproject.toml
├── setup.py                 ← C extension declaration
└── MANIFEST.in
```

---

## The C extension (`src/hseeker/_hdna.c`)

### Architecture

The C file has a `#ifdef STANDALONE` guard so it can be compiled either as:
- A **CPython extension** (`pip install`) — the default.
- A **standalone CLI binary** (`gcc -DSTANDALONE -O2 -o findHDNA _hdna.c -lm`).

### Key functions

| Function | Purpose |
|---|---|
| `findHDNA_core()` | Center-outward mirror scanner. Fills an `HDNA_HIT[]` array. |
| `remove_overlaps_core()` | Keeps longest non-overlapping hits (tiebreak: shorter spacer). |
| `py_scan_sequence()` | Python-callable wrapper; releases GIL during `findHDNA_core`. |
| `PyInit__hdna()` | Module initialization, registers `scan_sequence` method. |

### Algorithm invariants — never break these

1. **Center-outward scan**: for each `(ctr, sp)`, left pointer starts at `ctr` and moves LEFT; right pointer starts at `ctr + sp + 1` and moves RIGHT.
2. **Composition tracking uses the RIGHT arm** (`rb` base), not the left.
3. **`is_perfect`** is true iff `(best_ga == 1.0 || best_ct == 1.0) && best_mir == 1.0` — this is checked using `float` equality against `1.0f`. Do not change this to a threshold comparison.
4. **Coordinates are 1-based inclusive** — `start = left_start + 1 + (seq_offset - 1)`.
5. **N bases terminate arm extension** — `if (lb == 'n' || rb == 'n') break;`
6. **Overlap removal** keeps the longest arm; among ties it keeps the shorter spacer.

### When editing `_hdna.c`

- After any edit, run `pip install -e .` to recompile, then run `pytest -v tests/`.
- Do NOT change the `HDNA_HIT` struct field order without updating `py_scan_sequence()` which iterates those fields to build the Python dict.
- The `DEFAULT_MAX_HITS` constant (1,000,000) is a safety cap, not a target. Do not lower it without a compelling reason.
- Memory for `hits[]` is allocated in `py_scan_sequence()` with `PyMem_Malloc`. Ensure any added code does not leak on error paths (all error paths goto `cleanup:`).

---

## Python API (`src/hseeker/__init__.py`)

### Public surface

```python
scan_sequence(seq, *, minrep=6, maxrep=50, maxspacer=7,
              purity=0.80, mismatch=0.20,
              remove_overlaps=True, seq_offset=1) -> list[dict]

scan_fasta(path, *, minrep=6, ...) -> list[dict]

scan_fasta_iter(path, *, minrep=6, ...) -> Generator[dict, None, None]

scan_fasta_parallel(path, *, minrep=6, ..., workers=None) -> list[dict]
```

### Parameter naming contract

The Python parameter names (`minrep`, `maxrep`, `maxspacer`, `purity`, `mismatch`) are the **canonical names** used everywhere — in `__init__.py`, `__main__.py`, the C extension, the test file, and the README. Never rename them without updating all five locations.

---

## CLI (`src/hseeker/__main__.py`)

- Entry point: `hseeker` (declared in `pyproject.toml` under `[project.scripts]`).
- Also callable as `python -m hseeker`.
- Output: `<prefix>_HDNA.tsv` — a tab-separated file with the column order defined in `fieldnames` in `main()`.
- The CLI flag `-skipoverlap` **inverts** the Python `remove_overlaps` parameter (i.e. passing `-skipoverlap` sets `remove_overlaps=False`).
- All progress/diagnostic output goes to `stderr`; only the TSV header/rows go to the output file.

---

## Output dict / TSV columns

Every hit (from `scan_sequence` or `scan_fasta`) has these keys:

| Key | Type | Notes |
|---|---|---|
| `start` | int | 1-based, inclusive, left arm start |
| `end` | int | 1-based, inclusive, right arm end |
| `arm_length` | int | Length of each arm (both equal) |
| `spacer_length` | int | Hinge loop length |
| `total_length` | int | `arm_length * 2 + spacer_length` |
| `ga_pct` | float | % GA in right arm (0–100) |
| `ct_pct` | float | % CT in right arm (0–100) |
| `mirror_identity` | float | Mirror match % (0–100) |
| `is_perfect` | bool | 100 % pure AND 100 % mirror |
| `left_arm` | str | Left arm sequence (5′→3′) |
| `spacer` | str | Hinge loop sequence |
| `right_arm` | str | Right arm sequence (5′→3′) |
| `full_sequence` | str | `left_arm + spacer + right_arm` |

`scan_fasta` additionally adds:
- `seq_id` — FASTA record identifier.
- `source` — always `"findHDNA"`.

The invariant `start + total_length - 1 == end` always holds. Tests verify this.

---

## Test suite (`tests/test_hdna.py`)

120 tests across 25 sections. Before adding a new test:

1. Read the section it belongs to and follow the existing naming pattern.
2. Add reference sequences as **module-level constants** (e.g., `MY_SEQ = "AAAGGGAAAGGG"`).
3. Run the diagnostic snippet below if you are unsure what the algorithm actually returns for your sequence:

```python
import hseeker
print(hseeker.scan_sequence("YOUR_SEQ", minrep=6, remove_overlaps=False))
```

### Critical known values (do not change without re-running diagnostics)

| Constant | Sequence | Expected result |
|---|---|---|
| `GAMIR_SEQ` | `"GGGAAATTAAAGGG"` | arm=7, sp=0, ga=85.71%, ct=14.29%, mir=100%, is_perfect=False |
| `CTMIR_SEQ` | `"CCCTTTAATTTCCC"` | arm=7, sp=0, ct=85.71%, ga=14.29%, mir=100%, is_perfect=False |
| `PERFECT_SEQ` | `"A" * 14` | arm=7, sp=0, ga=100%, mir=100%, is_perfect=True |
| `SPACER8_SEQ` | `"AAAAAA" + "CTTTTTTT" + "AAAAAA"` | arm=6, sp=8 visible only with `remove_overlaps=False` and `maxspacer=8`; zero hits with `maxspacer=7` |

### Why GAMIR_SEQ returns arm=7 (not arm=6)

`"GGGAAATTAAAGGG"` — the algorithm at `(ctr=6, sp=0)` extends to arm=7 using `left_arm="gggaaat"`, `right_arm="taaaggg"`. The arm=7 hit survives overlap removal over any arm=6 hit at other `(ctr, sp)` pairs because it is longer. This is correct and expected behaviour.

### Why SPACER8_SEQ needs `remove_overlaps=False` to see sp=8

`"AAAAAA" + "CTTTTTTT" + "AAAAAA"` — at `(ctr=6, sp=6)` the algorithm finds arm=7, sp=6 (T at both sides → ct match, then A's). Overlap removal keeps arm=7 over arm=6 because it is longer, even when `maxspacer=8`. Use `remove_overlaps=False` to verify the arm=6, sp=8 hit exists before overlap removal.

---

## Build system

| File | Role |
|---|---|
| `pyproject.toml` | Project metadata, dependencies, setuptools config, cibuildwheel config, pytest config |
| `setup.py` | `Extension("hseeker._hdna", sources=["src/hseeker/_hdna.c"])` |
| `MANIFEST.in` | Ensures `_hdna.c` is included in the sdist |

Build backend: `setuptools.build_meta` (NOT the legacy backend). This is declared in `[build-system]` in `pyproject.toml` — never change it to `setuptools.build_meta:__legacy__`.

Compile and install in editable mode:

```bash
pip install -e ".[dev]"
```

---

## CI / CD (`.github/workflows/build_wheels.yml`)

The workflow:
1. Triggers on `v*` tags and `workflow_dispatch`.
2. Runs `cibuildwheel` on `ubuntu-latest`, `macos-latest`, `windows-latest`.
3. QEMU is set up on Linux runners to enable `aarch64` cross-compilation.
4. Builds CPython 3.9–3.13 wheels; skips 32-bit and musl targets.
5. Builds an sdist on `ubuntu-latest`.
6. Uploads all artifacts; on tag events, publishes to PyPI via OIDC Trusted Publishing.

To trigger a release: create and push an annotated tag `vX.Y.Z`.

---

## Benchmarks (`benchmarks/benchmark.py`)

Self-contained performance benchmark suite. Generates synthetic FASTA datasets and measures wall time, peak RAM, CPU utilisation, and parallelism speedup for all public API paths.

### Running

```bash
# Small tier only (30 MB, ~3–5 min) — recommended default
python benchmarks/benchmark.py --no-cli

# Include medium tier (300 MB)
python benchmarks/benchmark.py --medium --no-cli

# Benchmark a real FASTA file
python benchmarks/benchmark.py --real path/to/genome.fa --no-cli

# Save results as JSON
python benchmarks/benchmark.py --no-cli --json results.json

# Parallelism scaling table
python benchmarks/benchmark.py --scaling --no-cli
```

### Key facts

- Generated FASTA files are cached in `benchmarks/data/` (gitignored). Pass `--no-cache` to regenerate.
- Three sequence profiles: `uniform` (sparse hits, baseline), `ga_biased` (dense hits, stresses buffer), `realistic` (medium density, mimics human chromosomes).
- Three size tiers: `small` (30 MB, default), `medium` (300 MB, `--medium`), `large` (3 GB, `--large`).
- Do NOT commit anything under `benchmarks/data/`.
- When editing `benchmark.py`, do not change the `DatasetSpec` field names — they are used as JSON keys in `--json` output.

---

## Environment

The maintainers use the `biomni_e1` conda environment. Run tests with:

```bash
conda run -n biomni_e1 pytest -v tests/
```

or after activating the environment:

```bash
pytest -v tests/
```

---

## Code style

- Python: PEP 8, type annotations on all new public functions, `from __future__ import annotations`.
- C: K&R brace style, `/* block comments */`, no C99 VLAs.
- Commit messages: imperative mood, e.g. `Add maxrep validation`, `Fix overlap removal tiebreak`.
- Do not add docstrings, type annotations, or comments to code that was not changed in a PR.

---

## Common pitfalls

| Mistake | Correct approach |
|---|---|
| Expecting `scan_sequence` to return arm=6 for `GAMIR_SEQ` | It returns arm=7. The center-outward scan at sp=0 finds the longer arm. |
| Testing `hits == []` for sequences with a shorter valid arm | Test `[h for h in hits if h["spacer_length"] == N] == []` instead. |
| Forgetting to recompile after editing `_hdna.c` | Run `pip install -e .` before running tests. |
| Using `build-backend = "setuptools.build_meta:__legacy__"` | Always use `"setuptools.build_meta"` (no `:__legacy__` suffix). |
| Passing `seq_offset` as a positional argument | All `scan_sequence` parameters after `seq` are keyword-only (`*,`). |
