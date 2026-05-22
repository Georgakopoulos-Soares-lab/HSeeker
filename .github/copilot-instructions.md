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
├── webapp/
│   ├── main.py              ← FastAPI application
│   ├── requirements.txt     ← webapp-only dependencies
│   ├── templates/
│   │   ├── base.html        ← shared nav/footer layout
│   │   ├── index.html       ← upload / paste page
│   │   ├── job.html         ← job results dashboard
│   │   └── about.html       ← informational page
│   └── static/
│       ├── css/style.css    ← dark theme styles
│       └── hseeker_validation_motifs.fasta  ← 39-sequence example
├── Dockerfile               ← multi-stage build (builder → runtime)
├── railway.toml             ← Railway deployment config
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

The maintainers use the `hseeker_bench` conda environment. Run tests with:

```bash
conda run -n hseeker_bench pytest -v tests/
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

---

## Installation

### Developer / local install (recommended)

The project uses a `hseeker_bench` conda environment. All dev commands should be run inside it.

```bash
conda activate hseeker_bench

# Install hseeker + dev extras (compiles _hdna.c automatically)
pip install -e ".[dev]"

# Verify
python -c "import hseeker; print(hseeker.__version__)"
pytest -v tests/
```

The `[dev]` extra installs `pytest`, `cibuildwheel`, `build`, and `twine`.

### Webapp install (separate requirement set)

The webapp has its own `requirements.txt` because `fastapi`, `uvicorn`, `plotly`, `pandas`, and `numpy` are **not** listed as `hseeker` package dependencies. Install them separately inside the same conda environment:

```bash
pip install -r webapp/requirements.txt
```

### Running the webapp locally

```bash
conda run -n hseeker_bench uvicorn webapp.main:app --reload --port 8000
# or after activating:
uvicorn webapp.main:app --reload --port 8000
```

The `--reload` flag auto-reloads on Python source changes. Template and static file changes are picked up immediately without reload.

---

## Webapp (`webapp/`)

### Architecture

- **FastAPI** application in `webapp/main.py`.
- **Jinja2** templates in `webapp/templates/`.
- **Alpine.js 3.14.1** (CDN) for all client-side reactivity.
- **Plotly.js 2.35.2** (CDN) for chart rendering. Charts are generated server-side via `plotly.graph_objects` → `fig.to_json()`, embedded as `<script type="application/json">` in `job.html`, and rendered client-side lazily.
- **No database** — jobs are stored in a Python dict (`_jobs`) in process memory, protected by `threading.Lock`. Jobs expire after `JOB_TTL_SECONDS` (default 2 h). This means **all jobs are lost on restart**.

### Routes

| Method | Path | Purpose |
|---|---|---|
| `GET` | `/` | Upload / paste page (`index.html`) |
| `GET` | `/about` | Informational page (`about.html`) |
| `POST` | `/submit` | Accept upload or pasted FASTA; create job; redirect to `/job/{id}` |
| `GET` | `/job/{job_id}` | Job status + results page (`job.html`) |
| `GET` | `/api/job/{job_id}/status` | JSON status polling (used by `job.html` Alpine component) |
| `GET` | `/api/job/{job_id}/hits` | Paginated + filtered hit data (JSON) |
| `GET` | `/job/{job_id}/download` | Download results as TSV |
| `GET` | `/health` | Health check — returns `{"status":"ok","version":"..."}` |

### Job processing

1. `POST /submit` writes the upload to a temp file, creates a job dict, and spawns a `threading.Thread` running `_run_job()`.
2. `_run_job()` calls `hseeker.parse_fasta()` record by record, calls `hseeker.scan_sequence()` for each, accumulates hits, then calls `compute_stats()`, `generate_charts()`, and `_write_tsv()`.
3. Charts are generated **once** at job completion and stored in `job["charts"]` as Plotly figure dicts. They are serialised with `json.dumps` and injected into `job.html` at render time.
4. The temp FASTA file is deleted in the `finally` block of `_run_job()` regardless of success or failure.

### `/api/job/{job_id}/hits` filter parameters

| Param | Type | Default | Effect |
|---|---|---|---|
| `page` | int | 1 | Page number |
| `per_page` | int | 50 | Rows per page |
| `q` | str | `""` | Substring match against `seq_id` |
| `min_arm` | int | 0 | `arm_length >= min_arm` |
| `max_arm` | int | 9999 | `arm_length <= max_arm` |
| `min_mirror` | float | 0.0 | `mirror_identity >= min_mirror` |
| `perfect_only` | bool | false | Keep only `is_perfect == true` |
| `seq_type` | str | `""` | `"ga"` → `ga_pct >= 80`; `"ct"` → `ct_pct >= 80` |

### Charts

Eight charts are generated in `generate_charts()`, keyed by name:

| Key | Type | Notes |
|---|---|---|
| `hits_per_seq` | Horizontal bar | Top 25 sequences by hit count; height scales with row count |
| `perfect_donut` | Pie/donut | Perfect vs Imperfect; `plot_bgcolor="#000000"` required for label contrast |
| `arm_length_dist` | Bar | Arm length value counts |
| `spacer_dist` | Bar | Spacer length value counts |
| `mirror_dist` | Bar | 50-bin histogram; dashed vline at mismatch threshold |
| `composition_scatter` | Scatter | GA% vs CT% (right arm); up to 4,000 points sampled |
| `arm_spacer_heatmap` | Heatmap | 2-D density of arm × spacer lengths |
| `locus_map` | Scatter | Top 15 sequences; backbone line + marker per site; up to 600 sites sampled per sequence |

**Critical**: `_base_layout()` returns a dict that already contains `xaxis`, `yaxis`, and `margin` keys. Do **not** pass those keys again as kwargs to `fig.update_layout(**_base_layout(), xaxis=..., yaxis=...)` — it causes a `TypeError`. Use `fig.update_xaxes()` / `fig.update_yaxes()` separately.

### Templates

| File | Role |
|---|---|
| `base.html` | Nav (logo, About link, GitHub, PyPI, version badge), footer, CDN scripts (Alpine.js, Plotly.js, Inter font), global CSS link |
| `index.html` | Input form: file upload tab + paste textarea tab. "Load 39-sequence validation example" button in `.input-mode-bar` — always visible in both modes. Clicking loads the example into textarea and switches to Paste mode via `loadExample()`. Parameter accordion. |
| `job.html` | Status poller (Alpine `dataTable(jobId)` component). Stats bar. Four tabs: Data Table, Locus Map, Structure, Composition. Collapsible filter panel (`filtersOpen: false` default). Nucleotide coloring via `colorSeq(s)` (`nt-A`/`nt-C`/`nt-G`/`nt-T` CSS classes). TSV download button. |
| `about.html` | Static page: What is H-DNA → How HSeeker Works → Parameter Modes → Validation → Limitations → Links |

### Alpine.js constraints

- **Never** use `<template x-if>` inside an Alpine-managed scope (especially `<tbody>`). Use `<tr x-show>` instead — `<template x-if>` causes `_x_dataStack` null errors.
- `hseeker_version` is injected into all templates as a global Jinja2 variable via `templates.env.globals["hseeker_version"]`.

### CSS design tokens (`webapp/static/css/style.css`)

| Variable | Value | Usage |
|---|---|---|
| `--accent` | `#39d5ab` | Primary teal — buttons, highlights |
| `--bg` | `#0d1117` | Page background |
| `--surface` | `#161b22` | Card background |
| `--surface2` | `#1c2128` | Nested surfaces |
| `--border` | `#30363d` | Default borders |
| `--border2` | `#21262d` | Subtle borders |
| `--text` | `#c9d1d9` | Body text |
| `--text2` | `#e6edf3` | Headings |
| `--text-sub` | `#7d8590` | Secondary / muted text |

Nucleotide color classes: `.nt-A {#4caf50}` `.nt-C {#ff9800}` `.nt-G {#ef5350}` `.nt-T {#64b5f6}`.

### Environment variables

| Variable | Default | Effect |
|---|---|---|
| `MAX_UPLOAD_MB` | `200` | Max upload file size in MB |
| `JOB_TTL_SECONDS` | `7200` | How long job results are kept in memory (2 h) |
| `PORT` | `8000` | Uvicorn listen port (set by Railway automatically) |

---

## Deployment

### Docker (local)

```bash
# Build
docker build -t hseeker-webapp .

# Run (maps port 8000)
docker run -p 8000:8000 hseeker-webapp

# Custom limits
docker run -p 8000:8000 -e MAX_UPLOAD_MB=500 -e JOB_TTL_SECONDS=3600 hseeker-webapp
```

The `Dockerfile` is a two-stage build:
1. **builder** (`python:3.11-slim`): installs `gcc` + `python3-dev`, compiles the C extension via `pip install ".[app]"`, installs webapp requirements.
2. **runtime** (`python:3.11-slim`): copies installed packages from builder, copies `webapp/` only, runs as non-root `appuser`.

The entrypoint is:
```
CMD uvicorn webapp.main:app --host 0.0.0.0 --port "${PORT:-8000}"
```

### Railway

The `railway.toml` declares:
```toml
[build]
builder = "DOCKERFILE"
dockerfilePath = "Dockerfile"

[deploy]
healthcheckPath = "/health"
healthcheckTimeout = 300
restartPolicyType = "ON_FAILURE"
restartPolicyMaxRetries = 3
```

To deploy:
1. Push to `main` — Railway auto-deploys on any push.
2. Railway injects `PORT` automatically; no manual port configuration needed.
3. Health checks hit `/health` which returns `{"status":"ok","version":"..."}`.

**Important**: Railway does not mount persistent volumes by default. All in-memory jobs and generated TSV files are lost on redeploy/restart. This is by design — jobs are ephemeral.

---

## Package updating and releasing to PyPI

### Updating the version

The version is declared once in `pyproject.toml`:
```toml
[project]
version = "0.1.0"
```

Update it there, then also update `src/hseeker/__init__.py` if it defines `__version__` manually (check that file). Keep them in sync.

### Local build verification

```bash
conda activate hseeker_bench

# Recompile and run tests
pip install -e ".[dev]"
pytest -v tests/

# Build source distribution and wheel
python -m build

# Check the distribution
twine check dist/*
```

### Releasing to PyPI (automated via CI)

The GitHub Actions workflow (`.github/workflows/build_wheels.yml`) handles all PyPI publishing automatically via OIDC Trusted Publishing — no tokens or secrets needed.

**To trigger a release:**
```bash
git tag -a v0.2.0 -m "Release v0.2.0"
git push origin v0.2.0
```

The workflow will:
1. Build wheels for Linux (x86_64 + aarch64), macOS, Windows across Python 3.9–3.13 using `cibuildwheel`.
2. Build an sdist.
3. Publish everything to PyPI.

**Do not** push to PyPI manually with `twine upload` — use the tag workflow.

### Updating webapp dependencies

The webapp dependencies in `webapp/requirements.txt` are **not** part of the `hseeker` package distribution. To update them:

```bash
# Edit webapp/requirements.txt with new versions
# Then reinstall
pip install -r webapp/requirements.txt

# Test the webapp still works
uvicorn webapp.main:app --reload --port 8000
```

The Docker image will pick up the new versions automatically on the next `docker build`.

### Updating the `[app]` extra in `pyproject.toml`

The `[project.optional-dependencies]` section has an `app` extra (currently `streamlit` and `plotly`) that is used by the Dockerfile's first stage (`pip install ".[app]"`). This is separate from `webapp/requirements.txt`. If you add a new Python import to `webapp/main.py`, add the package to **both** `webapp/requirements.txt` (for local dev) and the `[app]` extra in `pyproject.toml` (for Docker builds).
