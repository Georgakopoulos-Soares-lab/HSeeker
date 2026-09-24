# WS2 pre-work — first HSeeker baseline on benchmark v3

**Status:** PROVISIONAL. This is **not** WS2-A2 and must not be cited in the
response letter. See "Why provisional" below.
**Date:** 2026-09-23
**Input:** `hdna_benchmark_balanced_v3.csv`, sha256 `f4602d6453d0926c…`

## Why this exists

The deposited results in `analysis/results/` were run on **2026-06-24** against the
**old** CSV (`hdna_experimental_sequences_final.csv`, 80 rows, 71 forming / 9
non-forming), on a different machine, with code predating the AT/homopolymer filter
commit (`6ecaa84`, 2026-07-02). **No HSeeker result had ever been produced on the v3
benchmark.** Establishing that baseline is a prerequisite for every comparator
comparison (WS3-A3 depends on WS2-A2), so it was run before going further with
NeSSie/nBMST.

## Result — HSeeker detection call, current script parameters

`minrep=8, maxspacer=10, purity=0.90, mismatch=0.10, purity_rmq=True` — the values
the direct script uses today, **not** the manuscript's (DEC-5 unresolved).

| subset | n | pos | neg | TP | FN | TN | FP | sens | spec | prec | MCC | acc | AUC |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| balanced (all) | 130 | 65 | 65 | 63 | 2 | 61 | 4 | 0.969 | 0.938 | 0.940 | 0.908 | 0.954 | 0.963 |
| experimental | 71 | 65 | 6 | 63 | 2 | 5 | 1 | 0.969 | 0.833 | 0.984 | 0.749 | 0.958 | 0.977 |
| synthetic | 59 | 0 | 59 | 0 | 0 | 56 | 3 | — | 0.949 | — | — | 0.949 | n/a |

The balanced set does what it was built to do: specificity now rests on **65**
negatives instead of 6, so it is no longer hostage to one or two records — the exact
objection in R3.1.

## The finding that matters for DEC-5 — every false positive is a homopolymer

All four FPs on the balanced set:

| Record | Type | Sequence |
|---|---|---|
| HDNA0053 `pRW1405` | experimental | `A`×20 |
| SYN0054 `poly-A30` | synthetic | `A`×30 |
| SYN0055 `poly-T25` | synthetic | `T`×25 |
| SYN0056 `poly-T35` | synthetic | `T`×35 |

These are precisely what the production CLI removes (`at_threshold=0.8`,
`filter_homopolymers=True`, `__main__.py:97-110`) and what the library defaults used
by the direct script do **not** (`at_threshold=None`, `filter_homopolymers=False`).
Register **D-4** is therefore not cosmetic — it decides the headline number:

### Balanced set (n=130)

| configuration | TP | FN | TN | FP | sens | spec | MCC |
|---|---|---|---|---|---|---|---|
| library defaults (**what the direct script uses**) | 63 | 2 | 61 | 4 | 0.969 | 0.938 | 0.908 |
| + AT filter only | 62 | 3 | 65 | 0 | 0.954 | **1.000** | 0.955 |
| + homopolymer filter only | 62 | 3 | 65 | 0 | 0.954 | **1.000** | 0.955 |
| **production CLI (both filters)** | 61 | 4 | 65 | 0 | 0.938 | **1.000** | 0.940 |

### Experimental only (n=71)

| configuration | TP | FN | TN | FP | sens | spec | MCC |
|---|---|---|---|---|---|---|---|
| library defaults | 63 | 2 | 5 | 1 | 0.969 | 0.833 | 0.749 |
| production CLI (both filters) | 61 | 4 | 6 | 0 | 0.938 | 1.000 | 0.750 |

**Reading:** either filter alone removes all four FPs and drives specificity to
1.000. The filters are not free — they cost 1–2 true positives (sens 0.969 → 0.938)
— but on the balanced set they raise MCC from 0.908 to 0.940. On the
experimental-only set MCC is flat (0.749 vs 0.750) because six negatives cannot
resolve the difference, which is itself an argument for dual reporting.

**Consequence:** the two arms must be pinned to one configuration before any number
is quoted (WS2-A1), and Methods must state which. Reporting the unfiltered direct
arm gives specificity 0.938; the production configuration gives 1.000. Both are
"HSeeker", and the difference is entirely the filter switch.

## The two false negatives are the same motif class

| Record | Sequence | Note |
|---|---|---|
| HDNA0058 `cMYC_NSE_ACCCTCCCC4` | `ACCCTCCCC` ×4 | pyrimidine-rich tandem **direct** repeat |
| HDNA0078 `gamma_globin_228_189` | `CTCCCAGC…` ×4 | pyrimidine-rich tandem **direct** repeat |

Both are ~36–37 nt tandem repeats of a ~9-nt unit, i.e. direct-repeat architecture
rather than the arm–spacer–arm mirror HSeeker searches for. Worth an explicit
sentence in the FN discussion (WS2-A5) rather than leaving them unexplained; a
reviewer who checks will notice they share a structure.

## Why provisional — four reasons it is not WS2-A2

1. **DEC-5 unresolved.** Run with the script's current `mismatch=0.10, minrep=8`, not
   the manuscript's `mismatch=0.15`. Numbers will move.
2. **Benchmark not frozen (WS0-A1).** Blocked on DEC-2/DEC-3. DEC-3 would add a
   positive (→132 rows, 66:66); the pCON terminal-C correction (D-9) changes one
   negative's sequence.
3. **HSeeker only — no Triplex comparator.** R/Rscript is not installed in this
   environment, so the direct script's Triplex arm could not run. Every
   tool-vs-tool claim is still unsupported.
4. **No confidence intervals, no family-aware evaluation.** Point estimates only;
   WS2-A2/WS2-A3 remain entirely outstanding, and R2.7/R3.1/R3.2/R3.3 require them.

## Reproduce

```
python3 -c "…"   # uses load_rows/run_hseeker_direct/metrics_from_binary
                  # from analysis/scripts/run_direct_sensitivity_analysis.py
                  # with --subset all|experimental|synthetic
```
