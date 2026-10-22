# DEC-8: restore the 1:1 balance with four replacement synthetic negatives

**Status:** done, verified (2026-10-07)
**Decision:** author instruction, 2026-10-07: "Add 4 non forming to balance."
**Follows:** DEC-7 (2026-10-02), which removed every pure single-base homopolymer from the
benchmark and left it at 124 records (64 forming / 60 non-forming).
**Result:** `hdna_benchmark_balanced_v3.csv` = **128 records, 64 forming / 64 non-forming**
(69 experimental + 59 synthetic), sha256 `7bb495444e91709fb0e322f08e5c9737e3a9419f94d9c881059073f725907cf9`.

## 1. Background: what DEC-7 removed

DEC-7 removed six records whose sequence is a single repeated base:

| Record | Set | Sequence | Label |
|---|---|---|---|
| HDNA0053 `pRW1405` | experimental | A×20 | non-forming |
| HDNA0054 `poly-(dC)30` | experimental | C×30 | **forming** |
| SYN0053 `poly-A15` | synthetic class E | A×15 | non-forming |
| SYN0054 `poly-A30` | synthetic class E | A×30 | non-forming |
| SYN0055 `poly-T25` | synthetic class E | T×25 | non-forming |
| SYN0056 `poly-T35` | synthetic class E | T×35 | non-forming |

In the experimental CSV, the two experimental records are marked `curation_decision = removed`,
with the justification "Removed (pure homopolymer): author decision (2026-10-02)…". The removal
took out one positive and five negatives, so the benchmark was no longer 1:1. DEC-8 restores
the balance by adding four negatives.

## 2. How the four replacements were generated

The replacements are made by `RTR/generate_synthetic_negatives.py`, using the generator's own
class definitions. Three choices keep the change minimal and auditable:

1. **Same recipe as the perturbation core.** The generator splits its non-fixed negatives
   roughly 50/25/25 across mirror-disrupted mutants (A), dinucleotide shuffles (B) and random
   GC/length-matched sequences (C). The four replacements follow that split:
   **2 × A, 1 × B, 1 × C** (`REPLACEMENT_COUNTS`).
2. **Every earlier record is unchanged.** The replacements are drawn *after* the full v3 pool,
   from an independent RNG stream (`REPLACEMENT_SEED = 20261002`). The first 124 data rows of
   the new file are byte-identical to the DEC-7 file (checked with `cmp`), and record IDs are
   unchanged. The new records continue the numbering at SYN0060–SYN0063, after SYN0057–0059
   (class F).
3. **No reuse, no homopolymer sources.** The A/B replacements take the next perturbation sources
   that the v3 pool did not use, so no source contributes twice. The C replacement is
   length/GC-matched only to forming records that remain in the benchmark, which excludes
   poly-dC HDNA0054.

To support this, the class A–C loop bodies were factored into helpers (`add_mutant`,
`add_shuffle`, `add_random`) and reused for the replacements. The v3 pool comes out identical,
which confirms the refactor did not change any RNG draw. Each new row records `seed = 20261002`,
and its `expected_label_rationale` ends with "Added 2026-10-07 to restore the 1:1 balance after
the pure-homopolymer exclusion."

## 3. The four new records

| Record | Name | Class | Derived from | Family | nt | Residual arm≥8 mirror (arm / identity) | Purine fraction |
|---|---|---|---|---|---|---|---|
| SYN0060 | pAA32_mirrordisrupted | A mirror-disrupted mutant | HDNA0001 pAA32 | SYN:fam:pAA32 | 32 | 9 / 0.778 | 0.938 |
| SYN0061 | GA32_mirrordisrupted | A mirror-disrupted mutant | HDNA0055 GA32 | SYN:fam:GA32 | 30 | 8 / 0.750 | 0.933 |
| SYN0062 | pGC32_dinucshuffle | B dinucleotide shuffle | HDNA0007 pGC32 | SYN:fam:pAA32 | 32 | 8 / 0.750 | 0.906 |
| SYN0063 | random_matched_to_HDNA0004 | C random GC/length-matched | HDNA0004 pTG32 | SYN:random | 32 | 8 / 0.750 | 0.125 |

Sequences (5'→3'):

- SYN0060 `AAGAGGAAAGAAAGTATAGGGGAAAGAGGGAA`
- SYN0061 `AAAAGGGGGGGGGTATAGGGAAAGAGGGAA`
- SYN0062 `AGGGAGGGTAAGAATAAGGGCAGAGGGGGGAA`
- SYN0063 `CCTTCCCCCAGCTTCCTCATCTATTTTCCCCT`

Three of the four are purine-rich (0.91–0.94) with their strong mirror broken to identity
≤ 0.78. They are therefore hard negatives of the kind reviewer R2.1 asked for, unlike the A/T
homopolymers they replace.

**Caveat for R3.2:** three of the four derive from pAA32-family constructs (SYN0063 is
length/GC-matched to pTG32, also in that family). Under the family-aware grouping, synthetic
records inherit their source's family. The SYN0060 and SYN0062 mutant/shuffle records join
fam:pAA32, already the largest group, which now has 37 of 128 records. This was not chosen: the
replacements take the next unused sources in the generator's fixed shuffled order, and pAA32
constructs make up 16 of the 37 eligible perturbation sources (43%). After DEC-8, only 4
eligible sources remain unused. The family-aware results in `WS2-A3_family_aware_validation.md` account for it,
because those records are held out together with their source family.

## 4. Validation (rerun 2026-10-07)

| Check | Result |
|---|---|
| Class balance | 64 forming / 64 non-forming |
| Exact or reverse-complement duplicates across all 128 sequences | none |
| H-DNA-competent synthetic rows (strong mirror arm ≥ 8 at identity ≥ 0.80 **and** purine ≥ 0.85; generator's own check) | `[]` |
| Maximum residual strong-mirror identity, classes A–C | 0.778 (< 0.80) |
| Pure homopolymers remaining | none |
| First 124 rows vs the DEC-7 file | byte-identical |
| HSeeker detection on the 4 new records (library defaults, direct script) | 4 / 4 true negatives |

**Synthetic composition now:** A 22, B 11, C 11, D 12 (4 G4 / 3 Z-DNA / 5 B-DNA),
E 0 (excluded), F 3 = **59**.

## 5. Effect on the reported analyses

All analyses were regenerated on the 128-record file. The rebalance moves headline numbers only
slightly, because the four new negatives are all called correctly:

| Analysis | 124 records (DEC-7) | 128 records (DEC-8) |
|---|---|---|
| Baseline, library defaults (balanced) | 62/2/60/0, MCC 0.968 | 62/2/64/0, MCC 0.969 |
| Nested CV, per-fold tuned (primary) | MCC 0.936 | MCC 0.938 |
| Nested CV, fixed script default / manuscript params | 0.968 / 0.968 | 0.969 / 0.969 |
| LOFO, per-fold tuned | 0.921 | 0.923 |
| Repeated splits: tuning beats the default | 0/200 | 0/200 |
| Filters selected across nested-CV folds | none, 100% | none, 100% |

The held-out EVAL split is class-balanced again (20 / 20, n = 40). Details are in
`WS2-PRE_hseeker_baseline_v3.md`, `WS4-A3_parameter_robustness.md` and
`WS2-A3_family_aware_validation.md`.

## Reproduce

```
python3 RTR/generate_synthetic_negatives.py   # writes hdna_benchmark_balanced_v3.csv (128 rows)
sha256sum hdna_benchmark_balanced_v3.csv      # 7bb49544…07cf9
```
