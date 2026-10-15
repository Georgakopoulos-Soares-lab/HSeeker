# Kouzine 2017: signal at HSeeker loci isolated from other sources of ssDNA (2026-10-11)

**Design:** `RTR/evidence/X8_kouzine_analysis_design_2026-10-09.md`, deviation 5 (written before any
of these numbers; added after the human pre-planned signal comparison had been seen; its paired
statistic corrected after a chr19-only smoke test, before any genome-wide number).
**Script:** `analysis/scripts/kouzine_hdna_isolated.py`; outputs
`analysis/results/kouzine_v1/published/<genome>.<config>.p4.json`.
**Data root:** `/scratch/10899/kimopro/kouzine` (TACC). HSeeker scans reproduce the committed locus
counts exactly (mm9 `cli` 801,282, `cli_minrep8` 852,240; hg19 508,864 / 555,753).

**Which configuration is primary.** The design's general rule makes the configuration frozen in
Phase 1.3 primary; deviation 5, written on 2026-10-09 before the freeze, named `cli` (minrep 10). The
freeze on 2026-10-11 chose minrep 8, and the author decided the same day to run **only minimum arm 8**
(`cli_minrep8`) for the signal comparisons; the `cli` (minrep 10) runs were stopped before they started.
Earlier `cli` results (recall, selectivity, and the human pre-planned comparison of 2026-10-09) stay on
record but are not part of the answer.

**Statistics.** Ratio of mean signal, loci / their own controls, with a bootstrap CI over loci;
P(sup) = mean over loci of the share of its own controls it beats (ties ½), null 0.5, with a bootstrap
CI and a one-sided t-test. Placebo: each locus's first control against its other controls (ratio
should be ≈ 1 and P(sup) ≈ 0.5). Positive control: the authors' ssDNA+ H-DNA that pass tier B.

## Mouse (mm9), activated and resting B cells, `cli_minrep8`

Funnel: 852,240 loci → 539,027 mappable (tier A) → 286,092 isolated (tier B). Positive control: 2,799
of the authors' 17,109 ssDNA+ H-DNA pass tier B. 20,000 loci sampled per tier.

**Activated B cells** (the cells the authors' calls come from):

| Loci | Controls | n | Ratio (95 % CI) | P(sup) (95 % CI), p | Placebo ratio / P(sup) |
|---|---|---|---|---|---|
| **isolated (B)** | **local, purine-matched (primary)** | 15,899 | **1.254 (1.179–1.332)** | **0.514 (0.511–0.517), 2×10⁻¹⁷** | 1.07 / 0.502 |
| isolated (B) | chromosome, purine-matched | 17,077 | 1.306 (1.236–1.374) | 0.517 (0.514–0.520), 1×10⁻²⁶ | 1.00 / 0.499 |
| isolated (B) | local, A/T-count only | 19,953 | 1.031 (0.982–1.083) | 0.511 (0.508–0.514), 8×10⁻¹⁶ | 0.96 / 0.501 |
| mappable (A) | local, purine-matched | 15,530 | 1.037 (0.967–1.113) | 0.515 (0.511–0.518), 2×10⁻¹⁶ | 1.05 / 0.501 |
| mappable (A) | local, A/T-count only | 19,722 | 0.443 (0.418–0.469) | 0.493 | 0.97 / 0.500 |
| positive control (B) | local, purine-matched | 1,839 | 7.44 (5.82–9.27) | 0.758 (0.744–0.770) | 0.86 / 0.505 |

Arm-length dependence (isolated, local, purine-matched):

| Arm (nt) | n | Ratio (95 % CI) | P(sup) |
|---|---|---|---|
| 8–9 | 1,724 | 1.12 (0.93–1.35) | 0.510 |
| 10–11 | 11,697 | 1.20 (1.11–1.30) | 0.510 |
| 12–14 | 2,049 | 1.58 (1.34–1.83) | 0.533 |
| 15–19 | 240 | 1.84 (1.02–3.28) | 0.559 |
| ≥ 20 | 189 | 1.63 (1.11–2.45) | 0.516 |

Spearman ρ (arm length vs per-locus P(sup)) = 0.075, p = 3×10⁻²¹. Stability-score quartiles: Q1 1.25,
Q2 1.04, Q3 1.28, Q4 1.40 (not monotonic).

**Resting B cells:** isolated, local, purine-matched ratio 0.858 (0.831–0.885), P(sup) 0.474; no
arm-length trend (ρ = −0.014, p = 0.08). The authors' own activated-B ssDNA+ H-DNA also show no excess
in resting cells (0.674, 0.599–0.759). Activated / resting ratio at HSeeker loci 1.46.

**Reading (against the criteria fixed in deviation 5).**
1. In activated B cells, isolated HSeeker loci carry 25 % more ssDNA signal than local controls of
   the same length, A/T content and purine asymmetry; the CI excludes 1; the placebo is calibrated
   (P(sup) 0.502) and the positive control is enriched 7.4-fold. The criterion "not explained by
   composition, mappability, local context, the other predicted non-B structures or promoter
   transcription" is met.
2. The excess grows with arm length, from no detectable excess at 8–9 nt to 1.6–1.8-fold at ≥ 12 nt,
   which is the H-DNA-specific signature the design names. G4, Z-DNA, melting and transcription do not
   depend on mirror symmetry once composition is matched.
3. The enrichment depends on cell state: it is absent in resting B cells, where even the authors'
   formed sites show no excess. This agrees with the supercoiling- and activation-dependent formation
   Kouzine et al. report.
4. The per-locus effect is small (P(sup) 0.514): most loci and controls have no signal (16 % and 13 %
   non-zero), and a minority of loci carry the excess. HSeeker loci are candidates, most of which are
   not single-stranded in a given cell state at a given time.
5. Without isolation (tier A), the mean-signal ratio is 1.04 (CI includes 1) although P(sup) is
   still above 0.5. Other sources of single-stranded DNA near promoters and other structures dominate
   the signal at both loci and controls until they are removed.
6. Controls matched on A/T count but not purine asymmetry carry far more signal (tier A ratio 0.44):
   purine/pyrimidine-asymmetric DNA has less signal than mixed DNA of the same A/T content, so matching
   on asymmetry is necessary, and the human pre-planned result matched on GC only (0.57) reflects this.
7. Association, not proof of triplex formation at a given locus; unannotated transcription, structures
   the authors' predictors miss, and the 100 bp tag extension are not controlled.

## Human (hg19), Raji cells, `cli_minrep8`

Funnel: 555,753 loci → 443,656 mappable (tier A) → 281,375 isolated (tier B). Positive control: 2,072
of the authors' 8,356 ssDNA+ H-DNA pass tier B.

| Loci | Controls | n | Ratio (95 % CI) | P(sup) (95 % CI), p | Placebo ratio / P(sup) |
|---|---|---|---|---|---|
| **isolated (B)** | **local, purine-matched (primary)** | 17,456 | **1.023 (0.988–1.061)** | **0.503 (0.499–0.507), 0.08** | 1.01 / 0.500 |
| isolated (B) | chromosome, purine-matched | 18,201 | 1.029 (0.993–1.069) | 0.504 (0.500–0.508), 0.017 | 1.02 / 0.498 |
| isolated (B) | local, A/T-count only | 19,989 | 0.905 (0.874–0.937) | 0.498 | 0.98 / 0.499 |
| mappable (A) | local, purine-matched | 17,103 | 0.932 (0.897–0.968) | 0.499, 0.6 | 1.02 / 0.502 |
| positive control (B) | local, purine-matched | 1,645 | 2.10 (1.92–2.28) | 0.724 (0.708–0.740) | 1.00 / 0.502 |
| positive control (B) | chromosome, purine-matched | 1,742 | 4.76 (4.32–5.26) | 0.813 | 1.10 / 0.516 |

Arm length (isolated, local, purine-matched): 8–9 nt 1.09 (0.98–1.20); 10–11 nt 1.00 (0.96–1.05);
12–14 nt 1.07 (0.96–1.20); 15–19 nt 1.52 (0.69–3.08, n = 141); ≥ 20 nt 1.24 (0.67–2.51, n = 39).
Spearman ρ = 0.033, p = 1.5×10⁻⁵. Score quartiles 1.02 / 1.01 / 0.98 / 1.07.

**Reading.** In Raji cells the isolated HSeeker loci carry no detectable excess of single-stranded signal
over matched controls (1.02, CI includes 1; P(sup) 0.503, p = 0.08). The pipeline can detect
enrichment in these data (positive control 2.1-fold locally), and the controls are calibrated (placebo
0.500). The weak arm-length correlation is statistically detectable at this sample size but no arm-length
class differs from 1. The criterion fixed in deviation 5 is **not met in human**.

Differences from mouse that are facts, not explanations: Raji has no input control; the signal is far
less sparse (38 % of isolated loci and 37 % of controls have non-zero signal, against 16 % and 13 % in
mouse activated B cells); and the authors' own formed sites are enriched 2.1-fold over matched local
controls in Raji against 7.4-fold in mouse, so the dynamic range of the human data is smaller.

## Pre-planned comparison (all loci, genome-wide pool controls), `cli_minrep8`

| | Purine-matched controls: ratio (95 % CI) | GC-matched only |
|---|---|---|
| Mouse, activated B (825,749 loci matched) | **1.137 (1.127–1.147)** | 0.405 |
| Mouse, resting B | 0.822 (0.818–0.827) | 0.540 |
| Human, Raji (550,935 loci matched) | 0.977 (0.969–0.984) | 0.585 |

Mouse activated / resting (purine-matched) 1.38. The pre-specified Mann–Whitney test (one-sided,
loci > controls) gives p ≈ 0 for mouse activated and p ≈ 1 for human. Control sampling: mouse 26,478
loci without an asymmetry-matched cell (20,374 widened by 1, 11,000 by 2); human 4,818 unmatched.

## Overall reading (minimum arm 8)

- **Mouse activated B cells: enrichment defended.** Both the pre-planned comparison (1.14) and the
  isolated comparison (1.25) are above 1 with CIs excluding 1; the excess grows with arm length; it is
  absent in resting B cells; placebo and positive control behave as required.
- **Human Raji cells: no enrichment beyond composition.** Pre-planned 0.98, isolated 1.02 (CI includes
  1). This is reported as such; recall (82.3 %) and selectivity (odds ratio 3.39; 4.16–5.70 on
  mappable motifs) remain strong in human, because they compare HSeeker's calls with the authors'
  other predicted mirror repeats rather than with non-mirror DNA.
- Recall and selectivity were re-computed on TACC and reproduce the committed values exactly (mouse
  79.0 %, odds ratio 1.65; human 82.3 %, 3.39).
