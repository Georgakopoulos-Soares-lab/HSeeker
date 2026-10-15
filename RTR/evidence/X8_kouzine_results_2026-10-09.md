# X8 results: HSeeker vs the authors' published Kouzine 2017 calls

**Design:** `RTR/evidence/X8_kouzine_analysis_design_2026-10-09.md`, deviation 4. Tests P1–P3 were
fixed before any result was computed. Everything marked *sensitivity* or *exploratory* below was
added afterwards and is not pre-specified.

**Scripts:**
- `analysis/scripts/run_kouzine_published.py`: P1–P3, plus the P2 mappability sensitivity.
- `analysis/scripts/kouzine_missed_hdna.py`: the missed-motif analysis and the AT stratification
  (exploratory).

**Outputs:** `analysis/results/kouzine_v1/published/<genome>.<config>.json` and
`<genome>.<config>.missed.json`.

**HSeeker configurations.** All three use maxrep 1000, maxspacer 10, purity 0.90, mismatch 0.10
and greedy overlap removal.

| Name | minrep | AT filter 0.80 | Homopolymer filter | Role |
|---|---|---|---|---|
| `cli` | 10 | on | on | shipped CLI defaults |
| `cli_minrep8` | 8 | on | on | candidate default (Phase 1.3 minrep decision) |
| `benchmark` | 8 | off | off | benchmark/script configuration; mm9 P1/P2 only |

**Status (2026-10-09):**
- mm9 (mouse activated B) and hg19 (human Raji): P1 and P2 are done for `cli` and `cli_minrep8`,
  and for `benchmark` on mm9.
- P3 waits on composition-matched control sampling. The `cli` sampling is running; the
  `cli_minrep8` sampling restarts once memory allows.
- The mm10 / hg38 reprocessing of the raw reads (secondary) is still downloading and aligning.

## P1: recall of the authors' ssDNA+ H-DNA

| | mm9, 17,109 ssDNA+ (= paper 16,876 + 233) | hg19, 8,356 ssDNA+ (= paper 8,008 + 348) |
|---|---|---|
| `cli` (minrep 10) | 12,676 = **74.1 %** (801,282 loci) | 5,679 = **68.0 %** (508,864 loci) |
| `cli_minrep8` | 13,517 = **79.0 %** (852,240 loci, +6.4 %) | 6,879 = **82.3 %** (555,753 loci, +9.2 %) |
| `benchmark` (minrep 8, no filters) | 16,375 = 95.7 % (1,285,045 loci) | — |

## P2: selectivity among the authors' predicted H-DNA motifs

**Universe.** The authors' predicted H-DNA motifs with no SINE within a 500 bp window centred on the
motif (their filter, re-applied from UCSC rmsk):
- mm9: 428,237 motifs, against the paper's 341,431;
- hg19: 206,902 motifs, against the paper's 158,370.

Their mappability filter is not reproduced in this primary universe.

| | HSeeker-called: ssDNA+ rate | not called: ssDNA+ rate | Odds ratio (95 % CI) | Score, ssDNA+ vs ssDNA− (median; one-sided MW p) |
|---|---|---|---|---|
| mm9 `cli` | 4.56 % | 2.95 % | 1.57 (1.52–1.63) | 161.5 vs 172.9; p = 0.93 |
| mm9 `cli_minrep8` | 4.51 % | 2.79 % | 1.65 (1.59–1.71) | 150.2 vs 156.0; p = 0.03 |
| mm9 `benchmark` | 3.92 % | 7.04 % | 0.54 (0.50–0.58) | 122.3 vs 109.8; p = 1e-87 |
| hg19 `cli` | 5.59 % | 2.54 % | 2.27 (2.16–2.38) | 101.6 vs 100.9; p = 3e-5 |
| hg19 `cli_minrep8` | 5.65 % | 1.73 % | 3.39 (3.20–3.59) | 97.1 vs 93.3; p = 4e-23 |

**Sensitivity: mappable motifs only.** Unmappable motifs cannot be ssDNA+, whatever their
structure. Motifs are kept if the mean ENCODE CRG 36-mer mappability over the start positions of
every 36 bp read overlapping the motif (the reads are 36 bp) is ≥ 0.5, or ≥ 0.999 ("unique").

| | mappability ≥ 0.5: odds ratio; score p | unique: odds ratio; score p |
|---|---|---|
| mm9 `cli` | 1.62 (1.55–1.70); 107.9 vs 101.6, p = 4e-25 | 1.87 (1.72–2.05); p = 6e-32 |
| mm9 `cli_minrep8` | 1.68 (1.60–1.76); p = 5e-43 | 2.07 (1.87–2.28); p = 3e-52 |
| mm9 `benchmark` | 0.54 (0.48–0.61) | 0.61 (0.42–0.89) |
| hg19 `cli` | 2.64 (2.51–2.77); p = 3e-44 | 2.75 (2.58–2.93); p = 5e-146 |
| hg19 `cli_minrep8` | 4.16 (3.91–4.43); p = 5e-94 | 5.70 (5.18–6.27); p = 6e-258 |

Universe sizes, with the authors' ssDNA+ motifs retained:

| | mm9 motifs | mm9 ssDNA+ | hg19 motifs | hg19 ssDNA+ |
|---|---|---|---|---|
| SINE filter only | 428,237 | 17,109 | 206,902 | 8,356 |
| mappability ≥ 0.5 | 238,998 | 8,147 | 169,963 | 7,188 |
| unique | 88,267 | 2,110 | 83,971 | 4,278 |

In hg19 the ≥ 0.5 universe (169,963) is close to the paper's 158,370. In mm9 the proxy is stricter
than the authors' filter: half of their ssDNA+ motifs fall below 0.5. This is plausibly because
their 100 bp tag extension lets reads from mappable flanks cover a poorly mappable motif.

**Reading P2.**
1. With the filters on, HSeeker is selective in both species: a motif it calls is 1.6–5.7 times as
   likely (odds) to be ssDNA+, across every universe.
2. **minrep 8 adds recall without diluting precision.** The per-call ssDNA+ rate is unchanged
   (mm9 4.56 → 4.51 %, hg19 5.59 → 5.65 %), and the odds ratio rises in every universe.
3. **Without the filters, HSeeker loses selectivity.** It calls 97.6 % of the authors' mm9 motifs,
   the per-call rate drops to 3.92 %, and the odds ratio falls below 1.
4. **Score, corrected.** An earlier version of this file said the score does not discriminate
   (mm9, primary universe, p = 0.93). That null result does not survive the mappability
   restriction: on mappable motifs, ssDNA+ motifs score higher in mm9 (p = 4e-25 to 3e-52), as they
   do in hg19 in every universe. A likely cause is that long, high-scoring repeats are also poorly
   mappable, so they read as ssDNA−. The effect is small (on mappable motifs, medians differ by about 5–9 %). The
   pre-specified primary result for mm9 stays non-significant and must be reported as such,
   together with this sensitivity analysis.

## Exploratory: which ssDNA+ H-DNA motifs does HSeeker miss, and why?

**Method.** Each missed motif is first re-scanned with ±50 bp flanks under the same configuration.
If that recovers it, the reason is overlap removal. Otherwise settings are relaxed one step at a
time, cumulatively, in a fixed order, and the first step that recovers the motif is recorded:
1. AT filter off;
2. homopolymer filter off;
3. minrep 8 (`cli` only);
4. minrep 6;
5. mismatch 0.20;
6. purity 0.80.

Attribution depends on this order. With `cli`, an A/T-rich motif with 8-bp arms is attributed to
the arm length; with `cli_minrep8`, the same motif is attributed to the AT filter.

| First relaxation that recovers the motif | mm9 `cli` | mm9 `cli_minrep8` | hg19 `cli` | hg19 `cli_minrep8` |
|---|---|---|---|---|
| AT filter (motif ≥ 80 % A/T) | 1,819 (41.0 %) | 1,984 (55.2 %) | 709 (26.5 %) | 778 (52.7 %) |
| homopolymer filter | 686 (15.5 %) | 879 (24.5 %) | 435 (16.2 %) | 553 (37.4 %) |
| arm 8–9 bp (below minrep 10) | 1,191 (26.9 %) | — | 1,383 (51.7 %) | — |
| overlap removal (a longer neighbouring locus kept) | 463 (10.4 %) | 456 (12.7 %) | 48 (1.8 %) | 44 (3.0 %) |
| mirror mismatch 10–20 % | 235 (5.3 %) | 234 (6.5 %) | 76 (2.8 %) | 76 (5.1 %) |
| arm 6–7 bp | 20 | 20 | 21 | 21 |
| not recovered by any relaxation | 18 | 18 | 3 | 3 |
| purity 80–90 % | 1 | 1 | 2 | 2 |
| **missed, total** | **4,433 (25.9 %)** | **3,592 (21.0 %)** | **2,677 (32.0 %)** | **1,477 (17.7 %)** |

**Overlap-removal and mismatch groups.** 97–98 % of the overlap-removal misses have an HSeeker
locus within 100 bp, as do 72–92 % of the mismatch misses. Five overlap-removal cases were checked
directly. With ±300 bp of context, the window scan reproduces the genome loci exactly: greedy
removal keeps a longer neighbouring candidate, for example an arm of 63 bp starting 79 bp away,
which overlaps the motif's own candidate but not the motif. These are boundary disagreements
inside tracts HSeeker does call.

**At minrep 8, the composition filters explain almost all misses:** 80 % in mm9 and 90 % in hg19.

Missed motifs at `cli_minrep8` vs caught motifs (medians unless stated):

| | mm9 caught | mm9 missed | hg19 caught | hg19 missed |
|---|---|---|---|---|
| length | 27 bp | 21 bp | 21 bp | 20 bp |
| A/T fraction | 0.42 | 0.88 | 0.36 | 0.92 |
| share ≥ 80 % A/T | 3 % | 76 % | 2 % | 89 % |
| 1-mer periodic | 0.9 % | 30 % | 0.6 % | 45 % |

## Exploratory: are the AT-rich motifs that HSeeker filters less likely to form?

This test uses the authors' analysed universe, independent of HSeeker. Their predicted H-DNA
motifs are split by composition, and the ssDNA+ (formed) rate is compared with that of the
A/T < 0.8 motifs. Single-base motifs are their own group; "A/T ≥ 0.8" excludes them.

| | A/T < 0.8: ssDNA+ rate | A/T ≥ 0.8: rate; OR (95 % CI) | single-base: rate; OR (95 % CI) |
|---|---|---|---|
| mm9, SINE filter only | 4.69 % | 2.59 %; 0.54 (0.52–0.56) | 2.24 %; 0.47 (0.44–0.50) |
| mm9, mappability ≥ 0.5 | 4.25 % | 2.49 %; 0.58 (0.55–0.61) | 2.08 %; 0.48 (0.44–0.52) |
| mm9, unique | 3.07 % | 1.41 %; 0.45 (0.40–0.51) | 1.66 %; 0.53 (0.47–0.61) |
| hg19, SINE filter only | 5.83 % | 1.53 %; 0.25 (0.23–0.27) | 1.85 %; 0.30 (0.28–0.33) |
| hg19, mappability ≥ 0.5 | 6.62 % | 1.49 %; 0.21 (0.20–0.23) | 1.61 %; 0.23 (0.21–0.25) |
| hg19, unique | 8.01 % | 1.26 %; 0.15 (0.13–0.17) | 1.42 %; 0.17 (0.14–0.19) |

**Share of the authors' ssDNA+ H-DNA that the filters target** (SINE-only universe):
- mm9: 13.6 % A/T ≥ 0.8 + 5.8 % single-base = 19.4 %;
- hg19: 10.2 % + 7.3 % = 17.5 %.

**Melting (SIDD) as a partial explanation, tested within composition.** SIDD means
stress-induced duplex destabilisation. The comparison is ssDNA+ odds for motifs that overlap one
of the authors' *predicted* SIDD motifs (sequence-based) vs those that do not, within the same
composition group:

| | mm9: SINE only / ≥ 0.5 / unique | hg19: SINE only / ≥ 0.5 / unique |
|---|---|---|
| A/T ≥ 0.8 | 1.35 (1.24–1.46) / 1.58 (1.44–1.74) / 2.00 (1.59–2.50) | 1.54 (1.34–1.77) / 1.65 (1.42–1.92) / 1.35 (1.04–1.74) |
| single-base | 0.85 (0.74–0.97) / 0.86 / 0.86 (0.63–1.18) | 1.80 (1.53–2.12) / 2.01 / 1.64 (1.18–2.26) |
| A/T < 0.8 | 0.57 (0.52–0.62) | 0.44 (0.38–0.51) |

**Correction.** An earlier version of this file compared the SIDD overlap of missed A/T-rich motifs
(44–50 %) with that of caught motifs (2.6 %). That comparison is confounded by composition: 41 % of
*all* A/T-rich predicted motifs in mm9 overlap a SIDD motif, whether formed or not. The
within-composition odds ratios above replace it.

**Reading.**
1. HSeeker's remaining misses are a sequence class: short, A/T-rich or single-base motifs that
   the AT and homopolymer filters exclude by design.
2. In the authors' own data, that class is ssDNA+ at about half the rate of other predicted H-DNA
   in mouse (OR 0.45–0.58) and a quarter or less in human (OR 0.15–0.30). This holds in all three
   universes, so mappability does not explain it.
3. This is the first data-based support for the filters. The benchmark cannot evaluate them
   (R2.6): they remove the predicted-H-DNA class least likely to show formation signal. Removing
   them raises recall to 95.7 % but costs all selectivity (P2, `benchmark`).
4. The filters are not free. A/T-rich and single-base motifs are 17–19 % of the authors' ssDNA+
   H-DNA, and HSeeker misses nearly all of them.
5. A/T-rich motifs that sit in predicted melting regions are ssDNA+ more often (OR 1.35–2.00).
   Melting under supercoiling probably accounts for part of the signal at these motifs. It is a
   modest effect and not the main explanation.

**Cannot claim:**
- that the A/T-rich ssDNA+ motifs are false positives in the authors' data, or never form
  triplexes;
- that HSeeker's score is a strong predictor of formation (the effect is small and the mm9 primary
  test is null);
- anything about P3 before the controls finish.
