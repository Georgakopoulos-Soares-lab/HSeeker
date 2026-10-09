# X8 results: HSeeker vs the authors' published Kouzine 2017 calls

**Design:** `RTR/evidence/X8_kouzine_analysis_design_2026-10-09.md`, deviation 4. Tests P1–P3 were
fixed before any result was computed.
**Scripts:**
- `analysis/scripts/run_kouzine_published.py` (P1–P3)
- `analysis/scripts/kouzine_missed_hdna.py` (missed-motif analysis, exploratory and **not**
  pre-specified)

**Outputs:** `analysis/results/kouzine_v1/published/{mm9.cli.json, mm9_missed_analysis.json}`.
**HSeeker configuration:** shipped CLI defaults (minrep 10, maxrep 1000, maxspacer 10, purity 0.90,
mismatch 0.10, AT filter 0.80, homopolymer filter on, greedy overlap removal). The mm9 scan gives
801,282 loci on chromosomes 1–19, X and Y.

| Item | Status (2026-10-09) |
|---|---|
| mm9 (mouse activated B), P1 + P2 | done, below |
| mm9 P3 | pending: composition-matched control sampling is running |
| hg19 (human Raji), P1–P3 | pending: control sampling is running |
| mm9 benchmark-configuration sensitivity scan (minrep 8, no filters) | scan running |
| mm10 / hg38 reprocessing of the raw reads (secondary) | download and alignment running |

## P1: recall of the authors' ssDNA+ H-DNA (mm9, activated B cells)

12,676 of the 17,109 H-DNA motifs the authors called as ssDNA+ (formed) are overlapped by an HSeeker
locus: **74.1 %**. The count of 17,109 equals the paper's 16,876 + 233 (Table 1).

## P2: selectivity among the authors' predicted H-DNA motifs (mm9)

- **Universe.** The authors' predicted H-DNA motifs, restricted to those with no SINE within a 500 bp
  window centred on the motif (their filter, re-applied from UCSC rmsk): 428,237 motifs. The paper
  analysed 341,431. Their mappability filter was not reproduced, so our universe is larger.

| | ssDNA+ | ssDNA− | ssDNA+ rate |
|---|---|---|---|
| HSeeker calls the motif | 12,676 | 265,141 | 4.56 % |
| HSeeker does not call it | 4,433 | 145,987 | 2.95 % |

- **Odds ratio:** 1.57 (95 % CI 1.52–1.63); one-sided Fisher p = 1.4 × 10⁻¹⁵³.
- **Stability score.** Within HSeeker-called motifs, the score does **not** separate ssDNA+ from
  ssDNA− motifs: median 161.45 vs 172.93, one-sided Mann–Whitney p = 0.93. This bears on R3.3 and
  R3.4, and must be reported. The score ranks thermodynamic stability of the predicted triplex; it
  is not validated as a predictor of formation in vivo.

## Exploratory: which ssDNA+ H-DNA motifs does HSeeker miss, and why? (4,433 motifs, 25.9 %)

**Method.** Each missed motif is first re-scanned with ±50 bp flanks under the CLI settings. If
that recovers it, the reason is overlap removal; if not, the setting is relaxed one step at a time
until a hit overlaps the motif. The steps are cumulative, applied in this fixed order: AT filter
off → homopolymer filter off → minrep 8 → minrep 6 → mismatch 0.20 → purity 0.80. The first step
that recovers the motif is recorded. The attribution depends on that order, and composition
filters come first.

| First relaxation that recovers the motif | n | Share of missed | HSeeker locus within 100 bp |
|---|---|---|---|
| AT filter (motif ≥ 80 % A/T) | 1,819 | 41.0 % | 15 % |
| arm 8–9 bp (below minrep 10) | 1,191 | 26.9 % | 19 % |
| homopolymer filter | 686 | 15.5 % | 7 % |
| overlap removal (a longer neighbouring locus was kept) | 463 | 10.4 % | 97 % |
| mirror mismatch 10–20 % | 235 | 5.3 % | 92 % |
| arm 6–7 bp | 20 | 0.5 % | 25 % |
| not recovered by any relaxation | 18 | 0.4 % | 100 % |
| purity 80–90 % | 1 | 0.0 % | 0 % |

**Overlap-removal group checked directly on 5 examples.** With ±300 bp of context, the window scan
reproduces the genome loci exactly. Greedy removal keeps a longer candidate next to the motif
(e.g. an arm of 63 bp starting 79 bp away). That candidate overlaps the motif's own candidate but
not the motif itself. These are boundary disagreements inside long purine/pyrimidine tracts that
HSeeker does call. The same holds for most of the mismatch and unrecovered groups (92 % and 100 %
have an HSeeker locus within 100 bp).

**Missed motifs differ in composition** (medians unless stated):

| | caught (12,676) | missed (4,433) |
|---|---|---|
| length | 28 bp | 20 bp |
| A/T fraction | 0.43 | 0.87 |
| share with A/T ≥ 0.8 | 3 % | 62 % |
| longest single-base run / length | 0.09 | 0.32 |
| 1-mer periodic (≥ 90 % self-match at shift 1) | 0.9 % | 24 % |
| 2-, 3-, 4-mer periodic | 13 %, 12 %, 7 % | 2 %, 2 %, 2 % |

**Overlap with the authors' own SIDD calls** (stress-induced duplex destabilisation, i.e.
melting-prone DNA), from their `nonB_DNA_predicted/SIDD` and `*_ssDNA_enriched_SIDD` files:

| Group | n | overlaps a predicted SIDD motif | overlaps an ssDNA+ SIDD call |
|---|---|---|---|
| caught, A/T < 0.8 | 12,215 | 3.5 % | 2.6 % |
| caught, A/T ≥ 0.8 | 360 | 36.7 % | 30.6 % |
| caught, homopolymer | 101 | 7.9 % | 5.9 % |
| missed, A/T < 0.8 | 1,578 | 4.9 % | 3.8 % |
| missed, A/T ≥ 0.8 | 1,966 | 50.3 % | 44.1 % |
| missed, homopolymer | 889 | 35.4 % | 30.5 % |

**Reading.**
1. The misses are a sequence class, not random. Most are short (≈ 20 bp), A/T-rich or
   single-base motifs that HSeeker's composition filters exclude by design.
2. About half of the A/T-rich misses and a third of the homopolymer misses coincide with an ssDNA+
   SIDD call, against 2.6 % of the typical caught motif. Their single-stranded signal therefore has
   a competing explanation: duplex melting of A/T-rich DNA under supercoiling. The KMnO4 chemistry
   also reports unpaired thymines, so T-rich DNA is the most readily detected. This supports the
   filters' rationale (R3.5) but does not prove that these motifs never form triplexes.
3. About 27 % of misses (1,191 motifs, 7.0 % of all ssDNA+ H-DNA) have 8–9 bp arms. This is direct
   input to the open reporting-configuration decision (shipped minrep 10 vs 8). The benchmark scan
   (minrep 8, no filters) is running and gives the recall without the cumulative-ladder caveat.
4. About 15 % of misses (684 motifs, almost all in the overlap-removal, mismatch and unrecovered
   groups) sit within 100 bp of an HSeeker locus in the same tract. They are boundary
   disagreements, not missed tracts.

**Cannot claim:**
- that the A/T-rich misses are false positives in the authors' data;
- that HSeeker's score predicts formation;
- anything about P3 or the human data before those runs finish.
