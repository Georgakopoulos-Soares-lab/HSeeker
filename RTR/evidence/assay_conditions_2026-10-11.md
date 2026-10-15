# Assay conditions per benchmark record (promised in the R3.4 answer) — 2026-10-11

**Done:** both benchmark CSVs now carry, for every experimental record, `assay_method`,
`assay_dna_context`, `assay_ph`, `assay_temperature`, `assay_ionic_conditions`, `assay_supercoiling`,
plus `assay_conditions_source`, `assay_conditions_source_type`, `assay_conditions_confidence`,
`assay_conditions_quote` and `assay_conditions_note`. Designed negatives read
"not applicable (designed negative)". No sequence, label or call changed.

**Sources.**
- 53 records (42 of them kept): the condition fields of the original supplementary curation
  (`hdna_experimental_sequences_final.csv`: method, DNA context, pH, temperature, ions, supercoiling,
  with its source location and paper URL).
- 28 records (all kept; HDNA0027–0028, HDNA0055–0079, HDNA0081), which the original curation left
  without pH/ions/supercoiling: curated on 2026-10-11 from the cited papers (full text where open access,
  otherwise the abstract, or the primary study a cited review relies on), with a verbatim quote, source
  location, source type and confidence. Raw curation: `assay_conditions_curation_2026-10-11.json`.

**Coverage among the 70 kept records:** pH reported for 64; "not reported" for HDNA0055 (cited paper is a
recombination assay), HDNA0058 (full text not accessible), HDNA0065 (no experiment exists), HDNA0069 and
HDNA0070 (cited studies are polymerase/replication assays without structural probing), HDNA0072 (full
text not accessible; abstract says "strongly pH-dependent").

**Source types of the 28 curated records:** full text 15, abstract only 6, primary study behind a cited
review 6, not found 1 (HDNA0065).

## Curation findings to decide on (labels and sequences were not changed)

| Record | Finding |
|---|---|
| HDNA0065 XDP (CCCTCT)10 | The cited review states that no in vitro or in vivo data exist; label inferred from sequence. Kept in the benchmark (author decision 2026-10-11) but marked `primary`. |
| HDNA0055 GA32, HDNA0056 AA32 | Cited paper (Rooney & Moore 1995) is a recombination assay; the H-form is cited from Mirkin 1987. GA32 forms H-DNA "with difficulty". |
| HDNA0058 cMYC NSE (ACCCTCCCC)4 | Full text not accessible; also a tandem direct repeat (one of the two false negatives). |
| HDNA0059 BCL6 | Triplex shown on oligonucleotides at pH 8 with 10 mM Mg²⁺, not on supercoiled DNA. |
| HDNA0060–0062 FXN (GAA)n | Class-level match; the study's no-flank control needed σ ≤ −0.09. |
| HDNA0063 FXN (GAA)66 | "Sticky DNA" (purine-motif triplex dimer), not classic H-DNA; abstract only. |
| HDNA0064 CANVAS (AAGGG)10 | Primary study (Hisey 2024 NAR) shows a purine-motif H-r3 triplex. |
| HDNA0066–0068 mouse S1-seq | S1 at pH 4.5 with Zn²⁺; the authors cannot exclude formation after lysis. |
| HDNA0069 AG11 | Cited primary study measures polymerase error rates; no structural assay. |
| HDNA0070 GA20 SV40 | Primary construct is (GA)27 with flanks; replication-fork pausing in vivo, no structural assay. |
| HDNA0071 PKD1 CT44 | No standalone 44-nt construct in the paper; 44 bp is an overlapping repeat inside a P1-sensitive region. |
| HDNA0072 GA16 sea urchin | Mirkin 1987 abstract does not mention (GA)16 or a sea-urchin histone gene; full text not accessible. |
| HDNA0075 SV40 GA·CT22 | The structure shown is *H-DNA (Pu·Pu·Py) induced by 0.4 mM Zn²⁺ at neutral pH. |
| HDNA0077 TTCCC avian | Construct is (TTCCC)48, not (TTCCC)20; supercoiling not required. |
| HDNA0078 γ-globin | Conditions from the primary study behind the cited paper (acidic pH, supercoiling). |

Most of these are secondary records already. The ones marked `primary` that rest on inference or on a
different construct (HDNA0065, HDNA0055/0056, HDNA0069, HDNA0070, HDNA0071, HDNA0072) are candidates for
moving to the secondary set; that is an author decision.

**New sha256 (metadata columns only):** experimental `bcaee4ef…e801a0`, balanced `c3f3b842…a9ed3e`.
Results computed earlier on 2026-10-11 used `02cdd769…985ea20b`; every call is identical.
