# DEC-2 — FXN records (HDNA0060–0063): intended citation

**Status:** **RESOLVED 2026-10-09** (author decision by Kimon; commit `8e781ca`); see "Resolution" below. Sections 1–4 and the options table are the evidence as it stood before the decision.
**Records:** HDNA0060 `FXN_GAA10`, HDNA0061 `FXN_GAA20`, HDNA0062 `FXN_GAA33`, HDNA0063 `FXN_GAA66`
(all `label=forming`, `confidence_tier=literature_curated`, `disposition=secondary`, `curation_decision=kept`)
**Sequences:** pure (GAA)n tracts, n = 10, 20, 33, 66 (30, 60, 99, 198 nt)
**Verified:** 2026-09-23, PubMed E-utilities + Crossref REST
*Counts updated 2026-10-07 for the 128-record (64:64) benchmark after DEC-8 rebalance
and DEC-9 score-informed overlap default; 124-record values in this file's 2026-10-02
version. The citation evidence itself is unchanged.*

## Resolution (2026-10-09)

Neither review candidate (options A and B) was adopted. The records now cite **primary experimental
studies**, which settles the open question in §3:

| Records | Source now recorded in `corrected_pmid` / `corrected_doi` | What it shows |
|---|---|---|
| HDNA0060–0062 ((GAA)10/20/33) | Potaman VN et al. 2004, *Nucleic Acids Res* 32(3):1224-31; PMID 14978261; DOI 10.1093/nar/gkh274 | in supercoiled plasmids at pH 7.4, (GAA)9 forms a stable intramolecular H-DNA, (GAA)23 a family of H-DNAs, and (GAA)42 a bi-triplex at higher supercoiling |
| HDNA0063 ((GAA)66) | Sakamoto N et al. 1999, *Mol Cell* 3(4):465-75; PMID 10230399; DOI 10.1016/s1097-2765(00)80474-8 | R·R·Y triplex / sticky DNA for more than 59 repeats |

- The records stay **secondary**: the evidence is class-level, and the exact constructs (pure
  (GAA)n tracts of n = 10, 20, 33, 66) were not tested in either study.
- The originally cited PMID 18024960 / DOI 10.1074/jbc.R700013200 stay in the `cited_*` fields
  as an audit trail. The recorded title matches only the review Rajeswari 2012 (PMID 22750988).
- Only `corrected_pmid`, `corrected_doi`, `resolved_origin`, `evidence` and
  `curation_justification` changed, in these 4 rows of both CSVs. Sequences, labels and
  `study_id` are unchanged, so no detection call changes. `verification_status=citation_corrected`
  is now truthful, which resolves the inconsistency in §4.
- New sha256: experimental `dfa31654559e366fd666e98c0476349d9e6cb51c4e0fe36813f8caddbae0e273`,
  balanced `c5177b2479a5da4414aefa3c6ac6251f5645868260d2b5cefe9189de182f1ae4` (previously
  `7197c6df…` and `7bb49544…`). `analysis/results/*/metadata.json` still record the previous
  balanced sha; they are refreshed at the final benchmark freeze (only metadata differs).
- The §3 open question (no primary support for the `forming` label of (GAA)10/20/33) is answered
  by Potaman et al. 2004.

## 1. The citation in the CSV before the resolution was wrong — confirmed two independent ways

The four records carry `cited_pmid=18024960` **and** `cited_doi=10.1074/jbc.R700013200`.
These do not point to the same paper, and neither is about H-DNA:

| Field in CSV | Resolves to (verified) | Related to FXN/H-DNA? |
|---|---|---|
| PMID 18024960 | Romney et al. 2008, *J Biol Chem* 283(2):716-25 — "An iron enhancer element in the FTN-1 gene directs iron-dependent expression in *Caenorhabditis elegans* intestine." DOI 10.1074/jbc.M707043200 | No |
| DOI 10.1074/jbc.R700013200 | Samuel C.E. 2007, *J Biol Chem* 282:15313-15314 — "Innate Immunity Minireview Series: Making Biochemical Sense of Nucleic Acid Sensors That Trigger Antiviral Innate Immunity" | No |

**Likely mechanism of the PMID error:** the wrong paper is about **FTN-1** (ferritin);
the records are **FXN** (frataxin). A gene-symbol confusion is the most economical
explanation. This is consistent with the systematic PMID defect in the
`literature_curated` batch (register D-10) rather than being an isolated slip.

## 2. The cited *title* points to a third paper — not the audit's candidate

CSV `cited_paper` = "DNA triplex structures in neurodegenerative and repeat expansion diseases".
No paper carries that exact title. The closest real publication:

- **Rajeswari M.R. 2012**, *J Biosci* 37:519-32, "DNA triplex structures in
  neurodegenerative disorder, Friedreich's ataxia." PMID 22750988,
  DOI 10.1007/s12038-012-9219-1 — title differs only in the tail.

The provenance audit proposed a different candidate:

- **Wells R.D. 2008**, *FASEB J* 22(6):1625-34, "DNA triplexes and Friedreich
  ataxia." PMID 18211957, DOI 10.1096/fj.07-097857 — verified real; correct topic;
  but its title is a poorer match to the recorded title than Rajeswari 2012.

**Both candidates are review articles, not primary experimental reports.**

## 3. No primary source located for these specific constructs

The canonical primary study of GAA·TTC triplex formation is Sakamoto et al. 1999,
*Mol Cell* 3:465-475 ("Sticky DNA…"), PMID 10230399. It does **not** cover these
records: it reports tracts of **75, 90, 115, 150 and 270** repeats and finds the
sticky-DNA conformation requires **>59 repeats**. The records here are n = 10, 20,
33, 66 — three of the four fall below that threshold.

**Open scientific question for the author (beyond the citation itself):** the
`forming` label for FXN_GAA10 (10 repeats), FXN_GAA20 and FXN_GAA33 is not
supported by any primary source found in this audit. Intramolecular R·R·Y triplex
at shorter tracts under negative supercoiling is plausible and distinct from the
dimeric "sticky DNA" association, but it needs its own citation. Reviewer 3.5 asks
precisely this kind of question, so it is better raised by us than by the reviewer.

## 4. Internal inconsistency in the current CSV

All four records have `verification_status=citation_corrected` while
`corrected_pmid` and `corrected_doi` are **empty** — marked as corrected, carrying
no correction. These are the only 4 records in the file in that state. Whatever
DEC-2 decides, this field must be made truthful.

## Options (as presented before the decision; the outcome is closest to C — primary sources were found — but the records stay secondary, see Resolution)

| # | Option | Consequence |
|---|---|---|
| A | Cite **Wells 2008** (PMID 18211957) | Audit's candidate; real, on-topic; review, so records stay `secondary` |
| B | Cite **Rajeswari 2012** (PMID 22750988) | Best match to the recorded title; review, so records stay `secondary` |
| C | Find/supply the primary source for (GAA)10/20/33/66 | Only route to promoting these to `primary`; may not exist for n<59 |
| D | Drop the 4 records | Removes an unverifiable citation; costs 4 forming records (64→60 forming, leaving 128→124 records at 60:64; was 65→61 before the 2026-10-02 homopolymer removal) and forces a rebalance |

## Evidence commands

```
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=18024960,18211957,2835375&retmode=json"
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=22750988&retmode=json"
curl -s "https://api.crossref.org/works/10.1074/jbc.R700013200"
```
