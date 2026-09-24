# DEC-2 — FXN records (HDNA0060–0063): intended citation

**Status:** evidence complete, awaiting author decision
**Records:** HDNA0060 `FXN_GAA10`, HDNA0061 `FXN_GAA20`, HDNA0062 `FXN_GAA33`, HDNA0063 `FXN_GAA66`
(all `label=forming`, `confidence_tier=literature_curated`, `disposition=secondary`, `curation_decision=kept`)
**Sequences:** pure (GAA)n tracts, n = 10, 20, 33, 66 (30, 60, 99, 198 nt)
**Verified:** 2026-09-23, PubMed E-utilities + Crossref REST

## 1. The citation currently in the CSV is wrong — confirmed two independent ways

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

## Options

| # | Option | Consequence |
|---|---|---|
| A | Cite **Wells 2008** (PMID 18211957) | Audit's candidate; real, on-topic; review, so records stay `secondary` |
| B | Cite **Rajeswari 2012** (PMID 22750988) | Best match to the recorded title; review, so records stay `secondary` |
| C | Find/supply the primary source for (GAA)10/20/33/66 | Only route to promoting these to `primary`; may not exist for n<59 |
| D | Drop the 4 records | Removes an unverifiable citation; costs 4 forming records (65→61) and forces a rebalance |

## Evidence commands

```
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=18024960,18211957,2835375&retmode=json"
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=pubmed&id=22750988&retmode=json"
curl -s "https://api.crossref.org/works/10.1074/jbc.R700013200"
```
