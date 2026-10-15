# DEC-3 — add the Hanvey 1988 (TTC)8 insert as a new primary record?

**Status:** evidence complete, awaiting author decision
**Source:** Hanvey J.C., Klysik J., Wells R.D. 1988, *J Biol Chem* 263(15):7386-96,
"Influence of DNA sequence on the formation of non-B right-handed helices in
oligopurine.oligopyrimidine inserts in plasmids." PMID **2835375** (verified)
**Verified:** 2026-09-23, PubMed E-utilities `efetch` (full abstract retrieved)
*Counts updated 2026-10-07 for the 128-record (64:64) benchmark after DEC-8 rebalance
and DEC-9 score-informed overlap default; 124-record values in this file's 2026-10-02
version. The primary-source evidence is unchanged.*

## 1. The primary abstract resolves both the sequence and the label

Verbatim from the abstract:

> "A systematic study was conducted on **seven** recombinant plasmids harboring
> synthetic inserts… The inserts ranged in G+C content from 100% [G19.C19] to 0%
> [A20.T20] with intermediate contents at 66% [(TCC)8.(GGA)8], 50% [(CT)12.(AG)12
> and (TTCC)6.(GGAA)6], **33% [(TTC)8.(GAA)8]**, and 25% [(GAAA)6.(TTTC)6]."

> "We conclude that the **five inserts with 66-25% G+C adopt a non-B right-handed
> conformation** which is stabilized by negative supercoiling… An **intramolecular
> triple-stranded model** for the unusual structure of the insert accounts most
> favorably for these observations."

> "Unexpectedly, the **A20.T20 insert seems to remain in an orthodox right-handed
> B-conformation** under all conditions tested. The **G19.C19 insert does adopt a
> non-B right-handed structure** as for the five inserts with 66-25% G+C, **but the
> pattern of reactivities and hence its conformation is different**."

The five inserts at 66–25% G+C are (TCC)8, (CT)12, (TTCC)6, **(TTC)8**, (GAAA)6.
So **(TTC)8·(GAA)8 is one of the five that form the intramolecular triplex**:
sequence and label are both established by the primary source.

## 2. Current coverage of this study in the benchmark — 4 of 7 inserts

*Updated 2026-10-02:* was 5 of 7; HDNA0053 `pRW1405` (A20·T20) was removed from the
benchmark with the other pure homopolymers (author decision 2026-10-02).

| Insert (abstract) | G+C | Abstract verdict | In benchmark? |
|---|---|---|---|
| (TCC)8·(GGA)8 | 66% | non-B / triplex | HDNA0050 `pRW1411` — forming ✓ |
| (CT)12·(AG)12 | 50% | non-B / triplex | HDNA0049 `pRW1404` — forming ✓ |
| (TTCC)6·(GGAA)6 | 50% | non-B / triplex | HDNA0051 `pRW1408` — forming ✓ |
| **(TTC)8·(GAA)8** | **33%** | **non-B / triplex** | **absent — this is DEC-3** |
| (GAAA)6·(TTTC)6 | 25% | non-B / triplex | HDNA0052 `pRW1410` — forming ✓ |
| A20·T20 | 0% | stays B-form | HDNA0053 `pRW1405` — non_forming; **removed 2026-10-02** (pure homopolymer) |
| G19·C19 | 100% | non-B but *different* conformation | absent — see §4 |

Every record's label (including the removed HDNA0053) agrees with the abstract. The benchmark is missing
exactly one unambiguous forming insert from this study.

## 3. The plasmid *name* is not verifiable; the sequence and label are

The abstract names only the vector (`pRW790`), never the individual construct
numbers. The `pRW####` → insert mapping comes from the same secondary
"H-DNA Thermostability document" that was **already shown to be unreliable for this
exact study** — it listed pRW1406 with the (TCC)8 sequence, the conflation the audit
resolved in favour of pRW1411=(TCC)8 (plan §A.2).

So adding the record **as `pRW1406`** would re-import the one mapping known to be
corrupted, to gain nothing: the sequence and label stand on the abstract alone.

**Recommended form** — name it from the verified content, not the unverified
plasmid number:

```
sequence_name : Hanvey1988_(TTC)8
sequence_5to3 : TTCTTCTTCTTCTTCTTCTTCTTC      (24 nt, 33% G+C)
label         : forming
label_scope   : sequence_specific
disposition   : primary          # construct-level, primary-source verified
confidence_tier: high
study_id      : 10.1016/S0021-9258(18)68654-5
cited_pmid    : 2835375
verification_status: verified_primary_abstract
```

Confirmed absent from the benchmark: no record contains `(TTC)8` or `(GAA)8` as an
exact sequence (checked across all 80 rows, exact and reverse-complement).

## 4. Second candidate the plan did not raise: G19·C19

G19·C19 is also missing and is **not** a clean addition: the abstract says it adopts a
non-B structure *with a different conformation*. Adding it as `forming` would repeat
the HDNA0080 class of error (label scope not matching the tool's target class).
Recommend **not** adding it, or adding it only to the secondary set with the caveat
recorded. Noted here so the decision is explicit rather than an omission.

## 5. Downstream consequence of adding one forming record

| | Now | If (TTC)8 added |
|---|---|---|
| Kept experimental | 69 (64 forming / 5 non-forming) | 70 (65 / 5) |
| Primary / secondary | 57 / 12 | 58 / 12 |
| Balanced benchmark | 128 rows, 64:64 | 129 rows, 65:64 if the synthetic set is left as is; restoring 1:1 needs 130 rows, 65:65 (generator output not recomputed) |
| Synthetic negatives | 59 | 59 if left as is; 60 needed for 1:1 (which records the generator would draw: not recomputed) |

*Updated 2026-10-02:* "Now" column was 71 (65/6), 59/12, 130 rows 65:65, 59
synthetic; changed after the 6 pure homopolymers were removed (HDNA0053/HDNA0054 and
SYN0053–SYN0056). The balanced file was then no longer 1:1, and the generator now drops
pure homopolymers after drawing the pool, so the "if added" balanced/synthetic counts
were not re-derived.

*Updated 2026-10-07:* "Now" column changed from 124 rows 64:60 / 55 synthetic to
128 rows 64:64 / 59 synthetic after DEC-8 added four replacement synthetic negatives
(SYN0060–SYN0063). The experimental counts (69 = 64/5; 57/12) are unchanged. The
"if added" column is now plain arithmetic: one extra forming record on the current
file gives 129 rows at 65:64; a 1:1 balance would need one more synthetic negative
(60, 130 rows at 65:65). The generator sizes its pool from the experimental class
counts, so re-running it with the new record would change the draw; the exact
records it would produce were not recomputed.

The generator is deterministic and was verified to reproduce the current balanced
file byte-for-byte, so the rebalance is mechanical: re-run
`RTR/generate_synthetic_negatives.py` and re-record the SHA-256.

## Options

| # | Option | Effect |
|---|---|---|
| A | Add as `Hanvey1988_(TTC)8`, disposition `primary` (**recommended**) | +1 primary forming record with full primary-source backing; avoids the corrupted plasmid mapping |
| B | Add as `pRW1406` | Same sequence/label, but re-imports the mapping the audit found conflated |
| C | Do not add | Benchmark stays at 69/128 (was 69/124 on 2026-10-02 and 71/130 before that); a verified primary positive is left on the table |

## Evidence command

```
curl -s "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&id=2835375&rettype=abstract&retmode=text"
```
