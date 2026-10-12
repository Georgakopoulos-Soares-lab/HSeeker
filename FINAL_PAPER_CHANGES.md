# HSeeker — final revision package

Everything needed to produce the revised manuscript and the response to reviewers. Text in
quoted blocks is final copy, written in the voice of the paper, and can be pasted directly.

Three things changed that touch almost every number in the paper, so they are stated once here
and not repeated in each section:

1. **The benchmark is now 130 sequences, 65 H-DNA-forming and 65 non-forming** (70 experimental
   records and 60 designed negatives), replacing the 79-sequence set at 70:9.
2. **The default minimum arm length is 8 bp and the default mismatch tolerance is 0.10**, and the
   AT-content and homopolymer filters are applied by the Python package, the command-line tool and
   the web application alike. The three interfaces now return identical results for identical input.
3. **The *E. coli* insertion experiment and the Youden-selected score threshold are removed.**
   Figure 3 is replaced by a comparison against permanganate/S1 single-stranded DNA maps.

### Reference values used throughout

| Quantity | Value |
|---|---|
| Benchmark | 130 records (65 forming / 65 non-forming); 70 experimental (65/5), 60 designed negatives |
| Experimental records, construct-level | 58 (53 forming / 5 non-forming); 12 further records are class-level |
| HSeeker on the benchmark | TP 63, FN 2, TN 65, FP 0 |
| Sensitivity | 0.969 (95% CI 0.895–0.992) |
| Specificity | 1.000 (0.944–1.000) |
| PPV / NPV | 1.000 (0.943–1.000) / 0.970 (0.898–0.992) |
| MCC | 0.970 (0.926–1.000) |
| Leave-one-family-out / leave-one-study-out MCC | 0.970 / 0.970 |
| Chromosome 1, default settings, 16 threads | 45,927 loci in 16.3 s, peak 1,334 MB |
| Complete hg38 primary assembly, 16 threads | 559,601 loci in 204.5 s, peak 2,296 MB |
| Chromosome 1, previous defaults (arm 10, no filters) | 96,731 loci |

The two sequences HSeeker misses are the *c-MYC* NSE construct (ACCCTCCCC)₄ and the γ-globin
−228/−189 element. Both are **tandem direct repeats**, not mirror repeats, and therefore fall
outside the structural model.

---

# Part 1 — Response to reviewers

Reuse the existing replies verbatim for Reviewer 2, comments 8, 9 and 10, and for the first half
of Reviewer 3, comment 4; they are accurate and need only the commit identifier updated to the
release tag. Everything below replaces or extends what is currently in the response document.

Where a comment maps onto a specific edit, the new manuscript wording is quoted so the editor can
see it without opening the manuscript. No such block is given where the change is only
consequential.

---

## Reviewer 1

> HSeeker, an algorithm for systematic H-DNA sequence identification, is well presented … nothing
> to change in this reviewer's opinion.

We thank the reviewer for the positive assessment. No changes were requested. We note that the
revisions made in response to the other two reviewers — an expanded and balanced benchmark,
family- and study-aware validation, a wider comparison with existing tools, a genome-scale
comparison against experimental data, and re-measured runtimes — also bear on the aspects this
reviewer highlighted, and we believe they strengthen the manuscript further.

---

## Reviewer 2

### Comment 1 — discriminating power and the size of the negative set

We agree that nine negative controls cannot establish discriminating power and that accuracy is
uninformative at a 7.8:1 imbalance. We have rebuilt the benchmark around this objection.

Every experimental record was traced back to its primary publication. Records that could not be
traced, and records whose annotation came only from a secondary compilation, were removed
(Comment 5 of Reviewer 3). This left 70 experimental records, 65 forming and five non-forming,
so the shortage of experimental negatives is real and we did not attempt to disguise it.

Instead, and following the reviewer's own suggestions, we added 60 designed negatives produced by
a documented procedure under a fixed random seed: 22 positives in which the mirror symmetry was
destroyed by purine-preserving substitutions, 11 dinucleotide-preserving shuffles of positives, 12
random sequences matched to positives in length and GC content, four G-quadruplex-forming
sequences, three Z-DNA-forming (CG)ₙ repeats, five random B-DNA sequences, and three perfect
mirror repeats of balanced base composition. The last group is the sharpest of these tests: the
sequences satisfy mirror symmetry but not purine richness, so a tool that keys only on symmetry
will call them. The benchmark is now balanced at 65:65.

On structure-generation tools, we note respectfully that BoltzGen and RFdiffusion are designed to
generate and predict macromolecular structures and cannot establish whether a DNA duplex adopts
H-DNA under a given set of conditions; they could not supply negative labels. The
perturbation-based controls above implement the alternative the reviewer proposed.

Poly-A and comparable cases are addressed directly. We report the six single-base homopolymers as
a separate panel rather than inside the benchmark metrics, because the composition filters decide
their calls by construction and because a run of one base is mirror-symmetric at every position
and so does not test mirror-repeat detection. At the default settings HSeeker rejects A₂₀ — the
insert shown by Hanvey and colleagues to remain B-form — as well as poly-A₁₅, poly-A₃₀, poly-T₂₅
and poly-T₃₅. The same panel exposes a genuine limitation that we now state in the manuscript:
the homopolymer filter also rejects poly-(dC)₃₀, which forms an Mg²⁺-dependent triplex.

On the balanced benchmark HSeeker detected 63 of 65 forming sequences and rejected all 65
non-forming sequences (sensitivity 0.969, 95% CI 0.895–0.992; specificity 1.000, 0.944–1.000;
MCC 0.970, 0.926–1.000). The two sequences it misses are tandem direct repeats rather than mirror
repeats and lie outside the structural model. Mirror-repeat detection on its own, before the
purity requirement is applied, calls seven of the 65 negatives, including the three balanced
perfect mirrors; the purity requirement removes all seven without losing a forming sequence.

> **The manuscript now reads** (Results, *Detection performance on the revised benchmark*):
>
> "HSeeker detected 63 of the 65 H-DNA-forming sequences and rejected all 65 non-forming
> sequences, corresponding to a sensitivity of 0.969 (95% CI 0.895–0.992), a specificity of 1.000
> (0.944–1.000) and a Matthews correlation coefficient of 0.970 (0.926–1.000) (Table 1). The two
> sequences that were not detected, the *c-MYC* nuclease-sensitive element construct (ACCCTCCCC)₄
> and the γ-globin −228/−189 element, are tandem direct repeats rather than mirror repeats and
> therefore fall outside the structural model that HSeeker implements. Restricting the evaluation
> to the 70 experimental records gives a sensitivity of 0.969 (0.895–0.992) and a specificity of
> 1.000 (0.566–1.000); the width of the latter interval reflects the five experimental
> non-forming constructs that survived curation, and conclusions about specificity rest largely on
> the designed negatives."

### Comment 2 — reproducibility of the sequences and a script to run the analysis

We agree. The benchmark is distributed as two machine-readable files. The first lists all curated
records with their provenance: the primary citation, the curation decision and its justification,
the family and study identifiers used for grouping, and the assay, DNA context, pH, temperature,
ionic conditions and supercoiling state under which each construct was tested, together with the
passage in the source that supports those values. The second is the balanced 130-record benchmark
used for the reported metrics. The designed negatives are regenerated from a single script under a
fixed seed.

The repository contains scripts that reproduce every analysis in the paper end to end — the
benchmark metrics, the family- and study-aware validation, the component ablation and parameter
sweeps, the comparison against the other tools, the homopolymer panel, the overlap-rule comparison
and the genome-wide comparison — with all settings pinned in a single configuration file.

> **The manuscript now reads** (Data Availability):
>
> "The benchmark sequences, their provenance and the assay conditions under which each construct
> was characterised are provided as machine-readable files in the repository and as Supplementary
> Table 1. Scripts that regenerate every analysis and figure in this work, together with the
> parameter settings used for each, are provided in the same repository. The human chromosome 1
> and hg38 primary assembly sequences used for performance benchmarking were obtained from the
> UCSC Genome Browser. The single-stranded DNA maps and non-B DNA motif calls used for the
> genome-wide comparison are those published by Kouzine and colleagues. Raw benchmark commands,
> timing output and result files are included in the repository and will be archived in a public
> repository with a permanent identifier on publication."

### Comment 3 — the *E. coli* insertion experiment and validation on a genome tested for non-B DNA

We agree, and we have removed the experiment. Inserting the benchmark sequences into the *E. coli*
genome tested whether sequences are recovered in a genomic background; it did not test whether
HSeeker's predictions correspond to structures that form. Figure 3 and the associated text are
gone, and the question of further hits in the modified genome no longer arises.

In its place we followed the reviewer's suggestion and compared HSeeker with the data of Kouzine
and colleagues, who mapped single-stranded DNA genome-wide by permanganate/S1 footprinting in
activated mouse B cells and in human Raji cells, predicted non-B DNA motifs, and called a motif
formed when it carried single-stranded signal. We used their published calls and signal tracks, on
the genome assemblies they used, and scanned both genomes with HSeeker at its default settings.

HSeeker loci overlap 79.0% of the H-DNA structures found to be formed in mouse (13,517 of 17,109)
and 82.3% in human Raji cells (6,879 of 8,356). Among their own predicted H-DNA motifs, those that
HSeeker calls are formed more often than those it does not: 4.5% against 2.8% in mouse (odds ratio
1.65, 95% CI 1.59–1.71) and 5.7% against 1.7% in human (3.39, 3.20–3.59). Restricting to motifs
that the 36-nt reads of that study can map to raises the odds ratios to 1.68–2.07 and 4.16–5.70.

On the subset criterion specifically, we would gently note that a strict subset is not attainable
by any sequence-based predictor on these data, because formation depends on cell state: only about
3–5% of the authors' own predicted H-DNA motifs carried single-stranded signal in a given cell
type. The questions that these data can answer are whether a tool's calls are enriched for
structures that form, which HSeeker's are, and what fraction of the formed structures a tool
recovers, which is about four in five.

Most of the formed H-DNA that HSeeker does not report are short A/T-rich or single-base motifs
excluded by the composition filters. In the authors' data this class is formed at roughly half the
rate of other predicted H-DNA in mouse and a quarter in human, which supports excluding it by
default; it is nonetheless 17–19% of the formed H-DNA, and we state this cost explicitly.

The enrichment of the measured signal at HSeeker loci is described in our response to Comment 8 of
Reviewer 3 and is the subject of the new Figure 3.

### Comment 4 — simplification of the thermodynamic approach, and a limitations section

We agree, and we have added a Limitations section. The score is now described throughout as an
experimentally informed heuristic rather than a free-energy calculation, and all kcal/mol labels
have been removed from the score itself.

> **The manuscript now reads** (Discussion, *Limitations*):
>
> "HSeeker has several limitations that bear on how its output should be read. The stability score
> is an experimentally informed heuristic, not a calculated free energy. It assigns fixed weights
> to matched G and A positions, a uniform penalty to every other arm pairing, a constant reward to
> adjacent matched positions and an escalating penalty to runs of disrupted ones, and it then
> selects the highest-scoring arm window. It does not model backbone geometry, superhelical
> stress, pH-dependent protonation, ionic composition, temperature, loop or spacer energetics, or
> competition with duplex and other non-canonical conformations, and it has not been calibrated
> against measured free energies, melting temperatures or B-to-H transition thresholds. It is
> therefore a relative ranking of candidates within the scheme defined here, and it carries no
> physical units. We report no score threshold, and a candidate is called on detection alone.
>
> The benchmark labels come from assays performed under different conditions, including acidic
> buffers, Mg²⁺-containing buffers and negatively supercoiled plasmids. We record those conditions
> for every record, but a sequence characterised as H-DNA-forming under one set of conditions does
> not necessarily form the structure under physiological conditions, and the pooled performance
> figures should be read as discrimination of reported sequence labels across conditions.
>
> The composition filters are fixed rules rather than fitted parameters, but they were chosen with
> knowledge of the H-DNA literature and of the composition of this benchmark, so the benchmark is
> not an independent test of them. They are not free of cost: the homopolymer filter rejects
> poly-(dC)·(dG) tracts, which form an Mg²⁺-dependent triplex, and the AT-content filter removes a
> class of motifs that constitutes 17–19% of the H-DNA found to be formed in published
> single-stranded DNA maps.
>
> Finally, HSeeker detects antiparallel, purine-motif H-DNA at mirror repeats. Tandem direct
> repeats, parallel triplexes requiring cytosine protonation, and intermolecular triplexes are
> outside its model, and the two forming sequences it misses in our benchmark are direct repeats.
> Its predictions identify sequences with the potential to adopt H-DNA; they are not evidence that
> the structure forms at a given locus in vivo."

### Comment 5 — complexity and sensitivity analyses belong in the Results

We agree. The complexity analysis and the sensitivity and specificity analyses have been moved out
of Methods and now appear in the Results, where the Methods sections retain only the design of each
analysis. The Results sections, in order, are: computational performance and parallel scalability;
algorithmic complexity; detection performance on the revised benchmark; family- and study-aware
validation; component ablation and parameter robustness; comparison with existing tools; and the
genome-wide comparison with single-stranded DNA maps.

### Comment 6 — the filters were defined on the same dataset

We thank the reviewer; this is a fair description of the original manuscript, and we address it in
three ways.

First, we state it plainly. The AT-content and homopolymer rules are fixed pipeline rules rather
than parameters fitted during evaluation, but they were chosen with knowledge of the H-DNA
literature and of this benchmark's composition, and the manuscript now says so.

Second, we removed the circularity from the reported numbers. In the original benchmark the only
false positives at the analysis settings were four homopolymers, and the filters removed them. In
the revised benchmark the homopolymers are reported as a separate panel, and consequently neither
filter changes any benchmark call: the reported metrics do not depend on the filters, and equally
the benchmark offers no evidence in their favour. We say both.

Third, we provide evidence for the filters from outside the benchmark. In the published
single-stranded DNA maps, predicted H-DNA motifs of the class the filters exclude are found formed
at roughly half the rate of other predicted H-DNA in mouse and a quarter in human. With the
filters disabled HSeeker calls almost all predicted motifs and its selectivity disappears.

Separately, every parameter is now selected out of sample. Across 576 combinations of purity,
mismatch tolerance, arm length, spacer length and filter state, with selection confined to the
training folds, per-fold tuning reached an MCC of 0.939 under nested cross-validation against 0.970
for the fixed default configuration, and tuning did not beat the fixed configuration in any of 200
random splits.

> **The manuscript now reads** (Methods, *Sequence composition filters*):
>
> "Two composition filters are applied to detected candidates. Motifs whose full sequence, counting
> both arms and the spacer, is at least 80% A/T are discarded before scoring, because very A/T-rich
> sequences destabilise triplex formation through the geometric incompatibility of Hoogsteen
> triplets and the reduced stability of A·A-T triads. Motifs whose sequence consists of a single
> repeated base are discarded, because a homopurine or homopyrimidine run satisfies the composition
> criterion without presenting a mirror repeat in any meaningful sense.
>
> Both rules are fixed and are not fitted during evaluation. They were, however, chosen with
> knowledge of the H-DNA literature and of the composition of our benchmark, and the benchmark is
> therefore not an independent test of them: because the single-base homopolymers are reported
> separately (see Results), neither filter changes any benchmark call, so the reported metrics
> neither depend on the filters nor support them. Evidence bearing on the filters is drawn instead
> from published single-stranded DNA maps and is presented in the Results."

### Comment 7 — confidence intervals in the tool-versus-tool comparison

We agree, and no comparative claim now rests on a point estimate. Sensitivity, specificity, PPV and
NPV are reported with Wilson intervals and MCC with bootstrap intervals, computed both over records
and over whole sequence families, since records within a family are not independent. For each
comparison with another tool we report the paired difference in sensitivity, specificity and MCC
with bootstrap intervals over records and over families, and an exact McNemar test on the records
where the two tools disagree.

We fixed the rule for reading these comparisons before running them: we describe HSeeker as better
than another tool only where the family-level interval of the MCC difference excludes zero and the
difference keeps its sign when the other tool's geometry is matched to HSeeker's and when its
parameters are tuned out of fold. By that rule HSeeker outperforms all four comparators on the
balanced benchmark, with MCC differences of +0.15 (95% CI 0.03–0.22) against triplex, +0.39
(0.18–0.56) against NeSSie, +0.52 (0.40–0.63) against nBMST and +0.95 (0.71–1.11) against
Triplexator. Against triplex the advantage comes entirely from the designed negatives; on the
experimental records alone the interval reaches zero (+0.06, 0.00–0.22) and we describe the two
tools as comparable there. The Key Points and the Results now agree on this wording.

### Comment 8 — `-march=native`

*(Keep the existing reply. Update the commit identifier to the release tag at acceptance.)*

### Comment 9 — unsynchronised profiling counters

*(Keep the existing reply. Update the commit identifier to the release tag at acceptance.)*

### Comment 10 — case inconsistency in `putative_triplex`

*(Keep the existing reply. Update the commit identifier to the release tag at acceptance.)*

### Comment 11 — two `app.py` files

We agree that two front-ends were confusing, and we thank the reviewer for the kind words about the
software engineering. Rather than relabel the demonstration interface, we removed it. HSeeker now
ships a single web application, the FastAPI service, which is also what the hosted deployment and
the container image serve; the Streamlit dependency is gone from the package. The command-line
tool, the Python package and the web application now share the same default settings and return
identical results for identical input.

> **The manuscript now reads** (Results, *HSeeker web interface*, first sentence):
>
> "HSeeker provides a single web interface, which uses the same scanning library and the same
> default settings as the command-line tool, so that the two return identical results for identical
> input."

---

## Reviewer 3

### Comment 1 — benchmark size, imbalance and confidence intervals

We agree, and we have both expanded the negative set and qualified every claim that depends on it.
The benchmark is now balanced at 65:65 (Comment 1 of Reviewer 2). Performance is reported with
confidence intervals on three sets, and we state the size of the experimental negative class
plainly rather than letting the balanced figure stand for it.

| Evaluation set | n | TP | FN | TN | FP | Sensitivity | Specificity | PPV | NPV | MCC |
|---|---|---|---|---|---|---|---|---|---|---|
| Balanced benchmark | 130 | 63 | 2 | 65 | 0 | 0.969 (0.895–0.992) | 1.000 (0.944–1.000) | 1.000 (0.943–1.000) | 0.970 (0.898–0.992) | 0.970 (0.926–1.000) |
| Experimental records | 70 | 63 | 2 | 5 | 0 | 0.969 (0.895–0.992) | 1.000 (0.566–1.000) | 1.000 (0.943–1.000) | 0.714 (0.359–0.918) | 0.832 (0.569–1.000) |
| Construct-level records | 58 | 51 | 2 | 5 | 0 | 0.962 (0.872–0.990) | 1.000 (0.566–1.000) | 1.000 (0.930–1.000) | 0.714 (0.359–0.918) | 0.829 (0.556–1.000) |

Only five experimental non-forming constructs survived source-level curation. The specificity
estimate therefore rests mainly on the designed negatives, and the manuscript says so wherever
specificity is discussed.

### Comment 2 — related sequence variants and family-aware validation

We agree that the records are not independent, and we are grateful for the suggestion, which turned
out to be straightforward to implement and informative. Records were grouped by sequence family and
by study using the curated provenance. The 70 experimental records span 48 families and 26 studies;
designed negatives inherit the family and study of the record they were derived from, giving 52
family groups and 27 study groups across the benchmark. The largest family, the pXY32 substitution
series and its derivatives, accounts for 37 of the 130 records.

| Validation scheme | Fixed default configuration | Per-fold tuned |
|---|---|---|
| Leave-one-family-out | 0.970 (family-cluster CI 0.907–1.000) | 0.924 |
| Leave-one-study-out | 0.970 (study-cluster CI 0.901–1.000) | 0.924 |
| One representative per family (n = 52) | 0.862 | — |
| One representative per family and class (n = 74) | 0.945 | — |

With the parameters fixed, each record's call depends only on its own sequence, so family-aware
validation does not move the point estimate. What it moves is the precision: the family-level
interval is wider than the record-level interval, and we report both, since they answer different
questions.

> **The manuscript now reads** (Results, *Family- and study-aware validation*):
>
> "Because many benchmark records are closely related variants from the same experimental series,
> we grouped records by sequence family and by study and repeated the evaluation holding out one
> whole group at a time. Under leave-one-family-out and leave-one-study-out validation the fixed
> default configuration reached an MCC of 0.970, against 0.924 when the parameters were tuned
> within each training fold. Collapsing the benchmark to one representative per family gave an MCC
> of 0.862 over 52 records, and to one representative per family and class 0.945 over 74 records.
> Because the parameters are fixed, a record's call depends only on its own sequence, and
> group-aware validation does not shift the point estimate; it widens the confidence interval, and
> both the record-level and the family-level intervals are reported throughout."

### Comment 3 — circularity in the stability-score threshold

We agree, and we are grateful for the suggestion, which proved decisive. We re-evaluated the
threshold exactly as proposed: within each training fold we selected the threshold maximising
Youden's J and applied it only to the held-out records, under repeated stratified five-fold
cross-validation, leave-one-family-out and leave-one-study-out validation.

A threshold chosen this way never improved performance. Out of fold it gave an MCC of 0.949 (95% CI
0.937–0.961) under cross-validation and 0.955 (0.898–1.000) under both group-aware schemes, against
0.970 (0.925–1.000) with no threshold at all. It removed no false positive, because after the
purity requirement no non-forming record in the benchmark retains any candidate, and out of fold it
lost forming records whose score fell below a threshold learned without them. For the same reason
the AUC on this benchmark, 0.985, only restates the detection result, since it equals
(1 + sensitivity)/2 when no negative carries a candidate.

We have therefore stopped using a score threshold. A sequence is called H-DNA-forming when it
contains a candidate that passes the detection criteria, and the stability score is reported as a
ranking of candidates. The Youden-selected threshold, the in-sample AUC and the threshold-based
confusion matrix are removed from the manuscript, and the web interface no longer labels sites
above a score of 50 as stable.

> **The manuscript now reads** (Results, *Detection performance on the revised benchmark*, closing
> paragraph):
>
> "We also asked whether the stability score could improve classification beyond detection. Within
> each training fold we selected the score threshold maximising Youden's J statistic and applied it
> only to the held-out records. Out of fold this produced an MCC of 0.949 (0.937–0.961) under
> repeated five-fold cross-validation and 0.955 (0.898–1.000) under leave-one-family-out and
> leave-one-study-out validation, below the 0.970 obtained with no threshold. The threshold removed
> no false positive, because no non-forming record retains a candidate once the purity requirement
> is applied, and it discarded forming records whose scores fell below a threshold learned without
> them. For the same reason the area under the ROC curve on this benchmark, 0.985, restates the
> detection result rather than adding to it. We therefore report no score threshold: a sequence is
> called H-DNA-forming when it contains a candidate satisfying the detection criteria, and the
> stability score is used to rank candidates."

### Comment 4 — the thermodynamic interpretation of the stability score

*(Keep the existing reply, which is accurate, and add the following.)*

We have also acted on the second half of the comment. The assay conditions are now recorded for
every experimental record — the assay itself, the DNA context, pH, temperature, ionic conditions
and supercoiling state, with the supporting passage from the source and an explicit marker where
the source does not report a value. Of the 70 experimental records, 64 have a reported pH. These
conditions are supplied with the benchmark and summarised in Supplementary Table 1, the pooled
figures are described as discrimination of sequence labels across assay conditions, and the new
Limitations section states that a label obtained under one set of conditions does not imply
formation under physiological conditions.

> **The manuscript now reads** (Methods, *Benchmark dataset*, closing paragraph):
>
> "For each experimental record we recorded the assay used, the DNA context, and the pH,
> temperature, ionic conditions and supercoiling state under which the structure was assessed,
> together with the passage in the source supporting those values; values that the source does not
> report are marked as such. Sixty-four of the 70 experimental records have a reported pH. Because
> the labels are outcomes under heterogeneous conditions, the pooled performance figures should be
> read as discrimination of reported sequence labels across assay conditions rather than as
> propensities under a single common condition."

### Comment 5 — benchmark curation

We agree, and we have audited every record at source level; the per-record audit is provided with
the benchmark.

Six records carried low-confidence annotations derived from a secondary compilation. Three had no
traceable primary source and were removed; three had traceable sources but were removed along with
all other low-confidence records. Two further records were excluded on validity grounds: an
inosine-substituted oligonucleotide that cannot be represented in A/C/G/T, and a record describing
an intermolecular triplex-forming oligonucleotide rather than an intramolecular structure. The two
single-base homopolymers were moved to the separate panel described above. Citation identifiers
were corrected for 25 records, 21 of which are retained in the benchmark, with the original
identifiers kept in the files so that each correction is auditable. For the four Friedreich ataxia (GAA)ₙ records, whose cited identifiers
pointed to unrelated papers, we now cite the primary experimental studies; because those studies
characterised closely related constructs rather than these exact sequences, the records are
retained in a clearly labelled class-level set that is reported separately from the construct-level
benchmark.

On pRW1406 and pRW1411 specifically, the primary source resolves the inconsistency. The abstract of
Hanvey and colleagues enumerates seven inserts by G+C content: the (TCC)₈ insert at 66% G+C and the
(TTC)₈ insert at 33% G+C are distinct constructs, and both are among the five inserts reported to
adopt the non-B conformation. The duplicated sequence in the secondary compilation was a
transcription error. We retained (TCC)₈ and added (TTC)₈ as a new record, named after its verified
insert rather than a plasmid number, since the plasmid-to-insert mapping is precisely what the
secondary compilation got wrong.

The construct-level benchmark now contains 58 records, 53 forming and five non-forming, with 12
further class-level records reported separately.

### Comment 6 — the comparator set and the strength of the claim

We agree on both points. We added the two mirror-repeat tools the reviewer names and ran all four
comparators on the revised benchmark. The settings, the rule converting each tool's output into a
call, the statistics and the wording rule were fixed before any comparator was run. Each tool was
run at its documented defaults, with its geometry matched to HSeeker's where it exposes the
corresponding parameters, and with a small parameter grid selected inside training folds only — the
same opportunity HSeeker's defaults were given. A record was called forming when a tool reported at
least one qualifying hit in it.

| Tool, default settings | Sensitivity | Specificity | MCC (95% CI) | Tuned out of fold |
|---|---|---|---|---|
| HSeeker | 63/65 | 65/65 | 0.970 (0.909–1.000) | 0.939 |
| triplex | 62/65 | 56/65 | 0.819 (0.743–0.935) | 0.819 |
| NeSSie | 33/65 | 65/65 | 0.583 (0.421–0.779) | 0.547 |
| nBMST | 24/65 | 64/65 | 0.449 (0.351–0.561) | 0.597 |
| Triplexator | 8/65 | 58/65 | 0.024 (−0.132–0.250) | 0.520 |

These differences follow from what each tool searches for, and the manuscript now describes this
rather than presenting the comparison as a like-for-like ranking. Triplexator searches for duplex
target sites for a third strand supplied separately and does not require a mirror repeat; at its
defaults it detects few of the forming records, and at the settings used in our original manuscript
it detects 64 of 65 but also calls 42 of the 65 negatives, including all 22 mirror-disrupted
mutants. NeSSie reports k-mers that are themselves mirror-symmetric and purine- or
pyrimidine-rich rather than arm–spacer–arm motifs; it made no false call but detected none of the
16 pXY32 constructs. nBMST's mirror-repeat module matches HSeeker's geometry, but its triplex
classification accepts only mirror repeats with more purines than pyrimidines on the reported
strand, so it detected none of the 19 forming records written on the pyrimidine-rich strand, whereas
its unclassified mirror-repeat output detected 15 of them. triplex, the only other tool designed
for intramolecular triplexes, detected nearly the same forming records as HSeeker but called nine
designed negatives forming, seven of them mirror-disrupted mutants.

We have harmonised the Key Points with the Results. The Key Point now reads: "On a balanced
benchmark of 130 experimentally characterised and designed sequences, HSeeker achieved higher
accuracy than four existing tools, matching triplex on experimentally characterised sequences while
rejecting the mirror-disrupted controls that triplex accepted." The corresponding Results sentence
makes the same distinction, and the earlier claim that HSeeker "outperformed existing tools in
benchmarking" without qualification has been removed.

### Comment 7 — the *E. coli* insertion experiment and the overlap threshold

We agree. The experiment tested computational recovery in a genomic background rather than H-DNA
formation, and Figure 3 gave it more weight than it could carry. We have removed the experiment and
its figure from the manuscript rather than merely toning them down. In their place the manuscript
compares HSeeker's genome-wide predictions with experimental single-stranded DNA maps (Comment 3 of
Reviewer 2 and Comment 8 below). Because the insertion experiment is no longer part of the paper,
the 80% overlap criterion it used no longer enters any reported result, and a sensitivity analysis
over that threshold is no longer required.

### Comment 8 — genome-scale biological validation

We thank the reviewer for this suggestion; it produced the most informative analysis in the revised
manuscript. We compared HSeeker's predictions with the structure-sensitive nuclease mapping of
Kouzine and colleagues: permanganate/S1 footprinting of single-stranded DNA in activated and
resting mouse B cells and in human Raji cells. Recall and selectivity against their published calls
are described in our response to Comment 3 of Reviewer 2. Here we describe two enrichment analyses
of the measured signal, which form the new Figure 3. The first was specified, with its controls and
statistics, before any result was computed; the second was added after the first and was likewise
specified before it was computed. We report both.

Single-stranded signal depends strongly on base composition, so control windows were matched to
each HSeeker locus on length, A/T count and purine/pyrimidine asymmetry, and were required not to
overlap any HSeeker locus or predicted H-DNA motif. Matching on asymmetry is essential: controls
matched on A/T content alone carry 1.7 to 2.5 times more signal than HSeeker loci, because
mixed-sequence DNA gives more signal than purine/pyrimidine-asymmetric DNA of the same composition.

In the first analysis, across all HSeeker loci with genome-wide matched controls, loci in activated
mouse B cells carried 14% more single-stranded signal than their controls (ratio of means 1.137,
95% CI 1.127–1.147; 825,749 loci). In resting B cells the ratio was 0.822 (0.818–0.827) and in
human Raji cells 0.977 (0.969–0.984). This analysis does not separate H-DNA from other sources of
single-stranded DNA and does not account for the poor mappability of repetitive purine tracts.

In the second analysis we removed from both loci and controls every interval that the reads cannot
map to, every interval within 100 bp of a motif predicted to form a G-quadruplex, Z-DNA, a
cruciform or a stress-induced melting site, and every interval within 2 kb of an annotated
transcription start site, and drew ten controls per locus within 50 kb. In activated mouse B cells
the isolated loci carried 25% more signal than their controls (1.254, 1.179–1.332; a locus exceeded
one of its own controls with probability 0.514, p = 2 × 10⁻¹⁷). In human Raji cells there was no
detectable excess (1.023, 0.988–1.061; p = 0.08). Two built-in checks behaved as required in both
species: comparing one control against the remaining controls of the same locus gave no difference
(0.500–0.502), and the motifs the authors called formed were enriched 7.4-fold in mouse and
2.1-fold in human under the identical procedure.

Because G-quadruplexes, Z-DNA, duplex melting and transcription do not depend on mirror symmetry
once composition is matched, whereas H-DNA stability increases with arm length, we also stratified
by arm length. In activated mouse B cells the excess rose from no detectable effect at 8–9 nt (1.12,
0.93–1.35) to 1.58 (1.34–1.83) at 12–14 nt and 1.6–1.8 at longer arms (Spearman p = 3 × 10⁻²¹). The
enrichment was absent in resting B cells (0.858, 0.831–0.885), in which even the authors' formed
motifs showed no excess, consistent with the activation- and supercoiling-dependent formation they
report. No arm-length class differed from unity in Raji cells.

We report the human result as it stands. The Raji data have no matched input control and a smaller
dynamic range — the authors' own formed motifs are enriched 2.1-fold there against 7.4-fold in
mouse — but we do not treat that as an explanation, and the manuscript states that we did not
observe the enrichment in human cells.

Finally, the chromosome 1 count quoted in the comment, 96,729 loci, was obtained with the previous
default minimum arm length of 10 bp and without the composition filters. With the revised defaults
HSeeker reports 45,927 loci on chromosome 1, or 106,481 with the composition filters disabled, and
the manuscript now states the settings alongside every count.

### Comment 9 — ablation and parameter robustness

We agree, and we have added both analyses on the revised benchmark; a notebook and machine-readable
tables are provided. A record is called forming when at least one candidate survives the stage in
question, and no score threshold is used at any stage.

Mirror-repeat detection alone detected 63 of 65 forming records and seven of 65 non-forming records
(sensitivity 0.969, specificity 0.892, MCC 0.864). The seven were the experimentally non-forming
*c-MYC* control B-33, one G-quadruplex sequence, two Z-DNA (CG)ₙ repeats and the three perfect
mirror repeats of balanced composition. Adding the purity requirement removed all seven without
losing a forming record (specificity 1.000, MCC 0.970). Adding the AT-content filter, pairing-only
scoring, pairing plus stacking, and the complete model with boundary optimisation and homopolymer
filtering changed no record-level call. Pairing-only and pairing-plus-stacking scores rank the same
candidates, and once purity is applied no non-forming record retains a candidate, so the score
AUROC of 0.985 restates the detection result. On the 70 experimental records alone the MCC is 0.708
for detection alone and 0.832 from the purity step onward. Retaining all candidates, longest-arm
selection and score-based selection gave identical record-level calls. The scorer exposes no
independent switch for the mismatch penalty, so we do not claim to have ablated that penalty
separately.

Varying one parameter at a time around the defaults: purity 0.80–0.90 gave an MCC of 0.970, while
0.95 and 1.00 lost nine forming records (0.843); mismatch tolerances of 0 and 0.05 detected 49 of 65
forming records (0.778), 0.10 detected 63 (0.970), and 0.15 and 0.20 detected 64 but added a false
positive (0.969); a minimum arm length of 6 added ten false positives (0.839), 7 and 8 gave 0.985
and 0.970, 9 and 10 lost eight forming records (0.884) and 12 lost sixteen (0.778); a maximum spacer
of 3 detected 41 forming records (0.679), 5 detected 55 (0.856), and 8 to 15 detected 63 (0.970).
The single record gained at a minimum arm length of 7, and at a mismatch tolerance of 0.15, is the
tandem direct repeat (ACCCTCCCC)₄, which lies outside the structural model.

Across all 576 combinations with selection confined to training folds, the fixed default
configuration performed as well as or better than per-fold tuning under nested cross-validation,
leave-one-family-out and leave-one-study-out validation (0.970 against 0.939, 0.924 and 0.924). A
minimum arm length of 8 was selected in 99% of nested folds and in every family and study fold, and
a mismatch tolerance of 0.10 in 48 of 52 family folds. On this basis we changed the default minimum
arm length from 10 to 8 bp. Performance therefore does depend on these choices, as the reviewer
anticipated, but the reported results do not rest on a single hand-selected parameterisation.

### Comment 10 — the overlap-removal strategy

*(Replace the current reply in full. The opening paragraph of the existing draft is garbled and the
comparison figures were computed on the superseded 79-sequence collection.)*

We thank the reviewer for this comment. We agree that arm length is not a validated proxy for
biological activity, and we now implement the alternative proposed and compare the two directly.

We retain the longest-arm rule as the default for three reasons, which the manuscript now states.
Among overlapping candidates covering the same mirror-repeat region, the longest arm describes the
full extent of the repeat, which is the region available for strand exchange. The criterion is
purely geometric and does not depend on the stability score, which as we acknowledge in response to
Comment 4 is a heuristic rather than a calibrated free energy; fine score differences between
overlapping candidates should not decide what a user sees by default. And it is faster at genome
scale, because overlaps are resolved before candidates are scored.

HSeeker now offers three modes in the Python package, the command-line tool and the web interface.
Retaining all candidates returns every candidate surviving the requested filters. Longest-arm
selection, the default, performs a start-ordered sweep in which a candidate with a longer arm
replaces the overlapping incumbent, with a shorter spacer breaking ties. Score-based selection
scores every candidate first and then considers them in descending order of score, retaining a
candidate when it does not overlap one already retained, with ties resolved by longer arm, then
shorter spacer, then earlier coordinate. Coordinates are 1-based and inclusive.

The order of operations matters and is now documented. AT-content filtering precedes scoring.
Longest-arm selection runs before scoring and before the scored homopolymer check, whereas
score-based selection runs after both, so a candidate removed by the scored filter leaves a
lower-scoring overlapping alternative eligible. We recommend the default for genome-wide
annotation, score-based selection for locus-level work such as choosing which variant to test
experimentally, and retaining all candidates for exhaustive or downstream statistical analysis. A
higher score helps prioritise candidates and is not experimental evidence of formation.

The two rules agree far more often than they differ. On the benchmark, 43 of the 70 experimental
records have overlapping candidates after filtering; the two rules return identical intervals for
38 of them and differ for five, and the record-level call is identical for all 130 records and for
all 74,880 record-by-configuration calls across our parameter grid. On human chromosome 1 at the
default settings, longest-arm selection reports 45,927 loci and score-based selection 47,217; every
longest-arm locus is overlapped by a score-selected locus, 41,867 of them (91.2%) have identical
coordinates, and score-based selection reports 991 additional loci. With the composition filters
disabled the counts are 106,481 and 108,224, with 82.6% identical. The two rules therefore locate
the same regions and differ mainly in which representation of a region they report, at a modest
cost in time: 23.9 s against 16.3 s on chromosome 1 with 16 threads.

The following worked examples, with scoring enabled and a minimum arm of 7, maximum arm 25, maximum
spacer 5, purity 0.9 and mismatch tolerance 0.2, show when the rules agree and when they diverge.

| Input | Candidate | Interval | Arm / spacer | Total score |
|---|---|---|---|---|
| `GAGAAAAGAGGGAGAAAAAGGAA` | gagaaaagagggagaaaaagg | 1–21 | 10 / 1 | 53.635 |
| `GAGAAAAGAGGGAGAAAAAGGAA` | gaaaagagggagaaaaa | 3–19 | 8 / 1 | 53.635 |
| `GGAAGAAGAAGGAGGAGAGAGGGAAGGAA` | ggaagaagaaggagg | 1–15 | 7 / 1 | 38.575 |
| `GGAAGAAGAAGGAGGAGAGAGGGAAGGAA` | aaggaggagagagggaaggaa | 9–29 | 10 / 1 | 20.360 |

In the first input both rules return 1–21: longest-arm selection because its arm is longer,
score-based selection because the two scores tie and the longer arm wins the tie. In the second,
longest-arm selection returns 9–29 because its arm is longer, while score-based selection returns
1–15 because its score is higher. Retaining overlaps returns both candidates in each case. An
overlap chain need not collapse to a single representative: scanning (GAA)₂₀ with a minimum arm of
10 produces 102 overlapping candidates, and both rules retain the two compatible hits at 3–32 and
35–60.

### Comment 11 — computational benchmarking and the whole-genome extrapolation

We agree on both counts. We re-measured runtime and memory on the final code at the default
settings, on the same hardware as the original benchmark, with three runs at each point and peak
memory recorded externally. On human chromosome 1 the complete pipeline took 173.8 s with one
worker and 16.3 s with 16, a 10.65-fold speedup, with a peak memory of 1,334 MB. Rather than
extrapolating, we ran the complete hg38 primary assembly: 204.5 s with 16 workers, peak memory
2,296 MB, reporting 559,601 loci, with the three runs producing identical output. Score-based
overlap selection took 23.9 s on chromosome 1.

We also re-timed the comparators on the same node on the same day, so that the comparison is
like for like. At 16 cores chromosome 1 took 45.0 s with Triplexator and 64.7 s with triplex run as
16 parallel processes; nBMST and NeSSie have no parallel mode and took 96.5 s and 295.7 s
respectively on one core. triplex and Triplexator reproduced the hit counts of our original
measurement exactly.

All benchmark commands, raw timing output, hardware descriptions and input checksums are included
in the repository and will be archived with a permanent identifier on publication. The Data
Availability statement no longer refers to materials available on request.

### Comment 12 — claims, terminology and figures

We agree and have moderated the Discussion. The three statements quoted have been removed or
rewritten, the terminology has been made consistent, and the figures have been simplified.

> **The manuscript now reads** (Discussion, replacing the three statements):
>
> "The stability score ranks candidates in a manner consistent with the biochemical and genetic
> evidence available for the sequences in our benchmark; it is an experimentally informed heuristic
> and has not been calibrated against measured thermodynamic quantities."
>
> "Among currently available tools for intramolecular H-DNA, HSeeker is the only one that reports a
> per-motif stability ranking alongside the detected coordinates."
>
> "This efficiency makes genome-wide profiling and cross-species comparison practical on commodity
> hardware."

On terminology, predictions are now described as putative H-DNA-forming sequences or predicted
H-DNA motifs throughout; "H-DNA-forming sequence" without qualification is reserved for sequences
with experimental support; and "mirror repeat" is used for the sequence feature rather than as a
synonym for the structure. The manuscript has been language-edited. Figure 1 has been simplified,
the algorithm panels now present the pipeline as a short ordering of stages followed by one panel
per stage, and the numerical content that was previously embedded in Figure 1 has moved to the
Results.

---

# Part 2 — Manuscript changes

Ordered as they appear in the paper. Every item is either required by a reviewer comment, required
for internal consistency once another change is made, or a factual correction found while preparing
the revision.

## 2.1 Title block, Abstract and Key Points

**Abstract.** Three changes: remove "state-of-the-art"; describe the score as a heuristic rather
than thermodynamic; describe all three overlap modes; and add one sentence of validation. Replace
from "Using a center-outward search algorithm approach" to the end of the scoring sentence, and add
the validation sentence before "HSeeker is also user-friendly".

> "To address this gap, we present "HSeeker", a computational tool that compiles a three-part
> algorithm to identify and score potential H-DNA-forming sequences. Using a center-outward search
> algorithm approach, HSeeker evaluates candidate hinge positions and spacer lengths while allowing
> configurable mirror mismatches and purine-pyrimidine composition thresholds. Overlapping
> candidates are resolved by one of three user-selectable rules: retaining the longest and most
> compact motif within each overlapping region, retaining the highest-scoring motif, or reporting
> every candidate. Finally, an experimentally informed stability scoring model evaluates the
> candidate motifs using Hoogsteen G·G and A·A pairing energies, consecutive-pair stacking, and
> penalties for mismatches and disrupted stacks; the resulting score ranks candidates relative to
> one another within this scheme and is not a calibrated free energy. The scoring procedure also
> optimises motif boundaries by trimming weak terminal positions and reassigning unstable arm
> positions to the spacer. HSeeker reports genomic coordinates, sequence information, pairing and
> stacking components, and an overall stability score. On a balanced benchmark of 130
> experimentally characterised and designed sequences, HSeeker detected 63 of 65 H-DNA-forming
> sequences and rejected all 65 non-forming sequences, and its genome-wide predictions were
> enriched for mirror repeats found single-stranded in cells by permanganate/S1 footprinting."

**Key Points.** Replace the third and fourth bullets; adjust the second.

> "· It combines a center-outward mirror-repeat search, three selectable overlap-resolution rules,
> and an experimentally informed stability-scoring model that incorporates Hoogsteen G·G and A·A
> pairing, base stacking, and mismatch penalties."
>
> "· On a balanced benchmark of 130 experimentally characterised and designed sequences, HSeeker
> achieved higher accuracy than four existing tools, matching triplex on experimentally
> characterised sequences while rejecting the mirror-disrupted controls that triplex accepted."
>
> "· HSeeker scales efficiently, analysing human chromosome 1 in 16 s and the complete human
> genome in 3.4 min on 16 cores, and is available as an open-source Python package and a
> user-friendly web application."

## 2.2 Introduction

Two factual corrections, both arising from running the tools rather than reading their
documentation.

**nBMST description.** The current text states that nBMST "only includes mirror repeats of at
least 90% purines or pyrimidines, without tolerating any mismatch in mirror symmetry". Its
mirror-repeat module in fact tolerates no mismatch but applies a far weaker composition test than
90%. Replace with:

> "For example, nBMST provides a general interface for searching multiple classes of non-B DNA
> motifs, but its mirror-repeat module tolerates no mismatch in mirror symmetry, and the subset it
> labels as triplex-forming is selected by a purine/pyrimidine count on the reported strand rather
> than by a purity threshold, so mirror repeats written on the pyrimidine-rich strand are not
> labelled as triplex-forming."

**NeSSie description.** The current text describes NeSSie as "optimized for detecting the mirror
repeats characteristic of H-DNA". It detects k-mers that are themselves mirror-symmetric, which is
a different object from an arm–spacer–arm motif. Replace with:

> "NeSSie is a C/C++ library and tool for the detection of approximate sequence symmetries,
> including mirror symmetry and triplex-forming potential, using an optimal global alignment based
> on the Needleman-Wunsch algorithm and allowing gaps and mismatches. Its search returns
> symmetric or purine-rich k-mers rather than arm–spacer–arm motifs, so a mirror repeat whose two
> arms flank a non-symmetric spacer is not returned as a single motif."

## 2.3 Methods

### New subsection — *Benchmark dataset* (insert before *Sensitivity analysis*, which is renamed)

> "**Benchmark dataset.** The benchmark contains 130 sequences, 65 reported to form H-DNA and 65
> reported or designed not to form it. Seventy are experimental records, 65 forming and five
> non-forming. Each was traced to its primary publication; records whose annotation derived only
> from a secondary compilation, and records with no traceable primary source, were removed.
> Fifty-eight of the retained records are construct-level, meaning the published experiment tested
> the sequence in the benchmark, and are used for the primary performance figures. The remaining
> twelve are class-level, meaning the published experiment tested a closely related construct of
> the same repeat family, and are reported separately. Two records were excluded on validity
> grounds: an inosine-substituted oligonucleotide that cannot be represented in A/C/G/T, and a
> record describing an intermolecular triplex-forming oligonucleotide rather than an intramolecular
> structure.
>
> Because only five experimental non-forming constructs survived curation, we added 60 designed
> negatives, generated by a documented procedure under a fixed random seed: 22 positives whose
> mirror symmetry was destroyed by purine-preserving substitutions, 11 dinucleotide-preserving
> shuffles of positives, 12 random sequences matched to positives in length and GC content, four
> G-quadruplex-forming sequences, three Z-DNA-forming (CG)ₙ repeats, five random mixed-sequence
> controls, and three perfect mirror repeats of balanced base composition. Apart from the last
> group, every designed negative was verified to contain no homopurine mirror repeat satisfying the
> detection criteria; the three balanced perfect mirrors are included precisely because they do
> contain such a repeat but are not purine-rich, and they therefore test whether mirror symmetry
> alone triggers a call.
>
> Single-base homopolymers are reported as a separate panel rather than within the benchmark,
> because the composition filters determine their calls by construction and because a run of one
> base is mirror-symmetric at every position and so does not test mirror-repeat detection.
>
> For each experimental record we recorded the assay used, the DNA context, and the pH,
> temperature, ionic conditions and supercoiling state under which the structure was assessed,
> together with the passage in the source supporting those values; values that the source does not
> report are marked as such. Sixty-four of the 70 experimental records have a reported pH. Because
> the labels are outcomes under heterogeneous conditions, the pooled performance figures should be
> read as discrimination of reported sequence labels across assay conditions rather than as
> propensities under a single common condition."

### New subsection — *Detection settings*

Required because the paper quotes a command line with superseded values, and because the three
interfaces now share defaults.

> "**Detection settings.** Unless stated otherwise, all analyses used the default settings: minimum
> arm length 8 bp, maximum arm length 1,000 bp, maximum spacer 10 bp, purity threshold 0.90,
> mismatch tolerance 0.10, AT-content threshold 0.80 with homopolymer exclusion enabled, and
> overlapping candidates resolved by the longest-arm rule. These are the defaults of the Python
> package, the command-line tool and the web application alike, so that the three interfaces return
> identical results for identical input. The corresponding command line is
>
>     hseeker -seq <fasta> -out <prefix>
>
> and parameters are given explicitly wherever an analysis departs from these values."

### Revise — *Sequence composition filters*

Replace in full with the text quoted under Reviewer 2, Comment 6 above.

### Revise — *Stability score*

Four corrections, of which the third is a genuine discrepancy between the published description and
the implementation.

1. **Remove the free-energy framing of the output.** Keep the attribution of the pairing weights to
   the published interaction energies, since that is where they come from, but stop describing the
   combined score as a thermodynamic quantity. Replace the opening sentence of the section:

   > "The pairing weights used by the stability score were taken from the interaction energies
   > reported for Hoogsteen base pairs within a DNA triad, computed by simulating triplex structures
   > without the constraints of the phosphate backbone. The reported interaction energy of a G·G
   > reverse Hoogsteen hydrogen bond in a G·G-C triplet is 7.63 kcal/mol and that of an A·A reverse
   > Hoogsteen hydrogen bond in an A·A-T triplet is 3.83 kcal/mol; the implementation uses 7.6 and
   > 3.83. The same work reports stacking interaction energies within the third strand ranging from
   > 2.15 to 6.65 kcal/mol, from which we adopt a single value of 5.0 as the stacking term and as
   > the first disrupted-stack penalty. The combination of these terms is a ranking score, not a
   > free energy: it has no validated physical units, and we do not report it in kcal/mol."

2. **State the clamp.** Add after the definition of the total score:

   > "Negative totals are clamped to zero, so the reported score is non-negative."

3. **Correct the spacer constraint.** The manuscript states the constraint as
   2(n − R\*) + m ≤ min(20, L\*). The implementation requires the spacer produced by trimming to be
   no longer than the selected window, that is 2(a − R\*) + ℓ_s ≤ (R\* − L\*)·v with v = 1, and the
   fixed 20 bp cap is not applied in the window search. Replace the two displayed constraints with:

   > "subject to R\* − L\* ≥ ℓ, where ℓ is the minimum arm window, and to
   >
   >     2(a − R*) + ℓ_s ≤ (R* − L*) · v,
   >
   > with a the detected arm length, ℓ_s the detected spacer length and v = 1, so that the spacer
   > formed by folding the trimmed 5′ and 3′ arm tails inward does not exceed the length of the
   > selected window. If no window satisfies the constraint, the detected arms are scored as they
   > stand."

4. **Resolve the Supplementary Table 1 collision.** The sentence "we consider the values in
   Supplementary Table 1" refers to the scoring weights, while the benchmark records are also cited
   as Supplementary Table 1. Renumber: the benchmark records remain Supplementary Table 1, and the
   scoring weights become **Supplementary Table 2**. Update both citations. See §4 for the full
   supplementary numbering.

### New subsection — *Validation design*

> "**Validation design.** A sequence is called H-DNA-forming when at least one candidate survives
> detection and filtering; no score threshold is applied. Sensitivity, specificity, positive and
> negative predictive values are reported with Wilson intervals and the Matthews correlation
> coefficient with percentile bootstrap intervals. Because many records are closely related
> variants from the same experimental series, intervals are computed both by resampling records and
> by resampling whole sequence families, and both are reported. Records were grouped by sequence
> family and by study from the curated provenance, and designed negatives inherit the group of the
> record from which they were derived. Parameter selection, where performed, was confined to the
> training portion of each split, under repeated stratified five-fold cross-validation,
> leave-one-family-out and leave-one-study-out schemes."

### New subsection — *Comparison with existing tools*

> "**Comparison with existing tools.** HSeeker was compared with triplex (v1.50.0), Triplexator
> (v1.3.2), NeSSie and nBMST on the same 130 records. Each tool was run at its documented defaults;
> with its geometry matched to HSeeker's where it exposes the corresponding parameters; and with a
> small parameter grid selected inside training folds only. A record was scored as positive when a
> tool reported at least one qualifying hit within it, the qualifying hit type being a triplex of
> any type for triplex, a duplex target site for Triplexator, a triplex-potential motif for NeSSie,
> and a mirror repeat classified as triplex-forming for nBMST; for NeSSie and nBMST the alternative
> hit types, a mirror motif and an unclassified mirror repeat respectively, are reported alongside.
> Tool versions, command lines and the hit rule for each tool are given in Supplementary Table 3,
> and the per-record calls in Supplementary Table 4. The settings and the statistical comparison
> were fixed before the tools were run."

### New subsection — *Genome-wide comparison with single-stranded DNA maps*

> "**Genome-wide comparison with single-stranded DNA maps.** We compared HSeeker's predictions with
> the permanganate/S1 footprinting data of Kouzine and colleagues, using their published non-B DNA
> motif calls and single-stranded DNA signal tracks for activated and resting mouse B cells and for
> human Raji cells, on the genome assemblies used in that work. Both genomes were scanned with
> HSeeker at its default settings.
>
> Recall was computed as the fraction of the motifs they called formed that are overlapped by an
> HSeeker locus. Selectivity was computed among their own predicted H-DNA motifs, as the odds that
> a motif is called formed when HSeeker reports it against when HSeeker does not, with the
> mappability of each motif to 36-nt reads used as a sensitivity analysis.
>
> Enrichment of the measured signal was assessed against composition-matched controls. For each
> HSeeker locus, control windows were drawn on the same chromosome with the same length, the same
> A/T count and purine/pyrimidine asymmetry within one count, containing no ambiguous base and
> overlapping no HSeeker locus or predicted H-DNA motif. Two analyses were performed. The first used
> all loci with controls drawn genome-wide. The second restricted both loci and controls to
> sequence that 36-nt reads can map to, excluded intervals within 100 bp of a predicted
> G-quadruplex, Z-DNA, cruciform or stress-induced melting motif and within 2 kb of an annotated
> transcription start site, and drew ten controls per locus within 50 kb. Enrichment is reported as
> the ratio of mean signal with a bootstrap interval over loci, together with the probability that
> a locus exceeds one of its own controls. Two checks accompany each comparison: one control per
> locus was compared against that locus's remaining controls, and the motifs called formed by the
> original authors were passed through the identical procedure."

### Revise — *Computational benchmarking*

> "**Computational benchmarking.** Runtime and peak resident set size were measured on a single
> node of the TACC Lonestar6 cluster with 16 virtual cores of an AMD EPYC 7763 processor, 32 GB of
> RAM and Rocky Linux 8.10, with no other workload on the node. Each configuration was run three
> times and wall-clock time and peak memory were recorded externally. Inputs were human chromosome 1
> and the hg38 primary assembly, both obtained from the UCSC Genome Browser. Comparators were timed
> on the same node under the same protocol."

### Remove — *Genomic context validation*

Delete the paragraph describing the *E. coli* insertion experiment and the 80% overlap criterion,
and the sentence in the preceding paragraph that introduces it.

## 2.4 Results

### Rewrite — *Computational performance and parallel scalability*

Replace the HSeeker paragraph's figures and the comparator paragraph in full.

> "We benchmarked the complete pipeline — scanning with the purity prefilter, overlap removal and
> stability scoring — on human chromosome 1 (hg38, 248.9 Mbp) using 1 to 16 worker threads, with
> all parameters at their defaults. All runs produced an identical set of 45,927 non-overlapping
> H-DNA loci, confirming deterministic output across every level of parallelism. With a single
> worker the pipeline completed in 173.8 s; with 2, 4, 8 and 16 workers it completed in 89.6 s,
> 47.4 s, 26.7 s and 16.3 s, corresponding to speedups of 1.94-, 3.67-, 6.51- and 10.65-fold
> (Figure 2). Parallel efficiency remained above 90% through four workers and was 81% at eight and
> 67% at 16, the decline at higher thread counts reflecting the serial phases that follow the scan.
> Peak memory was essentially flat, between 1,330 and 1,336 MB across all thread counts.
> Score-based overlap selection, which defers overlap resolution until after scoring, took 23.9 s
> at 16 workers and reported 47,217 loci.
>
> Rather than extrapolating from chromosome 1, we ran the complete hg38 primary assembly. At 16
> workers it completed in 204.5 s, that is 3.4 min, with a peak memory of 2,296 MB, reporting
> 559,601 non-overlapping H-DNA loci; three independent runs produced identical output.
>
> For comparison we timed four existing tools on the same chromosome, the same node and the same
> day (Figure 2C). Triplexator, which provides native multi-threading, completed in 45.0 s at 16
> threads and reported 275,476 duplex target sites. The triplex package is single-threaded; split
> into 1 Mbp windows and run as 16 parallel processes it completed in 64.7 s and reported 81,906
> hits. nBMST and NeSSie provide no parallel mode and completed in 96.5 s and 295.7 s on one core,
> reporting 147,222 mirror repeats, of which 33,519 were classified as triplex-forming, and 720,811
> motif occurrences respectively. Peak memory was 1,714 MB for Triplexator, 660 MB for each triplex
> worker process, 866 MB for nBMST and 1,056 MB for NeSSie, against 1,334 MB for HSeeker. The hit
> counts reported by triplex and Triplexator reproduce our earlier measurements exactly.
>
> The comparison should be read with the differences in target class in mind: Triplexator enumerates
> duplex target sites for a separately supplied third strand and NeSSie enumerates symmetric
> k-mers, neither of which is the same object as an intramolecular mirror repeat, and the loci
> counts are correspondingly not comparable."

### New subsection — *Algorithmic complexity*

Move the existing complexity analysis here from Methods unchanged, apart from the lead-in sentence.

### New subsection — *Detection performance on the revised benchmark*

Use the first manuscript block quoted under Reviewer 2, Comment 1, followed by the error-analysis
paragraph below and the score-threshold block quoted under Reviewer 3, Comment 3.

> "The errors are informative. Mirror-repeat detection alone, before the purity requirement is
> applied, calls seven of the 65 non-forming sequences: the experimentally non-forming *c-MYC*
> control B-33, one G-quadruplex-forming sequence, two Z-DNA-forming (CG)ₙ repeats and all three
> perfect mirror repeats of balanced composition. The purity requirement removes all seven without
> losing a forming sequence, which is the clearest demonstration in the benchmark that mirror
> symmetry alone is not sufficient and that the composition criterion is doing the discriminating
> work. No designed negative of any class is called at the default settings."

### New subsection — *Family- and study-aware validation*

Use the block quoted under Reviewer 3, Comment 2, preceded by the group counts.

### New subsection — *Component ablation and parameter robustness*

Convert the response text under Reviewer 3, Comment 9 into Results prose; the numbers are
unchanged. Suggested opening:

> "To establish which components of the pipeline account for its behaviour, we evaluated the stages
> cumulatively on the benchmark, calling a record forming when at least one candidate survives the
> stage in question."

### New subsection — *Comparison with existing tools*

Use the table and the target-class paragraph from Reviewer 3, Comment 6, then:

> "Using paired differences with confidence intervals, HSeeker outperformed all four tools on the
> balanced benchmark: the difference in MCC was +0.15 (95% CI 0.03–0.22) against triplex, with 10
> records called correctly by HSeeker and none by triplex alone (exact McNemar p = 0.002), +0.39
> (0.18–0.56) against NeSSie, +0.52 (0.40–0.63) against nBMST and +0.95 (0.71–1.11) against
> Triplexator. The advantage over triplex comes from the designed negatives: on the 70 experimental
> records alone the difference was +0.06 with an interval reaching zero (0.00–0.22), and on that
> set the two tools are comparable. Against the other three tools the advantage holds on the
> experimental records as well.
>
> In the homopolymer panel, triplex and NeSSie call four of the five A/T homopolymers forming, and
> both call poly-(dC)₃₀ forming; Triplexator calls none; nBMST calls poly-A₃₀ only; and HSeeker
> calls none, which is correct for the five A/T homopolymers and incorrect for poly-(dC)₃₀."

### New subsection — *Genome-wide comparison with single-stranded DNA maps*

Convert the response text under Reviewer 3, Comment 8 into Results prose and cite Figure 3. Close
with:

> "Taken together, these results associate HSeeker loci in activated mouse B cells with
> single-stranded DNA beyond what composition, mappability, local context, the other predicted non-B
> structures and promoter transcription account for, and the association strengthens with
> mirror-repeat arm length, as H-DNA formation should. We did not observe the same enrichment in
> human Raji cells, where HSeeker's calls are nonetheless strongly enriched among the motifs found
> formed. The per-locus effect is modest even in mouse: in a given cell state most HSeeker loci,
> like most control windows, carry no single-stranded signal, so the loci are candidates rather
> than sites demonstrated to be structured, and single-stranded signal is not by itself evidence of
> triplex formation at an individual locus."

### Remove

The paragraph reporting the *E. coli* insertion results, the sentence reporting MCC = 0.73 from it,
the sentences describing the ROC curve and the Youden-selected operating point, and the sentence
describing Figure 3C and 3D.

### Revise — *HSeeker web interface*

Three edits. Use the single-interface sentence quoted under Reviewer 2, Comment 11. Update the
preset values, since the defaults changed. Remove the "Stable sites (score ≥ 50)" tile from the
description of the results page, because the threshold is withdrawn.

## 2.5 Discussion

Apply the three rewritten sentences quoted under Reviewer 3, Comment 12. Then replace the
benchmarking paragraph's opening, which currently claims improved detection performance against two
tools, with a sentence that matches the new comparison; and add the Limitations section quoted
under Reviewer 2, Comment 4 as a numbered subsection before the closing paragraph.

> "We benchmarked HSeeker against four existing tools for triplex and mirror-repeat detection —
> triplex, Triplexator, NeSSie and nBMST — on a balanced set of experimentally characterised and
> designed sequences, and against genome-wide maps of single-stranded DNA. HSeeker achieved the
> highest accuracy of the five on the balanced benchmark and matched triplex, the only other tool
> designed for intramolecular triplexes, on the experimentally characterised sequences."

## 2.6 Data Availability and Code Availability

Replace Data Availability with the block quoted under Reviewer 2, Comment 2. In Code Availability,
keep the repository, PyPI and web-application links, and add the version tag of the release that
corresponds to the revised manuscript.

---

# Part 3 — Figures and tables

## Figure 1 — workflow

Simplification is requested, and two panels are now factually out of date.

- **Panel A** (candidate-identification flowchart) omits the two steps added since submission: the
  inward extension that absorbs matching spacer-end pairs into the arms, and the deduplication of
  candidates that normalise to the same partition. Add both as terminal steps before "Return
  candidate list", and reduce the decision-diamond chain, which is the densest part of the figure.
- **Panel B** contains a "Sensitivity Analysis / AUC" inset reporting the in-sample AUC comparison.
  That analysis is withdrawn; **remove the inset**. Removing it also relieves the crowding the
  reviewer noted.
- **Panel B**, middle row, labels the stability balance as "H-DNA Stability"; relabel as "Stability
  score (heuristic ranking)" so the figure does not reintroduce the free-energy reading.
- The "Detection → H-DNA-forming ✓ / non-forming ✗" decision tree in Panel B implies a
  score-based decision; relabel as detection only.

## Figure 2 — parallel scaling

All numbers change, and the comparator sweeps cannot be reproduced as drawn. The original four
panels plot 1–16 worker sweeps for all three tools; we have a full sweep for HSeeker and
single-point measurements at 16 cores for the four comparators, and two of the comparators have no
parallel mode at all, which makes a sweep meaningless for them.

**Recommended three-panel replacement:**

- **A** — HSeeker wall-clock time and speedup against worker count (1, 2, 4, 8, 16), with the ideal
  speedup line. Data: 173.8, 89.6, 47.4, 26.7, 16.3 s; speedups 1.00, 1.94, 3.67, 6.51, 10.65.
- **B** — HSeeker peak memory against worker count, flat between 1,330 and 1,336 MB, with the
  whole-genome run marked at 2,296 MB.
- **C** — grouped bars at 16 cores for all five tools: wall-clock time and peak memory, with the
  number of reported loci printed on each bar and a footnote that the loci counts are not
  comparable across tools because the target classes differ.

If the four-panel design is preferred, the comparators would need sweeps at 1, 2, 4, 8 and 16
cores; this is feasible for triplex and Triplexator, and not meaningful for nBMST and NeSSie.

**Caption:**

> "**Figure 2: Computational performance of HSeeker.** **A)** Wall-clock time and speedup for the
> complete pipeline on human chromosome 1 (hg38, 248.9 Mbp) across 1–16 worker threads, mean of
> three runs; the dashed line is ideal linear speedup. **B)** Peak resident set size across the same
> sweep, with the complete hg38 primary assembly shown for reference. **C)** Wall-clock time and
> peak memory at 16 cores for HSeeker, Triplexator, triplex, nBMST and NeSSie on chromosome 1,
> measured on the same node under the same protocol; nBMST and NeSSie provide no parallel mode and
> were run on one core. The number of loci reported by each tool is printed on the corresponding
> bar; these counts are not directly comparable, because the tools enumerate different target
> classes."

## Figure 3 — replaced in full

The *E. coli* figure is withdrawn. The replacement is generated and deposited. It has four panels:
recall against the published calls; selectivity as odds ratios; signal at HSeeker loci relative to
matched controls for both analyses, with the placebo and positive-control checks; and the
arm-length stratification.

**Caption:**

> "**Figure 3: HSeeker predictions compared with single-stranded DNA mapped in cells.** HSeeker was
> run at its default settings on the genome assemblies used by Kouzine and colleagues. **A)** Recall:
> the percentage of the H-DNA motifs they found single-stranded that are overlapped by an HSeeker
> locus, in activated mouse B cells and in human Raji cells. **B)** Selectivity: the odds ratio that
> one of their predicted H-DNA motifs is found single-stranded when HSeeker reports it against when
> it does not, for all analysed motifs and for motifs that 36-nt reads can map to. **C)**
> Single-stranded signal at HSeeker loci relative to control windows matched for length, A/T content
> and purine/pyrimidine asymmetry. "All loci" uses genome-wide controls; "isolated loci" restricts
> loci and controls to mappable sequence more than 100 bp from other predicted non-B motifs and more
> than 2 kb from annotated transcription start sites, with controls drawn within 50 kb; "authors'
> formed motifs" applies the identical procedure to the motifs they called formed; "placebo"
> compares one control against the remaining controls of the same locus. **D)** Isolated loci
> stratified by mirror-repeat arm length. Points are ratios of mean signal, bars are 95% bootstrap
> confidence intervals, and the dashed line marks no difference. Scales in B–D are logarithmic."

## Supplementary Figure 1 — regenerate

Currently three confusion matrices on the 79-sequence set. Regenerate as five confusion matrices —
HSeeker, triplex, NeSSie, nBMST, Triplexator — on the 130-record benchmark at each tool's default
settings, with MCC printed on each panel. The underlying per-record calls are already deposited.

## Supplementary tables — renumber

The current manuscript cites "Supplementary Table 1" for two different things. Proposed numbering:

| Table | Content |
|---|---|
| Supplementary Table 1 | Benchmark records: sequence, label, provenance, curation decision, family and study identifiers, and assay conditions |
| Supplementary Table 2 | Stability-score pairing and stacking values |
| Supplementary Table 3 | Comparator versions, command lines and hit rules |
| Supplementary Table 4 | Per-record calls for all five tools on the benchmark |
| Supplementary Table 5 | Homopolymer panel across tools and filter settings |

Update the two in-text citations in the Methods accordingly.

## New main-text table

**Table 1 — detection performance with confidence intervals.** The three-row table given under
Reviewer 3, Comment 1. This is the table the comparison claims refer back to, and the paper
currently has no numerical table in the main text.

---

# Part 4 — LaTeX and algorithm panels

Three `.tex` files are in the working tree and their relationship needs to be settled before the
revision is compiled.

| File | Status |
|---|---|
| `algorithms.tex` | The **submitted** panels: three algorithms, greedy overlap removal only. Superseded. |
| `main.tex` | The **current** panels: four algorithms — pipeline ordering, candidate enumeration with inward extension and deduplication, both overlap rules, and scoring. |
| `scoring.tex` | A byte-identical copy of the Phase 3 panel already inside `main.tex`, not included anywhere. |

**Action:** compile the revision from `main.tex`; delete `scoring.tex`, which is a duplicate; and
either delete `algorithms.tex` or move it out of the build so that the superseded panels cannot be
compiled by accident.

The differences between the submitted and current panels are substantive, and each needs a
corresponding sentence in the Methods so the text and the panels agree:

1. **A pipeline panel now precedes the stage panels** and fixes the order: enumeration, then
   longest-arm overlap removal if selected, then AT filtering, then scoring, then homopolymer
   filtering, then score-based overlap removal if selected. The submitted version had no such panel,
   and the Methods described filtering as preceding overlap removal in all cases, which is true only
   for the score-based rule.
2. **Inward extension and deduplication** are new steps in the enumeration panel. Matching base
   pairs at the two ends of the spacer are absorbed into the arms when purity permits, and
   candidates that normalise to the same (start, end, arm) partition are reported once. Add to the
   Methods paragraph on candidate discovery:

   > "After enumeration, matching base pairs at the inner ends of the spacer are absorbed into the
   > arms where the purity criterion permits, so that a motif is reported with its maximal arm
   > representation, and candidates arising from different centre and spacer combinations that
   > normalise to the same arm and spacer partition are reported once."

3. **The overlap panel now covers both rules.** The Methods sentence describing overlap resolution
   as a greedy sweep is replaced by the three-mode description under Reviewer 3, Comment 10.
4. **The scoring panel no longer returns nothing when no window fits**; it scores the detected arms
   as they stand. It also clamps the total to zero and reports lower-case sequence. Items 2 and 3 of
   §2.3 *Stability score* cover the text side.
5. **The purity prefilter is optional in the panel** (`rmq` flag) where the submitted version made
   it unconditional. The Methods should say that it is an exact prefilter that never changes the
   output and is enabled by default.

---

# Part 5 — Open items for the authors

These need a decision before submission; none blocks drafting.

1. **One benchmark record rests on inference rather than experiment.** The *XDP* (CCCTCT)₁₀ record
   is labelled forming on the basis of a review that states explicitly that no in vitro or in vivo
   structural data exist for that repeat. It is currently in the construct-level set. Recommend
   moving it to the class-level set, which keeps it in the benchmark but out of the primary figures
   and removes the inconsistency with the statement that construct-level records have experimental
   support.
2. **Six further records carry caveats** recorded with the benchmark and worth a glance before the
   supplementary table is finalised: two whose cited paper reports a recombination assay rather than
   a structural one, one whose triplex was shown on oligonucleotides at pH 8 rather than on
   supercoiled DNA, and three mouse records from footprinting at pH 4.5 with zinc where the original
   authors state they cannot exclude formation after lysis. No label needs changing; the conditions
   columns already record all of this.
3. **The published chromosome 1 count is internally inconsistent in the submitted manuscript.** The
   Results state default parameters with a minimum arm length of 8 bp and report 96,729 loci, but
   that count corresponds to a minimum arm length of 10 bp with the composition filters disabled.
   The revised text reports 45,927 at the stated defaults and 106,481 with filters disabled; the
   manuscript should not carry the old number forward in any form.
4. **Whether to extend the genome-wide comparison to the other tools.** Running the four comparators
   through the same recall and selectivity analysis against the published single-stranded DNA maps
   would be the strongest available evidence, because none of the tools was developed against those
   data. It is roughly half a day of computation and is not required by any reviewer comment.
5. **Release mechanics.** The default minimum arm length changed, which changes the output of
   existing user workflows. This warrants a minor version bump, a changelog entry naming the change,
   and a note to collaborators before the release is published.
