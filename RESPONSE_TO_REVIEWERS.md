# Response to Reviewers

**Manuscript:** HSeeker — an algorithm for systematic H-DNA sequence identification
**Journal:** *Briefings in Bioinformatics*
**To:** Dr Shuangge Ma (Editor in Chief) and Dr Zhenxia Chen (Deputy Editor)

We thank all three reviewers for their careful reading and for comments that have substantially
improved the work. The revision rests on three changes that affect most of the reported numbers, so
we summarise them once here rather than repeating them under every comment.

1. **The benchmark has been rebuilt and balanced.** Every experimental record was traced to its
   primary publication; records that could not be traced, and records annotated only from a
   secondary compilation, were removed. The benchmark now contains 130 sequences, 65 forming and 65
   non-forming, comprising 70 experimental records and 60 designed negatives.
2. **The defaults have changed and the interfaces have been unified.** The minimum arm length is now
   8 bp and the mismatch tolerance 0.10, and the composition filters are applied identically by the
   Python package, the command-line tool and the web application, which now return identical results
   for identical input.
3. **The Youden-selected score threshold has been withdrawn.** A sequence is called on detection
   alone, and the stability score is reported as a ranking. The genomic-insertion experiment is
   retained but narrowed to what it measures, and a new comparison against permanganate/S1
   single-stranded DNA maps has been added.

Reviewer comments are reproduced in full, followed by our response and, where a comment maps onto a
specific edit, the text as it now stands in the manuscript.

---

# Reviewer 1

> HSeeker, an algorithm for systematic H-DNA sequence identification, is well presented, its
> structural features recognized are properly made clear and tested, a suitable comparison with the
> two well known competitor tools is done and shows the strong performance, Results are nicely
> separated from an intelligent discussion; All data and scripts are made fully available; language
> is nice and well written manuscript.

We thank the reviewer for the positive assessment and interest in our work. No changes were
requested. We note that the revisions made in response to the other two reviewers — an expanded and
balanced benchmark, family- and study-aware validation, a wider comparison with existing tools, a
genome-scale comparison against experimental data, and re-measured runtimes — also bear on the
aspects this reviewer highlights, and we believe they strengthen the manuscript further.

---

# Reviewer 2

> I have read the manuscript by Provatas et al. subjected to my review. In general my opinion about
> the whole work is positive, although some major flaws are visible. Therefore I advise major review
> of the manuscript.

We thank the reviewer for the overall positive assessment and provide a point by point response
below.

---

## Comment 1 — discriminating power and the size of the negative set

> The authors show the time efficiency and scalability of their approach, but they fail to
> unequivocally show the discriminating power of the tool — the main reason is the size of the
> benchmark dataset — 79 cases with only 9 negative. The class imbalance makes any accuracy
> comparison unreliable — the model which always answers 'true' on this dataset would have a similar
> accuracy. … One way to deal with it is to try the tool on random sequences, or those, which have
> structures in PDB already … Another step is to test mirror sequences. If there are no more than 9
> such sequences, one might try to create artificial structures using BoltzGen, or RFDiffusion.
> Alternatively, one might artificially perturb the positive cases in such way to disrupt formation
> of H-DNA. Also, including poly-A, or similar cases, which are known not to form H-DNA should be
> included, as a hard test if the implementation really covers this case. One could also use
> sequences with other known structures, like B-DNA, or DNA quadruplexes.

We thank the reviewer and agree that the original benchmark had too few negative controls to
establish discriminating power, and that accuracy was misleading given the class imbalance. We have
therefore substantially revised the evaluation.

We traced all experimental records to their primary publications. Using the 70 records that
survived, we had only five non-forming sequences. To address this limitation, we followed the
reviewer's suggestion and added 60 designed negatives, generated through a documented, reproducible
procedure under a fixed random seed: 22 positives whose mirror symmetry was destroyed by
purine-preserving substitutions, 11 dinucleotide-preserving shuffles, 12 random sequences matched to
positives in length and GC content, four G-quadruplex-forming sequences, three Z-DNA-forming (CG)ₙ
repeats, five random B-DNA controls, and three perfect mirror repeats lacking purine richness. The
resulting benchmark is balanced at 65 positives and 65 negatives. The last group is the sharpest of
these tests: those sequences satisfy mirror symmetry but not purine richness, so a tool that keys
only on symmetry will call them.

Regarding BoltzGen and RFdiffusion, these tools cannot reliably determine whether a DNA duplex
adopts H-DNA under specific experimental conditions and therefore cannot provide validated negative
labels. Our perturbation-based controls address the reviewer's concern more directly.

We also evaluated homopolymers separately, since their classification is largely determined by
composition filters rather than mirror-repeat detection, and because a run of one base is
mirror-symmetric at every position and so does not test mirror-repeat detection. HSeeker rejects
poly-A₂₀, previously shown to remain B-form, along with the other poly-A and poly-T controls. We
explicitly acknowledge that it also rejects poly-(dC)₃₀, despite its ability to form an
Mg²⁺-dependent triplex.

On the balanced benchmark, HSeeker identified 63 of 65 forming sequences and rejected all 65
negatives, achieving 96.9% sensitivity, 100% specificity, and an MCC of 0.970. The two missed
positives are tandem direct repeats outside the model's intended scope. Importantly, mirror symmetry
alone produced seven false positives, all eliminated by the purine-richness requirement without
sacrificing sensitivity.

**The Results section of our manuscript now reflects the changes above:**

> HSeeker detected 63 of the 65 H-DNA-forming sequences and rejected all 65 non-forming sequences,
> corresponding to a sensitivity of 0.969 (95% CI 0.895–0.992), a specificity of 1.000 (0.944–1.000)
> and a Matthews correlation coefficient of 0.970 (0.925–1.000) (Table 1). The two sequences that
> were not detected, the *c-MYC* nuclease-sensitive element construct (ACCCTCCCC)₄ and the γ-globin
> −228/−189 element, are tandem direct repeats rather than mirror repeats and therefore fall outside
> the structural model that HSeeker implements. Restricting the evaluation to the 70 experimental
> records gives a sensitivity of 0.969 (0.895–0.992) and a specificity of 1.000 (0.566–1.000); the
> width of the latter interval reflects the five experimental non-forming constructs that survived
> curation, and conclusions about specificity rest largely on the designed negatives.
>
> The errors are informative. Mirror-repeat detection alone, before the purity requirement is
> applied, calls seven of the 65 non-forming sequences: the experimentally non-forming *c-MYC*
> control B-33, one G-quadruplex-forming sequence, two Z-DNA-forming (CG)ₙ repeats and all three
> perfect mirror repeats of balanced composition. The purity requirement removes all seven without
> losing a forming sequence, which is the clearest demonstration in the benchmark that mirror
> symmetry alone is not sufficient and that the composition criterion is doing the discriminating
> work. No designed negative of any class is called at the default settings.

**We also provide the updated Methods section below:**

> **Benchmark dataset.** The benchmark contains 130 sequences, 65 reported to form H-DNA and 65
> reported or designed not to form it. Seventy are experimental records, 65 forming and five
> non-forming. Each was traced to its primary publication; records whose annotation derived only
> from a secondary compilation, and records with no traceable primary source, were removed.
> Fifty-seven of the retained records are construct-level, meaning the published experiment tested
> the sequence in the benchmark, and are used for the primary performance figures. The remaining
> thirteen are class-level, meaning the published experiment tested a closely related construct of
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
> For each experimental record we recorded the assay used, the DNA context, and the pH, temperature,
> ionic conditions and supercoiling state under which the structure was assessed, together with the
> passage in the source supporting those values; values that the source does not report are marked
> as such. Sixty-four of the 70 experimental records have a reported pH. Because the labels are
> outcomes under heterogeneous conditions, the pooled performance figures should be read as
> discrimination of reported sequence labels across assay conditions rather than as propensities
> under a single common condition.

---

## Comment 2 — reproducibility of the sequences and a script to run the analysis

> The results are non-reproducible, as the sequences cannot be retrieved from the supplement. It
> would be much better if the sequences would be given as a csv file, which would be much easier to
> parse. Also, it would be good if in the repo there was a script to run the analysis on the 79 test
> cases.

We thank the reviewer for pointing this out, and we agree. The benchmark is distributed as two
machine-readable files. The first lists all curated records with their provenance: the primary
citation, the curation decision and its justification, the family and study identifiers used for
grouping, and the assay, DNA context, pH, temperature, ionic conditions and supercoiling state under
which each construct was tested, together with the passage in the source that supports those values.
The second is the balanced 130-record benchmark used for the reported metrics. The designed
negatives are regenerated from a single script under a fixed seed.

The repository contains scripts that reproduce every analysis in the paper end to end — the
benchmark metrics, the family- and study-aware validation, the component ablation and parameter
sweeps, the comparison against the other tools, the homopolymer panel, the overlap-rule comparison,
the genomic-insertion analysis and the genome-wide comparison — with all settings pinned in a single
configuration file.

**The manuscript now provides in Data Availability:**

> The benchmark sequences, their provenance and the assay conditions under which each construct was
> characterised are provided as machine-readable files in the repository and as Supplementary Table
> 1. Scripts that regenerate every analysis and figure in this work, together with the parameter
> settings used for each, are provided in the same repository. The human chromosome 1 and hg38
> primary assembly sequences used for performance benchmarking were obtained from the UCSC Genome
> Browser. The single-stranded DNA maps and non-B DNA motif calls used for the genome-wide
> comparison are those published by Kouzine and colleagues. Raw benchmark commands, timing output
> and result files are included in the repository and will be archived in a public repository with a
> permanent identifier on publication.

---

## Comment 3 — the E. coli insertion experiment and validation on a genome tested for non-B DNA

> The results on introducing the sequences into E.Coli genome are not convincing. In particular, the
> experiment does not prove the tool finds the H-DNA, but rather retrieves correctly the sequences on
> which it is fitted to work good. But the "genetic context" might disallow to create H-DNA entirely
> in those cases. Also, were there any other cases obtained in the modified E.Coli genome? In
> general, to strengthen the analysis, the authors might run their tool not on the modified E.Coli
> genome, but on a genome which was tested to contain non B-DNA sequences. For example in the paper
> of Kouzine (ref 4) such a analysis was done. If the set of sequences found by HSeek would be a
> subset of all non-B-DNA sequences, that would be a biologically sound result!

We thank the reviewer for this comment, and we agree with their reading of the experiment. We have
acted on it in two ways.

**First, on the insertion experiment itself.** We accept that inserting the benchmark sequences into
the *E. coli* genome does not test whether HSeeker's predictions correspond to structures that form.
We have therefore removed every claim of biological validation from it and narrowed it to the
question it can answer: whether a sequence is called the same way once embedded in a genome as when
scanned on its own. Presented that way it is informative. HSeeker returns the identical call for all
130 benchmark records in both settings, and its performance in the genomic background is unchanged
(63 of 65 forming sequences detected, no false positive, MCC 0.970). The comparator does not: triplex
returns a different call for 44 of the 130 records, in every case a detection that is present when
the sequence is scanned alone and absent once it is embedded.

On the reviewer's direct question about other cases in the modified genome: the unmodified *E. coli*
genome yields 83 loci across its 4.64 Mb, or 17.9 per Mb, and the hits in the modified genome are
accounted for by those background loci together with the recovered insertions. No inserted
non-forming sequence is called in either setting. This analysis is now Figure 3, and the Methods and
Results present it as a property of the software rather than as evidence about structure formation.

**Second, on the biological validation.** Following the reviewer's suggestion, we examined HSeeker
against the genome-wide experimental data of Kouzine et al. (PMID 28237796), who mapped
single-stranded DNA in mouse B cells and human Raji cells. Using their published annotations and
signal tracks, we scanned both genomes with HSeeker at default settings.

HSeeker recovered 79.0% of experimentally formed H-DNA loci in mouse cells and 82.3% in human cells.
Moreover, motifs identified by HSeeker were significantly more likely to show experimental evidence
of formation than those it excluded, with odds ratios of 1.65 in mouse B cells and 3.39 in Raji
cells. These associations became stronger when restricting the analysis to regions reliably mappable
by the study's sequencing reads.

We note that requiring predictions to be a strict subset of experimentally formed structures is not
realistic for a sequence-based predictor, since formation depends on cellular conditions. Indeed,
only around 3–5% of the original study's predicted motifs showed evidence of formation in a given
cell type. Additionally, permanganate/S1 nuclease footprinting identifies other non-B DNA structures
and single-stranded sequences, which are far more ubiquitous in the genome. The relevant measures
are therefore enrichment for experimentally formed structures and recovery of known formed loci,
both of which support HSeeker's performance.

Most missed loci were short A/T-rich or homopolymeric motifs excluded by our composition filters.
These motifs formed substantially less frequently in the experimental data, supporting their
exclusion by default. Nevertheless, they account for approximately 17–19% of formed loci, a
limitation we now explicitly acknowledge.

**The paper changes in the following passage reflect this new approach (Results, Comparison with
single-stranded DNA maps):**

> We compared HSeeker's predictions with the permanganate/S1 footprinting data of Kouzine and
> colleagues, who mapped single-stranded DNA in activated and resting mouse B cells and in human
> Raji cells, predicted non-B DNA motifs genome-wide, and called a motif formed when it carried
> single-stranded signal. HSeeker loci overlap 79.0% of the H-DNA structures found to be formed in
> mouse (13,517 of 17,109) and 82.3% in human Raji cells (6,879 of 8,356). Among their predicted
> H-DNA motifs, those HSeeker calls are formed more often than those it does not: 4.5% against 2.8%
> in mouse (odds ratio 1.65, 95% CI 1.59–1.71) and 5.7% against 1.7% in human (3.39, 3.20–3.59);
> restricting to motifs that the 36-nt reads can map to raises these to 1.68–2.07 and 4.16–5.70
> (Figure 4A, B).
>
> Taken together, these results associate HSeeker loci in activated mouse B cells with
> single-stranded DNA beyond what composition, mappability, local context, the other predicted non-B
> structures and promoter transcription account for, and the association strengthens with
> mirror-repeat arm length, as H-DNA formation should. We did not observe the same enrichment in
> human Raji cells, where HSeeker's calls are nonetheless strongly enriched among the motifs found
> formed. The per-locus effect is modest even in mouse: in a given cell state most HSeeker loci, like
> most control windows, carry no single-stranded signal, so the loci are candidates rather than sites
> demonstrated to be structured, and single-stranded signal is not by itself evidence of triplex
> formation at an individual locus.

Beyond recovery and selectivity, we also asked whether the measured single-stranded signal is
elevated at HSeeker loci. Because that signal depends strongly on base composition, control windows
were matched to each locus on length, A/T count and purine/pyrimidine asymmetry. Across all loci,
those in activated mouse B cells carried 14% more signal than their matched controls (ratio of means
1.137, 95% CI 1.127–1.147). Restricting loci and controls to mappable sequence away from other
predicted non-B motifs and from transcription start sites raised this to 25% (1.254, 1.179–1.332),
and the excess increased with mirror-repeat arm length, as H-DNA stability should but the other
structures should not. We did not observe the enrichment in human Raji cells. The full analysis,
with its controls and checks, is presented in the Results and in Figure 4.

---

## Comment 4 — simplification of the thermodynamic approach, and a limitations section

> The thermodynamical approach is interesting, but it is a major simplification. The geometry of the
> DNA backbone is not included, the energy of the structure is not just the algebraic sum, the
> mismatch have different values (as authors already note), not universally equal 5kcal/mol etc. The
> authors should be honest about the limitations of the tool. In general, there should be a section
> devoted to tool limitation.

We agree with the reviewer that our thermodynamic scoring approach has limitations. As the reviewer
noted, the energy of an entire DNA structure cannot be reduced to a simple algebraic sum of isolated
parts, and our calculation does not explicitly list backbone geometry, superhelical stress, or
complex environmental variables as separate dynamic parameters.

To clarify, the underlying interaction and stacking energy parameters are drawn from the *ab initio*
thermodynamical calculations by Park et al. (2021), which systematically evaluated nucleobase
interactions (Watson–Crick pairing, Hoogsteen pairing, and base stacking) under a water continuum
model. That measurement isolated nucleobase interactions by replacing backbones with hydrogen atoms
to compute fundamental interaction energies and left out other factors, while explicitly recognising
that factors like torsional stress, backbone constraints, and structural distortions must be
evaluated by comparing model energies against experimental structures.

Predicting nucleic acid structure involves multi-factor physics that goes beyond a fast
sequence-based tool and requires genome-wide experimental supporting data, and we agree that these
constraints must be stated transparently. We have added a Limitations section. The score is now
described throughout as an experimentally informed heuristic rather than a free-energy calculation,
and all kcal/mol labels have been removed from the score itself.

**The manuscript now reads (Discussion, Limitations):**

> HSeeker has several limitations that bear on how its output should be read. The stability score is
> an experimentally informed heuristic, not a calculated free energy. It assigns fixed weights to
> matched G and A positions, a uniform penalty to every other arm pairing, a constant reward to
> adjacent matched positions and an escalating penalty to runs of disrupted ones, and it then selects
> the highest-scoring arm window. The underlying interaction and stacking energy parameters are drawn
> from the *ab initio* thermodynamical calculations by Park et al., which systematically evaluated
> nucleobase interactions (Watson–Crick pairing, Hoogsteen pairing, and base stacking) under a water
> continuum model. That measurement isolated nucleobase interactions by replacing backbones with
> hydrogen atoms to compute fundamental interaction energies and left out other factors, such as
> torsional stress, backbone structural distortions, superhelical stress, pH-dependent protonation,
> ionic composition, temperature, loop or spacer energetics, or competition with duplex and other
> non-canonical conformations. It is therefore a relative ranking of candidates within the scheme
> defined here, and it carries no physical units. Predicting nucleic acid structure involves
> multi-factor measurements beyond what this fast sequence-based tool can include, and requires more
> detailed experimental supporting data in future. We report no score threshold, and a candidate is
> called on detection alone.
>
> The benchmark labels come from assays performed under different conditions, including acidic
> buffers, Mg²⁺-containing buffers and negatively supercoiled plasmids. We record those conditions
> for every record, but a sequence characterised as H-DNA-forming under one set of conditions does
> not necessarily form the structure under physiological conditions, and the pooled performance
> figures should be read as discrimination of reported sequence labels across conditions.
>
> The composition filters are fixed rules rather than fitted parameters, but they were chosen with
> knowledge of the H-DNA literature and of the composition of this benchmark, so the benchmark is not
> an independent test of them. They are not free of cost: the homopolymer filter rejects
> poly-(dC)·(dG) tracts, which form an Mg²⁺-dependent triplex, and the AT-content filter removes a
> class of motifs that constitutes 17–19% of the H-DNA found to be formed in published
> single-stranded DNA maps.
>
> Finally, HSeeker detects antiparallel, purine-motif H-DNA at mirror repeats. Tandem direct repeats,
> parallel triplexes requiring cytosine protonation, and intermolecular triplexes are outside its
> model, and the two forming sequences it misses in our benchmark are direct repeats. Its predictions
> identify sequences with the potential to adopt H-DNA; they are not evidence that the structure
> forms at a given locus in vivo.

---

## Comment 5 — complexity and sensitivity analyses belong in the Results

> The complexity and sensitivity analyses are in my opinion results, not methods, so should be moved
> there.

We agree. The complexity analysis and the sensitivity and specificity analyses have been moved out of
Methods and now appear in the Results, where the Methods sections retain only the design of each
analysis. The Results sections, in order, are: computational performance and parallel scalability;
algorithmic complexity; detection performance on the revised benchmark; family- and study-aware
validation; component ablation and parameter robustness; comparison with existing tools; genomic
context recovery; and the genome-wide comparison with single-stranded DNA maps.

---

## Comment 6 — the filters were defined on the same dataset

> The same 79-sequence dataset is used to define the empirical filters and used next to retrieve the
> sequences. So the tool is overfitted to these exact sequences. If not fixed by the authors, this
> have to be made clear to the reader.

We thank the reviewer; this is a fair description of the original manuscript, and we address it in
three ways.

First, the AT-content and homopolymer rules are fixed pipeline rules rather than parameters fitted
during evaluation, but they were chosen with knowledge of the H-DNA literature and of this
benchmark's composition, and the manuscript now says so.

Second, we removed the circularity from the reported numbers. In the original benchmark the only
false positives at the analysis settings were four homopolymers, and the filters removed them. In the
revised benchmark the homopolymers are reported as a separate panel, and consequently neither filter
changes any benchmark call: the reported metrics do not depend on the filters, and equally the
benchmark offers no evidence in their favour. We say both.

Third, we provide evidence for the filters from outside the benchmark. In the published
single-stranded DNA maps, predicted H-DNA motifs of the class the filters exclude are found formed at
roughly half the rate of other predicted H-DNA in mouse and a quarter in human. With the filters
disabled HSeeker calls almost all predicted motifs and its selectivity disappears.

Separately, every parameter is now selected out of sample. Across 576 combinations of purity,
mismatch tolerance, arm length, spacer length and filter state, with selection confined to the
training folds, per-fold tuning reached an MCC of 0.939 under nested cross-validation against 0.970
for the fixed default configuration, and tuning did not beat the fixed configuration in any of 200
random splits.

**The manuscript now reads (Methods, Sequence composition filters):**

> Two composition filters are applied to detected candidates. Motifs whose full sequence, counting
> both arms and the spacer, is at least 80% A/T are discarded before scoring, because very A/T-rich
> sequences destabilise triplex formation through the geometric incompatibility of Hoogsteen triplets
> and the reduced stability of A·A-T triads. Motifs whose sequence consists of a single repeated base
> are discarded, because a homopurine or homopyrimidine run satisfies the composition criterion
> without presenting a mirror repeat in any meaningful sense.
>
> Both rules are fixed and are not fitted during evaluation. They were, however, chosen with
> knowledge of the H-DNA literature and of the composition of our benchmark, and the benchmark is
> therefore not an independent test of them: because the single-base homopolymers are reported
> separately (see Results), neither filter changes any benchmark call, so the reported metrics neither
> depend on the filters nor support them. Evidence bearing on the filters is drawn instead from
> published single-stranded DNA maps and is presented in the Results.

---

## Comment 7 — confidence intervals in the tool-versus-tool comparison

> The results are not statistically rigorous, not only because of the small dataset, but also because
> of lack of confidence intervals in tool-vs-tool comparison. This should be added if the authors want
> draw conclusions about tool superiority.

We thank the reviewer for the helpful comment, and agree. We have revised the statistical analysis to
ensure that comparisons are supported by confidence intervals rather than point estimates alone. We
now report 95% confidence intervals for all performance metrics, accounting for both individual
records and sequence families to address potential dependence between related sequences. Tool
comparisons also include paired bootstrap intervals and exact McNemar tests.

We consider HSeeker's results high confidence only when the family-level confidence interval for the
MCC difference excludes zero and the advantage remains after matching detection criteria and tuning
comparator parameters. These criteria were fixed before any comparator was run. Under these
conditions, HSeeker outperforms all four tools on the balanced benchmark, with MCC improvements
ranging from +0.15 against triplex to +0.95 against Triplexator.

Importantly, HSeeker's advantage over triplex is driven by the designed negative controls. On
experimental records alone, their performance is statistically comparable (MCC difference +0.06, 95%
CI 0.00–0.23). We have made this distinction explicit and revised the Results and Key Points
accordingly.

**We provide the updated section below (Results, Comparison with existing tools):**

> HSeeker was compared with four existing tools on the same 130 records, each at its documented
> default settings.
>
> | Tool | Sensitivity | Specificity | MCC (95% CI, family) | Tuned out of fold |
> |---|---|---|---|---|
> | HSeeker | 63/65 | 65/65 | 0.970 (0.909–1.000) | 0.939 |
> | triplex | 62/65 | 56/65 | 0.819 (0.744–0.936) | 0.819 |
> | NeSSie | 33/65 | 65/65 | 0.583 (0.421–0.780) | 0.547 |
> | nBMST | 24/65 | 64/65 | 0.449 (0.346–0.560) | 0.597 |
> | Triplexator | 8/65 | 58/65 | 0.024 (−0.137–0.249) | 0.520 |
>
> Using paired differences with confidence intervals, HSeeker outperformed all four tools on the
> balanced benchmark: the difference in MCC was +0.15 (95% CI 0.03–0.22) against triplex, with 10
> records called correctly by HSeeker and none by triplex alone (exact McNemar p = 0.002), +0.39
> (0.18–0.55) against NeSSie, +0.52 (0.41–0.62) against nBMST and +0.95 (0.71–1.11) against
> Triplexator. The advantage over triplex comes from the designed negatives: on the 70 experimental
> records alone the difference was +0.06 with an interval reaching zero (0.00–0.23), and on that set
> the two tools are comparable. Against the other three tools the advantage holds on the experimental
> records as well.
>
> These differences follow from what each tool searches for. Triplexator enumerates duplex target
> sites for a separately supplied third strand and does not require a mirror repeat; at the settings
> used in our original manuscript it detects 64 of 65 forming records but also calls 42 of the 65
> negatives, including all 22 mirror-disrupted mutants. NeSSie reports k-mers that are themselves
> mirror-symmetric rather than arm–spacer–arm motifs; it made no false call but detected none of the
> 16 pXY32 constructs. nBMST's mirror-repeat module matches HSeeker's geometry, but its triplex
> classification accepts only mirror repeats with more purines than pyrimidines on the reported
> strand, so it detected none of the 19 forming records written on the pyrimidine-rich strand, whereas
> its unclassified mirror-repeat output detected 15 of them. triplex, the only other tool designed for
> intramolecular triplexes, detected nearly the same forming records as HSeeker but called nine
> designed negatives forming, seven of them mirror-disrupted mutants.
>
> In the homopolymer panel, triplex and NeSSie call four of the five A/T homopolymers forming, and
> both call poly-(dC)₃₀ forming; Triplexator calls none; nBMST calls poly-A₃₀ only; and HSeeker calls
> none, which is correct for the five A/T homopolymers and incorrect for poly-(dC)₃₀.

---

## Comment 8 — `-march=native`

> On most of the linux and mac distributions the package will be compiled with -march=native. Actually
> ARCHFLAGS is never set on Linux, so -march=native will always be set on Linux. Was that intentional?

We thank the reviewer for this comment. The old cross-compilation check only inspected macOS
`ARCHFLAGS`; on an ordinary Linux build it allowed `-march=native` by default. That was not
appropriate for redistributable binaries, whose users may have a different CPU from the build machine.
We changed `setup.py` so Linux and macOS builds now use portable compiler targets by default (`-O3
-Wall`, without `-march=native`). Native CPU tuning is available only when a user explicitly sets
`HSEEKER_NATIVE_BUILD=1` for a local build; the option rejects a different macOS target architecture.
Our wheel configuration does not enable that option. We compiled the extension in a clean default
build and confirmed that its compiler command contains no `-march=native`, then ran the test suite.

These changes are in the released version accompanying this revision (tag `v1.1.0`).

---

## Comment 9 — unsynchronised profiling counters

> There is no synchronization in mutating prof_inner_iters and prof_ctr_sp_pairs in
> src/hseeker/_hdna.c. It may lead to losing the counting.

We thank the reviewer for identifying this issue. We removed the shared profiling counters from the
scan loops. Each invocation now counts center/spacer pairs and inner-loop iterations in its own
`ScanProfile` while the GIL is released. After the scan reacquires the GIL, it publishes those counts
to per-module state. The existing `profiling_info()` fields describe the last completed scan. The new
cumulative fields and `reset_profiling()` support measurements across parallel scans. We also
corrected the standard scan path, which previously did not increment the inner-loop counter. A
regression test runs 40 standard and purity-prefiltered scans concurrently and verifies that the
cumulative counts equal the sum of independently measured counts.

These changes have been merged in the main branch of the HSeeker repository (commit `063e3b5`) and are
included in the released version accompanying this revision (tag `v1.1.0`).

---

## Comment 10 — case inconsistency in `putative_triplex`

> the left_arm/right_arm/full_sequence come in lowecase, while putative_triplex in uppercase. Am I
> right? If so, simple string comparison would fail for a user not knowing this.

We agree with the reviewer that this was an inconsistency in the output format. We now return
`putative_triplex` with lowercase DNA bases, matching `left_arm`, `spacer`, `right_arm`, and
`full_sequence` in the Python API and CLI output. We added regression tests covering scored sequence
hits, ablation scoring, and TSV output, and documented the casing convention. `putative_triplex` can
still differ from `full_sequence` after case normalisation because the scorer may optimise the arm
boundaries or represent a CT-rich motif on the reverse-complement strand; those differences are
inherent to the scoring result rather than a casing error.

These changes have been added to the HSeeker GitHub repository and web tool.

---

## Comment 11 — two `app.py` files

> It seems there is some Streamlit demo, while the app stands on FastAPI, giving rise to two app.py
> files. This might be misleading.
>
> Apart from this, I really acknowledge the IT part of the project — the care for the speed of the code
> as well as putting the package on PyPi and as a standalone on-line tool. This might really help the
> user!

We agree that two front-ends were confusing, and we thank the reviewer for the kind words about the
software engineering. Rather than relabel the demonstration interface, we removed it. HSeeker now
ships a single web application, the FastAPI service, which is also what the hosted deployment and the
container image serve; the Streamlit demo and its dependency have been removed from the package. The
command-line tool, the Python package and the web application now share the same default settings and
return identical results for identical input, which also removes a source of confusion the reviewer
did not see: in the submitted version the library and the command-line tool applied different
composition filters by default.

**The manuscript now reads (Results, HSeeker web interface, first sentence):**

> HSeeker provides a single web interface, which uses the same scanning library and the same default
> settings as the command-line tool, so that the two return identical results for identical input.

---

# Reviewer 3

> The manuscript introduces HSeeker, a computational framework for identifying putative antiparallel
> H-DNA-forming sequences and assigning a stability score. … Nevertheless, I have several major
> concerns regarding validation and interpretation that should be addressed before the work can
> support publication in *Briefings in Bioinformatics*.

We thank the reviewer for the helpful comments and overall positive assessment. We provide a point by
point response below.

---

## Comment 1 — benchmark size, imbalance and confidence intervals

> The benchmark is too small, imbalanced, and internally heterogeneous to establish predictive
> performance convincingly. … The reported raw specificity of 66.7% corresponds to only six correctly
> classified negative sequences out of nine … The authors should either substantially expand the
> negative-control dataset or provide a compelling justification for why this is not possible. At
> minimum, confidence intervals should be reported for sensitivity, specificity, MCC, PPV and NPV, and
> conclusions concerning specificity should be appropriately qualified.

We agree, and we have both expanded the negative set and been more careful about the claims
associated. The benchmark is now balanced at 65 forming and 65 non-forming sequences, and we report
sensitivity, specificity, PPV, NPV and MCC with confidence intervals for three overlapping sets rather
than one: the full balanced set, the experimental records alone, and the subset of those that are
construct-level rather than inferred from a related construct.

On the full set, HSeeker reaches a sensitivity around 0.97 and a specificity of 1.0, with an MCC of
0.97. Restricting to the experimental records alone, sensitivity stays the same but the specificity
interval widens considerably, and the same pattern holds for the construct-level subset. This is
because only five experimental sequences in the whole collection are confirmed non-forming after
curation, so any specificity estimate on the experimental data by itself is really being carried by
five records, and the tighter interval on the full benchmark reflects the designed negatives we added
rather than a larger pool of confirmed experimental ones.

**We provide the updated table below (Table 1 of the revised manuscript):**

| Evaluation set | n | TP | FN | TN | FP | Sensitivity | Specificity | PPV | NPV | MCC |
|---|---|---|---|---|---|---|---|---|---|---|
| Balanced benchmark | 130 | 63 | 2 | 65 | 0 | 0.969 (0.895–0.992) | 1.000 (0.944–1.000) | 1.000 (0.943–1.000) | 0.970 (0.898–0.992) | 0.970 (0.925–1.000) |
| Experimental records | 70 | 63 | 2 | 5 | 0 | 0.969 (0.895–0.992) | 1.000 (0.566–1.000) | 1.000 (0.943–1.000) | 0.714 (0.359–0.918) | 0.832 (0.560–1.000) |
| Construct-level records | 57 | 50 | 2 | 5 | 0 | 0.962 (0.870–0.989) | 1.000 (0.566–1.000) | 1.000 (0.929–1.000) | 0.714 (0.359–0.918) | 0.829 (0.556–1.000) |

Intervals are Wilson intervals for proportions and percentile bootstrap intervals for the MCC. The
manuscript qualifies every statement about specificity by noting that it rests mainly on the designed
negatives.

---

## Comment 2 — related sequence variants and family-aware validation

> More importantly, the 79 observations do not appear to represent 79 genuinely independent
> experimental examples. … I strongly recommend a family-aware validation strategy in which closely
> related constructs from the same experimental series are grouped. Ideally, performance should be
> evaluated using leave-one-family-out or leave-one-study-out validation. At minimum, results should be
> recalculated after collapsing highly related sequence families.

We agree that related benchmark records are not independent and have addressed this through family-
and study-aware validation. Using the curated provenance, we grouped the 130 records into 52 sequence
families and 27 studies, ensuring that related sequences were evaluated together. The largest family,
the pXY32 substitution series and its derivatives, accounts for 37 of the 130 records.

With the default parameters fixed, HSeeker achieved an MCC of 0.970 under both leave-one-family-out
and leave-one-study-out validation, compared with 0.924 when parameters were tuned within each fold.
To further assess the influence of related sequences, we reduced the benchmark to one representative
per family, obtaining an MCC of 0.862, or one representative per family and class, obtaining 0.945.

Because HSeeker's default parameters are fixed, holding out entire families or studies does not change
individual predictions. It does, however, affect the uncertainty of the estimates. The family- and
study-level confidence intervals are wider than those calculated over individual records, reflecting
their dependence. We now report these intervals alongside the original results and have added a
dedicated section describing the validation.

**The manuscript now reads (Results, Family- and study-aware validation):**

> Because many benchmark records are closely related variants from the same experimental series, we
> grouped records by sequence family and by study and repeated the evaluation holding out one whole
> group at a time. Under leave-one-family-out and leave-one-study-out validation the fixed default
> configuration reached an MCC of 0.970, against 0.924 when the parameters were tuned within each
> training fold. Collapsing the benchmark to one representative per family gave an MCC of 0.862 over
> 52 records, and to one representative per family and class 0.945 over 74 records. Because the
> parameters are fixed, a record's call depends only on its own sequence, and group-aware validation
> does not shift the point estimate; it widens the confidence interval, and both the record-level and
> the family-level intervals are reported throughout.

---

## Comment 3 — circularity in the stability-score threshold

> There is an important risk of circularity/overfitting in the stability-score validation. The authors
> determine an "optimal classification threshold" using Youden's J statistic and subsequently evaluate
> discrimination on the same benchmark. … The threshold must be determined using training data and
> evaluated on independent test data, or through appropriately nested cross-validation/bootstrapping.

We agree, and we are grateful for the suggestion, which proved decisive. We re-evaluated the threshold
exactly as proposed: within each training fold we selected the threshold maximising Youden's J and
applied it only to the held-out records, under repeated stratified five-fold cross-validation,
leave-one-family-out and leave-one-study-out validation.

A threshold chosen this way never improved performance. Out of fold it gave an MCC of 0.949 (95% CI
0.937–0.961) under cross-validation and 0.955 (0.898–1.000) under both group-aware schemes, against
0.970 (0.925–1.000) with no threshold at all. It removed no false positive, because after the purity
requirement no non-forming record in the benchmark retains any candidate, and out of fold it lost
forming records whose score fell below a threshold learned without them. For the same reason the AUC
on this benchmark, 0.985, only restates the detection result, since it equals (1 + sensitivity)/2 when
no negative carries a candidate.

We have therefore stopped using a score threshold. A sequence is called H-DNA-forming when it contains
a candidate that passes the detection criteria, and the stability score is reported as a ranking of
candidates. The Youden-selected threshold, the in-sample AUC and the threshold-based confusion matrix
are removed from the manuscript, and the web interface no longer labels sites above a score of 50 as
stable.

**The manuscript now reads (Results, Detection performance on the revised benchmark, closing
paragraph):**

> We also asked whether the stability score could improve classification beyond detection. Within each
> training fold we selected the score threshold maximising Youden's J statistic and applied it only to
> the held-out records. Out of fold this produced an MCC of 0.949 (0.937–0.961) under repeated
> five-fold cross-validation and 0.955 (0.898–1.000) under leave-one-family-out and
> leave-one-study-out validation, below the 0.970 obtained with no threshold. The threshold removed no
> false positive, because no non-forming record retains a candidate once the purity requirement is
> applied, and it discarded forming records whose scores fell below a threshold learned without them.
> For the same reason the area under the ROC curve on this benchmark, 0.985, restates the detection
> result rather than adding to it. We therefore report no score threshold: a sequence is called
> H-DNA-forming when it contains a candidate satisfying the detection criteria, and the stability score
> is used to rank candidates.

---

## Comment 4 — the thermodynamic interpretation of the stability score

> The thermodynamic interpretation of the "stability score" needs considerably more caution and
> validation. … The authors should therefore clearly distinguish an experimentally informed heuristic
> stability score from a quantitatively calibrated thermodynamic prediction unless they can demonstrate
> a direct relationship with experimentally measured ΔG, melting temperatures, transition superhelical
> densities, or equivalent quantitative measurements. This distinction is especially important because
> the benchmark itself combines assays performed under quite different conditions. … The manuscript
> should discuss this limitation explicitly.

We agree with this reviewer and have now described the output as an experimentally informed,
sequence-based *heuristic stability score*, not a calculated free energy or a quantitative prediction
of H-DNA formation. Inspection of the implementation clarifies the assumptions: matched G and A
positions receive weights of 7.6 and 3.83, respectively; every other arm-position pairing receives a
penalty of −5.715, the negative mean of those two weights. Adjacent matched positions receive a
stacking reward of 5.0, whereas disrupted adjacencies incur a penalty that grows with the length of
the run until two adjacent matches reset it. The algorithm selects the highest-scoring permissible arm
window and clips negative final scores to zero. Thus 5.0 is the fixed stacking term and the first
disrupted-adjacency penalty, not a universal mismatch energy.

Although some weights were motivated by reported triplex interactions, the combined score has no
validated physical units and has not been calibrated against measured H-DNA ΔG, melting temperature,
or transition superhelical density. We have removed kcal/mol labels and language that treated
`total_score` as a thermodynamic free energy or an absolute formation probability, and now present it
only as a relative ranking within the specified scoring scheme. The current benchmark evaluates
recovery and discrimination of reported sequence labels; it does not establish thermodynamic
calibration.

**To reflect the limitations of the stability score, we have revised the main text as follows:**

> The stability score is a heuristic model of the propensity of H-DNA formation within a given mirror
> repeat sequence. It does not model superhelical stress, pH-dependent protonation, ionic composition,
> temperature, loop or spacer free energy, or competition with duplex and other non-canonical DNA
> structures.

On the heterogeneity of assay conditions, we agree and have acted on it concretely. The benchmark
records include acidic assays, Mg²⁺-containing assays, and negatively supercoiled plasmids, so their
forming/non-forming labels are outcomes under the reported experimental conditions. They cannot be
interpreted as equivalent propensities under one common condition, and especially not as evidence of
formation in physiological chromatin. We have therefore recorded, for every experimental record, the
assay used, the DNA context, and the pH, temperature, ionic conditions and supercoiling state under
which the structure was assessed, together with the supporting passage from the source and an explicit
marker where the source does not report a value. Sixty-four of the 70 experimental records have a
reported pH. These conditions are supplied with the benchmark and summarised in Supplementary Table 1;
the pooled performance figures are now described as across-condition sequence-label discrimination;
and the Limitations section identifies matched, condition-controlled measurements of ΔG, melting
behaviour, or B-to-H transition thresholds as the validation required before any quantitative
thermodynamic claim could be made.

**The manuscript now reads (Methods, Benchmark dataset, closing paragraph):**

> For each experimental record we recorded the assay used, the DNA context, and the pH, temperature,
> ionic conditions and supercoiling state under which the structure was assessed, together with the
> passage in the source supporting those values; values that the source does not report are marked as
> such. Sixty-four of the 70 experimental records have a reported pH. Because the labels are outcomes
> under heterogeneous conditions, the pooled performance figures should be read as discrimination of
> reported sequence labels across assay conditions rather than as propensities under a single common
> condition.

---

## Comment 5 — benchmark curation

> The benchmark curation requires substantially greater scrutiny. … some sequences are explicitly
> described as having been added from a supplied "H-DNA Thermostability document," with no
> primary-paper citation attached to the specific construct and with the classification described as
> "low confidence pending source-level verification." … Low-confidence or secondary-source annotations
> should either be verified against primary literature or excluded from the primary benchmark.
> Relatedly, the supplementary table contains an acknowledged probable transcription inconsistency
> involving pRW1406/pRW1411.

We thank the reviewer for the comment. We agree, and we have rebuilt the benchmark around explicit
inclusion criteria. A record enters the primary benchmark only if (i) it cites a resolvable primary
experimental paper, (ii) its exact sequence, or an unambiguous construct definition, can be traced to a
figure, table or methods section of that paper, and (iii) its forming/non-forming label is stated in
that paper for that construct. Records with experimental support for the motif family but not the exact
sequence form a separate, clearly labelled secondary set. Records that meet neither standard were
removed. Detection results played no part in any of these decisions. The per-record audit, with the
evidence location for every record, is deposited with the benchmark
(`hdna_benchmark_experimental_v3.csv`). Of 81 audited records, 70 were kept, nine removed and two
excluded on validity grounds.

**Low-confidence records.** All six records carrying low-confidence annotations from the secondary
"H-DNA Thermostability" document were removed. Three (GG32-long, Mycm1, MycM2) could not be traced to
any primary source. Three (the 30-nt GG32 and AG32 versions and MYCAA) could be traced, but were
removed together with all other low-confidence-tier records. A further record whose construct could not
be traced (HDNA0039, a poly-dC representative) was also removed. No low-confidence-tier record remains
in the benchmark.

**Other records from the secondary document.** For the remaining 25 records originally transcribed
from that document, we checked each sequence and label against the primary paper. Twelve were
confirmed and remain in the primary set. Thirteen were moved to the secondary set: twelve because
the published experiment tested a closely related construct of the same repeat family rather than
the exact sequence recorded, and one for the reason given below.

**pRW1406/pRW1411.** The primary paper (Hanvey, Klysik & Wells 1988, *J Biol Chem* 263:7386) resolves
the inconsistency. The abstract enumerates seven inserts by G+C content: the (TCC)₈ insert at 66% G+C
is pRW1411, and the (TTC)₈ insert at 33% G+C is a distinct construct, both among the inserts the paper
lists as forming an intramolecular triplex. The duplicated sequence in the secondary compilation was a
transcription error. We kept pRW1411 and added (TTC)₈ as a separate primary record, named after its
verified insert rather than a plasmid number, because the number-to-insert mapping is precisely what
the secondary document got wrong. One matched negative was added to preserve the balance of the
benchmark.

**Invalid records.** Two records were excluded. One was an inosine-substituted oligonucleotide that
cannot be represented as A/C/G/T DNA. The other was an intermolecular triplex-forming oligonucleotide
mislabelled as intramolecular H-DNA. The two pure homopolymers were moved to a separate panel, because
the composition filters decide their calls by construction.

**Inferred labels.** One record whose label rested on a review's sequence-based inference, with no
experiment (XDP (CCCTCT)₁₀), was moved to the secondary set.

**Citations.** Corrected identifiers were applied to 25 records where the original citation was wrong
or unrelated, 21 of which are retained in the benchmark. These include the four Friedreich ataxia
(GAA)ₙ records, for which we now cite the primary experimental studies (Potaman et al. 2004; Sakamoto
et al. 1999); because those studies characterised closely related constructs rather than these exact
sequences, the records sit in the secondary set. The original citations are retained in the file so
that every correction is auditable.

The primary benchmark now contains 57 construct-level records (52 forming, 5 non-forming). We report
performance on this set separately from the full balanced benchmark: MCC 0.829
(95% CI 0.556–1.000), sensitivity 0.962 (0.870–0.989), specificity 1.000 (0.566–1.000). Only five
experimental non-forming constructs met the primary criteria, so specificity rests mainly on the
designed negatives. We state this explicitly in the text.

---

## Comment 6 — the comparator set and the strength of the claim

> The comparator analysis is useful but not fully equivalent. … The manuscript also discusses NeSSie
> and nBMST as relevant tools but does not include them in the main experimental benchmark. A more
> comprehensive comparison with algorithms designed specifically to identify mirror repeats/H-DNA would
> strengthen the manuscript considerably. … The statement that HSeeker "outperformed existing tools in
> benchmarking" in the Key Points is also stronger than the Results themselves … These claims should be
> harmonized.

We agree and have expanded the comparison to include both mirror-repeat tools suggested by the
reviewer, bringing the total to four comparators. All tools were evaluated on the revised balanced
benchmark using documented default settings, matched detection parameters where possible, and parameter
tuning restricted to training folds. The evaluation criteria and statistical analysis were established
before running the comparisons.

HSeeker achieved the highest MCC (0.970), followed by triplex (0.819), NeSSie (0.583), nBMST (0.449),
and Triplexator (0.024). These differences partly reflect the tools' distinct purposes and detection
criteria, which we now discuss explicitly. Triplexator targets intermolecular rather than intramolecular
triplexes, while NeSSie searches for shorter mirror-symmetric motifs. nBMST's classification is
sensitive to strand orientation, missing forming sequences represented on the pyrimidine-rich strand.

The closest comparator, triplex, detected nearly the same forming sequences as HSeeker but also
classified nine designed negatives as forming, primarily mirror-disrupted mutants. Importantly, the two
tools performed comparably on experimental records alone.

We have revised the Results and Key Points to reflect these distinctions, replacing the earlier general
claim of superiority with a more precise conclusion.

**The Key Point now reads:**

> On a balanced benchmark of 130 experimentally characterised and designed sequences, HSeeker achieved
> higher accuracy than four existing tools, matching triplex on experimentally characterised sequences
> while rejecting the mirror-disrupted controls that triplex accepted.

---

## Comment 7 — the E. coli insertion experiment and the overlap threshold

> The E. coli insertion experiment provides limited independent biological validation. … it does not
> constitute an independent validation of H-DNA formation. The authors appropriately call this a
> secondary recovery analysis, which I agree with. However, Figure 3 places considerable emphasis on
> this experiment. The manuscript should make clearer that this analysis tests computational recovery in
> genomic background, not biological prediction accuracy. Moreover, the >80% overlap criterion is
> somewhat arbitrary and potentially tool-dependent. Sensitivity analyses using several overlap
> thresholds would demonstrate whether the comparative conclusions are robust.

We agree that the experiment tested computational recovery in a genomic background rather than H-DNA
formation, and that the original presentation gave it more weight than it could carry. We have
rewritten it accordingly rather than leaving the framing to be inferred. It is now presented as a test
of whether embedding a sequence in a genome changes the call made on it; the biological claims attached
to it are gone; and the biological validation the reviewer asks for is supplied by the comparison
with the permanganate/S1 single-stranded DNA maps of Kouzine and colleagues, which we describe in
full under the reviewer's comment 8 below.

Measured that way the experiment earns its place: HSeeker's calls are identical for all 130 records
whether the sequences are scanned on their own or embedded, whereas triplex loses 44 detections once
the same sequences sit in a genomic background. That is a practical result for anyone applying either
tool to genomic sequence, and it is not visible from the benchmark alone.

We also ran the sensitivity analysis the reviewer asks for, since the criterion does still enter a
reported result. Varying the positional overlap required to match a hit to an inserted sequence from
0.50 to 0.95 leaves every conclusion unchanged: HSeeker's MCC is 0.970 at 0.80 and above and 0.954 at
0.70 and below, the single difference being one locus that partially overlaps a neighbouring insertion
at the loosest settings, while triplex is unchanged at 0.512 throughout. No ranking depends on the
choice. The sweep is reported in Figure 3C.

**The manuscript now reads (Results, Genomic context recovery):**

> Sequences are rarely scanned in isolation, so we asked whether embedding a sequence in a genome
> changes the call made on it. Inserting the 130 benchmark sequences into the *E. coli* K-12 genome
> left HSeeker's calls unchanged for every record: the same 63 forming sequences were detected, no
> non-forming sequence was called, and the resulting coefficient of 0.970 is identical to the value
> obtained on the sequences alone (Figure 3A, B). The matching hit covered the whole inserted sequence
> for the majority of records and at least half of it in every case (Figure 3D). triplex, run on the
> same modified genome, returned a different call for 44 of the 130 records, in each case a detection
> present when the sequence was scanned alone and absent once it was embedded, and its coefficient fell
> from 0.819 to 0.512. Varying the overlap required to match a hit to an insertion from 0.50 to 0.95
> changed neither result (Figure 3C). The unmodified genome yielded 83 loci across its 4.64 Mb, or 17.9
> per Mb, and the hits in the modified genome are accounted for by those background loci together with
> the recovered insertions.
>
> This experiment measures the stability of a call across sequence context, not whether H-DNA forms. We
> draw no biological conclusion from it; the comparison with measured single-stranded DNA follows.

---

## Comment 8 — genome-scale biological validation

> A genome-scale biological validation would substantially strengthen the manuscript. HSeeker identifies
> 96,729 non-overlapping H-DNA loci on chromosome 1 under default parameters. Yet there is little
> analysis demonstrating that these predictions correspond to independent experimental signatures of
> non-B DNA formation … the authors should ideally compare HSeeker predictions with an independent
> genome-wide experimental dataset, for example structure-sensitive nuclease mapping … while carefully
> accounting for sequence composition.

We thank the reviewer for this suggestion; it produced the most informative analysis in the revised
manuscript. We compared HSeeker's predictions with the structure-sensitive nuclease mapping of Kouzine
and colleagues: permanganate/S1 footprinting of single-stranded DNA in activated and resting mouse B
cells and in human Raji cells. We used their published calls and signal tracks, on the genome
assemblies they used, and scanned both genomes with HSeeker at its default settings.

**Recall and selectivity.** HSeeker loci overlap 79.0% of the H-DNA structures the authors found
formed in mouse (13,517 of 17,109) and 82.3% in human Raji cells (6,879 of 8,356). Among their own
predicted H-DNA motifs, those HSeeker calls are formed more often than those it does not: 4.5%
against 2.8% in mouse (odds ratio 1.65, 95% CI 1.59-1.71) and 5.7% against 1.7% in human (3.39,
3.20-3.59); restricting to motifs that the 36-nt reads can map to raises these to 1.68-2.07 and
4.16-5.70. We note that a strict subset of experimentally formed structures is not attainable by any
sequence-based predictor on these data, because formation depends on cell state: only about 3-5% of
the authors' own predicted H-DNA motifs carried single-stranded signal in a given cell type.

**Enrichment of the measured signal.** We then asked whether the signal itself is elevated at HSeeker
loci, which forms the new Figure 4. We report two analyses. The first was specified, with its
controls and statistics, before any result was computed; the second was added after the first and was
likewise specified before it was computed. We report both.

**Composition matching.** Single-stranded signal depends strongly on base composition, so control
windows were matched to each HSeeker locus on length, A/T count and purine/pyrimidine asymmetry, and
were required not to overlap any HSeeker locus or predicted H-DNA motif. Matching on asymmetry is
essential: controls matched on A/T content alone carry 1.7 to 2.5 times more signal than HSeeker loci,
because mixed-sequence DNA gives more signal than purine/pyrimidine-asymmetric DNA of the same
composition.

**First analysis: all HSeeker loci, genome-wide controls.** Loci in activated mouse B cells carried 14%
more single-stranded signal than their controls (ratio of means 1.137, 95% CI 1.127–1.147; 825,749
loci). In resting B cells the ratio was 0.822 (0.818–0.827) and in human Raji cells 0.977 (0.969–0.984).
This analysis does not separate H-DNA from other sources of single-stranded DNA and does not account
for the poor mappability of repetitive purine tracts.

**Second analysis: loci isolated from other sources of single-stranded DNA.** We removed from both loci
and controls every interval that the reads cannot map to, every interval within 100 bp of a motif
predicted to form a G-quadruplex, Z-DNA, a cruciform or a stress-induced melting site, and every
interval within 2 kb of an annotated transcription start site, and drew ten controls per locus within
50 kb. In activated mouse B cells the isolated loci carried 25% more signal than their controls (1.254,
1.179–1.332; a locus exceeded one of its own controls with probability 0.514, p = 2 × 10⁻¹⁷). In human
Raji cells there was no detectable excess (1.023, 0.988–1.061; p = 0.08). Two built-in checks behaved as
required in both species: comparing one control against the remaining controls of the same locus gave no
difference (0.500–0.502), and the motifs the authors called formed were enriched 7.4-fold in mouse and
2.1-fold in human under the identical procedure.

**Dependence on arm length.** Because G-quadruplexes, Z-DNA, duplex melting and transcription do not
depend on mirror symmetry once composition is matched, whereas H-DNA stability increases with arm
length, we also stratified by arm length. In activated mouse B cells the excess rose from no detectable
effect at 8–9 nt (1.12, 0.93–1.35) to 1.58 (1.34–1.83) at 12–14 nt and 1.6–1.8 at longer arms (Spearman
p = 3 × 10⁻²¹). The enrichment was absent in resting B cells (0.858, 0.831–0.885), in which even the
authors' formed motifs showed no excess, consistent with the activation- and supercoiling-dependent
formation they report. No arm-length class differed from unity in Raji cells.

We report the human result as it stands. The Raji data have no matched input control and a smaller
dynamic range — the authors' own formed motifs are enriched 2.1-fold there against 7.4-fold in mouse —
but we do not treat that as an explanation, and the manuscript states that we did not observe the
enrichment in human cells.

Finally, the chromosome 1 count quoted in the comment, 96,729 loci, was obtained with the previous
default minimum arm length of 10 bp and without the composition filters. With the revised defaults
HSeeker reports 45,927 loci on chromosome 1, or 106,481 with the composition filters disabled, and the
manuscript now states the settings alongside every count.

---

## Comment 9 — ablation and parameter robustness

> Parameter robustness and ablation analyses are needed. … The authors should provide an ablation
> analysis comparing, for example, sequence detection alone; detection plus composition filtering;
> pairing-only scoring; pairing plus stacking; and the complete model. Similarly, performance across
> reasonable purity, mismatch, minimum-arm and spacer thresholds should be reported. This would
> demonstrate that the results are not dependent on one hand-selected parameterization.

We agree, and we have added both analyses on the revised balanced benchmark; a notebook and
machine-readable tables are provided. A record is called forming when at least one candidate survives
the stage in question, and no score threshold is used at any stage.

**Cumulative ablation.** Mirror-repeat detection alone detected 63 of 65 forming records and seven of 65
non-forming records (sensitivity 0.969, specificity 0.892, MCC 0.864). The seven were the experimentally
non-forming *c-MYC* control B-33, one G-quadruplex sequence, two Z-DNA (CG)ₙ repeats and the three
perfect mirror repeats of balanced composition. Adding the purity requirement removed all seven without
losing a forming record (specificity 1.000, MCC 0.970). Adding the AT-content filter, pairing-only
scoring, pairing plus stacking, and the complete model with boundary optimisation and homopolymer
filtering changed no record-level call. Pairing-only and pairing-plus-stacking scores rank the same
candidates, and once purity is applied no non-forming record retains a candidate, so the score AUROC of
0.985 restates the detection result rather than adding to it. On the 70 experimental records alone the
MCC is 0.708 for detection alone and 0.832 from the purity step onward. Retaining all candidates,
longest-arm selection and score-based selection gave identical record-level calls.

**One parameter at a time.** Purity 0.80–0.90 gave an MCC of 0.970, while 0.95 and 1.00 lost nine forming
records (0.843). Mismatch tolerances of 0 and 0.05 detected 49 of 65 forming records (0.778), 0.10
detected 63 (0.970), and 0.15 and 0.20 detected 64 but added a false positive (0.969). A minimum arm
length of 6 added ten false positives (0.839), 7 and 8 gave 0.985 and 0.970, 9 and 10 lost eight forming
records (0.884) and 12 lost sixteen (0.778). A maximum spacer of 3 detected 41 forming records (0.679),
5 detected 55 (0.856), and 8 to 15 detected 63 (0.970). The single record gained at a minimum arm length
of 7, and at a mismatch tolerance of 0.15, is the tandem direct repeat (ACCCTCCCC)₄, which lies outside
the structural model.

**Out-of-sample selection.** Across all 576 combinations with selection confined to training folds, the
fixed default configuration performed as well as or better than per-fold tuning under nested
cross-validation, leave-one-family-out and leave-one-study-out validation (0.970 against 0.939, 0.924
and 0.924), and tuning did not beat the fixed configuration in any of 200 random splits. A minimum arm
length of 8 was selected in 99% of nested folds and in every family and study fold, and a mismatch
tolerance of 0.10 in 48 of 52 family folds. On this basis we changed the default minimum arm length from
10 to 8 bp and the default mismatch tolerance to 0.10.

Performance therefore does depend on these choices, as the reviewer anticipated, but the reported
results do not rest on a single hand-selected parameterisation. The notebook and machine-readable tables
are in `analysis/hseeker_ablation_comment9.ipynb` and `analysis/results/ablation_comment9/`. The scorer
does not currently expose an independent mismatch-penalty switch, so we do not claim to have ablated that
penalty separately.

---

## Comment 10 — the overlap-removal strategy

> The overlap-removal strategy may discard biologically relevant alternatives. The greedy procedure
> retains the candidate with the longest arm within overlapping groups, using shorter spacer length as a
> tie-breaker. It is not clear that longest-arm length is always the biologically optimal criterion,
> particularly when a shorter overlapping motif has a substantially higher predicted stability. The
> authors should justify this choice and compare it with alternatives such as selecting the
> highest-scoring candidate. Since the web interface apparently permits overlapping motifs to be
> retained, the manuscript should also explain which mode is recommended for biological analysis.

We thank the reviewer for the comment. We agree that arm length is not a validated proxy for biological
activity, and we now implement the alternative proposed and compare the two directly.

**Why the default is retained.** We keep the longest-arm rule as the default for three reasons, which
the manuscript now states. Among overlapping candidates covering the same mirror-repeat region, the
longest arm describes the full extent of the repeat, which is the region available for strand exchange.
The criterion is purely geometric and does not depend on the stability score, which is a heuristic
ranking rather than a calibrated free energy; fine score differences between overlapping candidates
should not decide what a user sees by default. And it is faster at genome scale,
because overlaps are resolved before candidates are scored.

**The three modes.** HSeeker now offers three overlap modes in the Python API, CLI and web interface:

- **Retain all candidates** (`remove_overlaps=False`; CLI `-skipoverlap`; web "Keep all raw hits")
  returns every candidate surviving the requested filters.
- **Longest-arm selection** (the default) performs a start-ordered sweep in which a candidate with a
  longer arm replaces the overlapping incumbent, with a shorter spacer breaking ties.
- **Score-based selection** (`overlap_strategy="score"`; CLI `-overlap-strategy score`; web "Highest
  scores first") scores every candidate first and then considers them in descending order of score,
  retaining a candidate when it does not overlap one already retained, with ties resolved by longer arm,
  then shorter spacer, then earlier coordinate.

Score-based selection produces a deterministic maximal set of non-overlapping intervals. An overlap
chain can retain several compatible candidates; the rule does not claim to maximise the number of hits
or their combined score. Coordinates are 1-based and inclusive.

**Order of operations.** This matters and is now documented. AT-content filtering precedes scoring.
Longest-arm selection runs before scoring and before the scored homopolymer check, whereas score-based
selection runs after both, so a candidate removed by the scored filter leaves a lower-scoring
overlapping alternative eligible.

**Comparison of the two rules.** The two rules agree far more often than they differ. On the benchmark,
43 of the 70 experimental records have overlapping candidates after filtering; the two rules return
identical intervals for 38 of them and differ for five, and the record-level call is identical for all
130 records and for all 74,880 record-by-configuration calls across our parameter grid. On human
chromosome 1 at the default settings, longest-arm selection reports 45,927 loci and score-based
selection 47,217; every longest-arm locus is overlapped by a score-selected locus, 41,867 of them
(91.2%) have identical coordinates, and score-based selection reports 991 additional loci. With the
composition filters disabled the counts are 106,481 and 108,224, with 82.6% identical. The two rules
therefore locate the same regions and differ mainly in which representation of a region they report, at
a modest cost in time: 23.9 s against 16.3 s on chromosome 1 with 16 threads.

The following worked examples, with scoring enabled and a minimum arm of 7, maximum arm 25, maximum
spacer 5, purity 0.9 and mismatch tolerance 0.2, show when the rules agree and when they diverge.

| Input | Candidate | Interval | Arm / spacer | Total score |
|---|---|---|---|---|
| `GAGAAAAGAGGGAGAAAAAGGAA` | gagaaaagagggagaaaaagg | 1–21 | 10 / 1 | 53.635 |
| `GAGAAAAGAGGGAGAAAAAGGAA` | gaaaagagggagaaaaa | 3–19 | 8 / 1 | 53.635 |
| `GGAAGAAGAAGGAGGAGAGAGGGAAGGAA` | ggaagaagaaggagg | 1–15 | 7 / 1 | 38.575 |
| `GGAAGAAGAAGGAGGAGAGAGGGAAGGAA` | aaggaggagagagggaaggaa | 9–29 | 10 / 1 | 20.360 |

In the first input both rules return 1–21: longest-arm selection because its arm is longer, score-based
selection because the two scores tie and the longer arm wins the tie. In the second, longest-arm
selection returns 9–29 because its arm is longer, while score-based selection returns 1–15 because its
score is higher. Retaining overlaps returns both candidates in each case. An overlap chain need not
collapse to a single representative: scanning (GAA)₂₀ with a minimum arm of 10 produces 102 overlapping
candidates, and both rules retain the two compatible hits at 3–32 and 35–60.

**Recommended mode.** We recommend the default longest-arm rule for genome-wide annotation, where the
full extent of each repeat is the quantity of interest and the run time matters; score-based selection
for locus-level work such as choosing which variant to test experimentally, where the relative ranking
of overlapping alternatives is the question being asked; and retaining all candidates for exhaustive or
downstream statistical analysis. A higher score helps prioritise candidates and is not experimental
evidence of formation. The manuscript and the README now state this recommendation explicitly.

---

## Comment 11 — computational benchmarking and the whole-genome extrapolation

> Computational benchmarking should be made more reproducible and somewhat less promotional. …
> extrapolating chromosome 1 runtime linearly to an approximately five-minute whole-genome runtime is not
> equivalent to actually benchmarking the whole genome. The authors should either run the complete hg38
> genome or explicitly label the estimate as an extrapolation throughout. Raw benchmark commands, output
> files and timing data should also be deposited rather than being available only "upon reasonable
> request."

We agree on both counts. We re-measured runtime and memory on the final code at the default settings, on
the same hardware as the original benchmark, with three runs at each point and peak memory recorded
externally. On human chromosome 1 the complete pipeline took 173.8 s with one worker and 16.3 s with 16,
a 10.65-fold speedup, with a peak memory of 1,334 MB. Rather than extrapolating, we ran the complete hg38
primary assembly: 204.5 s with 16 workers, peak memory 2,296 MB, reporting 559,601 loci, with the three
runs producing identical output. Score-based overlap selection took 23.9 s on chromosome 1.

We also re-timed the comparators on the same node on the same day, so that the comparison is like for
like. At 16 cores chromosome 1 took 45.0 s with Triplexator and 64.7 s with triplex run as 16 parallel
processes; nBMST and NeSSie have no parallel mode and took 96.5 s and 295.7 s respectively on one core.
triplex and Triplexator reproduced the hit counts of our original measurement exactly.

All benchmark commands, raw timing output, hardware descriptions and input checksums are included in the
repository and will be archived with a permanent identifier on publication. The Data Availability
statement no longer refers to materials available on request. Figure 2 has been regenerated from these
measurements.

---

## Comment 12 — claims, terminology and figures

> Several claims in the Discussion should be moderated. Statements such as the scoring system "works
> highly effectively and aligns seamlessly" with experimental findings, that HSeeker is "the only H-DNA
> identification algorithm with a proven scoring system built-in," and that broad evolutionary analyses
> were "previously computationally prohibitive" are insufficiently demonstrated by the presented
> analyses. … There are also several minor issues of terminology, grammar and presentation throughout the
> manuscript. … The manuscript would also benefit from careful language editing, and Figures 1 and the
> algorithm panels are visually dense and could be simplified substantially.

We appreciate the reviewer's detailed and thoughtful comments. We have revised the manuscript
accordingly. In detail:

**12.1 — "works highly effectively and aligns seamlessly."** We agree that we were over-enthusiastic
about the validation result. The score is an experimentally informed heuristic that has not been
calibrated against measured free energies, melting temperatures or transition superhelical
densities, and experimental data under comparable controlled conditions are limited. We have
therefore moderated the statement to:

> The stability score ranks candidates in a manner consistent with the biochemical and genetic evidence
> available for the sequences in our benchmark; it is an experimentally informed heuristic and has not
> been calibrated against measured thermodynamic quantities.

**12.2 — "the only H-DNA identification algorithm with a proven scoring system built-in."** We agree this
is an overstatement. It now reads:

> Among currently available tools for intramolecular H-DNA, HSeeker is the only one that reports a
> per-motif stability ranking alongside the detected coordinates.

**12.3 — "previously computationally prohibitive."** We agree. It has been rephrased as:

> This efficiency makes genome-wide profiling and cross-species comparison practical on commodity
> hardware.

**12.4 — terminology.** We have applied the four terms consistently throughout. Predictions are now
described as putative H-DNA-forming sequences or predicted H-DNA motifs; "H-DNA-forming sequence" without
qualification is reserved for sequences with experimental support; and "mirror repeat" is used for the
sequence feature rather than as a synonym for the structure.

**12.5 — language editing.** The manuscript has been language-edited throughout for grammar, consistency
of tense and terminology, and clarity of presentation.

**12.6 — figure density.** Figure 1 has been simplified: the sensitivity-analysis and AUC inset has been
removed, since that analysis is withdrawn (see comment 3); the candidate-identification flowchart has
been reduced and now shows the inward-extension and deduplication steps explicitly; and the stability
panel is relabelled as a heuristic ranking rather than an energy. The algorithm panels have been
restructured into a short pipeline overview that fixes the order of stages, followed by one compact panel
per stage, which removes the repetition that made the submitted panels dense.

---

## Note for the editor

Two items were settled during the revision and may be of interest:

1. The chromosome 1 locus count in the submitted manuscript (96,729) was obtained at a minimum arm
   length of 10 bp with the composition filters disabled, although the text described the defaults as
   8 bp with filters enabled. The revised manuscript states the settings alongside every count.
2. The submitted version applied different composition filters in the Python library and the
   command-line tool. These are now identical across the library, the command-line tool and the web
   application, which return identical results for identical input.
