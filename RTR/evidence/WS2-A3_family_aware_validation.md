# WS2-A3 / WS2-A2 — family- and study-aware validation, dual confidence intervals,
# and the evidence that HSeeker is not overfitted to the benchmark

**Regenerated 2026-10-07 on 128 records (64:64) after DEC-8 rebalance and DEC-9
score-informed overlap default; 124-record values in this file's 2026-10-02 version.**
Every number below is from the regenerated `analysis/results/grouped_validation_v3/`
and `analysis/results/robustness_v3/` (both 2026-10-07). The four DEC-8 negatives are
all true negatives under every configuration, so the confusion matrices change only in
TN (+4) and the point MCCs move by ~0.001–0.003; the conclusions of 2026-10-02 hold.
DEC-9 (score-based overlap removal by default) changed no detection call on the
benchmark — 0 binary-call differences greedy vs score on all 128 records, and 0 of
73,728 calls across the full WS4-A3 grid — so it has no effect on any number here. In
general, with the homopolymer filter on, score can report a near-homopolymer hit where
greedy reports none (see DEC-9 evidence). DEC-10 (2026-10-08) made greedy the overlap
default again; for the same reason it changes no number here.

**Label note (2026-10-08).** `script_default` / "script default" below is the
**benchmark/script configuration** (the analysis scripts' setting: purity 0.90, mismatch
0.10, minrep 8, maxspacer 10, no filters), **not** the software default. The shipped
software defaults use minrep 10 (library: no filters; CLI: both filters); on this
benchmark they predict exactly like `production_cli` (TP56 FN8 TN64 FP0, MCC 0.882).

**Regenerated 2026-10-02** after merging main (Nikol's scanner) and removing the 6
pure homopolymers; previous values in git history at `007288b`. Every number below
is from the regenerated `analysis/results/grouped_validation_v3/` and
`analysis/results/robustness_v3/`. Changed conclusions are flagged in place: the
composition filters are now no-ops (selection picks `none` in every fold), and the
script default and the manuscript configuration are tied rather than the manuscript
configuration winning outright.

**Status:** done for the detection call; **PROVISIONAL** pending WS0-A1 freeze and
DEC-5. The stability-score threshold half of WS2-A3 remains open, as does
tool-vs-tool CI comparison (no R/Triplex in this environment).
**Date:** 2026-09-23 (first run); regenerated 2026-10-02 and 2026-10-07
**Reviewer items:** R3.2 (primary), R2.7, R3.3, R2.6, R3.1
**Scripts:** `analysis/scripts/run_grouped_validation.py`,
`analysis/scripts/run_parameter_robustness.py`
**Outputs:** `analysis/results/grouped_validation_v3/`, `analysis/results/robustness_v3/`
**Input:** `hdna_benchmark_balanced_v3.csv` (128 records, 64 forming / 64 non-forming),
sha256 `7bb495444e91709f…`

---

## 1. What R3.2 asked, and what was built

R3.2, verbatim: *"I strongly recommend a family-aware validation strategy in which
closely related constructs from the same experimental series are grouped. Ideally,
performance should be evaluated using **leave-one-family-out or leave-one-study-out**
validation. At minimum, results should be recalculated after **collapsing highly
related sequence families**."*

All three were implemented, so the response meets the reviewer's *"ideally"* bar
rather than the minimum:

| scheme | groups | note |
|---|---|---|
| leave-one-family-out | 51 | 25 singletons; largest `fam:pAA32` = 37 records |
| leave-one-study-out | 26 | 4 singletons; largest study (`10.1093/nar/18.22.6621`) = 37 records |
| collapsed per family | 51 representatives | reviewer's literal "at minimum" |
| collapsed per (family, class) | 72 representatives | class-preserving variant |

**The grouping rule is the part that makes it honest.** A synthetic negative
generated from an experimental positive (e.g. a mirror-disrupted mutant of pGG32)
carries that positive's sequence. It therefore **inherits its source's family and
study**, so holding out the pGG32 family cannot leave its own mutant in the training
portion. 44 of the 59 synthetic records inherit a group this way (40 before DEC-8;
all four DEC-8 records have a source and inherit its group); the remaining 15
(G4, Z-DNA, B-DNA, perfect-mirror controls) form their own groups. The synthetic
homopolymer class is now empty (removed by author decision, 2026-10-02).

The reviewer's specific example checks out: `fam:pAA32` is the pXY32 substitution
series — all 16 X,Y combinations — and with its inherited synthetic derivatives it
is a single group of 37 of the 128 records.

## 2. Headline result — two fixed configurations tie at the top under every scheme

Configuration: **`mismatch 0.15, minrep 8, purity 0.90, maxspacer 10, both filters`**
— the parameters already stated in the submitted manuscript (`manuscript_stated`):

| validation scheme | sens | spec | **MCC** | MCC record-level CI | MCC **family-cluster** CI |
|---|---|---|---|---|---|
| **leave-one-family-out** | 0.984 | 0.984 | **0.969** | [0.922, 1.000] | **[0.901, 1.000]** |
| leave-one-study-out | 0.984 | 0.984 | 0.969 | [0.922, 1.000] | [0.892, 1.000] (study-cluster) |
| sequence-level 5-fold CV (nested, WS4-A3) | 0.984 | 0.984 | 0.969 | [0.922, 1.000] | — |
| collapsed per (family, class), n=72 | 0.977 | 0.964 | 0.942 | [0.851, 1.000] | [0.854, 1.000] |
| collapsed per family, n=51 | 0.977 | 0.857 | 0.834 | [0.547, 1.000] | [0.519, 1.000] |

The script's current configuration (`mismatch 0.10`, otherwise the same, no filters;
`script_default` — the benchmark/script configuration, not the software default) scores
the same or marginally higher in every scheme:

| validation scheme | sens | spec | **MCC** | MCC record-level CI | MCC cluster CI |
|---|---|---|---|---|---|
| leave-one-family-out | 0.969 | 1.000 | **0.969** | [0.924, 1.000] | [0.909, 1.000] |
| leave-one-study-out | 0.969 | 1.000 | 0.969 | [0.924, 1.000] | [0.899, 1.000] |
| sequence-level 5-fold CV (nested, WS4-A3) | 0.969 | 1.000 | 0.969 | [0.923, 1.000] | — |
| collapsed per (family, class), n=72 | 0.955 | 1.000 | 0.944 | [0.865, 1.000] | [0.868, 1.000] |
| collapsed per family, n=51 | 0.955 | 1.000 | 0.862 | [0.630, 1.000] | [0.630, 1.000] |

Comparators under leave-one-family-out: `production_cli` (= the shipped software
defaults, minrep 10) MCC 0.882, per-fold tuned 0.923 — both below the two fixed
configurations (0.969 each). (124-record values:
`manuscript_stated` 0.9677, `script_default` 0.9682, `production_cli` 0.8787, tuned
0.9205.)

**What changed.** Previously (130 records) `manuscript_stated` won outright (LOFO 0.939
vs `script_default` 0.908), because the script default lacked the filters and so
called the homopolymer negatives positive. With those records removed the filters
change nothing, and the two configurations differ only in mismatch: 0.15 recovers
HDNA0058 but admits HDNA0028 `B-33` as a false positive (63/1/63/1 vs 62/2/64/0).
MCC is 0.9688 vs 0.9692 — a tie, not a ranking. **The data no longer single out the
manuscript configuration**; they say either fixed configuration is better than
tuning.

**Caveat on the strictest collapsed variant.** Collapsing to one record per family
leaves 44 forming / 7 non-forming, because synthetic negatives inherit their source
positive's family and are collapsed away with it. That reinstates a ~6:1 imbalance
the balanced set exists to remove and its interval is nearly uninformative
(`manuscript_stated` [0.547, 1.000]; specificity rests on 7 negatives). The
class-preserving variant (one representative per family *per class*, n=72, 44:28)
removes within-family redundancy without discarding the design and is the figure we
recommend reporting, with the stricter one disclosed beside it.

## 3. Why family-aware validation does not move the point estimate — and what it does move

LOFO, LOSO and the pooled analysis return **numerically identical** results for a fixed
configuration. This was checked, not assumed (manuscript configuration, re-checked
2026-10-07):

```
full-data evaluation : TP63 FN1 TN63 FP1  MCC 0.9688
LOFO pooled out-of-fold: TP63 FN1 TN63 FP1  MCC 0.9688
```

**Reason:** with parameters fixed, HSeeker's call for a record depends only on that
record's sequence. Nothing is fitted, so there is no training set through which a
related sequence could leak. Grouping *cannot* bias the point estimate of a fixed
configuration.

The per-fold tuned row is the one place grouping could matter, and here LOFO and LOSO
also give the same confusion matrix (60/4/63/1, MCC 0.923), because selection lands on
`purity 0.85, minrep 8, filters none` in every fold of both schemes.

What family structure genuinely affects is **precision**, which is exactly what R3.2
says — *"can substantially inflate apparent statistical power"*:

| uncertainty on MCC 0.969 (LOFO, `manuscript_stated`) | interval | width |
|---|---|---|
| record-level bootstrap | [0.922, 1.000] | 0.078 |
| **family-cluster bootstrap** | [0.901, 1.000] | **0.099** |

The cluster interval is **1.26×** wider. Across all 14 reported scheme/config rows
the ratio ranges 0.91–1.75 (median 1.23). (124-record: 1.23× headline, range
0.92–1.71, median 1.20.) Both intervals are now truncated at 1.000
for the two top configurations, so the widening shows up at the lower bound only.

**Both intervals are reported everywhere, never one instead of the other**, because
they answer different questions: the record-level interval is the R2.7 interval on
the observed sample; the cluster interval is the R3.2 interval asking what happens if
the *families* had differed. Treating 16 pXY32 variants as 16 independent
observations overstates precision by roughly a quarter (median ratio 1.23;
1.26 for the headline row) — and we now say so with a number instead of a hedge.

## 4. The evidence that HSeeker is not overfitted to this benchmark

R2.6 asserts the tool "is overfitted to these exact sequences". Five independent
lines of evidence say otherwise:

1. **Tuning on the benchmark does not help.** Under LOFO the two fixed configurations
   (`manuscript_stated`, `script_default`) both score **0.969** while per-fold tuned
   selection scores **0.923**; under LOSO, 0.969 vs 0.923; under sequence-level nested
   CV, 0.969 vs 0.938. A tool fitted to this data would gain from re-fitting it. This
   one loses.
2. **Repeated-split selection never gains.** Over 200 random splits the median
   tuned-minus-`script_default` eval MCC is **+0.0000** (95% [−0.0955, +0.0000]):
   tuning beat the fixed script configuration in **0** splits, tied in 112 and lost in
   88. (124-record run: 0/109/91;
   original 130-record run 99/43/58 — a coin flip; now it is strictly no better.)
3. **The parameter surface is a plateau, not a peak.** Of 576 configurations, **104
   lie within 0.05 MCC of the maximum and 32 tie at it exactly** (the script default
   among them; all 32 use minrep 8). Performance is not balanced on a hand-picked point.
4. **Selection is unanimous where it matters.** Across 51 family folds, `minrep 8`
   and `purity 0.85` were chosen in **51/51**; mismatch 0.10 in 47/51 (0.15 in 4);
   maxspacer 10 in 49/51 (5 in 2). Across 26 study folds: 26/26, 26/26, mismatch 0.10
   23/26, maxspacer 10 25/26. No knife-edge. Filters `none` in 51/51 and 26/26 — but see
   below: that is a tie-break, not a preference.
5. **Nothing is fitted at all.** Section 3: held-out and full-data evaluation are
   identical for a fixed configuration, because the algorithm has no trained
   parameters. There is no mechanism by which benchmark-specific overfitting could
   occur for the reported configuration.

Point 5 is the strongest and should lead the R2.6 response. The honest residual — and
it must be stated — is that the **AT-content and homopolymer thresholds were chosen
with knowledge of the H-DNA literature and of this benchmark's composition**. That is
the one place R2.6's concern has genuine purchase, and **this benchmark can no longer
speak to it either way**. Previously the filter axis was swept and the three filter
states beat "none" because every false positive was a pure homopolymer. Those six
homopolymer records have now been removed by author decision, so all four filter
states give identical predictions on every record under every grid setting, and the
selector's choice of `none` is just the "fewer filters" tie-break. The filters are
neither supported nor refuted here; the mitigation left is disclosure — that the
thresholds were set from the literature, and that the benchmark excludes pure
homopolymers — rather than a data-derived claim.

## 5. Recommended reporting — best balance of reviewer parity and results

> **Primary:** a fixed, pre-specified configuration evaluated **leave-one-family-out**
> on the 128-record balanced benchmark (64 forming / 64 non-forming). With the
> manuscript's stated parameters (`mismatch 0.15, minrep 8, purity 0.90,
> maxspacer 10, both filters`): **sensitivity 0.984, specificity 0.984, MCC 0.969**,
> record-level 95% CI [0.922, 1.000], family-cluster 95% CI [0.901, 1.000].
> With the script's current parameters (`mismatch 0.10`, no filters; the
> benchmark/script configuration): **sensitivity 0.969, specificity 1.000, MCC 0.969**,
> record-level [0.924, 1.000], family-cluster [0.909, 1.000]. Supported by
> leave-one-study-out (identical point estimates) and by the class-collapsed set (MCC
> 0.942 / 0.944). With the shipped software defaults (minrep 10): **sensitivity 0.875,
> specificity 1.000, MCC 0.882**, record-level [0.803, 0.954], family-cluster [0.750,
> 0.961].

Which of the two to lead with is **DEC-5**, and the data do not decide it — they tie.
The case for the manuscript's values is consistency with the submitted paper; either
way the key sentence for the letter holds: **no selection was performed**, so there
is nothing to discount, and every tuned alternative scored lower (LOFO 0.923). If the
manuscript configuration is used, Methods should note that its "both filters" setting
has no effect on this benchmark. Whichever minrep-8 configuration is reported, it must
be labelled as the benchmark configuration, not as the software default; unless the
shipped minrep changes from 10 to 8, the software-default result (MCC 0.882) is reported
beside it (`TODO(authors)` in the plan, under WS2).

Reviewer parity achieved by this presentation:

| item | status |
|---|---|
| **R3.2** family-aware validation | **met at the "ideally" level** — LOFO *and* LOSO, plus both collapsed variants, with leakage-safe grouping of synthetic derivatives |
| **R3.2** "inflate apparent statistical power" | **met** — family-cluster CIs quantify it (1.26× wider on the headline row, median 1.23× across rows) |
| **R2.7** confidence intervals | **met** — Wilson + record bootstrap + cluster bootstrap on every reported figure |
| **R3.3** circularity | **met for the detection call** — the reported configuration involves no selection whatsoever. *Open:* the stability-score threshold is not yet treated this way |
| **R3.1** imbalance / specificity precision | **met** — 64 negatives (64:64, exactly 1:1), specificity 0.984–1.000 with a stated interval |
| **R2.6** filters defined on the same data | **open, and untestable on this benchmark** — the filter axis was swept but has no effect now that the homopolymer negatives are removed; disclosure is the remaining remedy (see section 4) |

## 6. Correctness checks

| check | result |
|---|---|
| cluster CI contains the point estimate (MCC, sens, spec), all 14 rows | PASS (re-checked 2026-10-07) |
| **every record its own group → cluster CI ≈ record CI** (`manuscript_stated`): record [0.922, 1.000], singleton-cluster [0.922, 1.000] — equal up to bootstrap noise (different RNG streams) | PASS — validates the estimator (re-run 2026-10-07) |
| cluster CI wider than record CI under real grouping (1.26× on the headline row) | PASS |
| cluster bootstrap deterministic under a fixed seed | PASS (re-run 2026-10-07: two calls, identical intervals; full script re-run into a scratch directory reproduced both deposited CSVs byte-for-byte) |
| filter states give identical per-record predictions for all 144 non-filter grid settings | PASS (re-run 2026-10-07, all 128 records) |
| held-out training portion never single-class (all 51 + 26 folds) | PASS (asserted in code; the regenerated run completed) |
| deposited CSV carries both interval families for all 14 rows | PASS |
| synthetic derivatives inherit source group (44 of 59) | PASS |
| full-data vs LOFO pooled confusion matrix identical for `manuscript_stated` (63/1/63/1) | PASS (re-run 2026-10-07) |

## 7. Still open

- Benchmark not frozen (**WS0-A1**, blocked on DEC-2/DEC-3); if adopted, DEC-3 adds a
  positive (→129 records, 65:64). The homopolymer removal (2026-10-02) and the DEC-8
  rebalance (2026-10-07) have already moved every number twice.
- **DEC-5** unresolved. After regeneration sections 2 and 4 no longer argue *for* the
  manuscript's values on performance — mismatch 0.10 and 0.15 tie — only that a fixed
  configuration beats tuning.
- **Reporting configuration / shipped minrep** unresolved (added 2026-10-08): both
  top configurations use minrep 8 (selected in 51/51 family folds and 26/26 study folds),
  the software ships minrep 10 (MCC 0.882 here).
- **AT/homopolymer filters** cannot be evaluated on this benchmark (no target
  negatives remain); any evidence for them has to come from elsewhere.
- **Stability-score threshold** not yet evaluated out-of-fold — this is the half of
  R3.3 that remains genuinely open, and it is where circularity originally arose
  (Youden on the full data).
- No Triplex/Triplexator comparison, so no tool-vs-tool CI (R2.7's comparative half).

## Reproduce

```
python3 analysis/scripts/run_grouped_validation.py
python3 analysis/scripts/run_parameter_robustness.py
```
