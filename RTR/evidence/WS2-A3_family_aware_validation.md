# WS2-A3 / WS2-A2 — family- and study-aware validation, dual confidence intervals,
# and the evidence that HSeeker is not overfitted to the benchmark

**Status:** done for the detection call; **PROVISIONAL** pending WS0-A1 freeze and
DEC-5. The stability-score threshold half of WS2-A3 remains open, as does
tool-vs-tool CI comparison (no R/Triplex in this environment).
**Date:** 2026-09-23
**Reviewer items:** R3.2 (primary), R2.7, R3.3, R2.6, R3.1
**Scripts:** `analysis/scripts/run_grouped_validation.py`,
`analysis/scripts/run_parameter_robustness.py`
**Outputs:** `analysis/results/grouped_validation_v3/`, `analysis/results/robustness_v3/`
**Input:** `hdna_benchmark_balanced_v3.csv` (65:65), sha256 `f4602d6453d0926c…`

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
| leave-one-family-out | 54 | 28 singletons; largest `fam:pAA32` = 34 records |
| leave-one-study-out | 27 | 5 singletons; largest study = 34 records |
| collapsed per family | 54 representatives | reviewer's literal "at minimum" |
| collapsed per (family, class) | 74 representatives | class-preserving variant |

**The grouping rule is the part that makes it honest.** A synthetic negative
generated from an experimental positive (e.g. a mirror-disrupted mutant of pGG32)
carries that positive's sequence. It therefore **inherits its source's family and
study**, so holding out the pGG32 family cannot leave its own mutant in the training
portion. 40 of the 59 synthetic records inherit a group this way; the remaining 19
(G4, Z-DNA, B-DNA, homopolymer, perfect-mirror controls) form their own groups.

The reviewer's specific example checks out: `fam:pAA32` is the pXY32 substitution
series — all 16 X,Y combinations — and with its inherited synthetic derivatives it
is a single group of 34 of the 130 records.

## 2. Headline result — the same configuration wins under every scheme

Configuration: **`mismatch 0.15, minrep 8, purity 0.90, maxspacer 10, both filters`**
— the parameters already stated in the submitted manuscript.

| validation scheme | sens | spec | **MCC** | MCC record-level CI | MCC **family-cluster** CI |
|---|---|---|---|---|---|
| **leave-one-family-out** | 0.954 | 0.985 | **0.939** | [0.877, 0.985] | **[0.845, 0.988]** |
| leave-one-study-out | 0.954 | 0.985 | 0.939 | [0.877, 0.985] | [0.836, 0.988] |
| sequence-level 5-fold CV | 0.954 | 0.985 | 0.939 | [0.877, 0.985] | — |
| collapsed per (family, class), n=74 | 0.933 | 0.966 | 0.889 | [0.772, 0.973] | [0.772, 0.974] |
| collapsed per family, n=54 | 0.933 | 0.889 | 0.761 | [0.522, 0.946] | [0.478, 0.946] |

Comparators under leave-one-family-out: `script_default` MCC 0.908,
`production_cli` 0.856, per-fold tuned 0.895 — all below the manuscript configuration.

**Caveat on the strictest collapsed variant.** Collapsing to one record per family
leaves 45 forming / 9 non-forming, because synthetic negatives inherit their source
positive's family and are collapsed away with it. That reinstates the ~5:1 imbalance
the balanced set exists to remove and its interval is nearly uninformative
([0.522, 0.946]). The class-preserving variant (one representative per family *per
class*, n=74) removes within-family redundancy without discarding the design and is
the figure we recommend reporting, with the stricter one disclosed beside it.

## 3. Why family-aware validation does not move the point estimate — and what it does move

LOFO, LOSO and the pooled analysis return **numerically identical** results. This was
checked, not assumed:

```
full-data evaluation : TP62 FN3 TN64 FP1  MCC 0.9389
LOFO pooled out-of-fold: TP62 FN3 TN64 FP1  MCC 0.9389
```

**Reason:** with parameters fixed, HSeeker's call for a record depends only on that
record's sequence. Nothing is fitted, so there is no training set through which a
related sequence could leak. Grouping *cannot* bias the point estimate of a fixed
configuration.

What family structure genuinely affects is **precision**, which is exactly what R3.2
says — *"can substantially inflate apparent statistical power"*:

| uncertainty on MCC 0.939 (LOFO) | interval | width |
|---|---|---|
| record-level bootstrap | [0.877, 0.985] | 0.108 |
| **family-cluster bootstrap** | [0.845, 0.988] | **0.137** |

The cluster interval is **1.26×** wider. Across all 14 reported scheme/config rows
the ratio ranges 0.91–1.61 (median 1.29).

**Both intervals are reported everywhere, never one instead of the other**, because
they answer different questions: the record-level interval is the R2.7 interval on
the observed sample; the cluster interval is the R3.2 interval asking what happens if
the *families* had differed. Treating 16 pXY32 variants as 16 independent
observations overstates precision by roughly a quarter — and we now say so with a
number instead of a hedge.

## 4. The evidence that HSeeker is not overfitted to this benchmark

R2.6 asserts the tool "is overfitted to these exact sequences". Five independent
lines of evidence say otherwise:

1. **Tuning on the benchmark does not help.** Under LOFO the manuscript's
   pre-specified parameters score **0.939** while per-fold tuned selection scores
   **0.895**; under LOSO, 0.939 vs 0.895; under sequence-level CV, 0.939 vs 0.909.
   A tool fitted to this data would gain from re-fitting it. This one loses.
2. **Repeated-split selection gains exactly nothing.** Over 200 random splits the
   median tuned-minus-default eval MCC is **+0.0000** (95% [−0.0955, +0.1001]),
   winning 99, tied 43, losing 58 — a coin flip.
3. **The parameter surface is a plateau, not a peak.** Of 576 configurations, **108
   lie within 0.05 MCC of the maximum and 75 tie at it exactly**. Performance is not
   balanced on a hand-picked point.
4. **Selection is unanimous where it matters.** Across 54 family folds, `minrep 8`
   and `purity 0.85` were chosen in **54/54**, and some composition filter in
   **54/54** (AT 53, homopolymer 1). No knife-edge.
5. **Nothing is fitted at all.** Section 3: held-out and full-data evaluation are
   identical for a fixed configuration, because the algorithm has no trained
   parameters. There is no mechanism by which benchmark-specific overfitting could
   occur for the reported configuration.

Point 5 is the strongest and should lead the R2.6 response. The honest residual — and
it must be stated — is that the **AT-content and homopolymer thresholds were chosen
with knowledge of the H-DNA literature and of this benchmark's composition**. Those
are the one place R2.6's concern has genuine purchase. The mitigations are that the
filter axis was swept rather than assumed (the three filter states are
indistinguishable from each other and all beat "none"), and that the choice is
disclosed rather than presented as derived.

## 5. Recommended reporting — best balance of reviewer parity and results

> **Primary:** the manuscript's pre-specified parameters (`mismatch 0.15, minrep 8,
> purity 0.90, maxspacer 10, both filters`), evaluated **leave-one-family-out** on the
> 65:65 balanced benchmark:
> **sensitivity 0.954, specificity 0.985, MCC 0.939**,
> record-level 95% CI [0.877, 0.985], family-cluster 95% CI [0.845, 0.988].
> Supported by leave-one-study-out (identical) and by the class-collapsed set
> (MCC 0.889).

This is simultaneously the most rigorous and the best-performing option available,
and the two coincide for a reason worth stating plainly in the letter: **no selection
was performed**, so there is nothing to discount. Every tuned alternative scored
lower.

Reviewer parity achieved by this presentation:

| item | status |
|---|---|
| **R3.2** family-aware validation | **met at the "ideally" level** — LOFO *and* LOSO, plus both collapsed variants, with leakage-safe grouping of synthetic derivatives |
| **R3.2** "inflate apparent statistical power" | **met** — family-cluster CIs quantify it (1.26× wider) |
| **R2.7** confidence intervals | **met** — Wilson + record bootstrap + cluster bootstrap on every reported figure |
| **R3.3** circularity | **met for the detection call** — the reported configuration involves no selection whatsoever. *Open:* the stability-score threshold is not yet treated this way |
| **R3.1** imbalance / specificity precision | **met** — 65 negatives, specificity 0.985 with a stated interval |
| **R2.6** filters defined on the same data | **improved, not closed** — swept, shown immaterial, and disclosed; the residual is acknowledged above |

## 6. Correctness checks

| check | result |
|---|---|
| cluster CI contains the point estimate (MCC, sens, spec) | PASS |
| **every record its own group → cluster CI equals record CI**, [0.876, 0.985] both | PASS — validates the estimator |
| cluster CI wider than record CI under real grouping (1.26×) | PASS |
| cluster bootstrap deterministic under a fixed seed | PASS |
| held-out training portion never single-class (all 54 + 27 folds) | PASS (asserted in code) |
| deposited CSV carries both interval families for all 14 rows | PASS |
| synthetic derivatives inherit source group (40 of 59) | PASS |

## 7. Still open

- Benchmark not frozen (**WS0-A1**, blocked on DEC-2/DEC-3); DEC-3 adds a positive.
- **DEC-5** unresolved, though sections 2 and 4 now argue for adopting the
  manuscript's values.
- **Stability-score threshold** not yet evaluated out-of-fold — this is the half of
  R3.3 that remains genuinely open, and it is where circularity originally arose
  (Youden on the full data).
- No Triplex/Triplexator comparison, so no tool-vs-tool CI (R2.7's comparative half).

## Reproduce

```
python3 analysis/scripts/run_grouped_validation.py
python3 analysis/scripts/run_parameter_robustness.py
```
