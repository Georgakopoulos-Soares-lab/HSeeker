# WS4-A3 — parameter robustness with out-of-sample selection (reviewer R3.9)

**Regenerated 2026-10-07 on 128 records (64:64) after DEC-8 rebalance and DEC-9
score-informed overlap default; 124-record values in this file's 2026-10-02 version.**
Every number below is from the regenerated `analysis/results/robustness_v3/`
(`metadata.json` dated 2026-10-07 14:31:56). DEC-9 does not affect any cell here: the
analysis scores with `score=True` and uses only "any hit", and greedy and score give the
same binary call on every record under every configuration of this grid (verified: 0 of
576 × 128 = 73,728 calls differ, all four filter states). In general a call can differ
when the homopolymer filter is on — score can report a near-homopolymer hit where greedy
reports none (see DEC-9 evidence) — but no benchmark record triggers it. The changes
come entirely from the four DEC-8 negatives, which shift the stratified split and the
folds. The 2026-10-02 conclusions all hold: the filters are still no-ops, tuning never
beats the fixed script configuration on a random split (0/200), and the script default and
the manuscript-stated configuration are still tied at the top and both beat per-fold
tuning. What changed: classes are 1:1 again, and the single held-out split now selects
maxspacer 5 and loses to the script default by more (0.905 vs 1.000). DEC-10 (2026-10-08)
made greedy the overlap default again; since the two rules agree on all 73,728 grid
calls, it does not affect any cell either.

**Label note (2026-10-08).** `script_default` / "script default" in this file is the
**benchmark/script configuration** — the analysis scripts' setting (`REFERENCES["script_default"]`:
purity 0.90, mismatch 0.10, minrep 8, maxspacer 10, no filters). It is **not** the software
default. The shipped software defaults use minrep 10 (library: no filters; CLI: both
filters); on this benchmark the filters change no call, so they predict exactly like
`production_cli` (nested CV MCC 0.882). Earlier wording that called `script_default` "the
default" or "the shipped defaults" has been corrected.

**Regenerated 2026-10-02** after merging main (Nikol's scanner) and removing the 6
pure homopolymers; previous values in git history at `007288b`.

**Status:** done for the detection-parameter half of R3.9; **PROVISIONAL** pending
WS0-A1 freeze and DEC-5. The scoring-component ablation half of R3.9 remains open
(WS4-A2).
**Date:** 2026-09-23 (first run); regenerated 2026-10-02 and 2026-10-07
**Script:** `analysis/scripts/run_parameter_robustness.py` (deposited, deterministic)
**Outputs:** `analysis/results/robustness_v3/`
**Input:** `hdna_benchmark_balanced_v3.csv` (128 records, 64 forming / 64 non-forming),
sha256 `7bb495444e91709f…`

## Design (as instructed)

- **Data:** the balanced set only, all 128 kept records (64:64 — exactly 1:1 again
  after DEC-8 replaced the excluded homopolymer negatives with four synthetic ones).
- **Split:** **sequence level**, stratified by (label, record_type), 70/30,
  seed 20260923 → **train n=88 (44:44)**, **eval n=40 (20:20)**. Stratification keeps
  the input class ratio, so eval is exactly balanced again.
  Explicitly *not* family- or study-aware; the cost of that is measured below.
- **Grid:** full factorial, **576 configurations** —
  purity {0.85, 0.90, 0.95, 1.0} × mismatch {0.0, 0.05, 0.10, 0.15} ×
  minrep {8, 10, 12} × maxspacer {5, 10, 15} × filter state
  {none, AT, homopolymer, both}.
- **Selection:** maximum MCC **on train only**; ties broken deterministically toward
  higher sensitivity, then the *simpler* model (fewer filters, then larger minrep,
  smaller maxspacer, lower mismatch), so a tie never silently buys complexity.
- **Evaluation:** the selected configuration scored once on the held-out eval split,
  against three reference configurations, with 95% CIs.

## PRIMARY ESTIMATE — repeated stratified 5-fold nested CV

Added 2026-09-23 after the single holdout proved too small to estimate anything
precisely (eval n=40 now gives the selected config an MCC CI of [0.775, 1.000]).
Design fixed **a priori**: 5 folds, 20 repeats, stratified on (label, record_type),
parameters selected on the training folds only, every record held out exactly once
per repeat, predictions pooled across folds so each estimate covers all 128 records.

The partition was **not** chosen by comparing results — selecting a split because it
flatters the metric would reintroduce, at the split level, precisely the circularity
R3.3 objects to.

Out-of-fold, pooled over all 128 records:

| configuration | TP | FN | TN | FP | sens | spec | **MCC** | MCC boot 95% |
|---|---|---|---|---|---|---|---|---|
| nested-CV (per-fold tuned) | 61 | 3 | 63 | 1 | 0.953 | 0.984 | 0.938 | [0.875, 0.985] |
| **fixed: `script_default`** | 62 | 2 | 64 | 0 | 0.969 | **1.000** | **0.969** | **[0.923, 1.000]** |
| fixed: `production_cli` (= software defaults, minrep 10) | 56 | 8 | 64 | 0 | 0.875 | 1.000 | 0.882 | [0.802, 0.954] |
| **fixed: `manuscript_stated`** | 63 | 1 | 63 | 1 | **0.984** | 0.984 | **0.969** | **[0.922, 1.000]** |

The per-fold tuned row is the median over the 20 repeats (0.938, 2.5–97.5% range
across repeats [0.923, 0.969]); the confusion counts, sens/spec and bootstrap
interval are those of the first repeat's pooled predictions, whose MCC (0.938) is
also the median. The fixed rows do not vary across repeats.

**Both fixed configurations beat per-fold tuning** (0.969 vs 0.938), and they are
**tied with each other** (script_default 0.9692, manuscript_stated 0.9688). Tuning
adds selection variance without adding signal — the classic small-fold failure mode.
The direction of the old claim ("the manuscript's parameters beat tuning", previously
0.939 vs 0.909) still holds, but it is no longer unique to the manuscript
configuration: the script default does exactly as well. (124-record values: tuned
0.936, script_default 0.968, production 0.879, manuscript 0.968.)

### Which configuration to report — the data no longer decide it

The two tied configurations differ only in mismatch (0.10 vs 0.15) and filter state
(none vs both). On this benchmark the filter state is immaterial (see below), so the
real difference is mismatch, and it comes down to one record each way:

| | script_default (mm 0.10) | manuscript_stated (mm 0.15) |
|---|---|---|
| out-of-fold | 62/2/64/0 | 63/1/63/1 |
| sensitivity / specificity | 0.969 / 1.000 | 0.984 / 0.984 |
| MCC [boot 95%] | 0.969 [0.923, 1.000] | 0.969 [0.922, 1.000] |
| what differs | misses HDNA0058 | recovers HDNA0058, admits HDNA0028 `B-33` as FP |

The argument that remains for `manuscript_stated` (`mismatch 0.15, minrep 8,
purity 0.90, maxspacer 10, both filters`) is **procedural, not performance**: its
parameters were fixed in the submitted manuscript before this benchmark existed in
its current form, so R3.3's circularity objection does not apply — there is no
selection step to be circular about. The same is true of `script_default`, which is
also fixed and not selected on these data. Either is defensible; the choice is DEC-5's
and should be made on consistency with the manuscript, not on these numbers. What the
data do say is that **either fixed configuration is better than a tuned one**, and a
tuned number would have been *lower* and would have needed defending.

### What selection prefers (the algorithm-informing part)

Across all 100 outer folds:

| axis | modal choice | agreement |
|---|---|---|
| purity | 0.85 | **100%** |
| minrep | 8 | **100%** |
| filters | none | **100%** |
| maxspacer | 10 | 80% (5: 20%) |
| mismatch | 0.10 | 54% (0.15: 46%) |

(124-record values: maxspacer 10 76%, mismatch 0.10 60%.)

Selection is unanimous on minimum arm length 8 and purity 0.85. **The filter result
is not evidence against the filters** — it is a tie-break artefact: on this benchmark
the four filter states produce identical predictions for every one of the 144
non-filter parameter combinations (checked per record on all 128 records, 2026-10-07),
so the "fewer filters" tie-break always picks `none`. Previously (130 records) selection
chose AT in 86% and homopolymer in 14% of folds; that preference existed only because
of the homopolymer negatives, which are now removed. Selection is least settled on
mismatch tolerance — exactly the axis DEC-5 asks about, now close to an even split —
and out-of-fold the fixed 0.10 and 0.15 configurations tie.

---

## Secondary — single 70/30 holdout (run first, retained)

Selected on train: **purity 0.85, mismatch 0.15, minrep 8, maxspacer 5,
filters = none** (train MCC 0.9775; 28 configurations tie at that train maximum).
(124-record run: maxspacer 10, train MCC 1.000.)

Held-out eval (n=40, 20:20):

| configuration | TP/FN/TN/FP | sens | spec | MCC | MCC 95% CI |
|---|---|---|---|---|---|
| SELECTED (train-tuned) | 18/2/20/0 | 0.900 | 1.000 | 0.905 | [0.775, 1.000] |
| **`script_default`** (purity .90, mm .10, minrep 8, spacer 10, no filters) | 20/0/20/0 | **1.000** | **1.000** | **1.000** | [1.000, 1.000] |
| `production_cli` (minrep 10, both filters; = software defaults) | 19/1/20/0 | 0.950 | 1.000 | 0.951 | [0.854, 1.000] |
| `manuscript_stated` (mm **0.15**, minrep 8, both filters) | 20/0/19/1 | 1.000 | 0.950 | 0.951 | [0.854, 1.000] |

On this draw **the tuned configuration does worse than the script default** (0.905
vs 1.000) — another instance of tuning not beating the fixed script configuration. Train selection moved
maxspacer from 10 to 5, and on eval that shorter spacer loses two positives, HDNA0047
`pRW1713` and HDNA0048 `pRW1714` (both detected by every other reference
configuration, all of which use maxspacer 10). Tuning bought nothing; on this split it
cost 0.095 MCC. The single eval false positive under `manuscript_stated` is again
HDNA0028 `B-33`, which mismatch 0.15 admits at minrep 8; the single `production_cli`
miss is HDNA0023 (minrep 10). The script default's perfect score is one draw of 40
records — its CI is degenerate at [1.000, 1.000] only because nothing was missed —
and should not be quoted as a performance estimate; the nested CV above is the
estimate. (124-record holdout: SELECTED and manuscript_stated both 0.894,
script_default 0.949. Old 130-record holdout: SELECTED and manuscript_stated both
0.951, script_default 0.900.)

## The robustness finding — tuning never beats the fixed script configuration

The single split above is one draw. Repeating the whole procedure (split → select on
train → score on eval) over **200 random splits**:

```
selected-config eval MCC : median 0.951   [0.900, 1.000]
tuned - default eval MCC : median +0.0000 [-0.0955, +0.0000]
tuning beat the default in 0/200 splits, tied in 112, lost in 88
```

("default" in this output is `script_default`, the benchmark/script configuration with
minrep 8 and no filters, not the software default.)

The median improvement from tuning is **exactly zero**, and the interval does not
straddle zero — its upper end *is* zero. **In none of 200 splits did re-selecting
the parameters on training data beat the script default on held-out data**; in 88 it
was worse. (124-record run: beat 0, tied 109, lost 91. Original 130-record run: beat
99, tied 43, lost 58.) Selection always chose purity 0.85, minrep 8 and filters none;
it split between mismatch 0.15 (109 splits) and 0.10 (91), and between maxspacer 10
(160) and 5 (40). Of the 112 ties, 91 are splits where it re-selected mismatch 0.10 /
maxspacer 10 (the script configuration's values at purity 0.85) and 21 are splits where
it moved to mismatch 0.15 at maxspacer 10 and still scored exactly what the script
configuration scored on that eval split; all 88 losses are
splits where it moved to mismatch 0.15 (40 of them also to maxspacer 5). This is the
direct answer to R3.9's stated purpose — *"demonstrate that the results are not
dependent on one hand-selected parameterization"* — and it is stronger evidence than
any single "optimal" cell would have been.

## The surface is broad, not peaked

Eval MCC across all 576 configurations: min 0.577, median 0.734, max 1.000.
**104 of 576 configurations fall within 0.05 MCC of the maximum; 32 tie at it
exactly** (purity 0.85 or 0.90 × minrep 8 × {mismatch 0.10 with maxspacer 10 or 15,
or mismatch 0.0 or 0.05 with maxspacer 15}, each under all four filter states —
`script_default` is one of them). The maximum of 1.000 is a perfect score on a
40-record eval split and says more about this draw than about the configurations.

Marginal effect of each axis (spread in median eval MCC):

| axis | effect | detail |
|---|---|---|
| **minrep** | **0.203** | 8 → 0.838, 10 → 0.734, 12 → 0.635 |
| **mismatch** | **0.183** | 0.10/0.15 → 0.838, 0.05/0.0 → 0.655 |
| maxspacer | 0.183 | 15 → 0.837, 10 → 0.746, 5 → 0.655 |
| purity | 0.144 | 0.85/0.90 → 0.838; 0.95/1.0 → 0.694 |
| filter state | **0.000** | all four states → 0.734 (predictions identical) |

(124-record spreads: minrep 0.238, mismatch 0.182, maxspacer 0.161, purity 0.143,
filters 0.000.)

Two honest readings, both true and not in conflict:

1. **Parameters do matter** — minimum arm length still dominates (0.203 spread in
   median MCC), and a poor choice on minrep, mismatch or maxspacer (minrep 12,
   mismatch 0, maxspacer 5) costs 0.18–0.20 median MCC on that axis alone.
2. **The analysis scripts' configuration already sits on the optimal plateau** —
   `script_default` (minrep 8, no filters) attains the eval maximum — which is why
   re-selecting parameters out of sample yields no gain. The shipped software defaults
   do **not** sit on it: they use minrep 10 (median eval MCC 0.734 vs 0.838 at minrep 8),
   and `production_cli`, which predicts exactly like them here, scores 0.951 on this
   split and 0.882 in nested CV.

**The filter axis cannot be evaluated on this benchmark.** Previously the three filter
states were indistinguishable from each other and all beat "none", because every
false positive was a homopolymer. Those records are now removed by author decision
(the DEC-8 replacements are mutants, a shuffle and a random control, not
homopolymers), so there is nothing left for either filter to catch and all four
states give identical predictions. This is not evidence that the filters are useless —
only that this benchmark no longer contains the negatives they target.

### Bearing on DEC-5

`mismatch = 0.15` (manuscript) and `mismatch = 0.10` (script) are **tied** on median
eval MCC (0.838 each), and tied out-of-fold in nested CV (0.969 each). The data do
not favour the manuscript's value over the script's; they say both are on the
plateau and 0.0/0.05 are clearly worse. Adopting the manuscript's parameters as
canonical is supported by consistency with the submitted paper, and is not
contradicted by the data — but it is not *preferred* by them.

### Bearing on the shipped default (added 2026-10-08)

Both top configurations use minrep 8; the software ships minrep 10. The data here are
clear on that axis: selection picks minrep 8 in **100/100** nested-CV folds, and minrep
is the dominant parameter (spread 0.203 in median eval MCC; 8 → 0.838, 10 → 0.734). Out
of fold the shipped defaults give 0.882 against 0.969 for either minrep-8 configuration.
Against changing the default: minrep 10 reproduces the published chr1 count (96,729), so
a change means re-running chr1 and the runtime benchmarks. Which configuration the
manuscript reports, and whether the shipped minrep changes, is a `TODO(authors)` in the
plan (under WS2).

## Correctness checks

Re-run on the regenerated outputs (2026-10-07):

| check | result |
|---|---|
| train/eval disjoint, union = all 128 | PASS (88 + 40, 128 unique IDs) |
| eval stratified 20:20; train 44:44 (class ratio of input kept) | PASS |
| grid complete: 576 unique configs, no duplicates (train and eval grids) | PASS |
| selected config **is** a train argmax (28 configs tie on train at MCC 0.9775) | PASS |
| selected config's eval MCC ≤ grid eval maximum (0.905 ≤ 1.000) | PASS |
| Wilson CI recomputed by hand for sens 18/20 → [0.6990, 0.9721] | PASS (matches CSV to 1e-4) |
| TP+FN+TN+FP = n across all 1152 grid rows | PASS |
| MCC recomputed from each confusion matrix | PASS (max abs diff 1.1e-16) |
| determinism: full re-run into a scratch directory, all 10 data files byte-identical to the deposited ones (metadata.json differs only in its timestamp) | PASS |
| prediction matrix vs direct `hseeker.scan_sequence` calls, 6 random eval-grid cells | PASS |
| fixed reference configs on all 128 records vs direct `scan_sequence` (62/2/64/0, 56/8/64/0, 63/1/63/1) | PASS |

Two checks failed on the first pass of the original run (2026-09-23) and were
**defects in the checking code, not the analysis**: one matched the label
`non_forming` while the loader normalises to `non-forming`; the other compared a 4-dp
rounded MCC against an unrounded one. Both re-verified after correction. Recorded here
rather than quietly fixed.

On this split the train-selected config does **not** reach the eval maximum (0.905 vs
1.000); in the original run it happened to tie it. Either outcome is luck on one draw
— selection saw only train, and the repeated-split analysis is the honest estimate of
the procedure.

## Reviewer-comment parity

| item | status | detail |
|---|---|---|
| **R3.9** parameter thresholds | **met** | full factorial over purity/mismatch/minrep/maxspacer, metrics per setting, marginals reported |
| **R3.9** component ablation | **partial** | filter states cover "detection alone" and "detection + composition filtering", but on the current benchmark they are indistinguishable (no homopolymer negatives remain). **Not covered:** pairing-only vs pairing+stacking scoring ablation — that is WS4-A2 and still open (Nikol's `hseeker.ablation` / `analysis/scripts/ablation_eval.py` exist but have not been run on v3) |
| **R3.3** circularity | **met for detection parameters** | selection strictly on train, evaluation on untouched held-out; procedure re-estimated over 200 splits. **Note:** R3.3 also asks for leave-one-study/family-out, which this design deliberately does not use (see WS2-A3) |
| **R3.2** family-aware validation | **NOT met here, by instruction** | see leakage audit below; addressed in WS2-A3 |
| **R2.7** confidence intervals | **met** | Wilson CIs for sensitivity/specificity, bootstrap CI for MCC (record-level resampling) |
| **R2.1 / R3.1** imbalance | **met** | 64:64 throughout; eval 20:20, so specificity rests on 20 held-out / 64 total negatives rather than 5 experimental ones |
| **R2.6** filters fitted on the same data | **open, and now untestable here** | the filter axis was swept, but on the current benchmark it has no effect, so the benchmark neither supports nor refutes the filter thresholds. The disclosure R2.6 asks for is still required, and should say the benchmark excludes pure homopolymers |

## Leakage audit — the cost of a sequence-level split

Measured, not assumed. Of the 40 eval records:

- **25/40** belong to a family that also appears in train
- **8/40** have a ≥90% identical sequence in train
- **7/40** are synthetic negatives whose **source record** is in train

So the held-out estimate is optimistic relative to a family-aware design. This is
exactly what R3.2 predicts, and it means these numbers answer R3.9 and R3.3 but
**must not be presented as answering R3.2**. The family-aware variant has since been
run (WS2-A3); the grouping metadata covers 25 studies and 47 families among the
experimental records — with the caveat that 42 of 47 experimental families are
singletons and `fam:pAA32` alone holds 16 records, so folds are very uneven.

## Why still provisional

1. Benchmark not frozen (**WS0-A1**, blocked on DEC-2/DEC-3); if adopted, DEC-3 would
   add a positive → 129 records (65:64), changing every cell. The 2026-10-02
   homopolymer removal and the 2026-10-07 DEC-8 rebalance have already changed every
   cell twice.
2. **DEC-5** unresolved — this analysis informs it but, after regeneration, does not
   prefer mismatch 0.15 over 0.10. Nor is it decided which configuration the manuscript
   reports or whether the shipped minrep changes from 10 to 8 (see above).
3. HSeeker only. No Triplex/Triplexator comparison: R is not installed in this
   environment, so no tool-vs-tool claim is supported.
4. Detection call only; the stability-score threshold (the actual subject of R3.3)
   is not swept here.

## Reproduce

```
python3 analysis/scripts/run_parameter_robustness.py \
    --csv hdna_benchmark_balanced_v3.csv --train-frac 0.7 --seed 20260923 --repeats 200
```
