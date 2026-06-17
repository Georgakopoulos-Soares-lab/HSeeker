# FastGatedHSeeker Next-Step Benchmark

## What Was Tested

- Real local FASTA chunks at 100 kb, 1 Mb, and 5 Mb when available.
- Synthetic adversarial sequences at 20 kb and 100 kb.
- Validation CSV gate sweep at fixed HSeeker biological parameters.

## Validation Preservation

- Original F1: 0.905
- Fast 2x F1: 0.905
- Fast 4x F1: 0.905
- Fast 8x F1: 0.905
- Fast 16x F1: 0.905
- Lost calls: No lost validation calls.

## Real-Genome Runtime

- chrY.hs1.fa:chrY:first_100000: 3.99x speedup, lost_raw=0, lost_final=0, extra_final=0
- chrY.hs1.fa:chrY:first_1000000: 3.32x speedup, lost_raw=0, lost_final=0, extra_final=0
- chrY.hs1.fa:chrY:first_5000000: 3.78x speedup, lost_raw=5, lost_final=0, extra_final=0

## Adversarial Runtime

Hardest classes for fast 4x by speedup:
- pure_ttc_repeat_100000: 1.20x speedup, lost_final=0
- pure_ttc_repeat_20000: 1.20x speedup, lost_final=0
- pure_gaa_repeat_100000: 1.22x speedup, lost_final=0
- noisy_gaa_repeat_5pct_100000: 1.22x speedup, lost_final=0
- noisy_ttc_repeat_5pct_100000: 1.23x speedup, lost_final=0

## Gate Limit Tradeoff

Validation data:
- fast_2x: lost_calls=0, changed_best_score=0, changed_arm_spacer=0
- fast_4x: lost_calls=0, changed_best_score=0, changed_arm_spacer=0
- fast_8x: lost_calls=0, changed_best_score=0, changed_arm_spacer=0
- fast_16x: lost_calls=0, changed_best_score=0, changed_arm_spacer=0
- fast_maxrep: lost_calls=0, changed_best_score=0, changed_arm_spacer=0

Adversarial final-call or best-call differences:
- noisy_ttc_repeat_5pct_100000 / fast_gated_fast_2x: lost_final=4, extra_final=4, changed_best_score=0, changed_arm_spacer=0

The 2x/4x/8x/16x gate settings are all reported in `gate_tradeoff_validation.csv`.
Lost calls are explicitly listed in `gate_lost_calls_validation.csv`.

## Recommended Default Gate

Use `gate_search_limit = 4 * minrep` as a practical default if validation-call preservation remains true for the target dataset. Use `8 * minrep` when prioritizing recall margin over speed.

## Caveats

- These benchmarks are local and parameter-specific.
- GA/CT-rich and repeat-heavy backgrounds can reduce speedup because many pairs pass purity and gate checks.
- Real-genome superiority should only be claimed for the tested local FASTA chunks, not genome-wide in general.
- Exact-gate mode is useful for equivalence checks but can be slower than original on adversarial inputs.
