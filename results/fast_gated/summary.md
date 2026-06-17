# FastGatedHSeeker Benchmark

## Dataset

- Total sequences: 54
- Forming: 45
- Non-forming: 9
- Score threshold: 60

## Primary Metrics

| Method | TP | FN | TN | FP | Sensitivity | Specificity | Precision | F1 | Accuracy | Runtime sec | Inner iterations |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| fast_gated_fast_8x | 38 | 7 | 8 | 1 | 0.844 | 0.889 | 0.974 | 0.905 | 0.852 | 0.000892 | 74863 |
| fast_gated_fast_4x | 38 | 7 | 8 | 1 | 0.844 | 0.889 | 0.974 | 0.905 | 0.852 | 0.000920 | 74863 |
| fast_gated_exact_gate | 38 | 7 | 8 | 1 | 0.844 | 0.889 | 0.974 | 0.905 | 0.852 | 0.000991 | 88117 |
| original | 38 | 7 | 8 | 1 | 0.844 | 0.889 | 0.974 | 0.905 | 0.852 | 0.000993 | 82162 |

## Speedup

- Fast 4x gate vs original wall-clock speedup on this dataset: 1.08x
- Fast 4x gate vs original on synthetic 20 kb balanced-DNA stress input: 8.17x
- Exact-gate F1: 0.905
- Fast 4x F1: 0.905

## Stress Runtime

| Method | Runtime sec | Inner iterations | Center/spacer pairs | Gate passed | Raw hits | Final hits |
|---|---:|---:|---:|---:|---:|---:|
| original | 0.179658 | 136238669 | 1017756 | 0 | 1 | 1 |
| fast_gated_fast_4x | 0.021998 | 12166734 | 1017756 | 21 | 1 | 1 |
| fast_gated_exact_gate | 1.467354 | 967870087 | 1017756 | 277 | 1 | 1 |


## Lost Calls In Fast Mode

- None at score threshold 60.

## Failure Cases

- HDNA0022 original label=forming score=0.0
- HDNA0023 original label=forming score=0.0
- HDNA0024 original label=forming score=0.0
- HDNA0025 original label=forming score=0.0
- HDNA0027 original label=forming score=0.0
- HDNA0030 original label=forming score=0.0
- HDNA0034 original label=forming score=0.0
- HDNA0053 original label=non-forming score=83.3
- HDNA0022 fast_gated_exact_gate label=forming score=0.0
- HDNA0023 fast_gated_exact_gate label=forming score=0.0
- HDNA0024 fast_gated_exact_gate label=forming score=0.0
- HDNA0025 fast_gated_exact_gate label=forming score=0.0
- HDNA0027 fast_gated_exact_gate label=forming score=0.0
- HDNA0030 fast_gated_exact_gate label=forming score=0.0
- HDNA0034 fast_gated_exact_gate label=forming score=0.0
- HDNA0053 fast_gated_exact_gate label=non-forming score=83.3
- HDNA0022 fast_gated_fast_4x label=forming score=0.0
- HDNA0023 fast_gated_fast_4x label=forming score=0.0
- HDNA0024 fast_gated_fast_4x label=forming score=0.0
- HDNA0025 fast_gated_fast_4x label=forming score=0.0
- HDNA0027 fast_gated_fast_4x label=forming score=0.0
- HDNA0030 fast_gated_fast_4x label=forming score=0.0
- HDNA0034 fast_gated_fast_4x label=forming score=0.0
- HDNA0053 fast_gated_fast_4x label=non-forming score=83.3
- HDNA0022 fast_gated_fast_8x label=forming score=0.0
- HDNA0023 fast_gated_fast_8x label=forming score=0.0
- HDNA0024 fast_gated_fast_8x label=forming score=0.0
- HDNA0025 fast_gated_fast_8x label=forming score=0.0
- HDNA0027 fast_gated_fast_8x label=forming score=0.0
- HDNA0030 fast_gated_fast_8x label=forming score=0.0
- HDNA0034 fast_gated_fast_8x label=forming score=0.0
- HDNA0053 fast_gated_fast_8x label=non-forming score=83.3

## Triplex

- Triplex comparison run: yes
- Output: `triplex_comparison.csv`

## Recommended Gate Defaults

- `gate_window = minrep`
- `gate_mirror_frac = 1 - mismatch`
- `gate_purity_frac = purity`
- Use `exact_gate` (`gate_search_limit=maxrep`, purity prefilter disabled) when recall/equivalence matters.
- Use `fast` with `gate_search_limit=4*minrep` for speed-sensitive scans, while reviewing `lost_calls` in this report.
