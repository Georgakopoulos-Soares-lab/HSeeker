# WS1-A3 — normalize public sequence-field case (reviewer R2.10)

**Status:** done, verified
**Date:** 2026-09-23
**Reviewer item R2.10 (verbatim):** "the left_arm/right_arm/full_sequence come in
lowecase, while putative_triplex in uppercase. Am I right? If so, simple string
comparison would fail for a user not knowing this."

## Reviewer's reading confirmed

Cause located, both halves:

- `src/hseeker/_hdna.c:629-634` — the C core builds a **lowercased** working copy of
  the input (`dna[i] = tolower(raw_seq[i])`) and extracts `left_arm`, `spacer`,
  `right_arm`, `full_sequence` from that buffer (`_hdna.c:699-752`), so all four
  come back lowercase.
- `src/hseeker/_scoring.py:77` — the scorer calls `s = s.upper()` before building
  `putative_triplex`, so that field comes back uppercase.

A user comparing `full_sequence` with `putative_triplex` got False for no visible
reason. The reviewer is correct.

## Fix

Normalized to **uppercase at the Python boundary**, per plan WS1-A3.

- `src/hseeker/__init__.py` — added `_SEQUENCE_FIELDS` and
  `_normalize_sequence_case()`, applied in `scan_sequence()` immediately after the
  C call and before any filtering or scoring.
- The C core is **unchanged**: it keeps its lowercase internal buffer, so the
  detection algorithm carries no risk from this change.
- All FASTA entry points (`scan_fasta`, `scan_fasta_iter`, `scan_fasta_parallel`)
  route through `scan_sequence`, so one normalization point covers every public
  path including the CLI's TSV output.
- Return-contract docstring updated to state the guarantee.

Filter safety checked before the change, not after:
`_filter_at_content` already matched `("a","t","A","T")` (case-safe);
`_filter_homopolymer_triplex` operates on the already-uppercase `putative_triplex`;
`_scoring` uppercases internally. No filter depended on lowercase input.

## Verification

**1. The reviewer's exact comparison now succeeds** (input `"gggaaattaaaggg"`):

```
left_arm           GGGAAA
right_arm          AAAGGG
full_sequence      GGGAAATTAAAGGG
spacer             TT
putative_triplex   GGGAAA[TT]AAAGGG
full_sequence == putative_triplex bases?  True     # was False before
```

**2. Only case changed — detection and scoring untouched.** Raw C extension output
compared field-by-field against the public API for the same input:

```
hit count raw vs api:            1 1
non-sequence field differences:  none
sequence fields equal ignoring case: True
```

Every numeric field (`start`, `end`, `arm_length`, `spacer_length`, `ga_pct`,
`ct_pct`, `mirror_identity`, `is_perfect`, scores) is byte-identical. This change
cannot move any benchmark number.

**3. Filters still behave identically:**

```
poly-A, at_threshold=0.8          -> []            (still filtered)
poly-G, filter_homopolymers=True  -> 1 hit -> 0    (still filtered)
```

**4. All public paths, including mixed-case input** (`cccTTTAaaATTTCCC`):

```
scan_fasta          t1 GGGAAATTAAAGGG   t2 CCCTTTAAAATTTCCC
scan_fasta_iter     t1 GGGAAATTAAAGGG   t2 CCCTTTAAAATTTCCC
scan_fasta_parallel t1 GGGAAATTAAAGGG   t2 CCCTTTAAAATTTCCC
CLI TSV             t1 GGGAAA TT        t2 CCCTTT AAAA
```

**5. Test suite: 149 passed, 3 skipped** (chr1 regression tests skip — they require
`benchmarks/data/chr1.fa`). Was 147 passed / 3 skipped before.

## Test changes — disclosed deliberately

The old suite **asserted the buggy behaviour**, so those assertions had to change.
This is a contract change, not tests bent to fit code:

| Test | Before | After |
|---|---|---|
| `test_output_sequences_are_always_lowercase` | asserted fields `== .lower()` | replaced by `test_output_sequences_are_always_uppercase`, now also checking `spacer` and looping over upper/lower/mixed input |
| `test_gamir_left_arm_sequence` etc. (4 sites) | `"gggaaat"`, `"taaaggg"`, `"cccttta"`, `"atttccc"` | uppercase literals |
| `test_at_filter_boundary_...` | counted `b in "at"` | counts `b in "at"` on `.lower()` (test-local, filter itself is case-safe) |

Two regression tests added for the reviewer's scenario:
- `test_sequence_fields_comparable_with_putative_triplex` — cross-field comparison
- `test_case_of_input_does_not_change_output_fields` — upper/lower input give
  byte-identical sequence fields

## Files changed

- `src/hseeker/__init__.py` (+23 normalization, docstring)
- `tests/test_hdna.py` (7 edits, +2 tests)

## Reproduce

```
pip install -e . && python3 -m pytest -q tests
```
