# WS1-A3 — normalize public sequence-field case (reviewer R2.10)

**Status (2026-10-02):** resolved — RTR implementation superseded by Nikol's fix on main (merged)
**Date:** 2026-09-23 (RTR fix); updated 2026-10-02 after merging main
**Reviewer item R2.10 (verbatim):** "the left_arm/right_arm/full_sequence come in
lowecase, while putative_triplex in uppercase. Am I right? If so, simple string
comparison would fail for a user not knowing this."

## Reviewer's reading confirmed

Cause located, both halves (line numbers as of the merged tree, 2026-10-02):

- `src/hseeker/_hdna.c:646-651` — the C core builds a **lowercased** working copy of
  the input (`dna[i] = tolower(raw_seq[i])`) and extracts `left_arm`, `spacer`,
  `right_arm`, `full_sequence` from that buffer (`_hdna.c:716-772`), so all four
  come back lowercase.
- `src/hseeker/_scoring.py` — the scorer calls `s = s.upper()` (`:84`) and, before the
  fix, built `putative_triplex` from that uppercase string, so that field came back
  uppercase.

A user comparing `full_sequence` with `putative_triplex` got False for no visible
reason. The reviewer is correct.

## Resolution (current tree) — everything lowercase

All public sequence fields — `left_arm`, `spacer`, `right_arm`, `full_sequence`
**and `putative_triplex`** — are lowercase, whatever the input case.

- The C core is unchanged: it already returned the four detector fields lowercase.
- `src/hseeker/_scoring.py:145-149` — the scorer still works in uppercase internally
  but lowercases `putative_triplex` on output; the `score_hdna()` docstring
  (`_scoring.py:183`) states it.
- `scan_sequence()` return-contract docstring (`src/hseeker/__init__.py:322-323`):
  "All sequence fields, including ``putative_triplex``, are lowercase regardless of
  input case."
- Because the fix sits in the scorer, every path that scores agrees:
  `scan_sequence`, `scan_fasta`, `scan_fasta_iter`, `scan_fasta_parallel`,
  `hseeker.ablation.score_candidates`, and the CLI TSV.

Filters are case-safe: `_filter_at_content` counts `("a","t","A","T")`
(`__init__.py:80`); the homopolymer check compares `base.lower()`
(`__init__.py:201-202`).

## Verification (merged tree, 2026-10-02)

**1. The reviewer's comparison succeeds** (same result for `"gggaaattaaaggg"` and
`"GGGAAATTAAAGGG"`):

```
left_arm gggaaa  spacer tt  right_arm aaaggg
full_sequence      gggaaattaaaggg
putative_triplex   gggaaa[tt]aaaggg
full_sequence == putative_triplex bases?  True
```

**2. All public paths, mixed-case input** (`t1 gggaaattaaaggg`, `t2 cccTTTAaaATTTCCC`,
minrep 3):

```
scan_fasta / scan_fasta_iter / scan_fasta_parallel / CLI TSV:
  t1  full_sequence gggaaattaaaggg    putative_triplex gggaaa[tt]aaaggg
  t2  full_sequence ccctttaaaatttccc  putative_triplex gggaaa[tttt]aaaggg
```

(For a CT-rich left arm the scorer reports `putative_triplex` on the GA-rich strand,
as before; only the case was at issue.)

**3. Regression tests** — `tests/test_output_case.py` (2 tests, both pass):

- `test_sequence_and_ablation_scores_match_detected_case` — all five sequence fields
  from `scan_sequence`, and `putative_triplex` from `score_candidates`, are lowercase.
- `test_cli_uses_same_case_for_detected_and_scored_motifs` — same for the CLI
  `_HDNA.tsv`.

**4. Full suite on the merged tree: 211 passed, 3 skipped** (chr1 tests skip without
`benchmarks/data/chr1.fa`).

## Draft response (R2.10)

The reviewer is right. All sequence fields, including `putative_triplex`, are now
consistently lowercase (in every API path and in the CLI output), so direct string
comparison between fields works; this is documented in the `scan_sequence`
docstring and covered by regression tests. The fix only lowercases
`putative_triplex` after it is built, so it changes no coordinate or score.

## History (superseded RTR implementation)

RTR (2026-09-23) normalised the other way: `_normalize_sequence_case()` in
`src/hseeker/__init__.py` uppercased `left_arm`/`spacer`/`right_arm`/`full_sequence`
at the Python boundary, matching the scorer's uppercase `putative_triplex`. It was
verified then (numeric fields byte-identical to the raw C output, filters unchanged,
all paths uppercase; suite 149 passed, 3 skipped) and changed the existing tests to
expect uppercase (`test_output_sequences_are_always_uppercase`, uppercase arm
literals, plus `test_sequence_fields_comparable_with_putative_triplex` and
`test_case_of_input_does_not_change_output_fields`). Nikol fixed the same defect
independently on main by lowercasing `putative_triplex`; on merging, the authors chose
Nikol's version (2026-10-02). It keeps the original lowercase contract of the
detector fields, so existing user code that relied on lowercase arms is unaffected.
`_normalize_sequence_case` and RTR's uppercase tests were removed; the reviewer's
bug (cross-field `==` failing) is fixed either way.

## Files changed (current tree, Nikol's fix)

- `src/hseeker/_scoring.py` (lowercase `putative_triplex`, docstring)
- `src/hseeker/__init__.py` (`scan_sequence` docstring)
- `tests/test_output_case.py` (new, 2 tests)

## Reproduce

```
pip install -e . && python3 -m pytest -q tests/test_output_case.py
```
