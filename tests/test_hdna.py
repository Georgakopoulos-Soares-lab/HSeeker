"""
tests/test_hdna.py — comprehensive pytest suite for hseeker.

Covers:
  1.  Empty / trivial inputs
  2.  Known GA motif — basic detection
  3.  Strict mode (purity=1.0, mismatch=0.0)
  4.  Relaxed mode — imperfect mirrors
  5.  Coordinate offset (seq_offset)
  6.  Multi-record FASTA via scan_fasta()
  7.  Overlap removal
  8.  Parameter validation
  9.  Known-sequence exact field values (GGGAAATTAAAGGG)
  10. CT arm detection
  11. Purity filter boundary
  12. maxspacer boundary
  13. Case insensitivity
  14. Output field self-consistency (total_length, full_sequence, arm bounds)
  15. Flanking N bases — N stops extension; coordinates remain correct
  16. Multiple non-overlapping motifs in one sequence
  17. is_perfect flag logic
  18. Mirror identity with known mismatch count
  19. minrep / maxrep arm-length boundaries
  20. Genomic coordinate propagation via FASTA offset header
  21. parse_fasta edge cases
  22. Determinism / reproducibility
  23. CLI integration (python -m hseeker)
  24. Parallel chunking edge cases (scan_fasta_parallel)
  25. Real-genome regression (hg38 chr1 — requires benchmarks/data/chr1.fa)

Reference values are derived by running the algorithm and recording observed output.

Run:
    pytest -v tests/
or:
    conda run -n biomni_e1 pytest -v tests/
"""

from __future__ import annotations

import csv
import subprocess
import sys
import tempfile
from pathlib import Path

import pytest

import hseeker


# ===========================================================================
# Helpers
# ===========================================================================

def write_fasta(records: list[tuple[str, str]], line_width: int = 80) -> str:
    lines = []
    for name, seq in records:
        lines.append(f">{name}")
        if line_width <= 0:
            lines.append(seq)
        else:
            for i in range(0, len(seq), line_width):
                lines.append(seq[i : i + line_width])
    return "\n".join(lines) + "\n"


def fasta_to_tmp(records: list[tuple[str, str]], line_width: int = 80) -> Path:
    text = write_fasta(records, line_width=line_width)
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    tmp.write(text)
    tmp.close()
    return Path(tmp.name)


# ---------------------------------------------------------------------------
# Reference sequences with verified expected values.
#
# GAMIR_SEQ = "GGGAAATTAAAGGG"  (14 bp)
#   The algorithm finds the best hit: arm=7, sp=0, 0 mismatches.
#   left_arm  = "gggaaat"  (dna[0:7])
#   spacer    = "."        (spacer_length=0)
#   right_arm = "taaaggg"  (dna[7:14])  — reverse("gggaaat") = "taaaggg" ✓
#   ga_pct    = 6/7 × 100 ≈ 85.71%   (right arm: t=CT, a,a,a,g,g,g = 6 GA)
#   ct_pct    = 1/7 × 100 ≈ 14.29%
#   mirror_id = 100.0%,  is_perfect = False  (ga_pct ≠ 100%)
#   start=1, end=14  (with seq_offset=1)
#
# CTMIR_SEQ = "CCCTTTAATTTCCC"  (14 bp)
#   Best hit: arm=7, sp=0, 0 mismatches.
#   left_arm  = "cccttta",  right_arm = "atttccc"
#   ct_pct    ≈ 85.71%  (right arm: a=GA, t,t,t,c,c,c = 6 CT)
#   ga_pct    ≈ 14.29%,  mirror_id = 100.0%,  is_perfect = False
#
# IMPMIR_SEQ = "AAAAAACAAAAA"  (12 bp)
#   Intentionally imperfect: arm=6, sp=0, 1 mismatch (dna[5]='a' vs dna[6]='c').
#   mirror_id = 5/6 × 100 ≈ 83.33%,  is_perfect = False.
#   Found with mismatch=0.20; NOT found with mismatch=0.0.
#
# IMPURE_SEQ = "GGTAAATTAAATGG"  (14 bp)
#   Perfect mirror (0 mismatches), arm=6, sp=2, but ga_pct=5/6≈83.3%.
#   Found with purity=0.80; NOT found with purity=0.90.
#   is_perfect = False (ga_pct ≠ 100%).
#
# PERFECT_SEQ = "A" * 14  (14 bp)
#   All-A sequence: arm=7, sp=0, ga_pct=100%, mirror_id=100%.
#   is_perfect = True.
#
# LOWPUR_SEQ = "GTTTAAAATTTG"  (12 bp)
#   arm=6, sp=0, mirror_id=100%, ga_pct=50%, ct_pct=50%.
#   Found with purity=0.50; NOT found with purity=0.80.
#
# SPACER8_SEQ = "GGGGGG" + "T"*8 + "GGGGGG"  (20 bp)
#   Perfect arm=6, sp=8, ga_pct=100%.  Requires maxspacer >= 8.
#   With mismatch=0.0: NOT found when maxspacer=7; found when maxspacer=8.
# ---------------------------------------------------------------------------

PURE_GA    = "GAGAGAGAGAGAGAGAGAGAGA"           # 22 bp — general GA mirror
GAMIR_SEQ  = "GGGAAATTAAAGGG"                   # 14 bp — exact-value fixture
CTMIR_SEQ  = "CCCTTTAATTTCCC"                   # 14 bp — CT arm fixture
IMPMIR_SEQ = "AAAAAACAAAAA"                      # 12 bp — imperfect mirror fixture
IMPURE_SEQ = "GGTAAATTAAATGG"                   # 14 bp — impure-arm fixture
PERFECT_SEQ = "A" * 14                          # 14 bp — is_perfect=True fixture
LOWPUR_SEQ = "GTTTAAAATTTG"                     # 12 bp — low-purity fixture
# SPACER8_SEQ: arms = "AAAAAA" (6 pure-A), spacer = "CTTTTTTT" (asymmetric 8-char)
# The 'C' at spacer[0] ≠ 'T' at spacer[7] prevents any arm from forming at sp ≤ 7.
# With mismatch=0.0: no hit when maxspacer=7; arm=6,sp=8,is_perfect=True when maxspacer=8.
SPACER8_SEQ = "AAAAAA" + "CTTTTTTT" + "AAAAAA"  # 20 bp — spacer boundary fixture


# ===========================================================================
# 1. Empty / trivial inputs
# ===========================================================================

def test_empty_sequence():
    assert hseeker.scan_sequence("") == []


def test_single_base():
    assert hseeker.scan_sequence("A") == []


def test_too_short_for_minrep():
    """Sequence shorter than 2 × minrep cannot contain any arm."""
    assert hseeker.scan_sequence("GAGAGA", minrep=10) == []


def test_all_n_bases():
    """'N' breaks arm extension — should produce no hits."""
    assert hseeker.scan_sequence("N" * 100, minrep=6) == []


def test_all_same_base_can_form_mirror():
    """A run of identical bases forms a trivial mirror; must be detected."""
    hits = hseeker.scan_sequence("A" * 30, minrep=6)
    assert len(hits) > 0


# ===========================================================================
# 2. Known GA motif — basic detection
# ===========================================================================

def test_pure_ga_mirror_detected():
    hits = hseeker.scan_sequence(PURE_GA, minrep=6)
    assert len(hits) > 0, "Expected at least one H-DNA hit on a pure GA mirror"


def test_pure_ga_all_required_keys_present():
    required = {
        "start", "end", "arm_length", "spacer_length", "total_length",
        "ga_pct", "ct_pct", "mirror_identity", "is_perfect",
        "left_arm", "spacer", "right_arm", "full_sequence",
    }
    hits = hseeker.scan_sequence(PURE_GA, minrep=6)
    assert len(hits) > 0
    for h in hits:
        missing = required - h.keys()
        assert not missing, f"Missing keys: {missing}"


def test_pure_ga_coords_are_valid():
    for h in hseeker.scan_sequence(PURE_GA, minrep=6):
        assert h["start"] >= 1
        assert h["end"] >= h["start"]
        assert h["arm_length"] >= 6
        assert h["spacer_length"] >= 0


def test_pure_ga_is_perfect_flag():
    """PURE_GA is all-GA — strict mode must return at least one is_perfect=True hit."""
    hits = hseeker.scan_sequence(PURE_GA, minrep=6, purity=1.0, mismatch=0.0)
    assert any(h["is_perfect"] for h in hits)


# ===========================================================================
# 3. Strict mode (purity=1.0, mismatch=0.0)
# ===========================================================================

PERFECT_MR = "GAGAGAG" + "GAGAGAG"   # pure GA direct adjacency


def test_strict_mode_detects_perfect():
    hits = hseeker.scan_sequence(PERFECT_MR, minrep=6, purity=1.0, mismatch=0.0)
    assert len(hits) > 0


def test_strict_mode_all_hits_have_100pct_mirror_identity():
    """With mismatch=0 every returned hit must have 100% mirror identity."""
    seq = "GAGAGAGCAGAGAG"
    for h in hseeker.scan_sequence(seq, minrep=6, purity=1.0, mismatch=0.0):
        assert h["mirror_identity"] == pytest.approx(100.0, abs=0.01), \
            "Strict mode returned a non-perfect mirror hit"


def test_strict_mode_all_hits_are_flagged_perfect():
    """GAMIR_SEQ with purity=1.0 finds the arm=6, sp=2, 100%-pure hit → is_perfect=True."""
    for h in hseeker.scan_sequence(GAMIR_SEQ, minrep=6, purity=1.0, mismatch=0.0):
        assert h["is_perfect"] is True


# ===========================================================================
# 4. Relaxed mode
# ===========================================================================

def test_relaxed_finds_at_least_as_many_as_strict():
    seq = "GAGAGAGCAGAGAG"
    strict  = hseeker.scan_sequence(seq, minrep=6, purity=1.0, mismatch=0.0)
    relaxed = hseeker.scan_sequence(seq, minrep=6, purity=0.8, mismatch=0.2)
    assert len(relaxed) >= len(strict)


def test_relaxed_can_find_imperfect_hit():
    """IMPMIR_SEQ has 1 mismatch — found under relaxed settings."""
    hits = hseeker.scan_sequence(IMPMIR_SEQ, minrep=6, purity=0.80, mismatch=0.20)
    assert len(hits) > 0, "Relaxed mode should find the 1-mismatch mirror"


def test_relaxed_mirror_identity_below_100():
    """IMPMIR_SEQ hit must have mirror_identity < 100% (1 mismatch out of 6)."""
    hits = hseeker.scan_sequence(IMPMIR_SEQ, minrep=6, purity=0.80, mismatch=0.20)
    imperfect = [h for h in hits if h["mirror_identity"] < 100.0]
    assert len(imperfect) > 0, "Expected at least one hit with mirror_identity < 100%"


# ===========================================================================
# 5. Coordinate offset (seq_offset)
# ===========================================================================

def test_seq_offset_shifts_start_and_end():
    offset   = 1_000_000
    hits_1   = hseeker.scan_sequence(PURE_GA, minrep=6, seq_offset=1)
    hits_off = hseeker.scan_sequence(PURE_GA, minrep=6, seq_offset=offset)
    assert len(hits_1) == len(hits_off)
    for h1, ho in zip(hits_1, hits_off):
        assert ho["start"] == h1["start"] + (offset - 1)
        assert ho["end"]   == h1["end"]   + (offset - 1)


def test_seq_offset_does_not_change_hit_count():
    for off in (1, 500, 43_585_222):
        hits = hseeker.scan_sequence(GAMIR_SEQ, minrep=6, seq_offset=off)
        assert len(hits) == len(
            hseeker.scan_sequence(GAMIR_SEQ, minrep=6, seq_offset=1)
        )


def test_seq_offset_does_not_change_arm_sequences():
    """Changing seq_offset must not alter the arm sequences — only coordinates shift."""
    for off in (1, 1_000_000):
        hits = hseeker.scan_sequence(GAMIR_SEQ, minrep=6, seq_offset=off)
        arm7 = [h for h in hits if h["arm_length"] == 7]
        if arm7:
            assert arm7[0]["left_arm"]  == "gggaaat"
            assert arm7[0]["right_arm"] == "taaaggg"


# ===========================================================================
# 6. Multi-record FASTA via scan_fasta()
# ===========================================================================

def test_scan_fasta_multi_record_seq_ids():
    path = fasta_to_tmp([("seq_A", PURE_GA * 3), ("seq_B", PURE_GA * 3)])
    try:
        hits = hseeker.scan_fasta(str(path), minrep=6)
        assert len(hits) > 0
        ids = {h["seq_id"] for h in hits}
        assert "seq_A" in ids
        assert "seq_B" in ids
    finally:
        path.unlink(missing_ok=True)


def test_scan_fasta_all_n_gives_no_hits():
    path = fasta_to_tmp([("empty", "N" * 80)])
    try:
        assert hseeker.scan_fasta(str(path), minrep=6) == []
    finally:
        path.unlink(missing_ok=True)


def test_scan_fasta_hit_count_matches_individual_scans():
    """scan_fasta result must equal the sum of individual scan_sequence calls."""
    path = fasta_to_tmp([("r1", GAMIR_SEQ), ("r2", CTMIR_SEQ)])
    try:
        fasta_hits = hseeker.scan_fasta(str(path), minrep=6)
        expected = (
            len(hseeker.scan_sequence(GAMIR_SEQ, minrep=6)) +
            len(hseeker.scan_sequence(CTMIR_SEQ, minrep=6))
        )
        assert len(fasta_hits) == expected
    finally:
        path.unlink(missing_ok=True)


def test_parse_fasta_offset_from_header():
    text = ">chr1:43585222-43586222\nGAGAGAGAGAGAGAGAGAGAGA\n"
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    tmp.write(text)
    tmp.close()
    path = Path(tmp.name)
    try:
        parsed = list(hseeker.parse_fasta(str(path)))
        assert len(parsed) == 1
        seq_id, seq, offset = parsed[0]
        assert seq_id == "chr1:43585222-43586222"
        assert offset == 43585222
        assert seq == "GAGAGAGAGAGAGAGAGAGAGA"
    finally:
        path.unlink(missing_ok=True)


def test_scan_fasta_offset_propagates_to_coordinates():
    """Coordinates from scan_fasta must be shifted by the FASTA header offset."""
    genomic_start = 1_000_000
    text = f">chr1:{genomic_start}-{genomic_start + 200}\n{GAMIR_SEQ}\n"
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    tmp.write(text)
    tmp.close()
    path = Path(tmp.name)
    try:
        fasta_hits  = hseeker.scan_fasta(str(path), minrep=6)
        direct_hits = hseeker.scan_sequence(GAMIR_SEQ, minrep=6, seq_offset=1)
        assert len(fasta_hits) == len(direct_hits)
        for fh, dh in zip(fasta_hits, direct_hits):
            assert fh["start"] == dh["start"] + (genomic_start - 1)
            assert fh["end"]   == dh["end"]   + (genomic_start - 1)
    finally:
        path.unlink(missing_ok=True)


# ===========================================================================
# 6b. scan_fasta_iter() — streaming generator
# ===========================================================================

def test_scan_fasta_iter_is_generator():
    """scan_fasta_iter must return a generator, not a list."""
    import types
    path = fasta_to_tmp([("s", GAMIR_SEQ)])
    try:
        result = hseeker.scan_fasta_iter(str(path), minrep=6)
        assert isinstance(result, types.GeneratorType), \
            f"Expected GeneratorType, got {type(result)}"
    finally:
        path.unlink(missing_ok=True)


def test_scan_fasta_iter_matches_scan_fasta():
    """Collecting scan_fasta_iter into a list must equal scan_fasta output."""
    path = fasta_to_tmp([("r1", GAMIR_SEQ), ("r2", CTMIR_SEQ)])
    try:
        iter_hits  = list(hseeker.scan_fasta_iter(str(path), minrep=6))
        fasta_hits = hseeker.scan_fasta(str(path), minrep=6)
        assert iter_hits == fasta_hits
    finally:
        path.unlink(missing_ok=True)


def test_scan_fasta_iter_seq_id_set_on_each_hit():
    """Every yielded hit must carry the correct seq_id from its FASTA record."""
    path = fasta_to_tmp([("seq_A", PURE_GA * 3), ("seq_B", CTMIR_SEQ)])
    try:
        hits = list(hseeker.scan_fasta_iter(str(path), minrep=6))
        assert len(hits) > 0
        ids = {h["seq_id"] for h in hits}
        assert "seq_A" in ids
        assert "seq_B" in ids
    finally:
        path.unlink(missing_ok=True)


def test_scan_fasta_iter_all_n_yields_nothing():
    """All-N FASTA record must yield zero hits."""
    path = fasta_to_tmp([("empty", "N" * 80)])
    try:
        hits = list(hseeker.scan_fasta_iter(str(path), minrep=6))
        assert hits == []
    finally:
        path.unlink(missing_ok=True)


def test_scan_fasta_iter_in_all():
    """scan_fasta_iter must be listed in hseeker.__all__."""
    assert "scan_fasta_iter" in hseeker.__all__


# ===========================================================================
# 6c. scan_fasta_parallel() — thread-parallel FASTA scan
# ===========================================================================

def test_scan_fasta_parallel_matches_scan_fasta():
    """scan_fasta_parallel must return the same hits as scan_fasta."""
    records = [("chr1", GAMIR_SEQ * 3), ("chr2", CTMIR_SEQ * 3)]
    path = fasta_to_tmp(records)
    try:
        expected = hseeker.scan_fasta(str(path), minrep=6)
        actual   = hseeker.scan_fasta_parallel(str(path), minrep=6)
        assert actual == expected
    finally:
        path.unlink(missing_ok=True)


def test_scan_fasta_parallel_preserves_record_order():
    """Hits from chr1 must appear before hits from chr2."""
    records = [("chr1", GAMIR_SEQ * 3), ("chr2", CTMIR_SEQ * 3)]
    path = fasta_to_tmp(records)
    try:
        hits = hseeker.scan_fasta_parallel(str(path), minrep=6)
        seq_ids = [h["seq_id"] for h in hits]
        chr1_indices = [i for i, s in enumerate(seq_ids) if s == "chr1"]
        chr2_indices = [i for i, s in enumerate(seq_ids) if s == "chr2"]
        if chr1_indices and chr2_indices:
            assert max(chr1_indices) < min(chr2_indices)
    finally:
        path.unlink(missing_ok=True)


def test_scan_fasta_parallel_explicit_workers():
    """workers=1 must produce the same result as the default."""
    records = [("seq1", PERFECT_SEQ * 5), ("seq2", PURE_GA)]
    path = fasta_to_tmp(records)
    try:
        default_hits  = hseeker.scan_fasta_parallel(str(path), minrep=6)
        workers1_hits = hseeker.scan_fasta_parallel(str(path), minrep=6, workers=1)
        assert workers1_hits == default_hits
    finally:
        path.unlink(missing_ok=True)


def test_scan_fasta_parallel_empty_records_yields_empty():
    """A FASTA file with only N-bases should produce no hits."""
    records = [("all_n", "N" * 50), ("also_n", "N" * 50)]
    path = fasta_to_tmp(records)
    try:
        hits = hseeker.scan_fasta_parallel(str(path), minrep=6)
        assert hits == []
    finally:
        path.unlink(missing_ok=True)


def test_scan_fasta_parallel_in_all():
    """scan_fasta_parallel must be listed in hseeker.__all__."""
    assert "scan_fasta_parallel" in hseeker.__all__


# ===========================================================================
# 7. Overlap removal
# ===========================================================================

def test_overlap_removal_reduces_hit_count():
    long_ga = "GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA"
    all_hits = hseeker.scan_sequence(long_ga, minrep=6, remove_overlaps=False)
    dedup    = hseeker.scan_sequence(long_ga, minrep=6, remove_overlaps=True)
    assert len(dedup) <= len(all_hits)


def test_no_overlapping_intervals_after_removal():
    long_ga = "GAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGAGA"
    hits = hseeker.scan_sequence(long_ga, minrep=6, remove_overlaps=True)
    for i in range(len(hits)):
        for j in range(i + 1, len(hits)):
            si, ei = hits[i]["start"], hits[i]["end"]
            sj, ej = hits[j]["start"], hits[j]["end"]
            assert not (si <= ej and sj <= ei), \
                f"Overlap not removed: [{si},{ei}] vs [{sj},{ej}]"


def test_skip_overlap_yields_more_or_equal_hits():
    seq = PURE_GA * 2
    with_removal    = hseeker.scan_sequence(seq, minrep=6, remove_overlaps=True)
    without_removal = hseeker.scan_sequence(seq, minrep=6, remove_overlaps=False)
    assert len(without_removal) >= len(with_removal)


@pytest.mark.parametrize("spacer,expected_spacer", [("AA", 0), ("AAA", 1)])
def test_matching_spacer_pairs_extend_arms_inward(spacer, expected_spacer):
    """Even and odd spacers absorb matching pairs from both ends."""
    seq = "A" * 6 + spacer + "A" * 6
    hits = hseeker.scan_sequence(
        seq, minrep=6, maxrep=7, maxspacer=3, purity=1.0,
        mismatch=0.0, remove_overlaps=False, score=False,
    )
    same_span = [h for h in hits if (h["start"], h["end"]) == (1, len(seq))]
    assert [(h["arm_length"], h["spacer_length"])
            for h in same_span] == [(7, expected_spacer)]


def test_inward_mismatch_preserves_both_representations():
    """A nonmatching spacer-end pair is not absorbed into the arms."""
    seq = "A" * 6 + "AG" + "A" * 6
    hits = hseeker.scan_sequence(
        seq, minrep=6, maxrep=7, maxspacer=2, purity=1.0,
        mismatch=0.2, remove_overlaps=False, score=False,
    )
    same_span = [h for h in hits if (h["start"], h["end"]) == (1, len(seq))]
    assert {(h["arm_length"], h["spacer_length"], round(h["mirror_identity"]))
            for h in same_span} == {(6, 2, 100), (7, 0, 86)}


def test_inward_extension_respects_purity_and_maxrep():
    """Matching spacer bases must still meet the existing arm limits."""
    seq = "G" * 6 + "CC" + "G" * 6
    hits = hseeker.scan_sequence(
        seq, minrep=6, maxrep=7, maxspacer=2, purity=1.0,
        mismatch=0.0, remove_overlaps=False, score=False,
    )
    same_span = [h for h in hits if (h["start"], h["end"]) == (1, len(seq))]
    assert [(h["arm_length"], h["spacer_length"]) for h in same_span] == [(6, 2)]

    capped = hseeker.scan_sequence(
        "A" * 14, minrep=6, maxrep=6, maxspacer=2, purity=1.0,
        mismatch=0.0, remove_overlaps=False, score=False,
    )
    same_span = [h for h in capped if (h["start"], h["end"]) == (1, 14)]
    assert [(h["arm_length"], h["spacer_length"]) for h in same_span] == [(6, 2)]


def test_inward_extension_without_preexisting_same_span_candidate():
    """Normalize a hit even if the inner-center scan grew farther outward."""
    seq = "GAAAAGGGAAGGGAAGGGAAGGGAAGGGAAGGGATGGGAAGAGA"
    hits = hseeker.scan_sequence(
        seq, minrep=10, maxspacer=10, purity=0.9, mismatch=0.1,
        remove_overlaps=False, seq_offset=112, score=False,
    )
    same_span = [h for h in hits if (h["start"], h["end"]) == (115, 151)]
    assert [(h["arm_length"], h["spacer_length"]) for h in same_span] == [(18, 1)]


# ===========================================================================
# 8. Parameter validation
# ===========================================================================

@pytest.mark.parametrize("kwargs", [
    {"minrep": 0},
    {"maxrep": 3, "minrep": 6},
    {"maxspacer": -1},
    {"purity": 1.5},
    {"purity": -0.1},
    {"mismatch": -0.1},
    {"mismatch": 1.1},
])
def test_bad_parameters_raise_value_error(kwargs):
    with pytest.raises(ValueError):
        hseeker.scan_sequence(PURE_GA, **kwargs)


# ===========================================================================
# 9. Known-sequence exact field values  (GAMIR_SEQ = "GGGAAATTAAAGGG")
#
# The center-outward algorithm keeps the longest valid arm at each (ctr, sp).
# After overlap removal the surviving hit has:
#   arm=7, sp=0  (direct adjacency: "gggaaat" | "taaaggg")
# rather than the shorter arm=6, sp=2 decomposition.
#
# Derived values (verified via diagnostic run):
#   left_arm  = "gggaaat"  →  right_arm = "taaaggg"  (= reverse left arm ✓)
#   ga_pct    = 6/7 × 100 ≈ 85.71%   (right arm: 1 T + 3 A + 3 G = 6 GA)
#   ct_pct    = 1/7 × 100 ≈ 14.29%
#   mirror_id = 100.0%
#   is_perfect = False  (requires ga_pct = 100% OR ct_pct = 100%)
# ===========================================================================

def _gamir_hit():
    """Return the surviving hit from GAMIR_SEQ (arm=7 after overlap removal).

    Uses purity=0.80 / mismatch=0.20 (the pre-v0.2.0 defaults) because
    GAMIR_SEQ's best arm has 85.71 % GA purity — detected by the relaxed
    threshold but not by the current default of 0.90.
    """
    hits = hseeker.scan_sequence(GAMIR_SEQ, minrep=6, purity=0.80, mismatch=0.20)
    arm7 = [h for h in hits if h["arm_length"] == 7]
    assert arm7, f"Expected arm=7 hit on GAMIR_SEQ, got: {hits}"
    return arm7[0]


def test_gamir_hit_exists():
    assert len(hseeker.scan_sequence(GAMIR_SEQ, minrep=6, purity=0.80, mismatch=0.20)) > 0


def test_gamir_arm_length():
    assert _gamir_hit()["arm_length"] == 7


def test_gamir_spacer_length():
    assert _gamir_hit()["spacer_length"] == 0


def test_gamir_spacer_sentinel():
    """Zero-spacer hit must use '.' as the spacer string."""
    assert _gamir_hit()["spacer"] == "."


def test_gamir_start_end_coordinates():
    h = _gamir_hit()
    assert h["start"] == 1
    assert h["end"]   == 14


def test_gamir_left_arm_sequence():
    assert _gamir_hit()["left_arm"] == "gggaaat"


def test_gamir_right_arm_sequence():
    assert _gamir_hit()["right_arm"] == "taaaggg"


def test_gamir_ga_pct():
    assert _gamir_hit()["ga_pct"] == pytest.approx(100.0 * 6 / 7, abs=0.1)


def test_gamir_ct_pct():
    assert _gamir_hit()["ct_pct"] == pytest.approx(100.0 * 1 / 7, abs=0.1)


def test_gamir_mirror_identity_is_100():
    assert _gamir_hit()["mirror_identity"] == pytest.approx(100.0, abs=0.1)


def test_gamir_is_perfect_false():
    """ga_pct ≈ 85.7% — not 100%, so is_perfect must be False."""
    assert _gamir_hit()["is_perfect"] is False


def test_gamir_total_length_formula():
    h = _gamir_hit()
    assert h["total_length"] == h["arm_length"] * 2 + h["spacer_length"]


def test_gamir_full_sequence_composition():
    h = _gamir_hit()
    spacer_str = h["spacer"] if h["spacer"] != "." else ""
    assert h["full_sequence"] == h["left_arm"] + spacer_str + h["right_arm"]


# ===========================================================================
# 10. CT arm detection  (CTMIR_SEQ = "CCCTTTAATTTCCC")
#
# Best hit (verified): arm=7, sp=0
#   left_arm = "cccttta",  right_arm = "atttccc"
#   ct_pct   = 6/7 × 100 ≈ 85.71%
#   ga_pct   = 1/7 × 100 ≈ 14.29%
#   is_perfect = False  (ct_pct ≠ 100%)
# ===========================================================================

def _ctmir_hit():
    """arm=7 hit on CTMIR_SEQ (ct_pct=85.71% — needs purity=0.80)."""
    hits = hseeker.scan_sequence(CTMIR_SEQ, minrep=6, purity=0.80, mismatch=0.20)
    arm7 = [h for h in hits if h["arm_length"] == 7]
    assert arm7, f"Expected arm=7 hit on CTMIR_SEQ, got: {hits}"
    return arm7[0]


def test_ct_arm_detected():
    assert len(hseeker.scan_sequence(CTMIR_SEQ, minrep=6, purity=0.80, mismatch=0.20)) > 0


def test_ct_arm_ct_pct():
    assert _ctmir_hit()["ct_pct"] == pytest.approx(100.0 * 6 / 7, abs=0.1)


def test_ct_arm_ga_pct():
    assert _ctmir_hit()["ga_pct"] == pytest.approx(100.0 * 1 / 7, abs=0.1)


def test_ct_arm_is_perfect_false():
    assert _ctmir_hit()["is_perfect"] is False


def test_ct_arm_ct_pct_exceeds_ga_pct():
    """CT motif must have ct_pct > ga_pct."""
    h = _ctmir_hit()
    assert h["ct_pct"] > h["ga_pct"]


def test_ct_arm_arm_sequences():
    h = _ctmir_hit()
    assert h["left_arm"]  == "cccttta"
    assert h["right_arm"] == "atttccc"


def test_ct_arm_coordinates():
    h = _ctmir_hit()
    assert h["start"] == 1
    assert h["end"]   == 14


# ===========================================================================
# 11. Purity filter boundary
# ===========================================================================

def test_low_purity_excluded_with_default_threshold():
    """LOWPUR_SEQ right arm is 50% GA/CT — below default purity=0.80."""
    hits = hseeker.scan_sequence(LOWPUR_SEQ, minrep=6, purity=0.80, mismatch=0.0)
    assert hits == [], "50% GA/CT arm should be excluded at purity=0.80"


def test_low_purity_included_with_relaxed_threshold():
    hits = hseeker.scan_sequence(LOWPUR_SEQ, minrep=6, purity=0.50, mismatch=0.0)
    assert len(hits) > 0, "50% GA/CT arm should be included at purity=0.50"


def test_purity_boundary_exact_50pct():
    """purity=0.50 must find it (>= comparison in C); purity=0.51 must not."""
    found     = hseeker.scan_sequence(LOWPUR_SEQ, minrep=6, purity=0.50, mismatch=0.0)
    not_found = hseeker.scan_sequence(LOWPUR_SEQ, minrep=6, purity=0.51, mismatch=0.0)
    assert len(found) > 0
    assert len(not_found) == 0


def test_purity_90_excludes_impure_arm():
    """IMPURE_SEQ has arm ga_pct≈83.3% → not found at purity=0.90."""
    assert hseeker.scan_sequence(IMPURE_SEQ, minrep=6, purity=0.90, mismatch=0.0) == []


def test_purity_80_finds_impure_arm():
    """IMPURE_SEQ arm ga_pct≈83.3% ≥ 80% → found at purity=0.80."""
    assert len(hseeker.scan_sequence(IMPURE_SEQ, minrep=6, purity=0.80, mismatch=0.0)) > 0


# ===========================================================================
# 12. maxspacer boundary
#
# SPACER8_SEQ = "GGGGGG" + "T"*8 + "GGGGGG"  (20 bp)
# With mismatch=0.0: arm=6, sp=8 is the only perfect hit.
# With maxspacer=7 this (ctr, sp=8) is never explored → no hit.
# With maxspacer=8 it is explored → hit found.
# ===========================================================================

def test_spacer8_not_found_with_maxspacer7_strict():
    """The only valid perfect arm in SPACER8_SEQ has sp=8.
    With maxspacer=7 that (ctr, sp=8) pair is never explored."""
    hits = hseeker.scan_sequence(SPACER8_SEQ, minrep=6, maxspacer=7, mismatch=0.0)
    assert hits == [], f"No hit should be found when maxspacer=7 (requires sp=8), got: {hits}"


def test_spacer8_found_with_maxspacer8_strict():
    hits = hseeker.scan_sequence(SPACER8_SEQ, minrep=6, maxspacer=8, mismatch=0.0)
    sp8  = [h for h in hits if h["spacer_length"] == 8 and h["arm_length"] == 6]
    assert len(sp8) > 0, f"arm=6, sp=8 hit must be found when maxspacer=8, got: {hits}"


def test_spacer8_hit_is_perfect():
    hits = hseeker.scan_sequence(SPACER8_SEQ, minrep=6, maxspacer=8, mismatch=0.0)
    sp8  = [h for h in hits if h["spacer_length"] == 8]
    assert sp8[0]["mirror_identity"] == pytest.approx(100.0, abs=0.01)
    assert sp8[0]["is_perfect"] is True


def test_spacer_zero_uses_dot_sentinel():
    """A spacer_length=0 hit must have spacer='.' in the output dict."""
    hits = hseeker.scan_sequence("A" * 20, minrep=6, maxspacer=0)
    sp0 = [h for h in hits if h["spacer_length"] == 0]
    assert len(sp0) > 0
    for h in sp0:
        assert h["spacer"] == ".", f"Expected spacer='.', got {h['spacer']!r}"


def test_maxspacer_0_finds_only_direct_adjacency():
    """maxspacer=0 → all returned hits must have spacer_length=0."""
    hits = hseeker.scan_sequence(PURE_GA * 2, minrep=6, maxspacer=0)
    for h in hits:
        assert h["spacer_length"] == 0


# ===========================================================================
# 13. Case insensitivity
# ===========================================================================

def test_uppercase_and_lowercase_give_same_hit_count():
    upper = hseeker.scan_sequence(GAMIR_SEQ.upper(), minrep=6)
    lower = hseeker.scan_sequence(GAMIR_SEQ.lower(), minrep=6)
    assert len(upper) == len(lower)


def test_mixed_case_same_as_lower():
    mixed = "gGgAaAtTaAaGgG"
    assert len(hseeker.scan_sequence(mixed, minrep=6)) == \
           len(hseeker.scan_sequence(mixed.lower(), minrep=6))


def test_output_sequences_are_always_lowercase():
    for h in hseeker.scan_sequence(GAMIR_SEQ.upper(), minrep=6):
        assert h["left_arm"]      == h["left_arm"].lower()
        assert h["right_arm"]     == h["right_arm"].lower()
        assert h["full_sequence"] == h["full_sequence"].lower()


# ===========================================================================
# 14. Output field self-consistency
# ===========================================================================

def test_total_length_formula_all_sequences():
    for seq in (GAMIR_SEQ, CTMIR_SEQ, PURE_GA, PURE_GA * 2):
        for h in hseeker.scan_sequence(seq, minrep=6):
            assert h["total_length"] == h["arm_length"] * 2 + h["spacer_length"]


def test_full_sequence_equals_left_spacer_right():
    for seq in (GAMIR_SEQ, CTMIR_SEQ, PURE_GA):
        for h in hseeker.scan_sequence(seq, minrep=6):
            spacer_str = h["spacer"] if h["spacer"] != "." else ""
            expected = h["left_arm"] + spacer_str + h["right_arm"]
            assert h["full_sequence"] == expected


def test_coordinate_span_equals_total_length():
    """end - start + 1 must equal total_length (1-based inclusive)."""
    for h in hseeker.scan_sequence(GAMIR_SEQ, minrep=6):
        assert h["end"] - h["start"] + 1 == h["total_length"]


def test_arm_string_lengths_match_arm_length():
    for h in hseeker.scan_sequence(PURE_GA, minrep=6):
        assert len(h["left_arm"])  == h["arm_length"]
        assert len(h["right_arm"]) == h["arm_length"]


def test_spacer_string_length_matches_spacer_length():
    for h in hseeker.scan_sequence(GAMIR_SEQ, minrep=6):
        if h["spacer_length"] == 0:
            assert h["spacer"] == "."
        else:
            assert len(h["spacer"]) == h["spacer_length"]


def test_pct_values_sum_at_most_100():
    """ga_pct + ct_pct <= 100 (no double-counting)."""
    for h in hseeker.scan_sequence(GAMIR_SEQ + CTMIR_SEQ, minrep=6):
        assert h["ga_pct"] + h["ct_pct"] <= 100.0 + 1e-4


def test_all_hits_are_dicts():
    for h in hseeker.scan_sequence(PURE_GA, minrep=6):
        assert isinstance(h, dict)


def test_scan_sequence_returns_list():
    assert isinstance(hseeker.scan_sequence(PURE_GA, minrep=6), list)


# ===========================================================================
# 15. Flanking N bases
#
# Embed GAMIR_SEQ at 0-based position N_PAD inside an N-padded sequence.
# N stops arm extension at the boundary.
# The surviving hit is arm=7, sp=0 (same as unpadded GAMIR_SEQ).
# ===========================================================================

N_PAD   = 30
FLANKED = "N" * N_PAD + GAMIR_SEQ + "N" * 30


def _flanked_hit():
    hits = hseeker.scan_sequence(FLANKED, minrep=6, purity=0.80, mismatch=0.20)
    arm7 = [h for h in hits if h["arm_length"] == 7]
    assert arm7, "Expected arm=7 hit inside N-flanked sequence"
    return arm7[0]


def test_flanking_n_does_not_prevent_detection():
    assert len(_flanked_hit()["left_arm"]) > 0


def test_flanking_n_correct_start_coordinate():
    """left_arm starts at 1-based position N_PAD + 1."""
    assert _flanked_hit()["start"] == N_PAD + 1


def test_flanking_n_correct_end_coordinate():
    """Right arm ends at N_PAD + len(GAMIR_SEQ) (= 44 for N_PAD=30)."""
    assert _flanked_hit()["end"] == N_PAD + len(GAMIR_SEQ)


def test_flanking_n_arm_sequences_identical_to_unflanked():
    base_arm7    = [h for h in hseeker.scan_sequence(GAMIR_SEQ, minrep=6, purity=0.80, mismatch=0.20) if h["arm_length"] == 7]
    flanked_arm7 = [h for h in hseeker.scan_sequence(FLANKED, minrep=6, purity=0.80, mismatch=0.20)   if h["arm_length"] == 7]
    assert base_arm7[0]["left_arm"]  == flanked_arm7[0]["left_arm"]
    assert base_arm7[0]["right_arm"] == flanked_arm7[0]["right_arm"]
    assert base_arm7[0]["spacer"]    == flanked_arm7[0]["spacer"]


# ===========================================================================
# 16. Multiple non-overlapping motifs in one sequence
# ===========================================================================

DUAL_SEQ = GAMIR_SEQ + "N" * 20 + CTMIR_SEQ


def test_two_motifs_both_detected():
    hits = hseeker.scan_sequence(DUAL_SEQ, minrep=6, purity=0.80, mismatch=0.20)
    assert len(hits) >= 2


def test_two_motifs_have_distinct_purity_profiles():
    """GA mirror has ga_pct > ct_pct; CT mirror has ct_pct > ga_pct."""
    arm7 = [h for h in hseeker.scan_sequence(DUAL_SEQ, minrep=6, purity=0.80, mismatch=0.20) if h["arm_length"] == 7]
    ga_dominant = [h for h in arm7 if h["ga_pct"] > h["ct_pct"]]
    ct_dominant = [h for h in arm7 if h["ct_pct"] > h["ga_pct"]]
    assert len(ga_dominant) >= 1, "Expected a GA-dominant hit"
    assert len(ct_dominant) >= 1, "Expected a CT-dominant hit"


def test_two_motifs_no_overlap_after_removal():
    hits = hseeker.scan_sequence(DUAL_SEQ, minrep=6, remove_overlaps=True)
    for i in range(len(hits)):
        for j in range(i + 1, len(hits)):
            si, ei = hits[i]["start"], hits[i]["end"]
            sj, ej = hits[j]["start"], hits[j]["end"]
            assert not (si <= ej and sj <= ei), \
                f"Unexpected overlap: [{si},{ei}] vs [{sj},{ej}]"


# ===========================================================================
# 17. is_perfect flag logic
#
# C definition:  is_perfect = (ga_pct==100% OR ct_pct==100%) AND mirror_id==100%
# ===========================================================================

def test_is_perfect_true_for_pure_sequence():
    """PERFECT_SEQ = 'A'*14: all-A arms → ga_pct=100%, mirror_id=100% → is_perfect=True."""
    hits = hseeker.scan_sequence(PERFECT_SEQ, minrep=6)
    assert any(h["is_perfect"] for h in hits), "All-A sequence must yield is_perfect=True"


def test_is_perfect_false_when_mirror_is_imperfect():
    """IMPMIR_SEQ: 1-mismatch arm → is_perfect=False."""
    hits = hseeker.scan_sequence(IMPMIR_SEQ, minrep=6, purity=0.80, mismatch=0.20)
    assert len(hits) > 0
    assert all(not h["is_perfect"] for h in hits)


def test_is_perfect_false_when_arm_not_100pct_pure():
    """IMPURE_SEQ: perfect mirror but ga_pct≈83.3% → is_perfect=False."""
    hits = hseeker.scan_sequence(IMPURE_SEQ, minrep=6, purity=0.80, mismatch=0.0)
    assert len(hits) > 0
    assert all(not h["is_perfect"] for h in hits)


def test_is_perfect_is_python_bool():
    for h in hseeker.scan_sequence(GAMIR_SEQ, minrep=6):
        assert isinstance(h["is_perfect"], bool)


def test_is_perfect_only_true_when_purity_and_mirror_are_100():
    """Every hit with is_perfect=True must satisfy the C formula."""
    for h in hseeker.scan_sequence(PURE_GA * 2, minrep=6, purity=0.80, mismatch=0.20):
        if h["is_perfect"]:
            pure_arm = (h["ga_pct"] >= 100.0 - 0.1 or h["ct_pct"] >= 100.0 - 0.1)
            exact_mir = h["mirror_identity"] >= 100.0 - 0.1
            assert pure_arm and exact_mir, \
                f"is_perfect=True but ga={h['ga_pct']:.2f} ct={h['ct_pct']:.2f} mir={h['mirror_identity']:.2f}"


# ===========================================================================
# 18. Mirror identity with known mismatch count
#
# IMPMIR_SEQ = "AAAAAACAAAAA"  (12 bp)
# At (ctr=5, sp=0): k=6 extension, 1 mismatch (dna[5]='a' vs dna[6]='c').
# mirror_id = (6-1)/6 × 100 = 5/6 × 100 ≈ 83.33%
# ===========================================================================

def test_impmir_mirror_identity_value():
    hits = hseeker.scan_sequence(IMPMIR_SEQ, minrep=6, purity=0.80, mismatch=0.20)
    arm6 = [h for h in hits if h["arm_length"] == 6]
    assert arm6, "Expected arm=6 hit on IMPMIR_SEQ"
    assert arm6[0]["mirror_identity"] == pytest.approx(100.0 * 5 / 6, abs=0.2)


def test_impmir_not_found_with_zero_mismatch():
    """With mismatch=0.0, the imperfect arm should not appear."""
    hits = hseeker.scan_sequence(IMPMIR_SEQ, minrep=6, purity=0.80, mismatch=0.0)
    assert hits == [], "Imperfect mirror must not be found with mismatch_tol=0.0"


def test_impmir_is_perfect_false():
    hits = hseeker.scan_sequence(IMPMIR_SEQ, minrep=6, purity=0.80, mismatch=0.20)
    assert all(not h["is_perfect"] for h in hits)


def test_mirror_identity_100_for_gamir():
    for h in hseeker.scan_sequence(GAMIR_SEQ, minrep=6):
        assert h["mirror_identity"] == pytest.approx(100.0, abs=0.01)


def test_mirror_identity_in_range_0_to_100():
    for h in hseeker.scan_sequence(PURE_GA * 2, minrep=6, purity=0.70, mismatch=0.30):
        assert 0.0 <= h["mirror_identity"] <= 100.0 + 1e-4


# ===========================================================================
# 19. minrep / maxrep arm-length boundaries
# ===========================================================================

def test_arm_length_always_at_least_minrep():
    for minrep in (6, 8, 10):
        for h in hseeker.scan_sequence(PURE_GA * 3, minrep=minrep):
            assert h["arm_length"] >= minrep


def test_arm_length_never_exceeds_maxrep():
    for maxrep in (8, 12, 20):
        for h in hseeker.scan_sequence(PURE_GA * 3, minrep=6, maxrep=maxrep):
            assert h["arm_length"] <= maxrep


def test_large_minrep_gives_no_hits():
    """minrep=8 on 14-bp GAMIR_SEQ: no room for arm(8)+sp+arm(8) → no hits."""
    assert hseeker.scan_sequence(GAMIR_SEQ, minrep=8) == []


def test_minrep_exactly_6_detects_gamir():
    """With minrep=6, at least one hit must exist on GAMIR_SEQ."""
    assert len(hseeker.scan_sequence(GAMIR_SEQ, minrep=6)) > 0


# ===========================================================================
# 20. Genomic coordinate propagation via FASTA offset header
# ===========================================================================

def test_embedded_motif_has_correct_genomic_coordinates():
    """Embed GAMIR_SEQ at 0-based position 99 inside a 200-bp genomic window."""
    pad     = 99
    seq     = "N" * pad + GAMIR_SEQ + "N" * 87   # total = 200 bp
    g_start = 5_000_000

    text = f">chr3:{g_start}-{g_start + 199}\n{seq}\n"
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    tmp.write(text)
    tmp.close()
    path = Path(tmp.name)
    try:
        hits = hseeker.scan_fasta(str(path), minrep=6, purity=0.80, mismatch=0.20)
        arm7 = [h for h in hits if h["arm_length"] == 7]
        assert arm7, "Expected arm=7 hit for GAMIR_SEQ embedded at offset 99"
        h = arm7[0]
        # start = (pad + 1) + (g_start - 1)
        expected_start = (pad + 1) + (g_start - 1)
        expected_end   = expected_start + len(GAMIR_SEQ) - 1
        assert h["start"] == expected_start, \
            f"Expected start={expected_start}, got {h['start']}"
        assert h["end"] == expected_end, \
            f"Expected end={expected_end}, got {h['end']}"
    finally:
        path.unlink(missing_ok=True)


# ===========================================================================
# 21. parse_fasta edge cases
# ===========================================================================

def test_parse_fasta_multiline_sequence_assembled_correctly():
    text = ">wrapped\nGGGAAA\nTTAAAG\nGG\n"
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    tmp.write(text)
    tmp.close()
    path = Path(tmp.name)
    try:
        _, seq, _ = list(hseeker.parse_fasta(str(path)))[0]
        assert seq == "GGGAAATTAAAGGG"
    finally:
        path.unlink(missing_ok=True)


def test_parse_fasta_header_description_does_not_leak():
    text = ">myseq some description text\nGAGAGAGAGAGAGA\n"
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    tmp.write(text)
    tmp.close()
    path = Path(tmp.name)
    try:
        seq_id, _, _ = list(hseeker.parse_fasta(str(path)))[0]
        assert seq_id == "myseq", f"Got seq_id={seq_id!r}"
    finally:
        path.unlink(missing_ok=True)


def test_parse_fasta_default_offset_is_1():
    text = ">plain_id\nGAGAGAGAGAGAGA\n"
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    tmp.write(text)
    tmp.close()
    path = Path(tmp.name)
    try:
        _, _, offset = list(hseeker.parse_fasta(str(path)))[0]
        assert offset == 1
    finally:
        path.unlink(missing_ok=True)


def test_parse_fasta_three_records_count_and_order():
    records = [("r1", "A" * 20), ("r2", GAMIR_SEQ), ("r3", "C" * 20)]
    path = fasta_to_tmp(records)
    try:
        parsed = list(hseeker.parse_fasta(str(path)))
        assert len(parsed) == 3
        assert [seq for _, seq, _ in parsed] == ["A" * 20, GAMIR_SEQ, "C" * 20]
    finally:
        path.unlink(missing_ok=True)


def test_parse_fasta_returns_uppercase_sequence():
    text = ">seq\ngggaaattaaaggg\n"
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    tmp.write(text)
    tmp.close()
    path = Path(tmp.name)
    try:
        _, seq, _ = list(hseeker.parse_fasta(str(path)))[0]
        assert seq == "GGGAAATTAAAGGG"
    finally:
        path.unlink(missing_ok=True)


def test_parse_fasta_single_line_sequence():
    text = write_fasta([("nowrap", GAMIR_SEQ)], line_width=0)
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    tmp.write(text)
    tmp.close()
    path = Path(tmp.name)
    try:
        _, seq, _ = list(hseeker.parse_fasta(str(path)))[0]
        assert seq == GAMIR_SEQ
    finally:
        path.unlink(missing_ok=True)


def test_parse_fasta_coordinate_offset_numeric_value():
    """chr1:43585222-... → offset=43585222 as an integer."""
    text = ">chr1:43585222-43586222\nAAAAAAAAAAAAAAAAAAAA\n"
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    tmp.write(text)
    tmp.close()
    path = Path(tmp.name)
    try:
        _, _, offset = list(hseeker.parse_fasta(str(path)))[0]
        assert isinstance(offset, int)
        assert offset == 43585222
    finally:
        path.unlink(missing_ok=True)


# ===========================================================================
# 22. Determinism / reproducibility
# ===========================================================================

def test_scan_sequence_is_deterministic():
    results = [hseeker.scan_sequence(PURE_GA * 2, minrep=6) for _ in range(5)]
    for r in results[1:]:
        assert r == results[0], "scan_sequence is not deterministic"


def test_scan_fasta_is_deterministic():
    path = fasta_to_tmp([("s", GAMIR_SEQ)])
    try:
        results = [hseeker.scan_fasta(str(path), minrep=6) for _ in range(3)]
        for r in results[1:]:
            assert r == results[0]
    finally:
        path.unlink(missing_ok=True)


# ===========================================================================
# 23. CLI integration  (python -m hseeker)
# ===========================================================================

def test_cli_zero_exit_code_and_tsv_created():
    path = fasta_to_tmp([("cli_test", GAMIR_SEQ)])
    out_prefix = str(path.parent / "cli_out")
    tsv_path = Path(out_prefix + "_HDNA.tsv")
    try:
        result = subprocess.run(
            [sys.executable, "-m", "hseeker",
             "-seq", str(path), "-out", out_prefix, "-minrep", "6"],
            capture_output=True, text=True,
        )
        assert result.returncode == 0, f"CLI exited non-zero:\n{result.stderr}"
        assert tsv_path.exists(), "Output TSV was not created"
    finally:
        path.unlink(missing_ok=True)
        tsv_path.unlink(missing_ok=True)


def test_cli_tsv_has_required_columns():
    path = fasta_to_tmp([("cli_test", GAMIR_SEQ)])
    out_prefix = str(path.parent / "cli_cols")
    tsv_path = Path(out_prefix + "_HDNA.tsv")
    try:
        subprocess.run(
            [sys.executable, "-m", "hseeker",
             "-seq", str(path), "-out", out_prefix, "-minrep", "6"],
            check=True, capture_output=True,
        )
        lines = tsv_path.read_text().splitlines()
        assert len(lines) >= 2, "TSV must have header + ≥1 data row"
        header = lines[0].split("\t")
        for col in ("seq_id", "source", "start", "end",
                    "arm_length", "mirror_identity", "is_perfect",
                    "left_arm", "right_arm"):
            assert col in header, f"Column missing from TSV: {col}"
    finally:
        path.unlink(missing_ok=True)
        tsv_path.unlink(missing_ok=True)


def test_cli_output_row_count_matches_python_api():
    path = fasta_to_tmp([("s1", GAMIR_SEQ), ("s2", CTMIR_SEQ)])
    out_prefix = str(path.parent / "cli_match")
    tsv_path = Path(out_prefix + "_HDNA.tsv")
    try:
        subprocess.run(
            [sys.executable, "-m", "hseeker",
             "-seq", str(path), "-out", out_prefix, "-minrep", "6"],
            check=True, capture_output=True,
        )
        with open(tsv_path) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        api_hits = hseeker.scan_fasta(str(path), minrep=6)
        assert len(rows) == len(api_hits), \
            f"CLI: {len(rows)} rows; API: {len(api_hits)} hits"
    finally:
        path.unlink(missing_ok=True)
        tsv_path.unlink(missing_ok=True)


def test_cli_skipoverlap_produces_more_or_equal_hits():
    path = fasta_to_tmp([("s", PURE_GA * 2)])
    out_a = str(path.parent / "cli_a")
    out_b = str(path.parent / "cli_b")
    try:
        subprocess.run(
            [sys.executable, "-m", "hseeker",
             "-seq", str(path), "-out", out_a, "-minrep", "6"],
            check=True, capture_output=True,
        )
        subprocess.run(
            [sys.executable, "-m", "hseeker",
             "-seq", str(path), "-out", out_b, "-minrep", "6", "-skipoverlap"],
            check=True, capture_output=True,
        )
        n_a = len(Path(out_a + "_HDNA.tsv").read_text().splitlines()) - 1
        n_b = len(Path(out_b + "_HDNA.tsv").read_text().splitlines()) - 1
        assert n_b >= n_a, f"-skipoverlap ({n_b}) should give >= hits vs default ({n_a})"
    finally:
        path.unlink(missing_ok=True)
        for p in (out_a, out_b):
            Path(p + "_HDNA.tsv").unlink(missing_ok=True)


def test_cli_source_column_is_findhdna():
    path = fasta_to_tmp([("s", GAMIR_SEQ)])
    out_prefix = str(path.parent / "cli_src")
    tsv_path = Path(out_prefix + "_HDNA.tsv")
    try:
        subprocess.run(
            [sys.executable, "-m", "hseeker",
             "-seq", str(path), "-out", out_prefix, "-minrep", "6"],
            check=True, capture_output=True,
        )
        with open(tsv_path) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        for row in rows:
            assert row["source"] == "findHDNA", \
                f"Expected source='findHDNA', got {row['source']!r}"
    finally:
        path.unlink(missing_ok=True)
        tsv_path.unlink(missing_ok=True)


# ===========================================================================
# 24. Parallel chunking edge cases  (scan_fasta_parallel)
#
# Tests that chunking, boundary overlap handling, and multi-worker
# dispatch produce results identical to the serial scan_fasta path.
# ===========================================================================

_PCHUNK_PARAMS = dict(minrep=10, maxrep=100, maxspacer=10, purity=0.90, mismatch=0.10)


def _tmp_fasta_parallel(*records):
    """(seq_id, seq) tuples → temp Path (caller must unlink)."""
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    for sid, seq in records:
        tmp.write(f">{sid}\n{seq}\n")
    tmp.close()
    return Path(tmp.name)


def _pkey(h):
    return (h["seq_id"], h["start"], h["end"], h["arm_length"])


def test_parallel_single_chunk_matches_serial():
    """Record < chunk_size → treated as single chunk; must match serial."""
    seq = "A" * 50 + "G" * 50 + "C" * 50 + "T" * 50 + "A" * 200 + "G" * 200
    p = _tmp_fasta_parallel(("test", seq))
    try:
        ser = hseeker.scan_fasta(str(p), **_PCHUNK_PARAMS)
        par = hseeker.scan_fasta_parallel(str(p), chunk_size=5000, workers=2, **_PCHUNK_PARAMS)
        assert sorted(_pkey(h) for h in ser) == sorted(_pkey(h) for h in par)
    finally:
        p.unlink(missing_ok=True)


def test_parallel_hit_across_chunk_boundary():
    """Pure poly-A tract that crosses the chunk boundary must still be found."""
    chunk_size = 500
    seq = "A" * chunk_size + "A" * 4000
    p = _tmp_fasta_parallel(("test", seq))
    try:
        ser = hseeker.scan_fasta(str(p), **_PCHUNK_PARAMS)
        par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=2, **_PCHUNK_PARAMS)
        ser_keys = sorted(_pkey(h) for h in ser)
        par_keys = sorted(_pkey(h) for h in par)
        assert ser_keys == par_keys, (
            f"Serial {len(ser)} hits vs parallel {len(par)} hits\n"
            f"Missing from par: {set(ser_keys)-set(par_keys)}\n"
            f"Extra in par: {set(par_keys)-set(ser_keys)}"
        )
    finally:
        p.unlink(missing_ok=True)


def test_parallel_record_exactly_chunk_size():
    """Record exactly chunk_size bases → single chunk; must match serial."""
    chunk_size = 1000
    seq = "G" * chunk_size
    p = _tmp_fasta_parallel(("exact", seq))
    try:
        ser = hseeker.scan_fasta(str(p), **_PCHUNK_PARAMS)
        par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=2, **_PCHUNK_PARAMS)
        assert sorted(_pkey(h) for h in ser) == sorted(_pkey(h) for h in par)
    finally:
        p.unlink(missing_ok=True)


def test_parallel_record_one_over_chunk():
    """Record chunk_size+1 bases → 2 chunks; must match serial."""
    chunk_size = 1000
    seq = "A" * (chunk_size + 1)
    p = _tmp_fasta_parallel(("one_over", seq))
    try:
        ser = hseeker.scan_fasta(str(p), **_PCHUNK_PARAMS)
        par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=2, **_PCHUNK_PARAMS)
        assert sorted(_pkey(h) for h in ser) == sorted(_pkey(h) for h in par)
    finally:
        p.unlink(missing_ok=True)


def test_parallel_empty_record():
    p = _tmp_fasta_parallel(("empty", ""))
    try:
        assert hseeker.scan_fasta_parallel(str(p), chunk_size=100, workers=2, **_PCHUNK_PARAMS) == []
    finally:
        p.unlink(missing_ok=True)


def test_parallel_too_short():
    p = _tmp_fasta_parallel(("short", "AAAAA"))
    try:
        assert hseeker.scan_fasta_parallel(str(p), chunk_size=10, workers=2, **_PCHUNK_PARAMS) == []
    finally:
        p.unlink(missing_ok=True)


def test_parallel_multi_record_mixed():
    """Multi-record FASTA with chunked, non-chunked, empty, and too-short records."""
    chunk_size = 500
    records = [
        ("chr1", "A" * 2000),
        ("chr2", "G" * 200),
        ("chr3", "C" * chunk_size),
        ("chr4", ""),
        ("chr5", "T" * 5),
        ("chr6", "A" * (chunk_size + 1)),
    ]
    p = _tmp_fasta_parallel(*records)
    try:
        ser = hseeker.scan_fasta(str(p), **_PCHUNK_PARAMS)
        par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=3, **_PCHUNK_PARAMS)
        assert sorted(_pkey(h) for h in ser) == sorted(_pkey(h) for h in par)
    finally:
        p.unlink(missing_ok=True)


@pytest.mark.parametrize("workers", [1, 2, 4])
def test_parallel_worker_independence(workers):
    """Results must be identical regardless of worker count."""
    chunk_size = 500
    seq = "A" * 2000 + "C" * 500 + "G" * 500
    p = _tmp_fasta_parallel(("wtest", seq))
    try:
        ser = hseeker.scan_fasta(str(p), **_PCHUNK_PARAMS)
        par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=workers, **_PCHUNK_PARAMS)
        assert sorted(_pkey(h) for h in ser) == sorted(_pkey(h) for h in par)
    finally:
        p.unlink(missing_ok=True)


def test_parallel_many_tiny_chunks():
    """Many small chunks must produce a valid non-overlapping hit set."""
    overlap = 2 * _PCHUNK_PARAMS["maxrep"] + _PCHUNK_PARAMS["maxspacer"] + 1
    chunk_size = overlap + 100
    seq = "A" * (chunk_size * 3)
    p = _tmp_fasta_parallel(("tiny", seq))
    try:
        ser = hseeker.scan_fasta(str(p), **_PCHUNK_PARAMS)
        par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=3, **_PCHUNK_PARAMS)

        def _no_overlap(hits):
            sh = sorted(hits, key=lambda h: h["start"])
            for i in range(len(sh) - 1):
                assert sh[i]["end"] < sh[i + 1]["start"], \
                    f"overlap: {sh[i]} vs {sh[i+1]}"
        _no_overlap(ser)
        _no_overlap(par)
        assert len(ser) == len(par), f"serial={len(ser)} parallel={len(par)}"
    finally:
        p.unlink(missing_ok=True)


def test_parallel_hit_at_exact_boundary_with_max_arm():
    """Max-arm hit whose center lands exactly at the chunk boundary exclusion point."""
    chunk_size = 500
    maxrep = _PCHUNK_PARAMS["maxrep"]
    seq = "A" * (chunk_size + 2 * maxrep + 100)
    p = _tmp_fasta_parallel(("boundary", seq))
    try:
        ser = hseeker.scan_fasta(str(p), **_PCHUNK_PARAMS)
        par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=2, **_PCHUNK_PARAMS)
        assert sorted(_pkey(h) for h in ser) == sorted(_pkey(h) for h in par)
    finally:
        p.unlink(missing_ok=True)


def test_parallel_n_bases_near_boundary():
    """N bases near chunk boundaries must be handled without duplicates or misses."""
    chunk_size = 500
    seq = "A" * 200 + "N" * 100 + "A" * 200 + "N" * 50 + "A" * 5000
    p = _tmp_fasta_parallel(("n_boundary", seq))
    try:
        ser = hseeker.scan_fasta(str(p), **_PCHUNK_PARAMS)
        par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=3, **_PCHUNK_PARAMS)
        assert sorted(_pkey(h) for h in ser) == sorted(_pkey(h) for h in par)
    finally:
        p.unlink(missing_ok=True)


def test_parallel_determinism():
    """Same input, same output across two independent parallel runs."""
    chunk_size = 800
    seq = "G" * 3000 + "C" * 2000 + "A" * 5000
    p = _tmp_fasta_parallel(("det", seq))
    try:
        a = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=3, **_PCHUNK_PARAMS)
        b = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=3, **_PCHUNK_PARAMS)
        assert sorted(_pkey(h) for h in a) == sorted(_pkey(h) for h in b)
    finally:
        p.unlink(missing_ok=True)


# ===========================================================================
# 25. Real-genome regression  (hg38 chr1 — requires benchmarks/data/chr1.fa)
#
# These tests are skipped automatically when the chr1 FASTA is not present
# (e.g. in CI without the benchmark data).  Run them locally after downloading:
#
#   python benchmarks/benchmark.py --chromosomes chr1
#
# Reference values were established on 2026-06-16 using hseeker v0.1.0 with
# default parameters (minrep=10, maxrep=1000, maxspacer=10, purity=0.90,
# mismatch=0.10, remove_overlaps=True, score=True).
#
#   Total hits : 96,729
#   Chromosome : chr1 (hg38/GRCh38, 249,698,492 bp)
#
# If the hit count changes after a code modification you MUST verify that
# the change is intentional (algorithm fix) and update this constant.
# ===========================================================================

_CHR1_FA = Path(__file__).parent.parent / "benchmarks" / "data" / "chr1.fa"

# Confirmed hit count — update only if the algorithm is intentionally changed.
_CHR1_EXPECTED_HITS = 96_729

_CHR1_DEFAULT_PARAMS = dict(
    minrep=10, maxrep=1000, maxspacer=10,
    purity=0.90, mismatch=0.10,
    remove_overlaps=True,
    score=False,   # scoring does not affect hit count; skip for speed
)


@pytest.mark.skipif(not _CHR1_FA.exists(), reason="benchmarks/data/chr1.fa not present")
def test_chr1_hit_count_regression():
    """Exact hit count on hg38 chr1 must not change between releases.

    Uses scan_fasta_parallel (all cores, default chunk_size=1_000_000) which
    is the same code path as the CLI.  score=False to keep the test fast.
    """
    hits = hseeker.scan_fasta_parallel(
        str(_CHR1_FA),
        **_CHR1_DEFAULT_PARAMS,
    )
    assert len(hits) == _CHR1_EXPECTED_HITS, (
        f"chr1 hit count changed: expected {_CHR1_EXPECTED_HITS}, "
        f"got {len(hits)}.  If this is intentional, update _CHR1_EXPECTED_HITS."
    )


@pytest.mark.skipif(not _CHR1_FA.exists(), reason="benchmarks/data/chr1.fa not present")
def test_chr1_hits_no_overlaps():
    """Every hit in the chr1 result must be non-overlapping."""
    hits = hseeker.scan_fasta_parallel(
        str(_CHR1_FA),
        **_CHR1_DEFAULT_PARAMS,
    )
    hits.sort(key=lambda h: h["start"])
    for i in range(len(hits) - 1):
        assert hits[i]["end"] < hits[i + 1]["start"], (
            f"Overlap found: [{hits[i]['start']},{hits[i]['end']}] "
            f"vs [{hits[i+1]['start']},{hits[i+1]['end']}]"
        )


@pytest.mark.skipif(not _CHR1_FA.exists(), reason="benchmarks/data/chr1.fa not present")
def test_chr1_coordinate_invariant():
    """start + total_length - 1 == end must hold for every hit on chr1."""
    hits = hseeker.scan_fasta_parallel(
        str(_CHR1_FA),
        **_CHR1_DEFAULT_PARAMS,
    )
    bad = [
        h for h in hits
        if h["start"] + h["total_length"] - 1 != h["end"]
    ]
    assert not bad, f"{len(bad)} hits violate coordinate invariant"


# ===========================================================================
# 26. AT-content pre-filter & homopolymer post-filter
#
# HSeeker can optionally drop two classes of hits that are unlikely to fold
# into real H-DNA: AT-rich arms (checked on left_arm, before scoring) and
# putative_triplex sequences that collapse to a single repeated base, e.g.
# poly-A / poly-G (checked after scoring). Both are opt-in (at_threshold
# defaults to None, filter_homopolymers defaults to False) so existing
# callers are unaffected; the CLI turns them on by default.
# ===========================================================================

POLY_A_SEQ = "A" * 30   # 100% AT-content left arm — dropped by AT filter
POLY_G_SEQ = "G" * 30   # 0% AT-content but a homopolymer putative_triplex


def test_at_filter_default_off_does_not_change_hit_count():
    """at_threshold=None (default) must not alter scan_sequence output."""
    default  = hseeker.scan_sequence(POLY_A_SEQ, minrep=6)
    explicit = hseeker.scan_sequence(POLY_A_SEQ, minrep=6, at_threshold=None)
    assert default == explicit


def test_at_filter_drops_pure_at_arm():
    """A pure-A arm has AT_content=1.0 and must be dropped at threshold=0.8."""
    assert len(hseeker.scan_sequence(POLY_A_SEQ, minrep=6)) > 0
    filtered = hseeker.scan_sequence(POLY_A_SEQ, minrep=6, at_threshold=0.8)
    assert filtered == []


def test_at_filter_keeps_low_at_arm():
    """PURE_GA has 0% AT content — must survive any AT threshold <= 1.0."""
    default  = hseeker.scan_sequence(PURE_GA, minrep=6)
    filtered = hseeker.scan_sequence(PURE_GA, minrep=6, at_threshold=0.8)
    assert len(filtered) == len(default)


def test_at_filter_boundary_is_exclusive_below_threshold():
    """AT_content >= threshold is dropped; AT_content < threshold is kept."""
    hits = hseeker.scan_sequence(POLY_A_SEQ, minrep=6)
    assert hits, "expected at least one hit on a pure-A run"
    at_content = sum(1 for b in hits[0]["left_arm"] if b in "at") / hits[0]["arm_length"]
    assert at_content == pytest.approx(1.0)
    assert hseeker.scan_sequence(POLY_A_SEQ, minrep=6, at_threshold=at_content) == []
    assert len(hseeker.scan_sequence(POLY_A_SEQ, minrep=6, at_threshold=at_content + 0.01)) > 0


def test_homopolymer_filter_default_off_does_not_change_hit_count():
    default  = hseeker.scan_sequence(POLY_G_SEQ, minrep=6)
    explicit = hseeker.scan_sequence(POLY_G_SEQ, minrep=6, filter_homopolymers=False)
    assert default == explicit


def test_homopolymer_filter_drops_poly_g_triplex():
    """POLY_G_SEQ scores to a putative_triplex of a single repeated 'g' — must be dropped."""
    hits = hseeker.scan_sequence(POLY_G_SEQ, minrep=6)
    assert hits and hits[0]["putative_triplex"], "expected a scored hit on a pure-G run"
    filtered = hseeker.scan_sequence(POLY_G_SEQ, minrep=6, filter_homopolymers=True)
    assert filtered == []


def test_homopolymer_filter_keeps_mixed_triplex():
    """GAMIR_SEQ's putative_triplex mixes G/A/T — must survive the homopolymer filter."""
    default  = hseeker.scan_sequence(GAMIR_SEQ, minrep=6, purity=0.80, mismatch=0.20)
    filtered = hseeker.scan_sequence(
        GAMIR_SEQ, minrep=6, purity=0.80, mismatch=0.20, filter_homopolymers=True
    )
    assert len(filtered) == len(default)


def test_homopolymer_filter_noop_when_score_false():
    """filter_homopolymers requires putative_triplex, so it must be a no-op when score=False."""
    unscored = hseeker.scan_sequence(POLY_G_SEQ, minrep=6, score=False, filter_homopolymers=True)
    assert len(unscored) > 0


def test_at_filter_applies_per_record_in_scan_fasta():
    """scan_fasta must apply at_threshold to every record, not just the first."""
    path = fasta_to_tmp([("poly_a", POLY_A_SEQ), ("ga", PURE_GA)])
    try:
        hits = hseeker.scan_fasta(str(path), minrep=6, at_threshold=0.8)
        ids = {h["seq_id"] for h in hits}
        assert "poly_a" not in ids
        assert "ga" in ids
    finally:
        path.unlink(missing_ok=True)


def test_homopolymer_filter_applies_in_scan_fasta():
    path = fasta_to_tmp([("poly_g", POLY_G_SEQ), ("gamir", GAMIR_SEQ)])
    try:
        hits = hseeker.scan_fasta(
            str(path), minrep=6, purity=0.80, mismatch=0.20, filter_homopolymers=True
        )
        ids = {h["seq_id"] for h in hits}
        assert "poly_g" not in ids
        assert "gamir" in ids
    finally:
        path.unlink(missing_ok=True)


def test_at_filter_and_homopolymer_filter_apply_in_scan_fasta_parallel():
    """Both filters must also work through the chunked/parallel entry point."""
    path = fasta_to_tmp([("poly_a", POLY_A_SEQ), ("poly_g", POLY_G_SEQ), ("ga", PURE_GA)])
    try:
        hits = hseeker.scan_fasta_parallel(
            str(path), minrep=6, at_threshold=0.8, filter_homopolymers=True
        )
        ids = {h["seq_id"] for h in hits}
        assert "poly_a" not in ids
        assert "poly_g" not in ids
        assert "ga" in ids
    finally:
        path.unlink(missing_ok=True)


def test_cli_default_at_threshold_drops_poly_a_record():
    """The CLI defaults to at_threshold=0.80 and filter_homopolymers=True."""
    path = fasta_to_tmp([("poly_a", POLY_A_SEQ), ("ga", PURE_GA)])
    out_prefix = str(path.parent / "cli_at_filter")
    tsv_path = Path(out_prefix + "_HDNA.tsv")
    try:
        subprocess.run(
            [sys.executable, "-m", "hseeker",
             "-seq", str(path), "-out", out_prefix, "-minrep", "6"],
            check=True, capture_output=True,
        )
        with open(tsv_path) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        seq_ids = {row["seq_id"] for row in rows}
        assert "poly_a" not in seq_ids
        assert "ga" in seq_ids
    finally:
        path.unlink(missing_ok=True)
        tsv_path.unlink(missing_ok=True)


def test_cli_at_threshold_flag_is_configurable():
    """Raising -at-threshold above the poly-A arm's AT content lets it through.

    Scoring is disabled here so the always-on homopolymer post-filter (which
    would also drop this poly-A hit) does not confound the AT-threshold check.
    """
    path = fasta_to_tmp([("poly_a", POLY_A_SEQ)])
    out_prefix = str(path.parent / "cli_at_relaxed")
    tsv_path = Path(out_prefix + "_HDNA.tsv")
    try:
        subprocess.run(
            [sys.executable, "-m", "hseeker",
             "-seq", str(path), "-out", out_prefix, "-minrep", "6",
             "-at-threshold", "1.5", "-no-score"],
            check=True, capture_output=True,
        )
        with open(tsv_path) as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert len(rows) > 0, "poly-A record should survive a >1.0 AT threshold"
    finally:
        path.unlink(missing_ok=True)
        tsv_path.unlink(missing_ok=True)
