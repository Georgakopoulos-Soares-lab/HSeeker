from __future__ import annotations

import pytest

import hseeker


def _fg(seq: str, **kwargs):
    params = dict(
        minrep=10,
        maxrep=50,
        maxspacer=5,
        purity=0.90,
        mismatch=0.10,
        score=False,
    )
    params.update(kwargs)
    return hseeker.search_fast_gated(seq, **params)


def test_fast_gated_detects_perfect_ga_mirror():
    left = "AAGAAGAAGAA"
    seq = left + "TTC" + left[::-1]
    hits = _fg(seq)
    assert hits
    assert any(h["arm_length"] >= 10 for h in hits)


def test_fast_gated_rejects_non_mirror_ga_rich():
    seq = "AAGAAGAAGAA" + "TTC" + "GGGGGGGGGGG"
    assert _fg(seq, mismatch=0.0) == []


def test_fast_gated_rejects_low_purity_mirror():
    left = "ACGTACGTAC"
    seq = left + "TTT" + left[::-1]
    assert _fg(seq, purity=0.80, mismatch=0.0) == []


def test_fast_gated_allows_one_mismatch_in_ten():
    left = "AAAAAAAAAA"
    right = list(left[::-1])
    right[3] = "C"
    seq = left + "".join(right)
    hits = _fg(seq, maxspacer=0, mismatch=0.10, purity=0.90)
    assert hits
    assert hits[0]["mirror_identity"] == pytest.approx(90.0, abs=0.01)


def test_fast_gated_rejects_one_mismatch_when_strict():
    left = "AAAAAAAAAA"
    right = list(left[::-1])
    right[3] = "C"
    seq = left + "".join(right)
    assert _fg(seq, maxspacer=0, mismatch=0.0, purity=0.90) == []


def test_fast_gated_right_arm_n_matches_original_behavior():
    seq = "AAAAAAAAAA" + "NN" + "AAAAAAAAAA"
    original = hseeker.scan_sequence(
        seq, minrep=10, maxrep=50, maxspacer=2, purity=0.90, mismatch=0.0, score=False
    )
    exact = hseeker.search_fast_gated(
        seq,
        minrep=10,
        maxrep=50,
        maxspacer=2,
        purity=0.90,
        mismatch=0.0,
        fast_mode=False,
        use_purity_prefilter=False,
        score=False,
    )
    assert exact == original


@pytest.mark.parametrize(
    "seq",
    [
        "AAGAAGAAGAA" + "TTC" + "AAGAAGAAGAA"[::-1],
        "AAAAAAAAAA" + "AAAAAAAAAA",
        "CCCCCCCCCC" + "GG" + "CCCCCCCCCC"[::-1],
    ],
)
def test_exact_gate_regression_matches_original_on_synthetic(seq):
    params = dict(minrep=10, maxrep=50, maxspacer=5, purity=0.90, mismatch=0.10, score=False)
    original = hseeker.scan_sequence(seq, **params)
    exact = hseeker.search_fast_gated(
        seq,
        **params,
        fast_mode=False,
        use_purity_prefilter=False,
    )
    assert exact == original
