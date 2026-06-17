from __future__ import annotations

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import hseeker  # noqa: E402


PARAMS = dict(minrep=10, maxrep=100, maxspacer=10, purity=0.90, mismatch=0.10, score=False)


def comparable(hits):
    return [
        (
            h["start"],
            h["end"],
            h["arm_length"],
            h["spacer_length"],
            round(float(h["ga_pct"]), 3),
            round(float(h["ct_pct"]), 3),
            round(float(h["mirror_identity"]), 3),
            h["left_arm"],
            h["spacer"],
            h["right_arm"],
        )
        for h in hits
    ]


def assert_safe_matches_original(seq: str, **overrides):
    params = dict(PARAMS)
    params.update(overrides)
    original = hseeker.search(seq, detector="original", **params)
    safe = hseeker.search(seq, detector="safe_pruned", safe_prune_purity=False, **params)
    safe_purity = hseeker.search(seq, detector="safe_pruned", safe_prune_purity=True, **params)
    assert comparable(safe) == comparable(original)
    assert comparable(safe_purity) == comparable(original)
    return original, safe


def test_safe_pruned_detects_perfect_ga_mirror():
    left = "AAGAAGAAGAA"
    seq = left + "TTC" + left[::-1]
    original, _safe = assert_safe_matches_original(seq)
    assert original


def test_safe_pruned_allows_one_mismatch_in_ten():
    left = "AAGAAGAAGAA"
    right = list(left[::-1])
    right[4] = "G" if right[4] != "G" else "A"
    seq = left + "C" + "".join(right)
    original, _safe = assert_safe_matches_original(seq, mismatch=0.10)
    assert original


def test_safe_pruned_rejects_one_mismatch_when_strict():
    left = "AAGAAGAAGAA"
    right = list(left[::-1])
    right[4] = "G" if right[4] != "G" else "A"
    seq = left + "C" + "".join(right)
    assert_safe_matches_original(seq, mismatch=0.0)


def test_safe_pruned_matches_nonmirror_ga_rich():
    seq = "AAGAAGAAGAA" + "TT" + "GGAGGAGGAGG"
    assert_safe_matches_original(seq)


def test_safe_pruned_matches_purity_failing_mirror():
    left = "ACGTACGTAC"
    seq = left + "A" + left[::-1]
    assert_safe_matches_original(seq)


def test_safe_pruned_matches_right_arm_invalid_n():
    left = "AAGGAAAAGA"
    right = left[::-1]
    seq = left + "T" + right[:5] + "N" + right[6:]
    assert_safe_matches_original(seq)


def test_safe_pruned_preserves_later_recovery_after_bad_prefix():
    left = "AAGGAAAAGA"
    right = list(left[::-1])
    right[0] = "T" if right[0] != "T" else "C"
    seq = left + "G" + "".join(right)
    original, _safe = assert_safe_matches_original(seq, mismatch=0.10)
    assert original


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok {name}")
