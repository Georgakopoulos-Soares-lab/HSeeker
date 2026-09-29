"""Boundary optimization respects its minimum arm size."""

from hseeker._scoring import _scorer, score_hit


def test_optimizer_does_not_select_one_base_arm():
    result = score_hit("GGTGAGAGGG", "", "GTATATAGGA", 10)
    assert result is not None
    left_arm = result["putative_triplex"].split("[", 1)[0]
    assert len(left_arm) >= 8


def test_short_detected_arm_can_still_be_scored():
    result = score_hit("G" * 6, "", "G" * 6, 6)
    assert result is not None
    assert result["putative_triplex"] == "g" * 6 + "[]" + "g" * 6


def test_stacking_penalty_grows_until_two_adjacent_matches_reset_it():
    stack = _scorer._calc_stacking
    assert sum(stack("1111111")) == 30.0
    assert sum(stack("1110111")) == 10.0
    assert sum(stack("11100111")) == -2.5
    assert sum(stack("111000111")) == -65.0
    # A lone matched position between mismatches is not a stacked pair.
    assert stack("111010111") == stack("111000111")
    # The next "11" resets the growing penalty before a later mismatch.
    assert stack("11101110111") == [0.0, 5.0, 5.0, -5, -5, 5.0, 5.0,
                                      -5, -5, 5.0, 5.0]
