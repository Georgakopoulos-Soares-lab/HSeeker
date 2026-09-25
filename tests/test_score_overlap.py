"""Behavior of the optional post-hoc overlap strategies."""

import csv
import subprocess
import sys

import pytest

from hseeker import (
    _filter_overlapping_hits,
    scan_fasta,
    scan_fasta_iter,
    scan_fasta_parallel,
    scan_sequence,
)


def hit(start, end, score, *, seq_id="chr1", arm=8, spacer=2):
    return {
        "seq_id": seq_id,
        "start": start,
        "end": end,
        "total_score": score,
        "arm_length": arm,
        "spacer_length": spacer,
    }


def test_score_priority_keeps_both_ends_of_overlap_chain():
    left = hit(1, 10, 90)
    bridge = hit(8, 18, 100)
    right = hit(16, 25, 80)
    # The bridge has the highest score, so neither end can be retained.
    assert _filter_overlapping_hits([right, left, bridge]) == [bridge]

    bridge["total_score"] = 70
    # Now the ends win independently even though all three are connected.
    assert _filter_overlapping_hits([bridge, right, left]) == [left, right]


def test_inclusive_coordinates_sequence_ids_and_ties():
    left = hit(1, 10, 40, arm=10)
    touching = hit(10, 20, 40, arm=9)
    other_seq = hit(1, 10, 20, seq_id="chr2")
    assert _filter_overlapping_hits([touching, other_seq, left]) == [left, other_seq]
    assert _filter_overlapping_hits([left, touching, other_seq]) == [left, other_seq]


def test_unscorable_hit_loses_and_unscored_input_is_rejected():
    scored = hit(1, 10, 0)
    unscorable = hit(2, 11, None)
    assert _filter_overlapping_hits([unscorable, scored]) == [scored]
    del unscorable["total_score"]
    with pytest.raises(ValueError, match="already scored"):
        _filter_overlapping_hits([scored, unscorable])


def test_greedy_remains_available_without_scores_and_can_choose_differently():
    high_score = hit(1, 20, 100, arm=8)
    long_arm = hit(2, 21, 10, arm=12)
    unscored = [{k: v for k, v in h.items() if k != "total_score"}
                for h in (high_score, long_arm)]
    assert _filter_overlapping_hits(unscored, strategy="greedy") == [unscored[1]]
    assert _filter_overlapping_hits([long_arm, high_score], strategy="score") == [high_score]
    assert _filter_overlapping_hits([long_arm, high_score], strategy="stability") == [high_score]


def test_dense_connected_cluster_keeps_several_separate_winners():
    # All candidates belong to one transitive overlap component. They are
    # schematic intervals with illustrative scores, not scanner output.
    candidates = [
        hit(1, 12, 90),    # A: blocked by B
        hit(4, 18, 120),   # B: kept first
        hit(13, 24, 110),  # C: blocked by B
        hit(22, 35, 95),   # D: blocked by E
        hit(30, 44, 100),  # E: kept
        hit(42, 55, 80),   # F: blocked by E and H
        hit(56, 65, 70),   # G: blocked by H
        hit(51, 60, 105),  # H: kept
    ]
    expected = [candidates[i] for i in (1, 4, 7)]
    assert _filter_overlapping_hits(candidates) == expected
    assert _filter_overlapping_hits(list(reversed(candidates))) == expected


@pytest.mark.parametrize("preferred,other", [
    (hit(2, 20, 50, arm=12, spacer=5), hit(1, 21, 50, arm=10, spacer=1)),
    (hit(2, 20, 50, arm=10, spacer=1), hit(1, 21, 50, arm=10, spacer=3)),
    (hit(1, 20, 50), hit(2, 21, 50)),
])
def test_equal_score_tiebreak_order(preferred, other):
    assert _filter_overlapping_hits([other, preferred]) == [preferred]
    assert _filter_overlapping_hits([preferred, other]) == [preferred]


def test_touching_base_overlaps_but_next_base_does_not():
    first = hit(1, 10, 30)
    touching = hit(10, 20, 20)
    next_base = hit(11, 20, 10)
    assert _filter_overlapping_hits([touching, next_base, first]) == [first, next_base]


def test_scores_are_compared_only_within_same_sequence():
    chr1 = hit(1, 20, 10)
    chr2 = hit(1, 20, 100, seq_id="chr2")
    chr2_loser = hit(5, 15, 20, seq_id="chr2")
    assert _filter_overlapping_hits([chr2_loser, chr2, chr1]) == [chr1, chr2]


def test_none_and_nan_scores_rank_after_real_scores():
    scored = hit(1, 20, 0)
    no_score = hit(2, 19, None)
    nan_score = hit(3, 18, float("nan"))
    assert _filter_overlapping_hits([nan_score, no_score, scored]) == [scored]


def test_score_priority_is_not_maximum_sum_of_scores():
    long_hit = hit(1, 30, 100)
    two_short_hits = [hit(1, 10, 60), hit(20, 30, 60)]
    # The two smaller scores sum to 120, but the rule gives priority to
    # the single highest individual score.
    assert _filter_overlapping_hits(two_short_hits + [long_hit]) == [long_hit]


def test_filter_does_not_change_input_order_or_hit_dicts():
    candidates = [hit(2, 20, 10), hit(1, 18, 20)]
    original = [candidate.copy() for candidate in candidates]
    assert _filter_overlapping_hits(candidates) == [candidates[1]]
    assert candidates == original


def test_invalid_strategy_is_rejected():
    with pytest.raises(ValueError, match="Invalid deduplication strategy"):
        _filter_overlapping_hits([], strategy="random")


def test_real_scanner_cluster_uses_real_stability_scores():
    sequence = "GA" * 21
    candidates = scan_sequence(
        sequence, minrep=8, maxrep=25, maxspacer=5,
        remove_overlaps=False, score=True,
    )
    assert len(candidates) > 10
    chosen = _filter_overlapping_hits(candidates)
    assert len(chosen) == 1
    assert (chosen[0]["start"], chosen[0]["end"]) == (1, 41)
    assert chosen[0]["total_score"] == max(h["total_score"] for h in candidates)


def test_score_strategy_is_applied_after_scoring_in_scan_sequence():
    sequence = "GA" * 21
    raw = scan_sequence(sequence, minrep=8, maxrep=25, maxspacer=5,
                        remove_overlaps=False, score=True)
    selected = scan_sequence(sequence, minrep=8, maxrep=25, maxspacer=5,
                             overlap_strategy="score")
    assert selected == _filter_overlapping_hits(raw, strategy="score")
    assert len(raw) > len(selected)


def test_score_strategy_agrees_across_fasta_apis_and_chunk_boundaries(tmp_path):
    path = tmp_path / "clusters.fa"
    path.write_text(">ga\n" + "GA" * 50 + "\n>ct\n" + "CT" * 30 + "\n")
    options = dict(minrep=8, maxrep=25, maxspacer=5,
                   remove_overlaps=True, overlap_strategy="score")
    sequential = scan_fasta(path, **options)
    streamed = list(scan_fasta_iter(path, **options))
    parallel = scan_fasta_parallel(path, workers=2, chunk_size=30, **options)
    key = lambda h: (h["seq_id"], h["start"], h["end"], h["total_score"])
    assert sorted(map(key, sequential)) == sorted(map(key, streamed))
    assert sorted(map(key, sequential)) == sorted(map(key, parallel))


def test_greedy_remains_the_default_in_public_api():
    sequence = "GA" * 21
    options = dict(minrep=8, maxrep=25, maxspacer=5)
    assert scan_sequence(sequence, **options) == scan_sequence(
        sequence, overlap_strategy="greedy", **options
    )


def test_score_overlap_requires_scoring_but_skipoverlap_does_not():
    sequence = "GA" * 21
    with pytest.raises(ValueError, match="score=True"):
        scan_sequence(sequence, score=False, overlap_strategy="score")
    assert scan_sequence(sequence, score=False, remove_overlaps=False,
                         overlap_strategy="score")


def test_cli_score_strategy_matches_api(tmp_path):
    path = tmp_path / "cluster.fa"
    path.write_text(">ga\n" + "GA" * 21 + "\n")
    prefix = tmp_path / "output"
    subprocess.run([
        sys.executable, "-m", "hseeker", "-seq", str(path), "-out", str(prefix),
        "-minrep", "8", "-maxrep", "25", "-maxspacer", "5",
        "-overlap-strategy", "score", "-at-threshold", "1.0",
    ], check=True, capture_output=True)
    with (tmp_path / "output_HDNA.tsv").open() as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    api = scan_fasta(path, minrep=8, maxrep=25, maxspacer=5,
                     overlap_strategy="score")
    assert [(r["seq_id"], int(r["start"]), int(r["end"])) for r in rows] == [
        (h["seq_id"], h["start"], h["end"]) for h in api
    ]
