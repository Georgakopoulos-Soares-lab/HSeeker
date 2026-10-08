"""Behavior of the optional post-hoc overlap strategies."""

import csv
import random
import subprocess
import sys

import pytest

from hseeker import (
    _apply_scoring,
    _filter_homopolymer_triplex,
    _filter_overlapping_hits,
    _hdna,
    _score_and_select,
    scan_fasta,
    scan_fasta_iter,
    scan_fasta_parallel,
    scan_sequence,
)
from hseeker._scoring import _scorer


def hit(start, end, score, *, seq_id="chr1", arm=8, spacer=2):
    return {
        "seq_id": seq_id,
        "start": start,
        "end": end,
        "total_score": score,
        "arm_length": arm,
        "spacer_length": spacer,
    }


def test_score_selection_compares_both_ends_with_bridge():
    left = hit(1, 10, 90)
    bridge = hit(8, 18, 100)
    right = hit(16, 25, 80)
    # The highest individual score has first choice, even though both ends fit.
    assert _filter_overlapping_hits([right, left, bridge]) == [bridge]

    bridge["total_score"] = 70
    assert _filter_overlapping_hits([bridge, right, left]) == [left, right]


def test_five_hit_chain_can_keep_three_or_two_motifs():
    chain = [
        hit(1, 10, 10),
        hit(8, 18, 1),
        hit(16, 26, 9),
        hit(24, 34, 1),
        hit(32, 42, 8),
    ]
    assert _filter_overlapping_hits(chain) == [chain[i] for i in (0, 2, 4)]
    assert _filter_overlapping_hits(list(reversed(chain))) == [
        chain[i] for i in (0, 2, 4)
    ]

    chain[1]["total_score"] = 13
    chain[3]["total_score"] = 12
    assert _filter_overlapping_hits(chain) == [chain[i] for i in (1, 3)]


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
        hit(1, 12, 90),
        hit(4, 18, 120),
        hit(13, 24, 110),
        hit(22, 35, 95),
        hit(30, 44, 100),
        hit(42, 55, 80),
        hit(56, 65, 70),
        hit(51, 60, 105),
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


def test_score_selection_prioritizes_individual_score_over_sum():
    long_hit = hit(1, 30, 100)
    two_short_hits = [hit(1, 10, 60), hit(20, 30, 60)]
    assert _filter_overlapping_hits(two_short_hits + [long_hit]) == [long_hit]


def test_filter_does_not_change_input_order_or_hit_dicts():
    candidates = [hit(2, 20, 10), hit(1, 18, 20)]
    original = [candidate.copy() for candidate in candidates]
    assert _filter_overlapping_hits(candidates) == [candidates[1]]
    assert candidates == original


def test_invalid_strategy_is_rejected():
    with pytest.raises(ValueError, match="Invalid deduplication strategy"):
        _filter_overlapping_hits([], strategy="random")


def test_score_selection_matches_literal_priority_rule_on_random_intervals():
    """Cross-check the overlap index against score-ordered acceptance."""
    rng = random.Random(2026)
    for _ in range(250):
        candidates = [
            hit(
                start := rng.randint(1, 100),
                start + rng.randint(0, 20),
                score=float(i),
                seq_id=rng.choice(("chr1", "chr2")),
            )
            for i in range(rng.randint(0, 40))
        ]
        rng.shuffle(candidates)
        expected = []
        for candidate in sorted(candidates, key=lambda h: -h["total_score"]):
            if all(
                candidate["seq_id"] != prior["seq_id"]
                or candidate["end"] < prior["start"]
                or candidate["start"] > prior["end"]
                for prior in expected
            ):
                expected.append(candidate)
        selected = _filter_overlapping_hits(candidates, strategy="score")
        assert selected == sorted(expected, key=lambda h: (
            h["seq_id"], h["start"], h["end"]
        ))


def test_real_sequence_score_selection_can_beat_greedy():
    sequence = "GAAGGAAAGAAAGAAAGGG"
    options = dict(minrep=8, maxrep=22, maxspacer=5,
                   purity=0.9, mismatch=0.15)
    greedy = scan_sequence(sequence, overlap_strategy="greedy", **options)
    scored = scan_sequence(sequence, overlap_strategy="score", **options)
    assert [(h["start"], h["end"]) for h in greedy] == [(1, 17)]
    assert [(h["start"], h["end"]) for h in scored] == [(3, 19)]
    assert scored[0]["total_score"] > greedy[0]["total_score"]


@pytest.mark.parametrize("sequence, expected_raw, greedy_winner, score_winner", [
    (
        "GAGAAAAGAGGGAGAAAAAGGAA",
        [(1, 21, 10, 1, 53.635), (3, 19, 7, 3, 41.035)],
        (1, 21), (1, 21),
    ),
    (
        "GGAAGAAGAAGGAGGAGAGAGGGAAGGAA",
        [(1, 15, 7, 1, 38.575), (9, 29, 10, 1, 20.360)],
        (9, 29), (1, 15),
    ),
])
def test_review_overlap_examples_across_apis(
    tmp_path, sequence, expected_raw, greedy_winner, score_winner
):
    """Real scanner candidates show both agreement and a shorter score winner."""
    options = dict(
        minrep=7, maxrep=25, maxspacer=5, purity=0.9, mismatch=0.2,
        at_threshold=0.8, filter_homopolymers=True,
    )
    path = tmp_path / "example.fa"
    path.write_text(f">example\n{sequence}\n")

    raw = scan_sequence(sequence, remove_overlaps=False, **options)
    assert [(h["start"], h["end"], h["arm_length"], h["spacer_length"])
            for h in sorted(raw, key=lambda h: h["start"])] == [
                row[:4] for row in expected_raw
            ]
    for hit, row in zip(sorted(raw, key=lambda h: h["start"]), expected_raw):
        assert hit["total_score"] == pytest.approx(row[4])

    # PRE (unscored full_sequence) and POST (scored putative_triplex) both
    # retain these mixed-base candidates.
    pre = scan_sequence(sequence, remove_overlaps=False, score=False, **options)
    assert {(h["start"], h["end"]) for h in pre} == {
        (h["start"], h["end"]) for h in raw
    }

    for strategy, winner in (("greedy", greedy_winner), ("score", score_winner)):
        kwargs = dict(overlap_strategy=strategy, **options)
        selected = scan_sequence(sequence, **kwargs)
        assert [(h["start"], h["end"]) for h in selected] == [winner]
        expected = [(h["start"], h["end"], h["arm_length"], h["total_score"])
                    for h in selected]
        for api_hits in (
            scan_fasta(path, **kwargs),
            list(scan_fasta_iter(path, **kwargs)),
            scan_fasta_parallel(path, workers=2, chunk_size=15, **kwargs),
        ):
            assert [(h["start"], h["end"], h["arm_length"], h["total_score"])
                    for h in api_hits] == expected


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


def test_real_scanner_keeps_compatible_hits_from_connected_cluster(tmp_path):
    sequence = "GAA" * 20
    options = dict(
        minrep=10, maxrep=15, maxspacer=10, purity=0.9, mismatch=0.1,
        at_threshold=0.8, filter_homopolymers=True,
    )
    raw = scan_sequence(sequence, remove_overlaps=False, **options)
    chosen = scan_sequence(sequence, overlap_strategy="score", **options)
    assert [(h["start"], h["end"], h["total_score"]) for h in chosen] == [
        (3, 32, 146.3), (35, 60, 124.87),
    ]
    assert any(h["start"] <= 32 and h["end"] >= 35 for h in raw)
    assert scan_sequence(sequence, overlap_strategy="greedy", **options) == chosen

    path = tmp_path / "chain.fa"
    path.write_text(f">chain\n{sequence}\n")
    expected = [(h["start"], h["end"], h["total_score"]) for h in chosen]
    for api_hits in (
        scan_fasta(path, overlap_strategy="score", **options),
        list(scan_fasta_iter(path, overlap_strategy="score", **options)),
        scan_fasta_parallel(path, workers=2, chunk_size=30,
                            overlap_strategy="score", **options),
    ):
        assert [(h["start"], h["end"], h["total_score"])
                for h in api_hits] == expected


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


def test_default_reports_the_longest_arm_where_score_would_differ():
    sequence = "GAAGGAAAGAAAGAAAGGG"
    options = dict(minrep=8, maxrep=22, maxspacer=5, purity=0.9, mismatch=0.15)
    assert [(h["start"], h["end"]) for h in scan_sequence(sequence, **options)] == [(1, 17)]
    assert [(h["start"], h["end"]) for h in scan_sequence(
        sequence, overlap_strategy="score", **options)] == [(3, 19)]


def test_cli_default_is_greedy(tmp_path):
    path = tmp_path / "motif.fa"
    path.write_text(">m\nGAAGGAAAGAAAGAAAGGG\n")
    options = ["-minrep", "8", "-maxrep", "22", "-maxspacer", "5",
               "-purity", "0.9", "-mismatch", "0.15", "-at-threshold", "1.0"]
    spans = {}
    for name, extra in (("default", []), ("score", ["-overlap-strategy", "score"])):
        prefix = tmp_path / name
        subprocess.run([sys.executable, "-m", "hseeker", "-seq", str(path),
                        "-out", str(prefix), *options, *extra],
                       check=True, capture_output=True)
        with (tmp_path / f"{name}_HDNA.tsv").open() as fh:
            spans[name] = [(int(r["start"]), int(r["end"]))
                           for r in csv.DictReader(fh, delimiter="\t")]
    assert spans == {"default": [(1, 17)], "score": [(3, 19)]}


def _score_everything_then_select(hits, filter_homopolymers):
    """The straightforward pipeline that lazy score selection must reproduce."""
    _apply_scoring(hits)
    if filter_homopolymers:
        hits = _filter_homopolymer_triplex(hits, mode="POST")
    return _filter_overlapping_hits(hits, strategy="score")


@pytest.mark.parametrize("sequence", [
    "G" * 120, "GA" * 60, "GAA" * 40, "CT" * 50 + "ACGT" * 5 + "G" * 70,
    "GAAGGAAAGAAAGAAAGGG", "AGGGAGGAGGGAGGAGGGAGGAGGGTTTAGGAGGGAGG",
])
@pytest.mark.parametrize("filter_homopolymers", [False, True])
def test_lazy_score_selection_matches_scoring_every_candidate(sequence,
                                                             filter_homopolymers):
    options = dict(minrep=6, maxrep=1000, maxspacer=10, purity=0.8,
                   mismatch=0.2, seq_offset=1, remove_overlaps=False)
    reference = _score_everything_then_select(
        _hdna.scan_sequence(sequence, **options), filter_homopolymers)
    lazy = _score_and_select(
        _hdna.scan_sequence(sequence, **options), filter_homopolymers)
    assert lazy == reference


def _naive_best_window(scorer, s1, s2, ll):
    """Rebuild every window score from scratch, as the scorer once did."""
    al = len(s1)
    scoring_array = "".join(
        "1" if (a == b == "G") or (a == b == "A") else "0" for a, b in zip(s1, s2))
    pairing = [scorer._pair_score(a, b) for a, b in zip(s1, s2)]
    best, best_pair = -float("inf"), None
    for L in range(al):
        for R in range(L + min(scorer.min_al, al), al + 1):
            if 2 * (al - R) + ll <= (R - L) * scorer.v:
                cur = round(sum(pairing[L:R])
                            + sum(scorer._calc_stacking(scoring_array[L:R])), 3)
                if cur > best:
                    best, best_pair = cur, (L, R)
    return best_pair


def test_incremental_window_search_matches_naive_recomputation():
    rng = random.Random(11)
    for _ in range(300):
        al = rng.randint(1, 40)
        ll = rng.randint(0, 12)
        s1 = "".join(rng.choice("GGGAAC") for _ in range(al))
        s2 = "".join(rng.choice("GGGAAT") for _ in range(al))
        full = s1 + "C" * ll + s2[::-1]
        result = _scorer.score(full, al)

        s1o, s2o, _, llo = _scorer._arms(full, al)
        L, R = _naive_best_window(_scorer, s1o, s2o, llo) or (0, al)
        scoring_array = "".join("1" if (a == b == "G") or (a == b == "A") else "0"
                                for a, b in zip(s1o, s2o))
        pairing = sum([_scorer._pair_score(a, b)
                       for a, b in zip(s1o[L:R], s2o[L:R])])
        stacking = sum(_scorer._calc_stacking(scoring_array[L:R]))
        assert result["pairing_score"] == pairing
        assert result["stacking_score"] == stacking
        assert result["total_score"] == max(round(pairing + stacking, 3), 0)
        assert result["total_score"] <= _scorer.upper_bound(full, al)



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
