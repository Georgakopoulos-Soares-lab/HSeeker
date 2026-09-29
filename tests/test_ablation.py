"""The ablation API must isolate model components without changing baseline scoring."""

import pytest

import hseeker
from hseeker.ablation import (
    MODELS,
    filter_at_content,
    filter_homopolymers,
    remove_overlaps,
    scan_fasta_ablation,
    scan_sequence_ablation,
    score_candidates,
)


SEQ = "GA" * 16
OPTIONS = dict(minrep=8, maxrep=20, maxspacer=4)


def test_complete_matches_public_api_with_same_options():
    expected = hseeker.scan_sequence(
        SEQ, remove_overlaps=False, at_threshold=0.80,
        filter_homopolymers=True, **OPTIONS,
    )
    assert scan_sequence_ablation(SEQ, **OPTIONS) == expected


def test_stage_scores_are_isolated_on_same_candidates():
    results = {model: scan_sequence_ablation(SEQ, model=model, **OPTIONS)
               for model in MODELS}
    assert len({len(results[m]) for m in MODELS}) == 1
    assert "total_score" not in results["detection"][0]
    assert "total_score" not in results["composition"][0]
    for pairing, stacked, full in zip(
        results["pairing"], results["pairing_stacking"], results["complete"]
    ):
        assert pairing["stacking_score"] == 0
        assert pairing["pairing_score"] == stacked["pairing_score"]
        assert stacked["stacking_score"] > 0
        assert full["total_score"] >= stacked["total_score"]


def test_detection_disables_purity_but_preserves_mismatch_threshold():
    seq = "ACGT" * 8
    raw = scan_sequence_ablation(seq, model="detection", mismatch=0,
                                 **OPTIONS)
    expected = hseeker.scan_sequence(seq, purity=0, mismatch=0, score=False,
                                     remove_overlaps=False, **OPTIONS)
    assert raw == expected


def test_public_stage_helpers_do_not_mutate_candidates():
    raw = scan_sequence_ablation(SEQ, model="detection", **OPTIONS)
    scored = score_candidates(raw, model="pairing")
    assert "total_score" not in raw[0]
    assert "total_score" in scored[0]
    assert filter_at_content(raw, 1.0) == raw
    assert filter_homopolymers(scored) == scored
    assert len(remove_overlaps(scored)) < len(scored)


def test_fasta_offsets_and_sequence_ids(tmp_path):
    path = tmp_path / "motifs.fa"
    path.write_text(">chr1:101-132\n" + SEQ + "\n>chr2\n" + SEQ + "\n")
    hits = scan_fasta_ablation(path, model="pairing", **OPTIONS)
    assert {h["seq_id"] for h in hits} == {"chr1:101-132", "chr2"}
    assert min(h["start"] for h in hits if h["seq_id"] == "chr1:101-132") == 101
    assert min(h["start"] for h in hits if h["seq_id"] == "chr2") == 1


def test_score_overlap_requires_score_model():
    with pytest.raises(ValueError, match="requires a scoring model"):
        scan_sequence_ablation(SEQ, model="detection",
                               remove_overlaps_enabled=True,
                               overlap_strategy="score")


def test_top_level_components_reproduce_scored_overlap_pipeline():
    sequence = "GGAAGAAGAAGGAGGAGAGAGGGAAGGAA"
    options = dict(minrep=7, maxrep=25, maxspacer=5,
                   purity=0.9, mismatch=0.2)
    raw = hseeker.detect_candidates(sequence, **options)
    assert raw == hseeker.scan_sequence(
        sequence, remove_overlaps=False, score=False, **options
    )
    composition = hseeker.filter_at_content(raw, threshold=0.8)
    scored = hseeker.score_candidates(composition, model="complete")
    filtered = hseeker.filter_homopolymers(scored, mode="POST")
    selected = hseeker.remove_overlaps(filtered, strategy="score")
    assert selected == hseeker.scan_sequence(
        sequence, overlap_strategy="score", at_threshold=0.8,
        filter_homopolymers=True, **options,
    )
    assert len(raw) == 2
    assert "total_score" not in raw[0]
    assert hseeker.ABLATION_MODELS == MODELS
    assert hseeker.scan_sequence_ablation is scan_sequence_ablation
    assert hseeker.scan_fasta_ablation is scan_fasta_ablation


def test_top_level_pre_and_post_homopolymer_filters():
    raw = hseeker.detect_candidates("G" * 18, minrep=6, maxrep=10,
                                    maxspacer=2)
    assert raw
    assert hseeker.filter_homopolymers(raw, mode="PRE") == []
    scored = hseeker.score_candidates(raw)
    assert hseeker.filter_homopolymers(scored, mode="POST") == []
    assert hseeker.score_hit_components(
        raw[0]["left_arm"], raw[0]["spacer"], raw[0]["right_arm"],
        raw[0]["arm_length"],
    ) == hseeker.score_hit(
        raw[0]["left_arm"], raw[0]["spacer"], raw[0]["right_arm"],
        raw[0]["arm_length"],
    )
