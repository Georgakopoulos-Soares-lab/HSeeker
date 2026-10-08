"""Evaluation checks for the manuscript ablation notebook."""

import json

import pytest

from analysis.scripts.ablation_eval import (
    confusion_summary,
    load_matched_cohort,
    load_verified_insertions,
    max_insertion_overlap_fraction,
    rank_auc,
)


def test_cohort_requires_matching_ids_and_sequences(tmp_path):
    table = tmp_path / "labels.csv"
    fasta = tmp_path / "sequences.fa"
    table.write_text(
        "record_id,sequence_5to3,label\n"
        "p,GAAGAA,forming\n"
        "n,CTCTCT,non-forming\n"
    )
    fasta.write_text(">p, description\ngaagaa\n>n, description\nctctct\n")
    records = load_matched_cohort(table, fasta)
    assert records == [
        {"record_id": "p", "sequence": "GAAGAA", "label": "forming"},
        {"record_id": "n", "sequence": "CTCTCT", "label": "non_forming"},
    ]
    fasta.write_text(">p, description\ngaagaa\n>n, description\nctctca\n")
    with pytest.raises(ValueError, match="sequence differs"):
        load_matched_cohort(table, fasta)
    fasta.write_text(">p, description\ngaagaa\n>n, description\nctctct\n"
                     ">n, duplicate\nctctct\n")
    with pytest.raises(ValueError, match="repeated record IDs"):
        load_matched_cohort(table, fasta)


def test_record_metrics_and_auc_use_fixed_labels():
    rows = [
        {"label": "forming", "has_hit": True, "max_score": 4},
        {"label": "forming", "has_hit": False, "max_score": 0},
        {"label": "non_forming", "has_hit": True, "max_score": 1},
        {"label": "non_forming", "has_hit": False, "max_score": 0},
    ]
    assert confusion_summary(rows) == {
        "TP": 1, "FN": 1, "TN": 1, "FP": 1,
        "sensitivity": 0.5, "specificity": 0.5,
    }
    # Four positive-negative pairs: one win twice, one loss, one tie.
    assert rank_auc(rows) == 0.625


def test_insertion_overlap_uses_inclusive_coordinates():
    hits = [{"start": 90, "end": 100}, {"start": 105, "end": 109}]
    assert max_insertion_overlap_fraction(hits, 100, 109) == 0.5
    assert max_insertion_overlap_fraction([], 100, 109) == 0
    with pytest.raises(ValueError, match="precedes"):
        max_insertion_overlap_fraction(hits, 109, 100)


def test_insertion_manifest_positions_are_verified_as_zero_based(tmp_path):
    manifest = tmp_path / "insertions.json"
    manifest.write_text(json.dumps([{
        "record_id": "p", "sequence": "GAAGAA", "label": "forming",
        "position": 4,
    }]))
    records = [{"record_id": "p", "sequence": "GAAGAA", "label": "forming"}]
    assert load_verified_insertions(manifest, "TTTTGAAGAACCCC", records) == [
        {"record_id": "p", "label": "forming", "start": 5, "end": 10}
    ]
    with pytest.raises(ValueError, match="not found"):
        load_verified_insertions(manifest, "TTTTTGAAGACCCC", records)
    manifest.write_text(json.dumps([{
        "record_id": "p", "sequence": "GAAGAA", "label": "forming",
        "position": -1,
    }]))
    with pytest.raises(ValueError, match="outside genome"):
        load_verified_insertions(manifest, "TTTTGAAGAACCCC", records)
