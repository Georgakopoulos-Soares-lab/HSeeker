"""All public motif sequence fields use the same lowercase convention."""

import csv
import subprocess
import sys

import hseeker
from hseeker.ablation import score_candidates


def test_sequence_and_ablation_scores_match_detected_case():
    sequence = "GA" * 16
    hits = hseeker.scan_sequence(sequence, minrep=8, remove_overlaps=False)
    assert hits
    for hit in hits:
        for key in ("left_arm", "spacer", "right_arm", "full_sequence",
                    "putative_triplex"):
            assert hit[key] == hit[key].lower()

    raw = hseeker.scan_sequence(sequence, minrep=8, remove_overlaps=False,
                                score=False)
    for hit in score_candidates(raw):
        assert hit["putative_triplex"] == hit["putative_triplex"].lower()


def test_cli_uses_same_case_for_detected_and_scored_motifs(tmp_path):
    path = tmp_path / "motif.fa"
    path.write_text(">motif\n" + "GA" * 16 + "\n")
    prefix = tmp_path / "out"
    subprocess.run([
        sys.executable, "-m", "hseeker", "-seq", str(path),
        "-out", str(prefix), "-minrep", "8", "-workers", "1",
    ], check=True, capture_output=True)
    with (tmp_path / "out_HDNA.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    for row in rows:
        for key in ("left_arm", "right_arm", "full_sequence", "putative_triplex"):
            assert row[key] == row[key].lower()
