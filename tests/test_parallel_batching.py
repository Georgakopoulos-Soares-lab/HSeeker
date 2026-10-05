"""scan_fasta_parallel scores and filters batches of whole records instead of the
whole assembly at once (bounds memory in score mode). Batching must not change
the hits, their scores, or their order."""

import random

import pytest

import hseeker


def _random_dna(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mirror_tract(rng, arm):
    half = "".join(rng.choice("GA") for _ in range(arm))
    return half + _random_dna(rng, rng.randint(0, 6)) + half[::-1]


@pytest.fixture
def multi_record_fasta(tmp_path):
    rng = random.Random(11)
    path = tmp_path / "genome.fa"
    with open(path, "w") as f:
        for i in range(8):
            parts = []
            for _ in range(rng.randint(2, 12)):
                parts.append(_random_dna(rng, rng.randint(100, 1500)))
                parts.append(_mirror_tract(rng, rng.randint(10, 30)))
            f.write(f">rec{7 - i} test record\n{''.join(parts)}\n")  # ids not in sorted order
    return path


@pytest.mark.parametrize("strategy", ["score", "greedy"])
@pytest.mark.parametrize("remove_overlaps", [True, False])
def test_batching_does_not_change_hits(multi_record_fasta, strategy, remove_overlaps):
    kwargs = dict(minrep=10, maxspacer=10, remove_overlaps=remove_overlaps, overlap_strategy=strategy,
                  at_threshold=0.8, filter_homopolymers=True)
    # Same chunk_size, so records are chunked identically; batches hold about
    # workers * chunk_size bases: ~8 kb (several records) vs ~1 kb (one record).
    single = hseeker.scan_fasta_parallel(str(multi_record_fasta), workers=8, chunk_size=1_000, **kwargs)
    batched = hseeker.scan_fasta_parallel(str(multi_record_fasta), workers=1, chunk_size=1_000, **kwargs)
    assert single, "fixture should produce hits"
    assert batched == single


def test_batched_output_sorted_by_seq_id_with_overlap_removal(multi_record_fasta):
    hits = hseeker.scan_fasta_parallel(str(multi_record_fasta), workers=2, chunk_size=1_000,
                                       remove_overlaps=True, overlap_strategy="score")
    seq_ids = [hit["seq_id"] for hit in hits]
    assert seq_ids == sorted(seq_ids)
