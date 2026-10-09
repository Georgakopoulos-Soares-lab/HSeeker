"""scan_fasta_parallel scores and filters batches of whole records instead of the
whole assembly at once (bounds memory in score mode). Batching must not change
the hits, their scores, or their order.

tests/data/batching_golden.json.gz was generated with the pre-batching
implementation (main before the batching change) on tests/data/batching_genome.fa
(10 synthetic records, 300 bp to 20 kb, with 10-30 bp purine/pyrimidine mirror
repeats, an N block and a lowercase record). Results differ between chunk sizes
in a few hits because of chunk-boundary handling, which batching does not touch.
After merging fix/hdna-maximal-representation-v3 (inward extension and removal
of duplicate arm/spacer representations in _hdna.c) the golden hits were
regenerated with that same pre-batching Python code (4784c4b) built against the
merged C detector: the 10 configurations with overlap removal were unchanged,
the 10 without it changed as the detector change intends.

Floats are compared up to 6 decimals: since Python 3.12 the built-in sum()
uses compensated summation, so unrounded scores (pairing_score, stacking_score)
differ in the last bits between Python versions.
"""

import gzip
import json
import math
import random
import shutil
from pathlib import Path

import pytest

import hseeker

DATA = Path(__file__).resolve().parent / "data"
GENOME = DATA / "batching_genome.fa"
with gzip.open(DATA / "batching_golden.json.gz", "rt") as _f:
    GOLDEN = json.load(_f)

FLAGS = [
    dict(overlap_strategy="score", remove_overlaps=True),
    dict(overlap_strategy="greedy", remove_overlaps=True),
    dict(overlap_strategy="score", remove_overlaps=False),
    dict(overlap_strategy="greedy", remove_overlaps=False),
]


def _assert_same_hits(got, expected):
    """Every field equal; floats equal up to 6 decimals (Python-version float summation)."""
    assert len(got) == len(expected)
    for g, e in zip(got, expected):
        assert g.keys() == e.keys()
        for key, value in e.items():
            if isinstance(value, float):
                assert math.isclose(g[key], value, rel_tol=0, abs_tol=1e-6), (key, g, e)
            else:
                assert g[key] == value, (key, g, e)


def _scan(path=GENOME, **kwargs):
    kwargs.setdefault("minrep", 10)
    kwargs.setdefault("maxspacer", 10)
    return hseeker.scan_fasta_parallel(str(path), **kwargs)


def _random_dna(rng, n):
    return "".join(rng.choice("ACGT") for _ in range(n))


def _mirror_tract(rng, arm):
    alphabet = rng.choice(["GA", "CT"])
    half = "".join(rng.choice(alphabet) for _ in range(arm))
    return half + _random_dna(rng, rng.randint(0, 6)) + half[::-1]


def _records(rng, sizes):
    """Random records of roughly the given sizes, each with mirror repeats; ids in reverse order."""
    records = []
    for i, size in enumerate(sizes):
        parts, n = [], 0
        while n < size:
            part = _random_dna(rng, rng.randint(80, 400)) + _mirror_tract(rng, rng.randint(10, 25))
            parts.append(part)
            n += len(part)
        records.append((f"r{len(sizes) - i:03d}", "".join(parts)))
    return records


def _write_fasta(path, records):
    with open(path, "w") as f:
        for name, seq in records:
            f.write(f">{name}\n{seq}\n")
    return path


@pytest.fixture
def scoring_calls(monkeypatch):
    """Record the number of hits in each batch's scoring pass (one entry per batch).

    Greedy mode scores a batch with one _apply_scoring call. Score mode hands
    the batch to _score_and_select, which scores hits lazily one at a time, so
    the batch is counted there and its nested _apply_scoring calls are not.
    """
    calls = []
    inside_select = [False]
    real_apply_scoring = hseeker._apply_scoring
    real_score_and_select = hseeker._score_and_select

    def recording_apply_scoring(hits):
        if not inside_select[0]:
            calls.append(len(hits))
        return real_apply_scoring(hits)

    def recording_score_and_select(hits, filter_homopolymers):
        calls.append(len(hits))
        inside_select[0] = True
        try:
            return real_score_and_select(hits, filter_homopolymers)
        finally:
            inside_select[0] = False

    monkeypatch.setattr(hseeker, "_apply_scoring", recording_apply_scoring)
    monkeypatch.setattr(hseeker, "_score_and_select", recording_score_and_select)
    return calls


# ---------------------------------------------------------------------------
# 1. Same output as the pre-batching implementation (golden file)
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", sorted(GOLDEN))
@pytest.mark.parametrize("workers", [1, 2, 4])
def test_matches_pre_batching_implementation(name, workers):
    expected = GOLDEN[name]
    hits = _scan(workers=workers, **expected["kwargs"])
    _assert_same_hits(hits, expected["hits"])


def test_golden_covers_every_mode():
    kwargs = [g["kwargs"] for g in GOLDEN.values()]
    assert {k.get("overlap_strategy", "greedy") for k in kwargs} == {"score", "greedy"}
    assert {k["remove_overlaps"] for k in kwargs} == {True, False}
    assert {k.get("score", True) for k in kwargs} == {True, False}
    assert {k["chunk_size"] for k in kwargs} == {1_000, 1_000_000}
    assert all(GOLDEN[name]["hits"] for name in GOLDEN), "every golden config should have hits"


# ---------------------------------------------------------------------------
# 2. Batch size never changes the result
#    Batches hold about workers * chunk_size bases, while chunk_size alone fixes
#    how records are chunked, so varying workers changes only the batching.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("flags", FLAGS, ids=lambda f: f"{f['overlap_strategy']}-ovl{int(f['remove_overlaps'])}")
@pytest.mark.parametrize("chunk_size", [500, 1_000, 5_000])
def test_batch_size_does_not_change_hits(flags, chunk_size):
    results = [_scan(workers=workers, chunk_size=chunk_size, at_threshold=0.8, filter_homopolymers=True, **flags)
               for workers in (1, 2, 3, 8)]
    assert results[0]
    for other in results[1:]:
        assert other == results[0]


@pytest.mark.parametrize("seed", range(6))
@pytest.mark.parametrize("flags", FLAGS, ids=lambda f: f"{f['overlap_strategy']}-ovl{int(f['remove_overlaps'])}")
def test_batch_size_invariance_random_genomes(tmp_path, seed, flags):
    rng = random.Random(seed)
    sizes = [rng.randint(200, 6_000) for _ in range(rng.randint(1, 12))]
    path = _write_fasta(tmp_path / "g.fa", _records(rng, sizes))
    one_batch = _scan(path, workers=8, chunk_size=20_000, **flags)
    many_batches = _scan(path, workers=1, chunk_size=20_000, **flags)
    assert one_batch == many_batches


# ---------------------------------------------------------------------------
# 3. Batches are really bounded (the point of the change)
# ---------------------------------------------------------------------------

def test_scoring_runs_per_batch_not_per_genome(scoring_calls):
    hits = _scan(workers=1, chunk_size=2_000, overlap_strategy="score", remove_overlaps=True)
    assert len(scoring_calls) > 1, "the genome should be split into several batches"
    assert max(scoring_calls) < sum(scoring_calls)
    _assert_same_hits(hits, GOLDEN["score-overlaps1-filters0-chunk1000000"]["hits"])


def test_one_batch_when_everything_fits(scoring_calls):
    _scan(workers=4, chunk_size=1_000_000, overlap_strategy="score")
    assert len(scoring_calls) == 1


def test_record_larger_than_batch_is_kept_whole(tmp_path, scoring_calls):
    path = _write_fasta(tmp_path / "g.fa", _records(random.Random(3), [30_000]))
    batched = _scan(path, workers=1, chunk_size=1_000, overlap_strategy="score")
    assert len(scoring_calls) == 1, "a record is never split across batches"
    assert batched == _scan(path, workers=8, chunk_size=1_000, overlap_strategy="score")


def test_batch_boundary_exactly_at_record_end(tmp_path, scoring_calls):
    records = _records(random.Random(5), [900, 900, 900, 900])
    path = _write_fasta(tmp_path / "g.fa", records)
    exact = len(records[0][1])  # the first batch closes exactly at the end of record 1
    a = _scan(path, workers=1, chunk_size=exact, overlap_strategy="score")
    assert len(scoring_calls) >= 2
    assert a == _scan(path, workers=4, chunk_size=exact, overlap_strategy="score")


# ---------------------------------------------------------------------------
# 4. Order and content of the output
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("strategy", ["score", "greedy"])
def test_overlap_removal_output_sorted_by_seq_id(strategy):
    hits = _scan(workers=1, chunk_size=1_000, remove_overlaps=True, overlap_strategy=strategy)
    seq_ids = [h["seq_id"] for h in hits]
    assert seq_ids == sorted(seq_ids)
    for seq_id in set(seq_ids):
        starts = [h["start"] for h in hits if h["seq_id"] == seq_id]
        assert starts == sorted(starts)


def test_without_overlap_removal_output_in_file_order():
    hits = _scan(workers=1, chunk_size=1_000, remove_overlaps=False, score=False)
    file_order = [line[1:].split()[0] for line in GENOME.read_text().splitlines() if line.startswith(">")]
    seen = list(dict.fromkeys(h["seq_id"] for h in hits))
    assert seen == [s for s in file_order if s in seen]


@pytest.mark.parametrize("strategy", ["score", "greedy"])
def test_overlap_removal_leaves_no_overlaps(strategy):
    hits = _scan(workers=1, chunk_size=1_000, remove_overlaps=True, overlap_strategy=strategy)
    for seq_id in {h["seq_id"] for h in hits}:
        spans = sorted((h["start"], h["end"]) for h in hits if h["seq_id"] == seq_id)
        for (_, end), (start, _) in zip(spans, spans[1:]):
            assert start > end


def test_every_hit_has_score_fields_in_score_mode():
    hits = _scan(workers=1, chunk_size=1_000, overlap_strategy="score")
    assert hits
    for key in ("stacking_score", "pairing_score", "total_score", "putative_triplex", "seq_id"):
        assert all(key in h for h in hits)


def test_unscored_hits_have_no_score_fields():
    hits = _scan(workers=1, chunk_size=1_000, overlap_strategy="greedy", score=False)
    assert hits and all("total_score" not in h for h in hits)


# ---------------------------------------------------------------------------
# 5. Agreement with the single-sequence API
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("strategy", ["score", "greedy"])
def test_matches_scan_sequence_per_record(tmp_path, strategy):
    records = _records(random.Random(9), [500, 2_000, 4_000, 1_200])
    path = _write_fasta(tmp_path / "g.fa", records)
    kwargs = dict(minrep=10, maxspacer=10, remove_overlaps=True, overlap_strategy=strategy)
    expected = []
    for name, seq in records:
        hits = hseeker.scan_sequence(seq, **kwargs)
        for h in hits:
            h["seq_id"] = name
        expected.extend(hits)
    expected.sort(key=lambda h: h["seq_id"])
    for workers in (1, 4):
        assert hseeker.scan_fasta_parallel(str(path), workers=workers, chunk_size=10_000, **kwargs) == expected


# ---------------------------------------------------------------------------
# 6. Inputs: parsers, gzip, empty and N-only files
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("strategy", ["score", "greedy"])
def test_inhouse_parser_matches_biopython(strategy):
    a = _scan(workers=1, chunk_size=1_000, overlap_strategy=strategy, parser="biopython")
    b = _scan(workers=1, chunk_size=1_000, overlap_strategy=strategy, parser="inhouse")
    assert a == b


def test_gzip_input_matches_plain(tmp_path):
    gz = tmp_path / "genome.fa.gz"
    with open(GENOME, "rb") as src, gzip.open(gz, "wb") as dst:
        shutil.copyfileobj(src, dst)
    assert _scan(gz, workers=1, chunk_size=1_000) == _scan(workers=1, chunk_size=1_000)


def test_empty_fasta(tmp_path):
    path = tmp_path / "empty.fa"
    path.write_text("")
    assert _scan(path, workers=2, chunk_size=1_000) == []


@pytest.mark.parametrize("strategy", ["score", "greedy"])
def test_n_only_records(tmp_path, strategy):
    path = _write_fasta(tmp_path / "n.fa", [("a", "N" * 5_000), ("b", "N" * 300)])
    assert _scan(path, workers=2, chunk_size=1_000, overlap_strategy=strategy) == []


def test_records_without_hits_between_records_with_hits(tmp_path):
    with_hits = _records(random.Random(21), [1_500, 1_500])
    records = [with_hits[0], ("empty1", "N" * 2_000), ("plain", "ACGT" * 500), with_hits[1]]
    path = _write_fasta(tmp_path / "g.fa", records)
    a = _scan(path, workers=1, chunk_size=1_000, overlap_strategy="score")
    assert a == _scan(path, workers=8, chunk_size=1_000, overlap_strategy="score")
    assert a and {h["seq_id"] for h in a} <= {with_hits[0][0], with_hits[1][0]}


def test_score_overlap_without_scoring_still_rejected():
    with pytest.raises(ValueError):
        _scan(workers=1, overlap_strategy="score", remove_overlaps=True, score=False)
