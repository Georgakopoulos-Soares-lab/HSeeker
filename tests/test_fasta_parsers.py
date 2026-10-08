"""FASTA parser selection must preserve IDs, genomic offsets, and scan results."""

import csv
import gzip
import subprocess
import sys

import pytest

import hseeker


FASTA = ">chr1:101-130 description\n" + "GA" * 15 + "\n>plain note\n" + "G" * 30 + "\n"


@pytest.mark.parametrize("compressed", [False, True])
def test_default_biopython_parser_handles_wrapped_mixed_case_and_empty_records(
    tmp_path, monkeypatch, compressed
):
    path = tmp_path / ("wrapped.fa.gz" if compressed else "wrapped.fa")
    fasta = (
        ">empty description\r\n"
        ">chr2:204-233 region description\r\n"
        + ("ga ga ga ga ga\r\n" * 3)
        + ">plain note\r\nnN n\r\n"
    )
    if compressed:
        with gzip.open(path, "wt") as handle:
            handle.write(fasta)
    else:
        path.write_bytes(fasta.encode())

    calls = []
    original = hseeker.SimpleFastaParser

    def tracked_parser(handle):
        calls.append(handle.name)
        yield from original(handle)

    monkeypatch.setattr(hseeker, "SimpleFastaParser", tracked_parser)
    expected = [("chr2:204-233", "GA" * 15, 204), ("plain", "NNN", 1)]
    assert list(hseeker.parse_fasta(path)) == expected
    assert len(calls) == 1

    hits = hseeker.scan_fasta(path, minrep=8, score=False,
                              remove_overlaps=False)
    assert hits
    assert {hit["seq_id"] for hit in hits} == {"chr2:204-233"}
    assert min(hit["start"] for hit in hits) == 204
    assert len(calls) == 2


@pytest.mark.parametrize("compressed", [False, True])
def test_parsers_return_same_records_for_plain_and_gzip(tmp_path, compressed):
    path = tmp_path / ("motifs.fa.gz" if compressed else "motifs.fa")
    if compressed:
        with gzip.open(path, "wt") as handle:
            handle.write(FASTA)
    else:
        path.write_text(FASTA)

    expected = [
        ("chr1:101-130", "GA" * 15, 101),
        ("plain", "G" * 30, 1),
    ]
    for parser in ("biopython", "inhouse"):
        assert list(hseeker.parse_fasta(path, parser=parser)) == expected
        hits = hseeker.scan_fasta(path, parser=parser, minrep=8,
                                  score=False, remove_overlaps=False)
        assert hits
        assert {h["seq_id"] for h in hits} == {"chr1:101-130", "plain"}
        assert min(h["start"] for h in hits if h["seq_id"] == "chr1:101-130") == 101


def test_parallel_and_iter_agree_for_both_parsers(tmp_path):
    path = tmp_path / "motifs.fa"
    path.write_text(FASTA)
    for parser in ("biopython", "inhouse"):
        options = dict(parser=parser, minrep=8, maxrep=20, maxspacer=4,
                       score=False, remove_overlaps=False)
        sequential = hseeker.scan_fasta(path, **options)
        streamed = list(hseeker.scan_fasta_iter(path, **options))
        parallel = hseeker.scan_fasta_parallel(path, workers=2, chunk_size=20,
                                               **options)
        key = lambda h: (h["seq_id"], h["start"], h["end"])
        assert sorted(map(key, sequential)) == sorted(map(key, streamed))
        assert sorted(map(key, sequential)) == sorted(map(key, parallel))


def test_invalid_parser_is_rejected(tmp_path):
    path = tmp_path / "motifs.fa"
    path.write_text(FASTA)
    with pytest.raises(ValueError, match="Invalid FASTA parser"):
        list(hseeker.parse_fasta(path, parser="unknown"))


def test_unscored_pre_filter_agrees_across_fasta_apis(tmp_path):
    path = tmp_path / "motifs.fa"
    path.write_text(FASTA)
    options = dict(minrep=8, score=False, filter_homopolymers=True)
    for hits in (
        hseeker.scan_fasta(path, **options),
        list(hseeker.scan_fasta_iter(path, **options)),
        hseeker.scan_fasta_parallel(path, workers=2, **options),
    ):
        assert hits
        assert {hit["seq_id"] for hit in hits} == {"chr1:101-130"}


@pytest.mark.parametrize("parser", ["biopython", "inhouse"])
def test_cli_parser_choice(tmp_path, parser):
    path = tmp_path / "motifs.fa.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(FASTA)
    prefix = tmp_path / parser
    subprocess.run([
        sys.executable, "-m", "hseeker", "-seq", str(path),
        "-out", str(prefix), "-parser", parser, "-minrep", "8",
    ], check=True, capture_output=True)
    with (tmp_path / f"{parser}_HDNA.tsv").open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert rows
    assert any(row["seq_id"] == "chr1:101-130" for row in rows)
