"""
hseeker — H-DNA / Triplex Mirror Repeat Detector
=================================================

Fast C core wrapped in a clean Python API.

Quick start
-----------
>>> import hseeker
>>> hits = hseeker.scan_sequence("GAGAGAGAGAGAGAGAGAGAGAGAGAGA", minrep=10)
>>> hits = hseeker.scan_fasta("genome.fa", minrep=10, purity=0.90)
"""

from __future__ import annotations

import os
import re
import gzip
from bisect import bisect_right
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Generator
from Bio.SeqIO.FastaIO import SimpleFastaParser

from hseeker import _hdna  # compiled C extension
from hseeker._scoring import score_hit, score_hit_components

__version__: str = "0.1.0"

# One-pass table: strips \n, \r, space AND uppercases a-z in a single .translate()
# call.  This avoids the extra string allocation from .upper() on genome-scale
# sequences (saves ~250 MB peak RSS when reading chr1).
_STRIP_UPPER_TABLE = str.maketrans(
    "abcdefghijklmnopqrstuvwxyz",
    "ABCDEFGHIJKLMNOPQRSTUVWXYZ",
    "\n\r ",
)
__all__ = [
    "scan_sequence",
    "scan_fasta",
    "scan_fasta_iter",
    "scan_fasta_parallel",
    "parse_fasta",
    "score_hit",
    "score_hit_components",
    "ABLATION_MODELS",
    "detect_candidates",
    "filter_at_content",
    "filter_homopolymers",
    "remove_overlaps",
    "score_candidates",
    "scan_sequence_ablation",
    "scan_fasta_ablation",
    "__version__",
]


# ---------------------------------------------------------------------------
# Scoring helper
# ---------------------------------------------------------------------------

_SCORE_KEYS = ("stacking_score", "pairing_score", "total_score", "putative_triplex")


def _filter_at_content(hits: list[dict], threshold: float) -> list[dict]:
    """Drop hits whose left-arm AT content is >= *threshold*.

    Applied before scoring: an arm dominated by A/T bases is unlikely to
    stack into a stable triplex (poly-purine/poly-pyrimidine stacking
    relies on G/A tracts), so these hits are pruned early to skip the
    O(n^2) scoring pass on unpromising candidates. Only the left arm is
    checked — for a mirror repeat the right arm is the same bases in
    reverse order, so its AT content is identical (or within one
    mismatch) and checking both would be redundant.
    """
    kept = []
    for h in hits:
        arm_length = h["arm_length"]
        if arm_length == 0:
            kept.append(h)
            continue
        at_count = sum(1 for b in h["left_arm"] if b in ("a", "t", "A", "T"))
        if (at_count / arm_length) < threshold:
            kept.append(h)
    return kept

class _AcceptedIntervals:
    """Check overlap with accepted inclusive intervals in O(log n) time."""

    def __init__(self, candidate_starts: list[int]):
        self.starts = sorted(set(candidate_starts))
        self.rank = {start: i + 1 for i, start in enumerate(self.starts)}
        self.max_end = [float("-inf")] * (len(self.starts) + 1)

    def overlaps(self, start: int, end: int) -> bool:
        i = bisect_right(self.starts, end)
        while i:
            if self.max_end[i] >= start:
                return True
            i -= i & -i
        return False

    def add(self, start: int, end: int) -> None:
        i = self.rank[start]
        while i < len(self.max_end):
            self.max_end[i] = max(self.max_end[i], end)
            i += i & -i


def _score_priority(hit: dict) -> tuple:
    """Visit stronger hits first, with deterministic tie breakers."""
    score = hit["total_score"]
    # None and NaN have no comparable stability evidence and come last.
    valid = score is not None and score == score
    return (
        0 if valid else 1,
        -score if valid else 0,
        -hit["arm_length"],
        hit["spacer_length"],
        str(hit.get("seq_id", "")),
        hit["start"],
        hit["end"],
        str(hit.get("full_sequence", "")),
    )


def _filter_overlapping_hits(hits: list[dict], strategy: str = "score") -> list[dict]:
    """Select nonoverlapping hits using ``score`` or ``greedy``.

    ``greedy`` preserves the original Python cross-chunk sweep: process hits
    by genomic start, replacing the last kept hit if the new one has a longer
    arm (or an equal arm and shorter spacer). It needs no scores.

    ``score`` (also called ``stability``) visits hits from highest individual
    total_score to lowest and keeps a hit if it overlaps none already kept on
    the same sequence. Ties favor longer arms, shorter spacers, then earlier
    coordinates. None/NaN scores follow numeric scores. This gives a
    deterministic maximal nonoverlapping subset; it does not maximize the
    number of hits or the sum of their scores. Coordinates are 1-based and
    inclusive in both strategies.
    """
    if strategy not in ("greedy", "score", "stability"):
        raise ValueError(f"Invalid deduplication strategy {strategy!r}.")
    if not hits:
        return []
    if strategy == "greedy":
        ordered = sorted(hits, key=lambda h: (
            str(h.get("seq_id", "")), h["start"],
            -h["arm_length"], h["spacer_length"],
        ))
        kept: list[dict] = []
        for hit in ordered:
            if not kept or hit.get("seq_id", "") != kept[-1].get("seq_id", ""):
                kept.append(hit)
            elif hit["start"] > kept[-1]["end"]:
                kept.append(hit)
            elif hit["arm_length"] > kept[-1]["arm_length"] or (
                hit["arm_length"] == kept[-1]["arm_length"]
                and hit["spacer_length"] < kept[-1]["spacer_length"]
            ):
                kept[-1] = hit
        return kept

    if any("total_score" not in hit for hit in hits):
        raise ValueError("Score overlap filtering requires already scored hits")

    by_seq: dict[str, list[dict]] = {}
    for hit in hits:
        by_seq.setdefault(hit.get("seq_id", ""), []).append(hit)

    kept: list[dict] = []
    for seq_hits in by_seq.values():
        accepted = _AcceptedIntervals([hit["start"] for hit in seq_hits])
        for hit in sorted(seq_hits, key=_score_priority):
            if accepted.overlaps(hit["start"], hit["end"]):
                continue
            kept.append(hit)
            accepted.add(hit["start"], hit["end"])

    return sorted(kept, key=lambda hit: (
        str(hit.get("seq_id", "")), hit["start"], hit["end"]
    ))


def _check_overlap_strategy(remove_overlaps: bool, score: bool, strategy: str) -> None:
    """Validate the overlap choice before scanning any sequence."""
    if strategy not in ("greedy", "score", "stability"):
        raise ValueError(f"Invalid deduplication strategy {strategy!r}.")
    if remove_overlaps and strategy in ("score", "stability") and not score:
        raise ValueError("Score overlap filtering requires score=True")


def _filter_homopolymer_triplex(hits: list[dict], mode: str = "POST") -> list[dict]:
    """Drop hits whose motif is a pure homopolymer run.

    ``PRE`` checks the detected ``full_sequence``. ``POST`` also checks the
    scored ``putative_triplex``. A failed/empty score still has its detected
    sequence checked; a mixed-base unscorable candidate is retained.
    """
    def is_homopolymer(seq: str) -> bool:
        if not seq:
            return False
        first = seq[0].lower()
        return all(base.lower() == first for base in seq)

    kept = []
    if mode == "POST":
        for h in hits:
            if is_homopolymer(h["full_sequence"]):
                continue
            seq = h["putative_triplex"].replace("[", "").replace("]", "")
            if is_homopolymer(seq):
                continue
            kept.append(h)
    elif mode == "PRE":
        for h in hits:
            if is_homopolymer(h["full_sequence"]):
                continue
            kept.append(h)
    else:
        raise ValueError(f"Invalid mode `{mode}`.")
    return kept


def _apply_scoring(hits: list[dict]) -> list[dict]:
    """Apply thermodynamic stability scoring to each hit in-place.

    Adds ``stacking_score``, ``pairing_score``, ``total_score``, and
    ``putative_triplex`` keys.  If scoring fails for a hit the fields
    are set to ``None`` / ``""``.
    """
    for h in hits:
        result = score_hit(
            h["left_arm"],  h["spacer"], h["right_arm"], h["arm_length"]
        )
        if result is not None:
            for k in _SCORE_KEYS:
                h[k] = result[k]
        else:
            h["stacking_score"] = None
            h["pairing_score"] = None
            h["total_score"] = None
            h["putative_triplex"] = ""
    return hits


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def scan_sequence(
    seq: str,
    *,
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    remove_overlaps: bool = True,
    overlap_strategy: str = "greedy",
    seq_offset: int = 1,
    score: bool = True,
    purity_rmq: bool = False,
    at_threshold: float | None = None,
    filter_homopolymers: bool = False,
) -> list[dict]:
    """Scan a raw DNA string for H-DNA / triplex mirror repeat motifs.

    The core computation is performed by a compiled C extension and the GIL
    is released during the scan, so other Python threads can continue running.

    Parameters
    ----------
    seq : str
        DNA sequence (ACGTN; case-insensitive).  'N' bases break arm extension.
    minrep : int
        Minimum arm length in bases (default 10).
    maxrep : int
        Maximum arm length in bases (default 1000).
    maxspacer : int
        Maximum spacer between arms in bases (default 10).
    purity : float
        Minimum fraction of GA or CT bases required in each arm (default 0.90).
        Set to 1.0 to reproduce strict non-B_gfa mirror-repeat results.
    mismatch : float
        Maximum fraction of mirror-position mismatches allowed (default 0.10).
        Set to 0.0 for exact mirror only.
    remove_overlaps : bool
        Remove overlapping hits (default True). Set False to keep all hits.
    overlap_strategy : {"greedy", "score", "stability"}
        ``greedy`` (default) keeps the longest arm before scoring, preserving
        historical behavior. ``score`` scores candidates first, then keeps
        each hit in descending individual ``total_score`` order if it does
        not overlap a previously kept hit. ``stability`` aliases ``score``.
        Score-based overlap removal requires score=True.
    seq_offset : int
        1-based genomic start coordinate of ``seq[0]`` (default 1).
        Pass the chromosomal start when ``seq`` is a genomic slice.
    score : bool
        Apply thermodynamic stability scoring to each hit (default True).
        Adds ``stacking_score``, ``pairing_score``, ``total_score``, and
        ``putative_triplex`` keys.
    purity_rmq : bool
        Use an exact right-arm purity feasibility prefilter before extension
        (default False). This preserves output semantics while skipping
        center/spacer pairs that cannot satisfy the purity rule.
    at_threshold : float | None
        Drop hits whose left-arm AT content is >= this value (default None,
        i.e. no filtering). Applied before scoring. AT-rich arms are
        unlikely to form stable H-DNA triplexes.
    filter_homopolymers : bool
        Drop homopolymer hits (default False). Always inspect the detected
        ``full_sequence``; with scoring, also inspect ``putative_triplex``.

    Returns
    -------
    list[dict]
        Each dict has keys: ``start``, ``end``, ``arm_length``,
        ``spacer_length``, ``total_length``, ``ga_pct``, ``ct_pct``,
        ``mirror_identity``, ``is_perfect``, ``left_arm``, ``spacer``,
        ``right_arm``, ``full_sequence``, and (if ``score=True``)
        ``stacking_score``, ``pairing_score``, ``total_score``,
        ``putative_triplex``.
        Coordinates are 1-based and inclusive.
    """
    _check_overlap_strategy(remove_overlaps, score, overlap_strategy)
    score_overlap = remove_overlaps and overlap_strategy in ("score", "stability")
    hits = _hdna.scan_sequence(
        seq,
        minrep=minrep,
        maxrep=maxrep,
        maxspacer=maxspacer,
        purity=purity,
        mismatch=mismatch,
        remove_overlaps=remove_overlaps and not score_overlap,
        seq_offset=seq_offset,
        purity_rmq=purity_rmq,
    )
    if at_threshold is not None and hits:
        hits = _filter_at_content(hits, at_threshold)
    if score and hits:
        _apply_scoring(hits)
        if filter_homopolymers:
            hits = _filter_homopolymer_triplex(hits, mode="POST")
    elif filter_homopolymers:
        hits = _filter_homopolymer_triplex(hits, mode="PRE")
    if score_overlap and hits:
        hits = _filter_overlapping_hits(hits, strategy=overlap_strategy)
    return hits


def parse_fasta(
    path: str | Path, *, parser: str = "biopython"
) -> Generator[tuple[str, str, int], None, None]:
    """Yield ``(seq_id, sequence, offset)`` tuples from a FASTA file.

    Handles multi-record FASTA files.  The genomic offset is parsed from
    headers of the form ``>seqid:start-end`` (UCSC / Ensembl region
    extracts). All other headers default to offset 1.

    Parameters
    ----------
    path : str | Path
        Path to a plain-text or gzip-compressed FASTA file.
    parser : {"biopython", "inhouse"}
        FASTA parser implementation (default "biopython"). Both yield the
        same record format, including offsets from region headers.

    Yields
    ------
    tuple[str, str, int]
        ``(seq_id, sequence, offset)`` where ``sequence`` is uppercase
        and ``offset`` is the 1-based genomic start of the first base.
    """
    if parser not in ("biopython", "inhouse"):
        raise ValueError(f"Invalid FASTA parser {parser!r}; choose 'biopython' or 'inhouse'")

    _header_re = re.compile(r"^(\S+)")
    _offset_re = re.compile(r":(\d+)[-–]")

    def record_info(header: str) -> tuple[str, int]:
        match = _header_re.match(header)
        seq_id = match.group(1) if match else "unknown"
        region = _offset_re.search(header)
        return seq_id, int(region.group(1)) if region else 1

    if parser == "biopython":
        opener = gzip.open if str(path).endswith(".gz") else open
        with opener(path, "rt", encoding="utf-8") as handle:
            for title, sequence in SimpleFastaParser(handle):
                if sequence:
                    seq_id, offset = record_info(title)
                    yield seq_id, sequence.upper(), offset
        return

    # Read the entire file in one shot.  For genome-scale FASTA files
    # (3–4 GB) this requires sufficient RAM, but avoids creating millions
    # of per-line Python string objects — the previous line-by-line
    # approach was a ~3× bottleneck on chr1.
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            raw = handle.read()
    else:
        raw = Path(path).read_text()

    # Split on '>' to get raw record blocks; the first (empty) block
    # before the initial '>' is discarded.  Free `raw` immediately after
    # the split so we don't hold both the original and the split copies.
    blocks = raw.split(">")
    del raw

    for block in blocks:
        if not block.strip():
            continue
        # First line is the header (up to first newline)
        newline_idx = block.find("\n")
        if newline_idx == -1:
            header = block
            seq_lines = ""
        else:
            header = block[:newline_idx]
            seq_lines = block[newline_idx + 1:]

        seq_id, offset = record_info(header)

        # Single-pass: strip \n/\r/space and uppercase in one .translate()
        # call.  Avoids the extra ~250 MB allocation that .upper() would create
        # on top of the already-translated string for chromosome-scale inputs.
        seq = seq_lines.translate(_STRIP_UPPER_TABLE)

        if seq:
            yield seq_id, seq, offset


def scan_fasta(
    path: str | Path,
    *,
    parser: str = "biopython",
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    remove_overlaps: bool = True,
    overlap_strategy: str = "greedy",
    score: bool = True,
    purity_rmq: bool = False,
    at_threshold: float | None = None,
    filter_homopolymers: bool = False,
) -> list[dict]:
    """Scan every record in a FASTA file for H-DNA motifs.

    This is the high-level convenience function.  It calls
    :func:`parse_fasta` to iterate records and :func:`scan_sequence` on
    each one, then returns all hits in a single flat list.

    Each returned dict includes an additional ``seq_id`` key identifying
    which FASTA record the hit belongs to.

    Parameters
    ----------
    path : str | Path
        Path to a FASTA file (may contain multiple records).
    parser : {"biopython", "inhouse"}
        FASTA parser implementation (default "biopython").
    minrep, maxrep, maxspacer, purity, mismatch, remove_overlaps,
    overlap_strategy, score,
    at_threshold, filter_homopolymers :
        Same as :func:`scan_sequence`.

    Returns
    -------
    list[dict]
        All hits from all records.  Each dict has the same keys as
        :func:`scan_sequence` plus ``seq_id``.
    """
    _check_overlap_strategy(remove_overlaps, score, overlap_strategy)
    score_overlap = remove_overlaps and overlap_strategy in ("score", "stability")
    results: list[dict] = []
    for seq_id, seq, offset in parse_fasta(path, parser=parser):
        hits = scan_sequence(
            seq,
            minrep=minrep,
            maxrep=maxrep,
            maxspacer=maxspacer,
            purity=purity,
            mismatch=mismatch,
            remove_overlaps=remove_overlaps and not score_overlap,
            seq_offset=offset,
            score=False,  # defer scoring to collect phase
            purity_rmq=purity_rmq,
            at_threshold=at_threshold,
        )
        for h in hits:
            h["seq_id"] = seq_id
        results.extend(hits)
    if score and results:
        _apply_scoring(results)
        if filter_homopolymers:
            results = _filter_homopolymer_triplex(results)
    elif results and filter_homopolymers:
        results = _filter_homopolymer_triplex(results, mode="PRE")
    if score_overlap and results:
        results = _filter_overlapping_hits(results, strategy=overlap_strategy)
    return results


def scan_fasta_iter(
    path: str | Path,
    *,
    parser: str = "biopython",
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    remove_overlaps: bool = True,
    overlap_strategy: str = "greedy",
    score: bool = True,
    purity_rmq: bool = False,
    at_threshold: float | None = None,
    filter_homopolymers: bool = False,
) -> Generator[dict, None, None]:
    """Scan a FASTA file and yield hits one at a time (streaming).

    Identical to :func:`scan_fasta` but returns a generator instead of a
    list, so peak RAM is proportional to the largest single record rather
    than the total number of hits.

    Yields
    ------
    dict
        Same keys as :func:`scan_fasta`.
    """
    _check_overlap_strategy(remove_overlaps, score, overlap_strategy)
    for seq_id, seq, offset in parse_fasta(path, parser=parser):
        hits = scan_sequence(
            seq,
            minrep=minrep,
            maxrep=maxrep,
            maxspacer=maxspacer,
            purity=purity,
            mismatch=mismatch,
            remove_overlaps=remove_overlaps,
            overlap_strategy=overlap_strategy,
            seq_offset=offset,
            score=score,  # score per-record for streaming
            purity_rmq=purity_rmq,
            at_threshold=at_threshold,
            filter_homopolymers=filter_homopolymers,
        )
        for h in hits:
            h["seq_id"] = seq_id
            yield h


def scan_fasta_parallel(
    path: str | Path,
    *,
    parser: str = "biopython",
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    remove_overlaps: bool = True,
    overlap_strategy: str = "greedy",
    workers: int | None = None,
    chunk_size: int = 1_000_000,
    score: bool = True,
    purity_rmq: bool = False,
    at_threshold: float | None = None,
    filter_homopolymers: bool = False,
) -> list[dict]:
    """Scan a FASTA file using a thread pool with intra-record chunk parallelism.

    Records shorter than *chunk_size* are scanned as a single task.  Longer
    records are split into overlapping chunks and scanned in parallel; hits
    that fall in the overlap zone are assigned exclusively to the chunk that
    owns them, so no hit is counted twice.

    The C extension releases the GIL during each scan, so true CPU parallelism
    is achieved across chunks and records.

    Parameters
    ----------
    workers : int | None
        Number of worker threads (default: ``os.cpu_count()``).
    chunk_size : int
        Target bases per chunk for large records (default 1,000,000).

    Returns
    -------
    list[dict]
        All hits from all records.  Each dict has the same keys as
        :func:`scan_sequence` plus ``seq_id``.
    """
    _check_overlap_strategy(remove_overlaps, score, overlap_strategy)
    score_overlap = remove_overlaps and overlap_strategy in ("score", "stability")
    n_workers = workers if workers is not None else (os.cpu_count() or 1)

    # Overlap must cover the largest possible hit that can straddle a chunk
    # boundary: left-arm start just before the boundary, right-arm end at
    # most maxrep bases into the next chunk.  Formula: 2*maxrep + maxspacer + 1.
    overlap = 2 * maxrep + maxspacer + 1

    def _build_tasks(
        seq_id: str, seq: str, genomic_offset: int
    ) -> list[tuple[str, str, int, int | None, bool]]:
        """Split one record into (seq_id, chunk_seq, chunk_offset, excl_end, defer_ovl) tasks.

        *excl_end* is the exclusive upper bound on the genomic start coordinate
        of hits that this chunk owns.  Hits with start >= excl_end are in the
        overlap region and will be found (and owned) by the next chunk.
        *defer_ovl* is True when overlap removal should be deferred to the
        global cross-chunk pass (all non-last chunks, plus the last chunk of
        any multi-chunk record).  Single-chunk records use the fast C-extension
        overlap removal."""
        seqlen = len(seq)
        if seqlen <= chunk_size:
            return [(seq_id, seq, genomic_offset, None, False)]

        tasks: list[tuple[str, str, int, int | None, bool]] = []
        pos = 0
        while pos < seqlen:
            is_last = (pos + chunk_size >= seqlen)
            if is_last:
                # Last chunk of a multi-chunk record: defer overlap removal
                # to the global pass so that secondary hits that survive
                # cross-chunk overlap resolution are not lost.
                tasks.append((seq_id, seq[pos:], genomic_offset + pos, None, True))
            else:
                chunk_seq = seq[pos : pos + chunk_size + overlap]
                excl_end = genomic_offset + pos + chunk_size
                tasks.append((seq_id, chunk_seq, genomic_offset + pos, excl_end, True))
            pos += chunk_size
        return tasks

    def _scan_chunk(
        task: tuple[str, str, int, int | None, bool],
    ) -> list[dict]:
        seq_id, chunk_seq, chunk_offset, excl_end, defer_ovl = task
        # Defer overlap removal to the global pass for all chunks except
        # single-chunk records (those with defer_ovl=False).  Also defer
        # scoring — each chunk may produce thousands of hits and scoring
        # them individually per chunk (especially with ProcessPoolExecutor
        # inside ThreadPoolExecutor threads) is very slow.  Both overlap
        # removal and scoring are applied once on the merged result set.
        hits = scan_sequence(
            chunk_seq,
            minrep=minrep,
            maxrep=maxrep,
            maxspacer=maxspacer,
            purity=purity,
            mismatch=mismatch,
            remove_overlaps=(remove_overlaps and not score_overlap and not defer_ovl),
            seq_offset=chunk_offset,
            score=False,
            purity_rmq=purity_rmq,
            at_threshold=at_threshold,
        )
        for h in hits:
            h["seq_id"] = seq_id
        if excl_end is not None:
            hits = [h for h in hits if h["start"] < excl_end]
        return hits

    # Flatten all records into a list of chunk tasks.
    # Each chunk_seq is a new str copy of a slice of the full sequence.
    # Explicitly delete the full-sequence loop variable once chunking is done
    # so the original full-sequence string can be freed before the scan starts.
    all_tasks: list[tuple[str, str, int, int | None, bool]] = []
    seq: str = ""  # pre-init so del below is safe even when FASTA has no records
    for seq_id, seq, offset in parse_fasta(path, parser=parser):
        all_tasks.extend(_build_tasks(seq_id, seq, offset))
    del seq  # free the last full-sequence string (or the "" sentinel)

    results: list[dict] = []
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        for hits in executor.map(_scan_chunk, all_tasks):
            results.extend(hits)
    del all_tasks  # free chunk strings now that all tasks have completed

    # In greedy mode, preserve the original pre-scoring overlap pass.
    # Score mode keeps raw candidates until their scores are available.
    if remove_overlaps and not score_overlap and results:
        results = _filter_overlapping_hits(results, strategy="greedy")

    if score and results:
        _apply_scoring(results)
        if filter_homopolymers:
            results = _filter_homopolymer_triplex(results, mode="POST")
    elif results and filter_homopolymers:
        results = _filter_homopolymer_triplex(results, mode="PRE")
    if score_overlap and results:
        results = _filter_overlapping_hits(results, strategy=overlap_strategy)
    return results


# Import these after the scanning functions: ablation composes them without
# duplicating scanner or scoring logic.
from hseeker.ablation import (  # noqa: E402
    MODELS as ABLATION_MODELS,
    detect_candidates,
    filter_at_content,
    filter_homopolymers,
    remove_overlaps,
    score_candidates,
    scan_sequence_ablation,
    scan_fasta_ablation,
)
