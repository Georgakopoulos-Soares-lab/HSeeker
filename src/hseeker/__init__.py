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
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Generator

from hseeker import _hdna  # compiled C extension
from hseeker._scoring import score_hit

__version__: str = "0.1.0"
__all__ = [
    "search",
    "search_fast_gated",
    "search_safe_pruned",
    "scan_sequence",
    "scan_fasta",
    "scan_fasta_iter",
    "scan_fasta_parallel",
    "parse_fasta",
    "profiling_info",
    "score_hit",
    "__version__",
]


# ---------------------------------------------------------------------------
# Scoring helper
# ---------------------------------------------------------------------------

_SCORE_KEYS = ("stacking_score", "pairing_score", "total_score", "putative_triplex")


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
    seq_offset: int = 1,
    score: bool = True,
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
        Remove overlapping hits, keeping the longest arm (default True).
    seq_offset : int
        1-based genomic start coordinate of ``seq[0]`` (default 1).
        Pass the chromosomal start when ``seq`` is a genomic slice.
    score : bool
        Apply thermodynamic stability scoring to each hit (default True).
        Adds ``stacking_score``, ``pairing_score``, ``total_score``, and
        ``putative_triplex`` keys.

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
    hits = _hdna.scan_sequence(
        seq,
        minrep=minrep,
        maxrep=maxrep,
        maxspacer=maxspacer,
        purity=purity,
        mismatch=mismatch,
        remove_overlaps=remove_overlaps,
        seq_offset=seq_offset,
    )
    if score and hits:
        _apply_scoring(hits)
    return hits


def search_fast_gated(
    seq: str,
    *,
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    gate_window: int | None = None,
    gate_search_limit: int | None = None,
    gate_mirror_frac: float | None = None,
    gate_purity_frac: float | None = None,
    fast_mode: bool = True,
    remove_overlaps: bool = True,
    use_purity_prefilter: bool = True,
    seq_offset: int = 1,
    score: bool = True,
) -> list[dict]:
    """Scan a DNA string with the FastGatedHSeeker detector.

    FastGatedHSeeker first applies a local density gate to each candidate
    center/spacer pair, then runs the original HSeeker-compatible extension
    only for pairs that pass the gate. ``fast_mode=False`` makes the default
    ``gate_search_limit`` equal to ``maxrep`` for exact-gate style runs.
    """
    hits = _hdna.scan_sequence_fast_gated(
        seq,
        minrep=minrep,
        maxrep=maxrep,
        maxspacer=maxspacer,
        purity=purity,
        mismatch=mismatch,
        gate_window=gate_window or 0,
        gate_search_limit=gate_search_limit or 0,
        gate_mirror_frac=-1.0 if gate_mirror_frac is None else gate_mirror_frac,
        gate_purity_frac=-1.0 if gate_purity_frac is None else gate_purity_frac,
        fast_mode=bool(fast_mode),
        remove_overlaps=remove_overlaps,
        use_purity_prefilter=bool(use_purity_prefilter),
        seq_offset=seq_offset,
    )
    if score and hits:
        _apply_scoring(hits)
    return hits


def search_safe_pruned(
    seq: str,
    *,
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    remove_overlaps: bool = True,
    safe_prune_purity: bool = False,
    seq_offset: int = 1,
    score: bool = True,
) -> list[dict]:
    """Scan a DNA string with the SafePrunedHSeeker detector.

    SafePrunedHSeeker keeps the original prefix-validity checks, but stops
    extending a center/spacer pair when the current mismatch count, and
    optionally right-arm purity, cannot mathematically recover before Kmax.
    The conservative default is mismatch-only pruning.
    """
    hits = _hdna.scan_sequence_safe_pruned(
        seq,
        minrep=minrep,
        maxrep=maxrep,
        maxspacer=maxspacer,
        purity=purity,
        mismatch=mismatch,
        remove_overlaps=remove_overlaps,
        safe_prune_purity=bool(safe_prune_purity),
        seq_offset=seq_offset,
    )
    if score and hits:
        _apply_scoring(hits)
    return hits


def search(
    seq: str,
    *,
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    detector: str = "original",
    gate_window: int | None = None,
    gate_search_limit: int | None = None,
    gate_mirror_frac: float | None = None,
    gate_purity_frac: float | None = None,
    fast_mode: bool = True,
    remove_overlaps: bool = True,
    use_purity_prefilter: bool = True,
    safe_prune_purity: bool = False,
    seq_offset: int = 1,
    score: bool = True,
) -> list[dict]:
    """Unified sequence search API for original and optimized detectors."""
    if detector == "original":
        return scan_sequence(
            seq,
            minrep=minrep,
            maxrep=maxrep,
            maxspacer=maxspacer,
            purity=purity,
            mismatch=mismatch,
            remove_overlaps=remove_overlaps,
            seq_offset=seq_offset,
            score=score,
        )
    if detector in {"fast_gated", "FastGatedHSeeker"}:
        return search_fast_gated(
            seq,
            minrep=minrep,
            maxrep=maxrep,
            maxspacer=maxspacer,
            purity=purity,
            mismatch=mismatch,
            gate_window=gate_window,
            gate_search_limit=gate_search_limit,
            gate_mirror_frac=gate_mirror_frac,
            gate_purity_frac=gate_purity_frac,
            fast_mode=fast_mode,
            remove_overlaps=remove_overlaps,
            use_purity_prefilter=use_purity_prefilter,
            seq_offset=seq_offset,
            score=score,
        )
    if detector in {"safe_pruned", "SafePrunedHSeeker"}:
        return search_safe_pruned(
            seq,
            minrep=minrep,
            maxrep=maxrep,
            maxspacer=maxspacer,
            purity=purity,
            mismatch=mismatch,
            remove_overlaps=remove_overlaps,
            safe_prune_purity=safe_prune_purity,
            seq_offset=seq_offset,
            score=score,
        )
    raise ValueError("detector must be 'original', 'fast_gated', or 'safe_pruned'")


def profiling_info() -> dict:
    """Return low-level detector profiling counters for the last scan."""
    return _hdna.profiling_info()


def parse_fasta(path: str | Path) -> Generator[tuple[str, str, int], None, None]:
    """Yield ``(seq_id, sequence, offset)`` tuples from a FASTA file.

    Handles multi-record FASTA files.  The genomic offset is parsed from
    headers of the form ``>seqid:start-end`` (UCSC / Ensembl region
    extracts). All other headers default to offset 1.

    Parameters
    ----------
    path : str | Path
        Path to a FASTA file (plain-text, not gzipped).

    Yields
    ------
    tuple[str, str, int]
        ``(seq_id, sequence, offset)`` where ``sequence`` is uppercase
        and ``offset`` is the 1-based genomic start of the first base.
    """
    _header_re = re.compile(r"^>(\S+)")
    _offset_re = re.compile(r":(\d+)[-–]")

    # Read the entire file in one shot.  For genome-scale FASTA files
    # (3–4 GB) this requires sufficient RAM, but avoids creating millions
    # of per-line Python string objects — the previous line-by-line
    # approach was a ~3× bottleneck on chr1.
    raw = Path(path).read_text()

    # Split on '>' to get raw record blocks; the first (empty) block
    # before the initial '>' is discarded.
    blocks = raw.split(">")
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

        m = _header_re.match(">" + header)
        seq_id = m.group(1) if m else "unknown"
        om = _offset_re.search(">" + header)
        offset = int(om.group(1)) if om else 1

        # Collapse all whitespace/newlines from the sequence block,
        # then uppercase.  This is a single string operation instead
        # of per-line allocations.
        seq = seq_lines.translate(
            {10: None, 13: None, 32: None}  # remove \n, \r, space
        ).upper()

        if seq:
            yield seq_id, seq, offset


def scan_fasta(
    path: str | Path,
    *,
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    remove_overlaps: bool = True,
    score: bool = True,
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
    minrep, maxrep, maxspacer, purity, mismatch, remove_overlaps, score :
        Same as :func:`scan_sequence`.

    Returns
    -------
    list[dict]
        All hits from all records.  Each dict has the same keys as
        :func:`scan_sequence` plus ``seq_id``.
    """
    results: list[dict] = []
    for seq_id, seq, offset in parse_fasta(path):
        hits = scan_sequence(
            seq,
            minrep=minrep,
            maxrep=maxrep,
            maxspacer=maxspacer,
            purity=purity,
            mismatch=mismatch,
            remove_overlaps=remove_overlaps,
            seq_offset=offset,
            score=False,  # defer scoring to collect phase
        )
        for h in hits:
            h["seq_id"] = seq_id
        results.extend(hits)
    if score and results:
        _apply_scoring(results)
    return results


def scan_fasta_iter(
    path: str | Path,
    *,
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    remove_overlaps: bool = True,
    score: bool = True,
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
    for seq_id, seq, offset in parse_fasta(path):
        hits = scan_sequence(
            seq,
            minrep=minrep,
            maxrep=maxrep,
            maxspacer=maxspacer,
            purity=purity,
            mismatch=mismatch,
            remove_overlaps=remove_overlaps,
            seq_offset=offset,
            score=score,  # score per-record for streaming
        )
        for h in hits:
            h["seq_id"] = seq_id
            yield h


def scan_fasta_parallel(
    path: str | Path,
    *,
    minrep: int = 10,
    maxrep: int = 1000,
    maxspacer: int = 10,
    purity: float = 0.90,
    mismatch: float = 0.10,
    remove_overlaps: bool = True,
    workers: int | None = None,
    chunk_size: int = 1_000_000,
    score: bool = True,
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
            remove_overlaps=(remove_overlaps and not defer_ovl),
            seq_offset=chunk_offset,
            score=False,
        )
        for h in hits:
            h["seq_id"] = seq_id
        if excl_end is not None:
            hits = [h for h in hits if h["start"] < excl_end]
        return hits

    # Flatten all records into a list of chunk tasks
    all_tasks: list[tuple[str, str, int, int | None, bool]] = []
    for seq_id, seq, offset in parse_fasta(path):
        all_tasks.extend(_build_tasks(seq_id, seq, offset))

    results: list[dict] = []
    with ThreadPoolExecutor(max_workers=n_workers) as executor:
        for hits in executor.map(_scan_chunk, all_tasks):
            results.extend(hits)

    # Cross-chunk overlap removal: each chunk's overlap removal is
    # independent, so hits from different chunks of the same record
    # can still overlap.  Use a sweep-line pass (O(n log n) sort +
    # O(n) linear scan) instead of the previous O(n²) nested-loop
    # approach which could not complete on chromosome-scale data.
    if remove_overlaps and results:
        # Sort by seq_id, then by start coordinate ascending.
        # Among hits that start at the same position, keep the
        # longest arm first (shortest spacer as tiebreak).
        results.sort(key=lambda h: (
            h["seq_id"],
            h["start"],
            -h["arm_length"],
            h["spacer_length"],
        ))
        kept: list[dict] = []
        for h in results:
            if not kept or h["seq_id"] != kept[-1]["seq_id"]:
                # First hit in this seq_id — always keep
                kept.append(h)
            elif h["start"] > kept[-1]["end"]:
                # No overlap with the last kept hit — keep
                kept.append(h)
            elif h["arm_length"] > kept[-1]["arm_length"] or (
                h["arm_length"] == kept[-1]["arm_length"]
                and h["spacer_length"] < kept[-1]["spacer_length"]
            ):
                # Overlaps but is better — replace
                kept[-1] = h
            # else: overlaps and is not better — discard
        results = kept

    if score and results:
        _apply_scoring(results)

    return results
