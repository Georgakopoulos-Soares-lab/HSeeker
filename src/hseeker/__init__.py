"""
hseeker — H-DNA / Triplex Mirror Repeat Detector
=================================================

Fast C core wrapped in a clean Python API.

Quick start
-----------
>>> import hseeker
>>> hits = hseeker.scan_sequence("GAGAGAGAGAGAGAGAGAGAGAGAGAGA", minrep=6)
>>> hits = hseeker.scan_fasta("genome.fa", minrep=10, purity=0.85)
>>> # memory-efficient streaming alternative:
>>> for hit in hseeker.scan_fasta_iter("genome.fa", minrep=10):
...     print(hit["start"], hit["end"])
"""

from __future__ import annotations

import os
import re
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from typing import Generator

from hseeker import _hdna  # compiled C extension

__version__: str = "0.1.0"
__all__ = ["scan_sequence", "scan_fasta", "scan_fasta_iter", "scan_fasta_parallel", "parse_fasta", "__version__"]


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def scan_sequence(
    seq: str,
    *,
    minrep: int = 6,
    maxrep: int = 50,
    maxspacer: int = 7,
    purity: float = 0.80,
    mismatch: float = 0.20,
    remove_overlaps: bool = True,
    seq_offset: int = 1,
) -> list[dict]:
    """Scan a raw DNA string for H-DNA / triplex mirror repeat motifs.

    The core computation is performed by a compiled C extension and the GIL
    is released during the scan, so other Python threads can continue running.

    Parameters
    ----------
    seq : str
        DNA sequence (ACGTN; case-insensitive).  'N' bases break arm extension.
    minrep : int
        Minimum arm length in bases (default 6).
    maxrep : int
        Maximum arm length in bases (default 50).
    maxspacer : int
        Maximum spacer between arms in bases (default 7).
    purity : float
        Minimum fraction of GA or CT bases required in each arm (default 0.80).
        Set to 1.0 to reproduce strict non-B_gfa mirror-repeat results.
    mismatch : float
        Maximum fraction of mirror-position mismatches allowed (default 0.20).
        Set to 0.0 for exact mirror only.
    remove_overlaps : bool
        Remove overlapping hits, keeping the longest arm (default True).
    seq_offset : int
        1-based genomic start coordinate of ``seq[0]`` (default 1).
        Pass the chromosomal start when ``seq`` is a genomic slice.

    Returns
    -------
    list[dict]
        Each dict has keys: ``start``, ``end``, ``arm_length``,
        ``spacer_length``, ``total_length``, ``ga_pct``, ``ct_pct``,
        ``mirror_identity``, ``is_perfect``, ``left_arm``, ``spacer``,
        ``right_arm``, ``full_sequence``.
        Coordinates are 1-based and inclusive.
    """
    return _hdna.scan_sequence(
        seq,
        minrep=minrep,
        maxrep=maxrep,
        maxspacer=maxspacer,
        purity=purity,
        mismatch=mismatch,
        remove_overlaps=remove_overlaps,
        seq_offset=seq_offset,
    )


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

    seq_id: str | None = None
    offset: int = 1
    parts: list[str] = []

    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n\r")
            if line.startswith(">"):
                if seq_id is not None:
                    yield seq_id, "".join(parts), offset
                m = _header_re.match(line)
                seq_id = m.group(1) if m else "unknown"
                om = _offset_re.search(line)
                offset = int(om.group(1)) if om else 1
                parts = []
            elif line:
                parts.append(line.upper().replace(" ", ""))

    if seq_id is not None:
        yield seq_id, "".join(parts), offset


def scan_fasta_iter(
    path: str | Path,
    *,
    minrep: int = 6,
    maxrep: int = 50,
    maxspacer: int = 7,
    purity: float = 0.80,
    mismatch: float = 0.20,
    remove_overlaps: bool = True,
) -> Generator[dict, None, None]:
    """Stream H-DNA hits from every record in a FASTA file one at a time.

    This is the memory-efficient alternative to :func:`scan_fasta`.  Only
    a single record's worth of hits is held in RAM at any moment.  For a
    24-chromosome genome this reduces peak memory by up to 24× compared
    to :func:`scan_fasta`.

    Parameters
    ----------
    path : str | Path
        Path to a FASTA file (may contain multiple records).
    minrep, maxrep, maxspacer, purity, mismatch, remove_overlaps :
        Same as :func:`scan_sequence`.

    Yields
    ------
    dict
        One hit dict per yield, with the same keys as
        :func:`scan_sequence` plus ``seq_id``.
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
        )
        for h in hits:
            h["seq_id"] = seq_id
            yield h


def scan_fasta(
    path: str | Path,
    *,
    minrep: int = 6,
    maxrep: int = 50,
    maxspacer: int = 7,
    purity: float = 0.80,
    mismatch: float = 0.20,
    remove_overlaps: bool = True,
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
    minrep, maxrep, maxspacer, purity, mismatch, remove_overlaps :
        Same as :func:`scan_sequence`.

    Returns
    -------
    list[dict]
        All hits from all records.  Each dict has the same keys as
        :func:`scan_sequence` plus ``seq_id``.
    """
    return list(
        scan_fasta_iter(
            path,
            minrep=minrep,
            maxrep=maxrep,
            maxspacer=maxspacer,
            purity=purity,
            mismatch=mismatch,
            remove_overlaps=remove_overlaps,
        )
    )


def scan_fasta_parallel(
    path: str | Path,
    *,
    minrep: int = 6,
    maxrep: int = 50,
    maxspacer: int = 7,
    purity: float = 0.80,
    mismatch: float = 0.20,
    remove_overlaps: bool = True,
    workers: int | None = None,
) -> list[dict]:
    """Scan all FASTA records in parallel using threads.

    The C scan releases the GIL (``Py_BEGIN_ALLOW_THREADS`` in ``_hdna.c``),
    so threads achieve true parallelism on the heavy C work.  For a
    24-chromosome genome on a 12-core machine this is roughly 10× faster
    than the sequential :func:`scan_fasta`.

    Parameters
    ----------
    path : str | Path
        Path to a FASTA file (may contain multiple records).
    minrep, maxrep, maxspacer, purity, mismatch, remove_overlaps :
        Same as :func:`scan_sequence`.
    workers : int | None
        Number of worker threads.  ``None`` (default) uses
        ``os.cpu_count()``.

    Returns
    -------
    list[dict]
        All hits from all records in record order.  Each dict has the same
        keys as :func:`scan_fasta`.
    """
    records = list(parse_fasta(path))
    if not records:
        return []

    n_workers = workers if workers is not None else (os.cpu_count() or 1)

    def _scan_record(
        record: tuple[str, str, int],
    ) -> list[dict]:
        seq_id, seq, offset = record
        hits = scan_sequence(
            seq,
            minrep=minrep,
            maxrep=maxrep,
            maxspacer=maxspacer,
            purity=purity,
            mismatch=mismatch,
            remove_overlaps=remove_overlaps,
            seq_offset=offset,
        )
        for h in hits:
            h["seq_id"] = seq_id
        return hits

    results: list[dict] = []
    with ThreadPoolExecutor(max_workers=n_workers) as pool:
        # submit all records; collect futures in submission order to
        # preserve deterministic record ordering in the output.
        futures = [pool.submit(_scan_record, r) for r in records]
        for fut in futures:
            results.extend(fut.result())

    return results
