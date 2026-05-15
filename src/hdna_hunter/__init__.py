"""
hdna_hunter — H-DNA / Triplex Mirror Repeat Detector
=====================================================

Fast C core wrapped in a clean Python API.

Quick start
-----------
>>> import hdna_hunter
>>> hits = hdna_hunter.scan_sequence("GAGAGAGAGAGAGAGAGAGAGAGAGAGA", minrep=6)
>>> hits = hdna_hunter.scan_fasta("genome.fa", minrep=10, purity=0.85)
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Generator

from hdna_hunter import _hdna  # compiled C extension

__version__: str = "0.1.0"
__all__ = ["scan_sequence", "scan_fasta", "parse_fasta", "__version__"]


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
        )
        for h in hits:
            h["seq_id"] = seq_id
        results.extend(hits)
    return results
