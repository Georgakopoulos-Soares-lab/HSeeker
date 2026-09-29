"""Explicit HSeeker component ablations for model evaluation.

Use the same input sequences and evaluation labels for every model. A scan
returns candidates; detection metrics and score-based ranking metrics are
different quantities and should be reported separately.
"""

from __future__ import annotations

from pathlib import Path

from hseeker import (
    _filter_at_content,
    _filter_homopolymer_triplex,
    _filter_overlapping_hits,
    parse_fasta,
    scan_sequence,
)
from hseeker._scoring import score_hit_components

MODELS = ("detection", "composition", "pairing", "pairing_stacking", "complete")

__all__ = (
    "MODELS",
    "detect_candidates",
    "filter_at_content",
    "filter_homopolymers",
    "remove_overlaps",
    "score_candidates",
    "scan_sequence_ablation",
    "scan_fasta_ablation",
)


def detect_candidates(
    seq: str, *, minrep: int = 10, maxrep: int = 1000,
    maxspacer: int = 10, purity: float = 0.90,
    mismatch: float = 0.10, seq_offset: int = 1,
    purity_rmq: bool = False,
) -> list[dict]:
    """Run the C detector without scoring or optional Python filters.

    All overlapping candidates are retained. Detection still applies the
    requested C purity, mismatch, arm-length, and spacer constraints.
    """
    return scan_sequence(
        seq, minrep=minrep, maxrep=maxrep, maxspacer=maxspacer,
        purity=purity, mismatch=mismatch, seq_offset=seq_offset,
        purity_rmq=purity_rmq, remove_overlaps=False, score=False,
        at_threshold=None, filter_homopolymers=False,
    )


def filter_at_content(hits: list[dict], threshold: float) -> list[dict]:
    """Return candidates below the left-arm AT fraction threshold."""
    return _filter_at_content(hits, threshold)


def filter_homopolymers(hits: list[dict], *, mode: str = "POST") -> list[dict]:
    """Exclude pure homopolymers using scored (POST) or raw (PRE) motifs.

    ``POST`` checks both ``full_sequence`` and ``putative_triplex``;
    ``PRE`` checks ``full_sequence``. Use PRE before candidates are scored.
    """
    return _filter_homopolymer_triplex(hits, mode=mode)


def remove_overlaps(hits: list[dict], strategy: str = "greedy") -> list[dict]:
    """Select nonoverlapping candidates; ``score`` prioritizes individual hits."""
    return _filter_overlapping_hits(hits, strategy=strategy)


def score_candidates(hits: list[dict], model: str = "complete") -> list[dict]:
    """Return copied hits scored with pairing, stacking, or the full model.

    ``pairing`` and ``pairing_stacking`` use the detected full-arm boundaries.
    ``complete`` also optimizes the arm boundaries. Invalid or unscorable
    motifs receive the same empty score fields as the main API.
    """
    if model not in ("pairing", "pairing_stacking", "complete"):
        raise ValueError(f"Invalid scoring model {model!r}")
    result = []
    for hit in hits:
        scored = score_hit_components(
            hit["left_arm"], hit["spacer"], hit["right_arm"], hit["arm_length"],
            include_stacking=model != "pairing",
            optimize_boundaries=model == "complete",
        )
        copy = hit.copy()
        if scored is None:
            copy.update(stacking_score=None, pairing_score=None,
                        total_score=None, putative_triplex="")
        else:
            copy.update(scored)
        result.append(copy)
    return result


def scan_sequence_ablation(
    seq: str, *, model: str = "complete", minrep: int = 10,
    maxrep: int = 1000, maxspacer: int = 10, purity: float = 0.90,
    mismatch: float = 0.10, at_threshold: float | None = 0.80,
    filter_homopolymers_enabled: bool = True,
    remove_overlaps_enabled: bool = False,
    overlap_strategy: str = "greedy", seq_offset: int = 1,
) -> list[dict]:
    """Run one stage on a DNA sequence with explicit parameter controls.

    ``detection`` uses the C mirror scanner with purity=0 and no composition
    or score filter. ``composition`` adds the requested C purity threshold
    and AT filter. ``pairing`` adds full-arm pairing energies including
    mismatch penalties. ``pairing_stacking`` also adds stacking energies.
    ``complete`` enables boundary optimization and optional homopolymer
    filtering. Overlap removal is controlled separately for every model;
    it is off by default so candidate recall can be compared fairly.

    Mismatch and arm/spacer thresholds apply to every stage. The C scanner
    retains only the longest valid arm for each center and spacer, even when
    purity=0. Thus this is an ablation of the available detector, not an
    exhaustive enumeration of all possible mirror-repeat arm lengths.
    """
    if model not in MODELS:
        raise ValueError(f"Invalid ablation model {model!r}; choose from {MODELS}")
    if overlap_strategy not in ("greedy", "score", "stability"):
        raise ValueError(f"Invalid overlap strategy {overlap_strategy!r}")
    if remove_overlaps_enabled and overlap_strategy != "greedy" and model in (
        "detection", "composition"
    ):
        raise ValueError("Score overlap filtering requires a scoring model")

    hits = detect_candidates(
        seq, minrep=minrep, maxrep=maxrep, maxspacer=maxspacer,
        purity=0.0 if model == "detection" else purity,
        mismatch=mismatch, seq_offset=seq_offset,
    )
    if model != "detection" and at_threshold is not None:
        hits = filter_at_content(hits, at_threshold)
    if model in ("pairing", "pairing_stacking", "complete"):
        hits = score_candidates(hits, model=model)
    if model == "complete" and filter_homopolymers_enabled:
        hits = filter_homopolymers(hits)
    if remove_overlaps_enabled:
        hits = remove_overlaps(hits, strategy=overlap_strategy)
    return hits


def scan_fasta_ablation(
    path: str | Path, *, parser: str = "biopython", **kwargs
) -> list[dict]:
    """Apply :func:`scan_sequence_ablation` to all FASTA records."""
    results = []
    for seq_id, seq, offset in parse_fasta(path, parser=parser):
        hits = scan_sequence_ablation(seq, seq_offset=offset, **kwargs)
        for hit in hits:
            hit["seq_id"] = seq_id
        results.extend(hits)
    return results
