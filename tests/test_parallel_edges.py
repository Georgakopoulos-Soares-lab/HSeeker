"""
Edge-case tests for scan_fasta_parallel chunking logic.
Run: conda run -n hseeker_test python -m pytest tests/test_parallel_edges.py -v
"""
import tempfile
from pathlib import Path

import pytest
import hseeker


# ── helpers ──

def _tmp_fasta(*records):
    """(seq_id, seq) tuples → temp file path."""
    tmp = tempfile.NamedTemporaryFile(mode="w", suffix=".fa", delete=False)
    for sid, seq in records:
        tmp.write(f">{sid}\n{seq}\n")
    tmp.close()
    return Path(tmp.name)


def _key(h):
    return (h["seq_id"], h["start"], h["end"], h["arm_length"])


PARAMS = dict(minrep=10, maxrep=1000, maxspacer=10, purity=0.90, mismatch=0.10)


# ═══════════════════════════════════════════════════════════════════
# 1. single chunk (record < chunk_size) — matches serial exactly
# ═══════════════════════════════════════════════════════════════════

def test_single_chunk_matches_serial():
    """1000 bp random-ish sequence: serial == parallel with chunk_size=5000."""
    seq = "A" * 50 + "G" * 50 + "C" * 50 + "T" * 50 + "A" * 200 + "G" * 200
    p = _tmp_fasta(("test", seq))
    ser = hseeker.scan_fasta(str(p), **PARAMS)
    par = hseeker.scan_fasta_parallel(str(p), chunk_size=5000, workers=2, **PARAMS)
    assert sorted(_key(h) for h in ser) == sorted(_key(h) for h in par)
    assert len(ser) == len(par)


# ═══════════════════════════════════════════════════════════════════
# 2. multi-chunk long homopolymer — hit near boundary
# ═══════════════════════════════════════════════════════════════════

def test_hit_across_chunk_boundary():
    """Place a pure poly-A tract that crosses chunk_size boundary."""
    chunk_size = 500
    # 300-A pad so a center near the boundary can extend left 3000+
    seq = "A" * (chunk_size) + "A" * 4000  # ensures hits straddle the 500 bp boundary
    p = _tmp_fasta(("test", seq))
    ser = hseeker.scan_fasta(str(p), **PARAMS)
    par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=2, **PARAMS)
    ser_keys = sorted(_key(h) for h in ser)
    par_keys = sorted(_key(h) for h in par)
    assert ser_keys == par_keys, (
        f"Serial {len(ser)} hits vs parallel {len(par)} hits\n"
        + f"Missing from par: {set(ser_keys)-set(par_keys)}\n"
        + f"Extra in par: {set(par_keys)-set(ser_keys)}"
    )


# ═══════════════════════════════════════════════════════════════════
# 3. record exactly at chunk_size — edge of single vs multi chunk
# ═══════════════════════════════════════════════════════════════════

def test_record_exactly_chunk_size():
    """Record exactly chunk_size bases → single chunk. Verify correctness."""
    chunk_size = 1000
    seq = "G" * chunk_size  # pure G — every center produces a hit
    p = _tmp_fasta(("exact", seq))
    ser = hseeker.scan_fasta(str(p), **PARAMS)
    par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=2, **PARAMS)
    assert sorted(_key(h) for h in ser) == sorted(_key(h) for h in par)


# ═══════════════════════════════════════════════════════════════════
# 4. record 1bp longer than chunk_size → 2 chunks
# ═══════════════════════════════════════════════════════════════════

def test_record_one_over_chunk():
    chunk_size = 1000
    seq = "A" * (chunk_size + 1)
    p = _tmp_fasta(("one_over", seq))
    ser = hseeker.scan_fasta(str(p), **PARAMS)
    par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=2, **PARAMS)
    assert sorted(_key(h) for h in ser) == sorted(_key(h) for h in par)


# ═══════════════════════════════════════════════════════════════════
# 5. empty sequence
# ═══════════════════════════════════════════════════════════════════

def test_empty_record():
    p = _tmp_fasta(("empty", ""))
    assert hseeker.scan_fasta_parallel(str(p), chunk_size=100, workers=2, **PARAMS) == []


# ═══════════════════════════════════════════════════════════════════
# 6. sequence shorter than minrep
# ═══════════════════════════════════════════════════════════════════

def test_too_short():
    p = _tmp_fasta(("short", "AAAAA"))
    assert hseeker.scan_fasta_parallel(str(p), chunk_size=10, workers=2, **PARAMS) == []


# ═══════════════════════════════════════════════════════════════════
# 7. multiple records — mixed chunked and non-chunked
# ═══════════════════════════════════════════════════════════════════

def test_multi_record_mixed():
    chunk_size = 500
    records = [
        ("chr1", "A" * 2000),           # 4 chunks
        ("chr2", "G" * 200),            # 1 chunk (below chunk_size)
        ("chr3", "C" * chunk_size),     # exactly chunk_size → 1 chunk
        ("chr4", ""),                   # empty
        ("chr5", "T" * 5),              # too short
        ("chr6", "A" * (chunk_size + 1)),  # 2 chunks
    ]
    p = _tmp_fasta(*records)
    ser = hseeker.scan_fasta(str(p), **PARAMS)
    par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=3, **PARAMS)
    assert sorted(_key(h) for h in ser) == sorted(_key(h) for h in par)


# ═══════════════════════════════════════════════════════════════════
# 8. set of workers (1, 2, 4) — determinism
# ═══════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("workers", [1, 2, 4])
def test_worker_independence(workers):
    chunk_size = 500
    seq = "A" * 2000 + "C" * 500 + "G" * 500
    p = _tmp_fasta(("wtest", seq))
    ser = hseeker.scan_fasta(str(p), **PARAMS)
    par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=workers, **PARAMS)
    assert sorted(_key(h) for h in ser) == sorted(_key(h) for h in par)


# ═══════════════════════════════════════════════════════════════════
# 9. minimal chunk_size (triggers many small chunks)
# ═══════════════════════════════════════════════════════════════════

def test_many_tiny_chunks():
    """Force many chunks (chunk_size just above overlap size).
    Overlap removal is order-dependent; serial and parallel may produce
    different but equally valid non-overlapping decompositions."""
    overlap = 2 * PARAMS["maxrep"] + PARAMS["maxspacer"] + 1  # 6021
    chunk_size = overlap + 100  # 6121
    seq = "A" * (chunk_size * 3)  # 3 chunks for a pure-A tract
    p = _tmp_fasta(("tiny", seq))
    ser = hseeker.scan_fasta(str(p), **PARAMS)
    par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=3, **PARAMS)

    # both must produce non-overlapping hit sets
    def _check_no_overlap(hits):
        sh = sorted(hits, key=lambda h: h["start"])
        for i in range(len(sh) - 1):
            assert sh[i]["end"] < sh[i + 1]["start"], \
                f"overlap: {sh[i]} vs {sh[i+1]}"
    _check_no_overlap(ser)
    _check_no_overlap(par)
    # both should have the same number of hits (4 for this pure-A tract)
    assert len(ser) == len(par) == 4, f"serial={len(ser)} parallel={len(par)}"


# ═══════════════════════════════════════════════════════════════════
# 10. hit exactly at chunk boundary with max arm length
# ═══════════════════════════════════════════════════════════════════

def test_hit_at_exact_boundary_with_max_arm():
    """Construct a sequence where the center producing a max-arm hit lands
    exactly at the chunk boundary exclusion point."""
    chunk_size = 500
    maxrep = PARAMS["maxrep"]
    # The exclusion condition: start = ctr - best_k + 1 >= chunk_size
    # We want ctr such that start = chunk_size (excluded from chunk N).
    # start = ctr - maxrep + 1 = chunk_size → ctr = chunk_size + maxrep - 1
    # This ctr is in the overlap of chunk N (since ctr >= chunk_size).
    # Chunk N+1 (local) should find it entirely.
    
    # seq: pad + boundary region where a pure-A tract starts at chunk_size - maxrep
    # This way the center at chunk_size + maxrep - 1 has left_i go back into the pad.
    seq = "A" * (chunk_size + 2 * maxrep + 100)
    p = _tmp_fasta(("boundary", seq))
    ser = hseeker.scan_fasta(str(p), **PARAMS)
    par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=2, **PARAMS)
    ser_keys = sorted(_key(h) for h in ser)
    par_keys = sorted(_key(h) for h in par)
    assert ser_keys == par_keys, f"diff: serial={len(ser_keys)} parallel={len(par_keys)}"


# ═══════════════════════════════════════════════════════════════════
# 11. N bases at chunk boundary
# ═══════════════════════════════════════════════════════════════════

def test_n_bases_near_boundary():
    """N bases at positions that would be centers — should be skipped."""
    chunk_size = 500
    seq = "A" * 200 + "N" * 100 + "A" * 200 + "N" * 50 + "A" * 5000
    p = _tmp_fasta(("n_boundary", seq))
    ser = hseeker.scan_fasta(str(p), **PARAMS)
    par = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=3, **PARAMS)
    assert sorted(_key(h) for h in ser) == sorted(_key(h) for h in par)


# ═══════════════════════════════════════════════════════════════════
# 12. determinism — same input, same output
# ═══════════════════════════════════════════════════════════════════

def test_determinism():
    chunk_size = 800
    seq = "G" * 3000 + "C" * 2000 + "A" * 5000
    p = _tmp_fasta(("det", seq))
    a = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=3, **PARAMS)
    b = hseeker.scan_fasta_parallel(str(p), chunk_size=chunk_size, workers=3, **PARAMS)
    assert sorted(_key(h) for h in a) == sorted(_key(h) for h in b)
