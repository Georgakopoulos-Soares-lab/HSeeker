"""Profiling counts remain correct when the C scanner runs concurrently."""

from concurrent.futures import ThreadPoolExecutor

from hseeker import _hdna


SEQUENCE = "GA" * 40
OPTIONS = dict(minrep=8, maxrep=20, maxspacer=4, purity=0.9,
               mismatch=0.1, remove_overlaps=False)


def scan(purity_rmq):
    return _hdna.scan_sequence(SEQUENCE, purity_rmq=purity_rmq, **OPTIONS)


def test_parallel_scans_preserve_every_profile_count():
    per_scan = {}
    for purity_rmq in (False, True):
        _hdna.reset_profiling()
        scan(purity_rmq)
        info = _hdna.profiling_info()
        assert info["scans_completed"] == 1
        assert info["inner_iters"] > 0
        assert info["ctr_sp_pairs"] > 0
        per_scan[purity_rmq] = info

    modes = [False, True] * 20
    _hdna.reset_profiling()
    with ThreadPoolExecutor(max_workers=8) as pool:
        list(pool.map(scan, modes))
    info = _hdna.profiling_info()
    assert info["scans_completed"] == len(modes)
    assert info["total_inner_iters"] == sum(
        per_scan[mode]["inner_iters"] for mode in modes
    )
    assert info["total_ctr_sp_pairs"] == sum(
        per_scan[mode]["ctr_sp_pairs"] for mode in modes
    )
    assert info["inner_iters"] in {per_scan[mode]["inner_iters"] for mode in modes}


def test_empty_scan_has_zero_last_counts_and_keeps_totals():
    _hdna.reset_profiling()
    scan(False)
    before = _hdna.profiling_info()
    assert _hdna.scan_sequence("") == []
    after = _hdna.profiling_info()
    assert after["inner_iters"] == after["ctr_sp_pairs"] == 0
    assert after["total_inner_iters"] == before["total_inner_iters"]
    assert after["total_ctr_sp_pairs"] == before["total_ctr_sp_pairs"]
    assert after["scans_completed"] == 2
