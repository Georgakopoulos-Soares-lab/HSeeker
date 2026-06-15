from __future__ import annotations

import datetime as dt
import json
import os
import platform
import time
from pathlib import Path
from typing import Any

import hseeker

try:
    import numpy as np
except Exception:  # pragma: no cover
    np = None

try:
    import psutil
except Exception:  # pragma: no cover
    psutil = None

WORKERS = [1, 4, 7, 10, 14]
MINREP = 10
ROOT = Path(__file__).resolve().parent
DATA_DIR = ROOT / "data"
OUT_JSON = ROOT / "parallel_workers_sweep.json"
OUT_MD = ROOT / "final_parallel_workers_report.md"


def get_system_info() -> dict[str, Any]:
    info: dict[str, Any] = {
        "timestamp": dt.datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "os": platform.platform(),
        "python": platform.python_version(),
        "hseeker": getattr(hseeker, "__version__", "unknown"),
        "cpu_cores_logical": os.cpu_count(),
    }
    if np is not None:
        info["numpy"] = np.__version__
    if psutil is not None:
        vm = psutil.virtual_memory()
        info["ram_total_gb"] = round(vm.total / 1e9, 1)
        info["ram_available_gb"] = round(vm.available / 1e9, 1)
        info["cpu_cores_physical"] = psutil.cpu_count(logical=False)
    return info


def run_one(path: Path, workers: int) -> dict[str, Any]:
    file_mb = path.stat().st_size / 1e6
    t0 = time.perf_counter()
    c0 = time.process_time()
    hits = hseeker.scan_fasta_parallel(str(path), workers=workers, minrep=MINREP)
    c1 = time.process_time()
    t1 = time.perf_counter()
    wall = t1 - t0
    cpu = c1 - c0
    return {
        "workers": workers,
        "wall_s": wall,
        "cpu_s": cpu,
        "hits": len(hits),
        "throughput_mbps": (file_mb / wall) if wall > 0 else 0.0,
    }


def sweep() -> dict[str, Any]:
    fasta_files = sorted(DATA_DIR.glob("*.fa"))
    out: dict[str, Any] = {
        "system_info": get_system_info(),
        "workers": WORKERS,
        "minrep": MINREP,
        "files": [],
    }

    for path in fasta_files:
        file_mb = path.stat().st_size / 1e6
        print(f"Running: {path.name} ({file_mb:.1f} MB)")
        rows: list[dict[str, Any]] = []
        baseline = None
        base_hits = None
        for w in WORKERS:
            print(f"  workers={w} ... ", end="", flush=True)
            r = run_one(path, w)
            if baseline is None:
                baseline = r["wall_s"]
                base_hits = r["hits"]
            r["speedup_vs_1w"] = (baseline / r["wall_s"]) if r["wall_s"] > 0 else 0.0
            r["efficiency_pct"] = (r["speedup_vs_1w"] / w) * 100.0
            rows.append(r)
            print(f"{r['wall_s']:.2f}s, hits={r['hits']:,}")

        hits_consistent = all(row["hits"] == base_hits for row in rows)
        out["files"].append(
            {
                "file": path.name,
                "path": str(path),
                "size_mb": file_mb,
                "hits_consistent": hits_consistent,
                "results": rows,
            }
        )

    return out


def write_markdown(payload: dict[str, Any]) -> str:
    si = payload["system_info"]
    lines: list[str] = []
    lines.append("# Final Parallel Worker Sweep Report")
    lines.append("")
    lines.append(f"**Date**: {si['timestamp']}  ")
    lines.append(f"**minrep**: {payload['minrep']}  ")
    lines.append(f"**Workers swept**: {', '.join(str(w) for w in payload['workers'])}  ")
    lines.append("")

    lines.append("## System")
    lines.append("")
    lines.append("| Property | Value |")
    lines.append("|---|---|")
    lines.append(f"| CPU logical cores | {si.get('cpu_cores_logical', 'n/a')} |")
    lines.append(f"| CPU physical cores | {si.get('cpu_cores_physical', 'n/a')} |")
    lines.append(f"| RAM total (GB) | {si.get('ram_total_gb', 'n/a')} |")
    lines.append(f"| RAM available (GB) | {si.get('ram_available_gb', 'n/a')} |")
    lines.append(f"| OS | {si.get('os', 'n/a')} |")
    lines.append(f"| Python | {si.get('python', 'n/a')} |")
    lines.append(f"| hseeker | {si.get('hseeker', 'n/a')} |")
    lines.append(f"| NumPy | {si.get('numpy', 'n/a')} |")
    lines.append("")

    lines.append("## Per-file Scaling")
    lines.append("")
    for f in payload["files"]:
        lines.append(f"### {f['file']} ({f['size_mb']:.1f} MB)")
        lines.append("")
        lines.append(f"Hits consistent across workers: **{f['hits_consistent']}**")
        lines.append("")
        lines.append("| Workers | Wall (s) | CPU (s) | Hits | Throughput (MB/s) | Speedup vs 1w | Efficiency |")
        lines.append("|---|---|---|---|---|---|---|")
        for r in f["results"]:
            lines.append(
                f"| {r['workers']} | {r['wall_s']:.3f} | {r['cpu_s']:.3f} | {r['hits']:,} | "
                f"{r['throughput_mbps']:.2f} | {r['speedup_vs_1w']:.2f}x | {r['efficiency_pct']:.1f}% |"
            )
        lines.append("")

    lines.append("## 14-worker Cross-file Summary")
    lines.append("")
    lines.append("| File | Size (MB) | 1w wall (s) | 14w wall (s) | Speedup | Efficiency | Throughput 14w (MB/s) |")
    lines.append("|---|---|---|---|---|---|---|")
    for f in payload["files"]:
        by_w = {row["workers"]: row for row in f["results"]}
        r1 = by_w[1]
        r14 = by_w[14]
        speedup = r1["wall_s"] / r14["wall_s"] if r14["wall_s"] > 0 else 0.0
        eff = speedup / 14.0 * 100.0
        lines.append(
            f"| {f['file']} | {f['size_mb']:.1f} | {r1['wall_s']:.3f} | {r14['wall_s']:.3f} | "
            f"{speedup:.2f}x | {eff:.1f}% | {r14['throughput_mbps']:.2f} |"
        )
    lines.append("")

    lines.append("## Findings")
    lines.append("")
    best = None
    worst = None
    for f in payload["files"]:
        by_w = {row["workers"]: row for row in f["results"]}
        speedup = by_w[1]["wall_s"] / by_w[14]["wall_s"] if by_w[14]["wall_s"] > 0 else 0.0
        if best is None or speedup > best[1]:
            best = (f["file"], speedup)
        if worst is None or speedup < worst[1]:
            worst = (f["file"], speedup)
    if best and worst:
        lines.append(f"- Best 14-worker speedup: **{best[0]}** at **{best[1]:.2f}x**.")
        lines.append(f"- Lowest 14-worker speedup: **{worst[0]}** at **{worst[1]:.2f}x**.")
    inconsistent = [f["file"] for f in payload["files"] if not f["hits_consistent"]]
    if inconsistent:
        lines.append(f"- Hit counts differed by worker for: {', '.join(inconsistent)}")
    else:
        lines.append("- Hit counts were identical across all worker settings for every file.")

    lines.append("")
    return "\n".join(lines)


def main() -> None:
    payload = sweep()
    OUT_JSON.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(f"\nSaved JSON: {OUT_JSON}")
    print(f"Saved report: {OUT_MD}")


if __name__ == "__main__":
    main()
