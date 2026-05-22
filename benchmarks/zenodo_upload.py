#!/usr/bin/env python3
"""
benchmarks/zenodo_upload.py — Upload HSeeker benchmark data to Zenodo
======================================================================

Generates benchmark FASTA files (if not already cached in benchmarks/data/),
compresses them with gzip, and uploads them to a Zenodo deposit.  After
publishing the deposit, the record ID, DOI, and a per-file manifest are
written to benchmarks/zenodo_record.json so that any user can later run

    python benchmarks/benchmark.py --from-zenodo

to download the files instead of regenerating them locally.

Prerequisites
-------------
1. Create a Zenodo personal access token at:
     https://zenodo.org/account/settings/applications/tokens/new/
   Required scopes:  deposit:write  deposit:actions
   (For the sandbox:  https://sandbox.zenodo.org/account/settings/…)

2. Export the token before running this script:
     export ZENODO_TOKEN=your_token_here

Usage
-----
    # Upload small tier only (default; ~30 MB compressed, ≈10 files)
    python benchmarks/zenodo_upload.py

    # Also upload the medium tier (~300 MB compressed)
    python benchmarks/zenodo_upload.py --medium

    # Test against the Zenodo sandbox first (no real DOI issued)
    python benchmarks/zenodo_upload.py --sandbox

    # Continue uploading files to an existing unpublished draft
    python benchmarks/zenodo_upload.py --record-id 12345678

    # Create a new version of an already-published record
    python benchmarks/zenodo_upload.py --new-version --record-id 12345678

    # Skip FASTA generation (benchmarks/data/ files must already exist)
    python benchmarks/zenodo_upload.py --no-generate

    # Upload without publishing (leaves deposit as draft for review)
    python benchmarks/zenodo_upload.py --no-publish

Notes
-----
- Files are stored on Zenodo as <name>.fa.gz  (gzip-compressed FASTA).
- The SEED_MANIFEST.json is uploaded alongside the FASTA files so that
  seeds are documented inside the Zenodo record.
- benchmark.py --from-zenodo downloads .fa.gz and decompresses to .fa
  transparently; the MD5 checksum is verified after download.
- Re-running this script after a published deposit will fail unless you
  pass --new-version or --record-id pointing to a fresh draft.
"""
from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import os
import shutil
import sys
import time
import urllib.error
import urllib.request
from pathlib import Path

# ---------------------------------------------------------------------------
# Bootstrap: resolve repo root so benchmark helpers can be imported
# ---------------------------------------------------------------------------
_BENCH_DIR  = Path(__file__).resolve().parent
_REPO_ROOT  = _BENCH_DIR.parent
sys.path.insert(0, str(_REPO_ROOT / "src"))   # for hseeker (if needed)
sys.path.insert(0, str(_BENCH_DIR))            # for benchmark.py helpers

from benchmark import (  # noqa: E402
    DATA_DIR,
    SMALL_DATASETS,
    MEDIUM_DATASETS,
    DatasetSpec,
    _ensure_dataset,
    _download_chr16,
    CHR16_LOCAL_NAME,
)

RECORD_FILE = _BENCH_DIR / "zenodo_record.json"

_ZENODO_API  = "https://zenodo.org/api"
_SANDBOX_API = "https://sandbox.zenodo.org/api"

_ZENODO_META = {
    "metadata": {
        "title": "HSeeker Benchmark FASTA Datasets",
        "upload_type": "dataset",
        "description": (
            "Synthetic FASTA files used by the HSeeker benchmark suite "
            "(benchmarks/benchmark.py). Generated deterministically with "
            "fixed NumPy/Python random seeds documented in SEED_MANIFEST.json. "
            "Four sequence profiles: uniform (25 % each ACGT), ga_biased "
            "(90 % purine), ct_biased (90 % pyrimidine), realistic (~41 % GC). "
            "Two size tiers: small (~30 MB) and medium (~300 MB). "
            "Download via: python benchmarks/benchmark.py --from-zenodo"
        ),
        "creators": [{"name": "HSeeker Contributors"}],
        "keywords": [
            "bioinformatics", "hseeker", "H-DNA", "triplex",
            "benchmark", "FASTA", "synthetic",
        ],
        "license": "cc-by-4.0",
        "access_right": "open",
    }
}


# ---------------------------------------------------------------------------
# HTTP helpers (stdlib only — no requests dependency)
# ---------------------------------------------------------------------------

def _api_request(
    method: str,
    url: str,
    token: str,
    *,
    json_data: dict | None = None,
) -> dict:
    """Execute a JSON API request and return the parsed response."""
    body = json.dumps(json_data).encode() if json_data is not None else None
    req = urllib.request.Request(
        url,
        data=body,
        method=method,
        headers={
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/json",
        },
    )
    try:
        with urllib.request.urlopen(req) as resp:
            return json.loads(resp.read().decode())
    except urllib.error.HTTPError as exc:
        body_text = exc.read().decode(errors="replace")
        raise RuntimeError(
            f"Zenodo API {method} {url} → HTTP {exc.code}: {body_text}"
        ) from exc


def _upload_file(bucket_url: str, token: str, local_path: Path, remote_name: str) -> dict:
    """Upload a file to a Zenodo S3-style bucket via PUT."""
    url = f"{bucket_url}/{remote_name}"
    size_mb = local_path.stat().st_size / 1e6
    print(f"    uploading {remote_name:45s} ({size_mb:6.1f} MB) … ", end="", flush=True)
    t0 = time.perf_counter()
    with open(local_path, "rb") as fh:
        data = fh.read()
    req = urllib.request.Request(
        url, data=data, method="PUT",
        headers={
            "Authorization": f"Bearer {token}",
            "Content-Type": "application/octet-stream",
        },
    )
    try:
        with urllib.request.urlopen(req) as resp:
            result = json.loads(resp.read().decode())
            elapsed = time.perf_counter() - t0
            speed = size_mb / elapsed if elapsed > 0 else 0
            print(f"done ({elapsed:.1f}s, {speed:.1f} MB/s)")
            return result
    except urllib.error.HTTPError as exc:
        body_text = exc.read().decode(errors="replace")
        raise RuntimeError(
            f"Upload failed for {remote_name}: HTTP {exc.code}: {body_text}"
        ) from exc


# ---------------------------------------------------------------------------
# Compression and checksum
# ---------------------------------------------------------------------------

def _gzip_file(src: Path, dst: Path) -> None:
    """Compress src → dst with gzip level 6."""
    with open(src, "rb") as f_in, gzip.open(dst, "wb", compresslevel=6) as f_out:
        shutil.copyfileobj(f_in, f_out)


def _md5(path: Path) -> str:
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------

def _parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument(
        "--medium", action="store_true",
        help="Also upload the medium-tier datasets (~300 MB compressed)",
    )
    p.add_argument(
        "--sandbox", action="store_true",
        help="Use the Zenodo sandbox (no real DOI; good for testing)",
    )
    p.add_argument(
        "--record-id", metavar="ID",
        help="Add files to an existing unpublished draft deposit",
    )
    p.add_argument(
        "--new-version", action="store_true",
        help="Create a new version of the record given by --record-id",
    )
    p.add_argument(
        "--no-generate", action="store_true",
        help="Skip FASTA generation; benchmarks/data/ files must exist",
    )
    p.add_argument(
        "--no-publish", action="store_true",
        help="Upload files but leave the deposit as a draft (do not publish)",
    )
    p.add_argument(
        "--chr16", action="store_true",
        help=(
            "Download hg38 chr16 from NCBI (~90 MB) and include it in the "
            "Zenodo upload as a real-genome benchmark dataset."
        ),
    )
    return p.parse_args()


def main() -> None:
    args = _parse_args()

    token = os.getenv("ZENODO_TOKEN")
    if not token:
        sys.exit(
            "ERROR: ZENODO_TOKEN environment variable is not set.\n"
            "Create a token at:\n"
            "  https://zenodo.org/account/settings/applications/tokens/new/\n"
            "  (sandbox: https://sandbox.zenodo.org/account/settings/"
            "applications/tokens/new/)\n"
            "Required scopes: deposit:write  deposit:actions\n\n"
            "Then:  export ZENODO_TOKEN=your_token_here"
        )

    api_base = _SANDBOX_API if args.sandbox else _ZENODO_API
    sandbox_note = "  [SANDBOX — no real DOI]" if args.sandbox else ""
    print(f"Zenodo endpoint: {api_base}{sandbox_note}\n")

    # ── 1. Select dataset specs ────────────────────────────────────────────
    specs: list[DatasetSpec] = list(SMALL_DATASETS)
    if args.medium:
        specs += MEDIUM_DATASETS

    # ── 2. Generate FASTA files if needed ─────────────────────────────────
    if not args.no_generate:
        print("Generating/verifying benchmark FASTA files…")
        for spec in specs:
            _ensure_dataset(
                spec.name, spec.n_records, spec.record_length,
                spec.profile, spec.seed, force=False,
            )
        print("Generation complete.\n")
    # ── 2b. Download hg38 chr16 if requested ─────────────────────────────
    include_chr16 = False
    if args.chr16:
        print("Downloading hg38 chr16 from NCBI…")
        chr16_path = _download_chr16()
        if chr16_path is None:
            sys.exit("ERROR: chr16 download failed. Aborting upload.")
        include_chr16 = True
        print()
    # ── 3. Compress files into a staging directory ─────────────────────────
    gz_dir = DATA_DIR / "_gz"
    gz_dir.mkdir(parents=True, exist_ok=True)

    print("Compressing FASTA files (gzip level 6)…")
    gz_files: list[tuple[Path, str]] = []   # (local_gz_path, remote_name)

    for spec in specs:
        src = DATA_DIR / f"{spec.name}.fa"
        if not src.exists():
            sys.exit(
                f"\nERROR: {src} does not exist.\n"
                "Run without --no-generate first, or re-run benchmark.py to "
                "regenerate benchmark data."
            )
        dst = gz_dir / f"{spec.name}.fa.gz"
        if dst.exists():
            print(f"  {spec.name}.fa.gz  already compressed ({dst.stat().st_size / 1e6:.1f} MB)")
        else:
            print(f"  compressing {spec.name}.fa … ", end="", flush=True)
            _gzip_file(src, dst)
            ratio = src.stat().st_size / dst.stat().st_size
            print(f"{dst.stat().st_size / 1e6:.1f} MB  ({ratio:.1f}× compression)")
        gz_files.append((dst, dst.name))

    # Also include SEED_MANIFEST.json (not compressed — already tiny JSON)
    manifest_src = DATA_DIR / "SEED_MANIFEST.json"
    if manifest_src.exists():
        manifest_dst = gz_dir / "SEED_MANIFEST.json"
        shutil.copy2(manifest_src, manifest_dst)
        gz_files.append((manifest_dst, "SEED_MANIFEST.json"))

    # Include chr16 if downloaded
    if include_chr16:
        chr16_src = DATA_DIR / CHR16_LOCAL_NAME
        chr16_gz  = gz_dir / (CHR16_LOCAL_NAME + ".gz")
        if chr16_gz.exists():
            print(f"  {chr16_gz.name}  already compressed ({chr16_gz.stat().st_size / 1e6:.1f} MB)")
        else:
            print(f"  compressing {CHR16_LOCAL_NAME} … ", end="", flush=True)
            _gzip_file(chr16_src, chr16_gz)
            ratio = chr16_src.stat().st_size / chr16_gz.stat().st_size
            print(f"{chr16_gz.stat().st_size / 1e6:.1f} MB  ({ratio:.1f}× compression)")
        gz_files.append((chr16_gz, chr16_gz.name))

    print(f"\n{len(gz_files)} file(s) ready to upload.\n")

    # ── 4. Create or retrieve Zenodo deposit ──────────────────────────────
    # Resolve which record to target.
    # Priority: --record-id CLI arg > zenodo_record.json > create brand-new deposit.
    record_id_to_use = args.record_id
    if not record_id_to_use and RECORD_FILE.exists():
        try:
            existing = json.loads(RECORD_FILE.read_text())
            if existing.get("record_id"):
                record_id_to_use = str(existing["record_id"])
                print(
                    f"Found existing record {record_id_to_use} in zenodo_record.json — "
                    f"will create a new version.\n"
                    f"  (Pass --record-id to override, or delete zenodo_record.json "
                    f"to force a brand-new deposit.)"
                )
        except Exception:
            pass

    if record_id_to_use and args.new_version:
        print(f"Creating new version of record {record_id_to_use}…")
        dep = _api_request(
            "POST",
            f"{api_base}/deposit/depositions/{record_id_to_use}/actions/newversion",
            token,
        )
        # The new-version response links to the latest draft
        new_url = dep["links"]["latest_draft"]
        dep = _api_request("GET", new_url, token)
        deposition_id = dep["id"]
        bucket_url    = dep["links"]["bucket"]
        print(f"  new-version draft ID: {deposition_id}")

    elif record_id_to_use:
        # Auto-detect: if the deposit is already published, create a new version.
        dep = _api_request(
            "GET", f"{api_base}/deposit/depositions/{record_id_to_use}", token
        )
        if dep.get("submitted"):
            print(f"Record {record_id_to_use} is published — creating new version…")
            dep = _api_request(
                "POST",
                f"{api_base}/deposit/depositions/{record_id_to_use}/actions/newversion",
                token,
            )
            new_url = dep["links"]["latest_draft"]
            dep = _api_request("GET", new_url, token)
            print(f"  new-version draft ID: {dep['id']}")
        else:
            print(f"Using existing draft deposit {record_id_to_use}…")
        deposition_id = dep["id"]
        bucket_url    = dep["links"]["bucket"]

    else:
        print("Creating new Zenodo deposit…")
        dep = _api_request(
            "POST", f"{api_base}/deposit/depositions", token, json_data={}
        )
        deposition_id = dep["id"]
        bucket_url    = dep["links"]["bucket"]
        print(f"  deposit ID: {deposition_id}")

        # Set metadata
        _api_request(
            "PUT",
            f"{api_base}/deposit/depositions/{deposition_id}",
            token,
            json_data=_ZENODO_META,
        )
        print("  metadata saved.")

    # ── 5. Upload files ────────────────────────────────────────────────────
    print(f"\nUploading {len(gz_files)} file(s) to deposit {deposition_id}…")
    file_manifest: list[dict] = []
    for local_path, remote_name in gz_files:
        _upload_file(bucket_url, token, local_path, remote_name)
        digest = _md5(local_path)
        file_manifest.append({
            "name":       remote_name,
            "dataset":    remote_name.replace(".fa.gz", ""),
            "size_bytes": local_path.stat().st_size,
            "md5":        digest,
            "checksum":   f"md5:{digest}",
        })

    # ── 6. Publish ─────────────────────────────────────────────────────────
    doi        = None
    record_url = None
    if not args.no_publish:
        print(f"\nPublishing deposit {deposition_id}…")
        pub = _api_request(
            "POST",
            f"{api_base}/deposit/depositions/{deposition_id}/actions/publish",
            token,
        )
        doi        = pub.get("doi") or pub.get("metadata", {}).get("doi")
        record_url = pub.get("links", {}).get("record_html")
        print(f"  ✓ Published!")
        print(f"  DOI : {doi}")
        print(f"  URL : {record_url}")
    else:
        web = api_base.replace("/api", "")
        print(f"\nDeposit saved as draft (not published).")
        print(f"  Review at: {web}/deposit/{deposition_id}")

    # ── 7. Save record info to zenodo_record.json ──────────────────────────
    record_data = {
        "_comment": (
            "Populated by benchmarks/zenodo_upload.py. "
            "Commit this file so contributors can run: "
            "python benchmarks/benchmark.py --from-zenodo"
        ),
        "record_id":   str(deposition_id),
        "doi":         doi,
        "record_url":  record_url,
        "base_url":    api_base.replace("/api", ""),
        "sandbox":     args.sandbox,
        "uploaded_at": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "files":       file_manifest,
    }
    RECORD_FILE.write_text(json.dumps(record_data, indent=2))
    print(f"\nRecord manifest saved → {RECORD_FILE.relative_to(_REPO_ROOT)}")
    print("\nNext steps:")
    print("  1. git add benchmarks/zenodo_record.json && git commit -m 'Add Zenodo record manifest'")
    print("  2. Users can now run: python benchmarks/benchmark.py --from-zenodo")


if __name__ == "__main__":
    main()
