"""
CLI entry point for hseeker.

Invoked as:
    hseeker -seq genome.fa -out results [options]

or equivalently:
    python -m hseeker -seq genome.fa -out results [options]

The output is a TSV file at <out>_HDNA.tsv with the same column layout as
the original findHDNA binary, making it a drop-in replacement.
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import hseeker


def main() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        prog="hseeker",
        description="HSeeker — H-DNA / Triplex Mirror Repeat Detector",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  hseeker -seq test.fa -out test\n"
            "  hseeker -seq genome.fa -out genome -purity 1.0 -mismatch 0.0\n"
            "  hseeker -seq genome.fa -out genome -minrep 10 -purity 0.85 -v\n"
        ),
    )

    # mirroring the original findHDNA CLI flags exactly
    parser.add_argument("-seq",         required=True,       help="Input FASTA file")
    parser.add_argument("-out",         required=True,       help="Output prefix → <prefix>_HDNA.tsv")
    parser.add_argument("-minrep",      type=int,   default=10,    metavar="INT",
                        help="Minimum arm length (default: 10)")
    parser.add_argument("-maxrep",      type=int,   default=1000,  metavar="INT",
                        help="Maximum arm length (default: 1000)")
    parser.add_argument("-maxspacer",   type=int,   default=10,   metavar="INT",
                        help="Maximum spacer length (default: 10)")
    parser.add_argument("-purity",      type=float, default=0.90, metavar="FLOAT",
                        help="Min GA or CT fraction in arm (default: 0.90)")
    parser.add_argument("-mismatch",    type=float, default=0.10, metavar="FLOAT",
                        help="Max mismatch fraction in mirror (default: 0.10)")
    parser.add_argument("-skipoverlap", action="store_true",
                        help="Skip overlap removal (keep all raw hits)")
    parser.add_argument("-score", action="store_true", default=True,
                        help="Apply stability scoring (default: on)")
    parser.add_argument("-no-score", action="store_false", dest="score",
                        help="Disable stability scoring")
    parser.add_argument("-v",           action="store_true",
                        help="Verbose: print per-sequence stats to stderr")

    args = parser.parse_args()

    out_path = Path(args.out + "_HDNA.tsv")

    # TSV column order matches original findHDNA output, plus scoring columns
    fieldnames = [
        "seq_id", "source", "start", "end",
        "arm_length", "spacer_length", "total_length",
        "ga_pct", "ct_pct", "mirror_identity", "is_perfect",
        "left_arm", "spacer", "right_arm", "full_sequence",
        "stacking_score", "pairing_score", "total_score",
        "putative_triplex",
    ]

    print(
        f"hseeker v{hseeker.__version__} — H-DNA / Triplex Mirror Repeat Detector\n"
        f"  Input : {args.seq}\n"
        f"  Output: {out_path}\n"
        f"  minrep={args.minrep}  maxrep={args.maxrep}  maxspacer={args.maxspacer}\n"
        f"  purity={args.purity:.2f}  mismatch={args.mismatch:.2f}\n",
        file=sys.stderr,
    )

    total_records = 0
    total_hits = 0

    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore",  # ignore seq_id added by scan_fasta
        )
        writer.writeheader()

        for seq_id, seq, offset in hseeker.parse_fasta(args.seq):
            total_records += 1
            if args.v:
                print(
                    f"Processing {seq_id} ({len(seq):,} bases, offset {offset})...",
                    file=sys.stderr,
                )

            hits = hseeker.scan_sequence(
                seq,
                minrep=args.minrep,
                maxrep=args.maxrep,
                maxspacer=args.maxspacer,
                purity=args.purity,
                mismatch=args.mismatch,
                remove_overlaps=not args.skipoverlap,
                seq_offset=offset,
                score=args.score,
            )

            if args.v:
                print(f"  {len(hits)} hits", file=sys.stderr)

            for h in hits:
                h["seq_id"] = seq_id
                h["source"] = "findHDNA"
            writer.writerows(hits)
            total_hits += len(hits)

    print(
        f"\nDone.\n"
        f"  Records processed : {total_records}\n"
        f"  Total hits written: {total_hits}\n"
        f"  Output: {out_path}",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
