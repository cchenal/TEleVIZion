#!/usr/bin/env python3
"""
Convert a FASTA file into a chroms.tsv, optionally filtering by sequence length.

chr    start   end     name    gieStain
<id>   1       <len>   <id>    chalk
"""

import argparse
from pathlib import Path

from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser(
        description="Create chroms.tsv from a FASTA input (using Biopython)."
    )
    parser.add_argument("fasta", help="Input FASTA file")
    parser.add_argument(
        "-o", "--out",
        help="Output TSV path. Default: <input_dir>/chroms.tsv"
    )
    parser.add_argument(
        "--gieStain", default="chalk",
        help='Value for gieStain column (default: "chalk")'
    )
    parser.add_argument(
        "--min-length", type=int, default=0,
        help="Only include sequences with length >= this value (bp). Default: 0"
    )
    parser.add_argument(
        "--no-header", action="store_true",
        help="Do not write the header line"
    )
    args = parser.parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise FileNotFoundError("Input FASTA not found: {}".format(fasta_path))

    out_path = Path(args.out) if args.out else (fasta_path.parent / "chroms.tsv")

    kept = 0
    skipped = 0

    with out_path.open("w", encoding="utf-8", newline="\n") as out:
        if not args.no_header:
            out.write("chr\tstart\tend\tname\tgieStain\n")

        for record in SeqIO.parse(str(fasta_path), "fasta"):
            length = len(record.seq)
            if length < args.min_length:
                skipped += 1
                continue

            seq_id = record.id
            out.write("{}\t1\t{}\t{}\t{}\n".format(
                seq_id, length, seq_id, args.gieStain
            ))
            kept += 1

    print("Wrote: {}".format(out_path))
    print("Sequences kept: {}, skipped (< {} bp): {}".format(
        kept, args.min_length, skipped
    ))


if __name__ == "__main__":
    main()
