#!/usr/bin/env python3
"""
Create a TEleVIZion genome metadata table from a FASTA file.

The output is a tab-separated table consumed by --genome:

chr    start   end     name    gieStain
<id>   1       <len>   <id>    chalk

Use this helper when FASTA sequence IDs match the chromosome/scaffold IDs in
your RepeatMasker or EDTA annotation. If annotation IDs differ, generate this
starter file and edit the chr/name columns as needed.
"""

import argparse
from pathlib import Path

from Bio import SeqIO


def parse_args():
    epilog = (
        "Output columns:\n"
        "  chr       FASTA record ID; must match IDs in the annotation file.\n"
        "  start     Always 1.\n"
        "  end       FASTA sequence length in bp.\n"
        "  name      Display label used in karyoplots; defaults to the FASTA ID.\n"
        "  gieStain  Cytoband-style label consumed by karyoploteR.\n"
        "\n"
        "Examples:\n"
        "  python3 scripts/televizion/create_chroms.py -f genome.fasta\n"
        "  python3 scripts/televizion/create_chroms.py -f genome.fasta \\\n"
        "    -o data/my_genome/chroms.tsv --min-length 1000000\n"
        "  python3 scripts/televizion/create_chroms.py -f genome.fasta \\\n"
        "    --gieStain gneg --no-header\n"
    )
    parser = argparse.ArgumentParser(
        description=(
            "Create a TEleVIZion genome metadata TSV from FASTA sequence lengths."
        ),
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument(
        "-f", "--fasta",
        required=True,
        metavar="FASTA",
        help=(
            "Input genome FASTA. Record IDs become the chr and name values in "
            "the output table."
        ),
    )
    parser.add_argument(
        "-o", "--out",
        metavar="TSV",
        help="Output TSV path. Default: <input_dir>/chroms.tsv."
    )
    parser.add_argument(
        "--gieStain", default="chalk",
        metavar="LABEL",
        help='Value written to the gieStain column. Default: "chalk".'
    )
    parser.add_argument(
        "--min-length", type=int, default=0,
        metavar="BP",
        help=(
            "Only include FASTA records with length >= BP. Use this to drop "
            "short unplaced scaffolds. Default: 0."
        ),
    )
    parser.add_argument(
        "--no-header", action="store_true",
        help="Do not write the chr/start/end/name/gieStain header line."
    )
    return parser.parse_args()


def main():
    args = parse_args()

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

    print("Wrote TEleVIZion genome metadata: {}".format(out_path))
    print("Sequences kept: {}, skipped (< {} bp): {}".format(
        kept, args.min_length, skipped
    ))
    print("Output columns: chr, start, end, name, gieStain")


if __name__ == "__main__":
    main()
