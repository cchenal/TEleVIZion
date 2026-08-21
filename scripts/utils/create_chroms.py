#!/usr/bin/env python3
"""
Create a TEleVIZion genome metadata table from a FASTA file.

The output is a tab-separated table consumed by --genome:

chr    start   end     name    gieStain
<id>   1       <len>   <id>    chalk

If an NCBI sequence_report.tsv is provided, FASTA IDs matching the
"GenBank seq accession" column will use the corresponding "Sequence name"
as the display name.
"""

import argparse
import csv
from pathlib import Path

from Bio import SeqIO


def parse_args():
    epilog = (
        "Output columns:\n"
        "  chr       FASTA record ID; must match IDs in the annotation file.\n"
        "  start     Always 1.\n"
        "  end       FASTA sequence length in bp.\n"
        "  name      Display label used in karyoplots. If --sequence-report is\n"
        "            provided, uses NCBI Sequence name when available.\n"
        "  gieStain  Cytoband-style label consumed by karyoploteR.\n"
        "\n"
        "Examples:\n"
        "  python3 scripts/utils/create_chroms.py -f genome.fasta\n"
        "  python3 scripts/utils/create_chroms.py -f genome.fasta \\\n"
        "    -o data/my_genome/chroms.tsv --min-length 1000000\n"
        "  python3 scripts/utils/create_chroms.py -f genome.fasta \\\n"
        "    --sequence-report sequence_report.tsv\n"
        "  python3 scripts/utils/create_chroms.py -f genome.fasta \\\n"
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
            "Input genome FASTA. Record IDs become the chr values in "
            "the output table."
        ),
    )
    parser.add_argument(
        "-o", "--out",
        metavar="TSV",
        help="Output TSV path. Default: <input_dir>/chroms.tsv."
    )
    parser.add_argument(
        "--sequence-report",
        metavar="TSV",
        help=(
            "Optional NCBI sequence_report.tsv. FASTA IDs matching the "
            "'GenBank seq accession' column will use the corresponding "
            "'Sequence name' as the output name."
        ),
    )
    parser.add_argument(
        "--gieStain", default="gneg",
        metavar="LABEL",
        help='Value written to the gieStain column. Default: "gneg".'
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


def load_sequence_names(report_path):
    """Return {GenBank accession: Sequence name} from an NCBI sequence report."""
    names = {}

    with report_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"GenBank seq accession", "Sequence name"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(
                "Sequence report is missing required column(s): {}".format(
                    ", ".join(sorted(missing))
                )
            )

        for row in reader:
            accession = row["GenBank seq accession"].strip()
            sequence_name = row["Sequence name"].strip()

            if accession and sequence_name:
                names[accession] = sequence_name

    return names


def main():
    args = parse_args()

    fasta_path = Path(args.fasta)
    if not fasta_path.exists():
        raise FileNotFoundError(
            "Input FASTA not found: {}".format(fasta_path)
        )

    out_path = (
        Path(args.out)
        if args.out
        else fasta_path.parent / "chroms.tsv"
    )

    sequence_names = {}

    if args.sequence_report:
        report_path = Path(args.sequence_report)

        if not report_path.exists():
            raise FileNotFoundError(
                "NCBI sequence report not found: {}".format(report_path)
            )

        sequence_names = load_sequence_names(report_path)

    kept = 0
    skipped = 0
    renamed = 0

    with out_path.open("w", encoding="utf-8", newline="\n") as out:
        if not args.no_header:
            out.write("chr\tstart\tend\tname\tgieStain\n")

        for record in SeqIO.parse(str(fasta_path), "fasta"):
            length = len(record.seq)

            if length < args.min_length:
                skipped += 1
                continue

            seq_id = record.id

            # Keep chr as the FASTA ID, but use the NCBI Sequence name
            # as the display label when a match is available.
            name = sequence_names.get(seq_id, seq_id)

            if name != seq_id:
                renamed += 1

            out.write(
                "{}\t1\t{}\t{}\t{}\n".format(
                    seq_id,
                    length,
                    name,
                    args.gieStain,
                )
            )

            kept += 1

    print("Wrote TEleVIZion genome metadata: {}".format(out_path))
    print(
        "Sequences kept: {}, skipped (< {} bp): {}".format(
            kept, args.min_length, skipped
        )
    )

    if args.sequence_report:
        print(
            "Sequence names replaced from NCBI report: {}".format(
                renamed
            )
        )

    print("Output columns: chr, start, end, name, gieStain")


if __name__ == "__main__":
    main()