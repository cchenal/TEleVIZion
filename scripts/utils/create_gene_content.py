#!/usr/bin/env python3
"""
Create a BED-like gene-density track colored by a continuous colormap.

Inputs:
1. Genome metadata TSV produced by create_chroms.py
2. GFF/GFF3 annotation
3. Optional NCBI sequence_report.tsv for RefSeq -> GenBank accession mapping

Genome metadata columns expected:
    chr    start    end    name    gieStain

Output columns:
    chr    start    end    name    itemRgb

Behavior:
- Genome windows are non-overlapping and 0-based half-open.
- Only GFF features with type "gene" are counted.
- Each gene is counted exactly once, according to its START coordinate.
- Gene counts are normalized using ONE global maximum across all windows.
- Windows with zero genes are retained.
- If --sequence-report is provided, GFF RefSeq accessions are translated
  to GenBank accessions before matching them to the genome metadata.
"""

import argparse
import csv
import sys

import seaborn as sns
from matplotlib.colors import to_hex


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Create a gene-density BED track from a TEleVIZion genome metadata "
            "TSV and a GFF/GFF3 annotation."
        )
    )

    parser.add_argument(
        "genome",
        help=(
            "Genome metadata TSV produced by create_chroms.py. "
            "Must contain chr and end columns."
        ),
    )

    parser.add_argument(
        "gff",
        help="Input GFF/GFF3 annotation file. Use '-' for stdin.",
    )

    parser.add_argument(
        "--sequence-report",
        metavar="TSV",
        help=(
            "Optional NCBI sequence_report.tsv used to map RefSeq sequence "
            "accessions in the GFF to GenBank sequence accessions in the "
            "genome metadata."
        ),
    )

    parser.add_argument(
        "--window",
        type=int,
        default=100000,
        help="Window size in bp. Default: 100000.",
    )

    parser.add_argument(
        "--out",
        default="-",
        help="Output BED path. Default: stdout.",
    )

    parser.add_argument(
        "--cmap",
        default="Blues",
        help='Seaborn/matplotlib colormap name. Default: "Blues".',
    )

    return parser.parse_args()


def load_genome(path):
    """
    Load chromosome lengths and order from create_chroms.py output.

    Returns:
        chrom_order: list of chromosome IDs in input order
        chrom_lengths: dict mapping chromosome ID -> length
    """
    chrom_order = []
    chrom_lengths = {}

    with open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {"chr", "end"}
        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                "Genome metadata is missing required column(s): {}".format(
                    ", ".join(sorted(missing))
                )
            )

        for row in reader:
            chrom = row["chr"].strip()

            if not chrom:
                continue

            try:
                length = int(row["end"])
            except ValueError:
                raise ValueError(
                    "Invalid chromosome length for '{}': {!r}".format(
                        chrom,
                        row["end"],
                    )
                )

            if length <= 0:
                raise ValueError(
                    "Chromosome length must be > 0 for '{}': {}".format(
                        chrom,
                        length,
                    )
                )

            if chrom in chrom_lengths:
                raise ValueError(
                    "Duplicate chromosome ID in genome metadata: {}".format(
                        chrom
                    )
                )

            chrom_order.append(chrom)
            chrom_lengths[chrom] = length

    return chrom_order, chrom_lengths


def load_sequence_report(path):
    """
    Load RefSeq -> GenBank sequence accession mappings from an
    NCBI sequence_report.tsv.

    Expected columns:
        RefSeq seq accession
        GenBank seq accession

    Returns:
        dict mapping RefSeq accession -> GenBank accession
    """
    refseq_to_genbank = {}

    with open(path, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")

        required = {
            "RefSeq seq accession",
            "GenBank seq accession",
        }

        missing = required - set(reader.fieldnames or [])

        if missing:
            raise ValueError(
                "NCBI sequence report is missing required column(s): {}".format(
                    ", ".join(sorted(missing))
                )
            )

        for row in reader:
            refseq = row["RefSeq seq accession"].strip()
            genbank = row["GenBank seq accession"].strip()

            # Some NCBI reports can contain empty values or "na".
            if not refseq or not genbank:
                continue

            if refseq.lower() == "na" or genbank.lower() == "na":
                continue

            refseq_to_genbank[refseq] = genbank

    return refseq_to_genbank


def initialize_windows(chrom_order, chrom_lengths, window_size):
    """
    Create a count array for every chromosome.
    """
    counts = {}

    for chrom in chrom_order:
        length = chrom_lengths[chrom]

        n_windows = (length + window_size - 1) // window_size

        counts[chrom] = [0] * n_windows

    return counts


def count_genes(
    gff_path,
    chrom_lengths,
    counts,
    window_size,
    refseq_to_genbank=None,
):
    """
    Count GFF 'gene' features according to their START position.

    If refseq_to_genbank is supplied, sequence IDs in the GFF are first
    translated from RefSeq accessions to GenBank accessions.

    GFF coordinates are 1-based inclusive. They are converted to 0-based
    coordinates before assigning genes to windows.

    Each gene contributes to exactly one window.
    """
    genes_counted = 0
    genes_mapped = 0
    genes_unknown_chrom = 0
    genes_out_of_bounds = 0

    handle = (
        sys.stdin
        if gff_path == "-"
        else open(gff_path, "rt", encoding="utf-8")
    )

    try:
        for line_number, line in enumerate(handle, start=1):
            line = line.rstrip("\n")

            if not line or line.startswith("#"):
                continue

            fields = line.split("\t")

            if len(fields) < 9:
                raise ValueError(
                    "Invalid GFF line {}: expected at least 9 columns".format(
                        line_number
                    )
                )

            gff_chrom = fields[0]
            feature_type = fields[2]

            if feature_type != "gene":
                continue

            chrom = gff_chrom

            # If the GFF uses RefSeq accessions but the genome metadata
            # uses GenBank accessions, translate here.
            if refseq_to_genbank is not None:
                mapped_chrom = refseq_to_genbank.get(gff_chrom)

                if mapped_chrom is not None:
                    chrom = mapped_chrom
                    genes_mapped += 1

            if chrom not in chrom_lengths:
                genes_unknown_chrom += 1
                continue

            try:
                gff_start = int(fields[3])
            except ValueError:
                raise ValueError(
                    "Invalid gene start coordinate on GFF line {}: {!r}".format(
                        line_number,
                        fields[3],
                    )
                )

            # GFF: 1-based inclusive
            # BED/windows: 0-based half-open
            start_0 = gff_start - 1

            if start_0 < 0 or start_0 >= chrom_lengths[chrom]:
                genes_out_of_bounds += 1
                continue

            window_index = start_0 // window_size

            counts[chrom][window_index] += 1
            genes_counted += 1

    finally:
        if handle is not sys.stdin:
            handle.close()

    return (
        genes_counted,
        genes_mapped,
        genes_unknown_chrom,
        genes_out_of_bounds,
    )


def get_global_max(counts):
    """
    Return the maximum gene count across ALL windows and chromosomes.
    """
    global_max = 0

    for chrom_counts in counts.values():
        if chrom_counts:
            global_max = max(global_max, max(chrom_counts))

    return global_max


def count_to_hex(count, global_max, cmap):
    """
    Map a raw gene count from 0..global_max onto colormap position 0..1.
    """
    if global_max == 0:
        norm = 0.0
    else:
        norm = count / global_max

    rgba = cmap(norm)

    return to_hex(rgba[:3])


def write_output(
    out_path,
    chrom_order,
    chrom_lengths,
    counts,
    window_size,
    cmap,
    global_max,
):
    """
    Write the BED-like colored gene-density track.
    """
    out_handle = (
        sys.stdout
        if out_path == "-"
        else open(out_path, "wt", encoding="utf-8", newline="")
    )

    try:
        out_handle.write("chr\tstart\tend\tname\titemRgb\n")

        uid = 0

        for chrom in chrom_order:
            chrom_length = chrom_lengths[chrom]

            for window_index, count in enumerate(counts[chrom]):
                start = window_index * window_size
                end = min(start + window_size, chrom_length)

                uid += 1

                color = count_to_hex(
                    count=count,
                    global_max=global_max,
                    cmap=cmap,
                )

                name = "genes{}-{}".format(uid, count)

                out_handle.write(
                    "{}\t{}\t{}\t{}\t{}\n".format(
                        chrom,
                        start,
                        end,
                        name,
                        color,
                    )
                )

    finally:
        if out_handle is not sys.stdout:
            out_handle.close()


def main():
    args = parse_args()

    if args.window <= 0:
        raise ValueError("--window must be greater than 0")

    chrom_order, chrom_lengths = load_genome(args.genome)

    if not chrom_order:
        raise ValueError(
            "No chromosomes found in genome metadata."
        )

    refseq_to_genbank = None

    if args.sequence_report:
        refseq_to_genbank = load_sequence_report(
            args.sequence_report
        )

        print(
            "Loaded {} RefSeq -> GenBank sequence mappings.".format(
                len(refseq_to_genbank)
            ),
            file=sys.stderr,
        )

    counts = initialize_windows(
        chrom_order=chrom_order,
        chrom_lengths=chrom_lengths,
        window_size=args.window,
    )

    (
        genes_counted,
        genes_mapped,
        genes_unknown_chrom,
        genes_out_of_bounds,
    ) = count_genes(
        gff_path=args.gff,
        chrom_lengths=chrom_lengths,
        counts=counts,
        window_size=args.window,
        refseq_to_genbank=refseq_to_genbank,
    )

    global_max = get_global_max(counts)

    cmap = sns.color_palette(
        args.cmap,
        as_cmap=True,
    )

    write_output(
        out_path=args.out,
        chrom_order=chrom_order,
        chrom_lengths=chrom_lengths,
        counts=counts,
        window_size=args.window,
        cmap=cmap,
        global_max=global_max,
    )

    print(
        "Genes counted: {}".format(genes_counted),
        file=sys.stderr,
    )

    if refseq_to_genbank is not None:
        print(
            "Genes whose sequence ID was mapped RefSeq -> GenBank: {}".format(
                genes_mapped
            ),
            file=sys.stderr,
        )

    print(
        "Global maximum genes per window: {}".format(
            global_max
        ),
        file=sys.stderr,
    )

    if genes_unknown_chrom:
        print(
            "Genes skipped because chromosome was absent from genome metadata: "
            "{}".format(genes_unknown_chrom),
            file=sys.stderr,
        )

    if genes_out_of_bounds:
        print(
            "Genes skipped because start coordinate was outside chromosome "
            "bounds: {}".format(genes_out_of_bounds),
            file=sys.stderr,
        )


if __name__ == "__main__":
    main()