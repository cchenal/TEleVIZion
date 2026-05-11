"""
CLI entrypoint for TEleVIZion.
python3 scripts/televizion_cli.py --name MyGenome --genome data/Acol_lib_GCA_943734845.1/chroms.tsv --repeatmasker data/Acol_lib_GCA_943734845.1/GCA_943734845.1.out --kimura data/Acol_lib_GCA_943734845.1/GCA_943734845.1.kimura  --windowsize 5000000
"""

import argparse

from televizion import aggregation as televizion_aggregation
from televizion import io as televizion_io
from televizion import plotting as televizion_plotting


def parse_figsize(value):
    try:
        width_s, height_s = value.split(",")
        return int(width_s), int(height_s)
    except (ValueError, AttributeError):
        raise argparse.ArgumentTypeError("figsize must be in the format W,H (e.g., 10,8).")


def parse_args():
    epilog = (
        "Examples:\n"
        "  RepeatMasker + Kimura:\n"
        "    python3 scripts/televizion_cli.py --name MyGenome --genome data/genome.tsv \\\n"
        "      --repeatmasker data/repeats.out --kimura data/repeats.divsum --windowsize 500000 \\\n"
        "      --chromtoplot Chr_2R,Chr_2L,Chr_3R,Chr_3L,Chr_X --perchromosome\n"
        "  EDTA:\n"
        "    python3 scripts/televizion_cli.py --name MyGenome --genome data/genome.tsv \\\n"
        "      --edta data/edta.gff3 --windowsize 500000 --chromtoplot 2RL,3RL,X\n"
    )
    parser = argparse.ArgumentParser(
        description="Generate TEleVIZion plots and tables from RepeatMasker or EDTA annotations.",
        epilog=epilog,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--name", required=False, type=str, default="output", help="Output prefix.")
    parser.add_argument(
        "--genome",
        required=True,
        type=str,
        help="Genome metadata TSV (chrom, start, end, name).",
    )
    parser.add_argument("--fasta", required=False, type=str, default=None, help="Genome FASTA.")
    parser.add_argument(
        "--kimura",
        required=False,
        type=str,
        default=None,
        help="RepeatMasker Kimura .divsum file.",
    )
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument("--repeatmasker", type=str, default=None, help="RepeatMasker .out file.")
    input_group.add_argument("--edta", type=str, default=None, help="EDTA GFF3 annotation file.")
    parser.add_argument(
        "--windowsize",
        required=False,
        type=int,
        default=10000,
        help="Window size in bp for visualisation.",
    )
    parser.add_argument(
        "--chromtoplot",
        required=False,
        type=str,
        default="all",
        help="Comma-separated chromosome list or 'all'.",
    )
    parser.add_argument("--perchromosome", action="store_true", help="Generate per-chromosome plots.")
    parser.add_argument(
        "--classesorder",
        required=False,
        type=str,
        default=None,
        help="Comma-separated class order override.",
    )
    parser.add_argument("--perclass", action="store_true", help="Generate per-class plots.")
    parser.add_argument(
        "--accessibility",
        required=False,
        type=str,
        default=None,
        help="Accessibility track path for karyotype plotting.",
    )
    parser.add_argument(
        "--figsize",
        required=False,
        type=parse_figsize,
        default=None,
        help="Figure size as W,H (e.g., 10,8).",
    )
    parser.add_argument(
        "--palette",
        required=False,
        type=str,
        default=None,
        help="TSV palette file to override class colors.",
    )
    parser.add_argument(
        "--layout",
        required=False,
        type=str,
        default="horizontal",
        help="Karyotype layout: horizontal or vertical.",
    )
    parser.add_argument(
        "--zoom",
        required=False,
        type=str,
        default=None,
        help="Optional karyotype zoom region as chrom:start-end, e.g. Chr1:1000000-2000000.",
    )
    parser.add_argument(
        "--zoom-chromosome",
        required=False,
        type=str,
        default=None,
        help="Optional chromosome for karyotype zoom. Use with --zoom-start and --zoom-end.",
    )
    parser.add_argument(
        "--zoom-start",
        required=False,
        type=int,
        default=None,
        help="Optional start coordinate for karyotype zoom. Use with --zoom-chromosome and --zoom-end.",
    )
    parser.add_argument(
        "--zoom-end",
        required=False,
        type=int,
        default=None,
        help="Optional end coordinate for karyotype zoom. Use with --zoom-chromosome and --zoom-start.",
    )
    args = parser.parse_args()
    if args.kimura is not None and args.repeatmasker is None:
        parser.error("--kimura requires --repeatmasker.")
    zoom_fields = [
        args.zoom_chromosome is not None,
        args.zoom_start is not None,
        args.zoom_end is not None,
    ]
    if args.zoom is not None and any(zoom_fields):
        parser.error("Use either --zoom or --zoom-chromosome/--zoom-start/--zoom-end, not both.")
    if any(zoom_fields) and not all(zoom_fields):
        parser.error("--zoom-chromosome, --zoom-start, and --zoom-end must be used together.")
    if args.zoom_start is not None and (args.zoom_start < 1 or args.zoom_end < args.zoom_start):
        parser.error("Zoom coordinates must satisfy 1 <= start <= end.")
    return args


def main():
    print("\nTEleVIZion is on!\n")

    args = parse_args()

    name = args.name
    genome = args.genome
    fasta = args.fasta

    if args.repeatmasker is not None:
        input_path = args.repeatmasker
        file_type = "RepeatMasker"
    elif args.edta is not None:
        input_path = args.edta
        file_type = "EDTA"
    else:
        parser = argparse.ArgumentParser()
        parser.error("Provide an input file using --repeatmasker or --edta.")

    kimura_file = args.kimura
    window_size = args.windowsize
    chrom_to_plot = args.chromtoplot if args.chromtoplot != "all" else "all"
    classes_order = args.classesorder.split(",") if args.classesorder is not None else None
    accessibility = args.accessibility
    per_chrom = args.perchromosome
    per_class = args.perclass

    if args.figsize is not None:
        width, height = args.figsize
    else:
        width, height = 10, 8

    palette = args.palette
    layout = args.layout
    zoom = args.zoom
    zoom_chromosome = args.zoom_chromosome
    zoom_start = args.zoom_start
    zoom_end = args.zoom_end

    windows, chromosomes_to_plot, chrom_names = televizion_aggregation.build_windows(
        genome_file=genome,
        window_size=window_size,
        chrom_to_plot=chrom_to_plot,
    )

    gc_path = televizion_io.compute_gc_content(
        fasta_file=fasta,
        gc_windows=10000,
    )

    print("\nInputs")
    print(f"- name: {name}")
    print(f"- genome: {genome}")
    print(f"- fasta: {fasta}")
    print(f"- input ({file_type}): {input_path}")
    print(f"- kimura: {kimura_file}")

    print("\nOptions")
    print(f"- window size: {window_size}")
    print(f"- chromosomes to plot: {chrom_to_plot}")
    print(f"- classes order: {classes_order}")
    print(f"- per chromosome: {per_chrom}")
    print(f"- per class: {per_class}")
    print(f"- figure size: {width}, {height}")
    print(f"- palette: {palette if palette is not None else 'default'}")
    print(f"- layout: {layout}")
    if zoom is not None:
        print(f"- karyotype zoom: {zoom}")
    elif zoom_chromosome is not None:
        print(f"- karyotype zoom: {zoom_chromosome}:{zoom_start}-{zoom_end}")
    else:
        print("- karyotype zoom: full selected chromosomes")
    print(f"- accessibility: {accessibility}")
    print(f"- gc content: {gc_path}")
    print("\n")

    kimura_div = None
    if args.repeatmasker is not None:
        if kimura_file is not None:
            kimura_dict = televizion_io.parse_repeatmasker_kimura_bins(
                repeatmasker_kimura_file=kimura_file
            )
            repeats, insertions, kimura_div = televizion_io.parse_repeatmasker_annotations(
                repeatmasker_file=input_path,
                windows_dict=windows,
                kimura_dict=kimura_dict,
            )
        else:
            repeats, insertions = televizion_io.parse_repeatmasker_annotations(
                repeatmasker_file=input_path,
                windows_dict=windows,
            )
    else:
        repeats, insertions = televizion_io.parse_edta_annotations(
            edta_file=input_path,
            windows_dict=windows,
        )

    agg_global, class_order, class_colors_hex, fam_colors = televizion_plotting.build_color_maps(
        annotations_by_window=insertions,
        classes_order=classes_order,
        palette=palette,
    )

    televizion_plotting.plot_repeat_family_bars(
        name=name,
        annotations_by_window=insertions,
        per_chromosome=per_chrom,
        chrom_names=chrom_names,
        agg_global=agg_global,
        class_order=class_order,
        class_colors_hex=class_colors_hex,
        fam_colors=fam_colors,
        width=width,
        height=height,
    )

    reversed_classes, reversed_colors = televizion_aggregation.export_window_class_table(
        name=name,
        window_size=window_size,
        windows=windows,
        annotations_by_window=insertions,
        class_colors=class_colors_hex,
        class_order=class_order,
    )

    kimura_bed = None
    if kimura_file is not None and kimura_div is not None:
        kimura_bed = televizion_aggregation.export_window_kimura_table(
            name=name,
            window_size=window_size,
            kimura_bins=kimura_div,
            windows=windows,
            class_order=class_order,
        )

    identity_bed = None
    if args.edta is not None:
        identity_bed = televizion_aggregation.export_window_identity_table(
            name=name,
            window_size=window_size,
            annotations_by_window=insertions,
            windows=windows,
            class_order=class_order,
        )

    televizion_plotting.plot_karyotype_tracks(
        name=name,
        window_size=window_size,
        genome_file=genome,
        chromosomes=chromosomes_to_plot,
        accessibility=accessibility,
        gc_content=gc_path,
        classes=reversed_classes,
        colors=reversed_colors,
        plot_per_class=per_class,
        kimura_bed=kimura_bed,
        identity_bed=identity_bed,
        layout=layout,
        zoom=zoom,
        zoom_chromosome=zoom_chromosome,
        zoom_start=zoom_start,
        zoom_end=zoom_end,
    )


if __name__ == "__main__":
    main()
