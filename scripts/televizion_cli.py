"""
CLI entrypoint for TEleVIZion.
python3 scripts/televizion_cli.py --name MyGenome --genome data/Acol_lib_GCA_943734845.1/chroms.tsv --repeatmasker data/Acol_lib_GCA_943734845.1/GCA_943734845.1.out --kimura data/Acol_lib_GCA_943734845.1/GCA_943734845.1.kimura  --windowsize 5000000
"""

import argparse

from televizion import aggregation as televizion_aggregation
from televizion import io as televizion_io
from televizion import plotting_stats as televizion_plotting_stats


def parse_figsize(value):
    try:
        width_s, height_s = value.split(",")
        return int(width_s), int(height_s)
    except (ValueError, AttributeError):
        raise argparse.ArgumentTypeError("figsize must be in the format W,H (e.g., 10,8).")


def parse_output_formats(value):
    allowed_formats = {"pdf", "png", "jpg"}
    output_formats = []
    for raw_format in value.split(","):
        fmt = raw_format.strip().lower()
        if not fmt:
            continue
        if fmt not in allowed_formats:
            raise argparse.ArgumentTypeError(
                "output formats must be one or more of: pdf,png,jpg"
            )
        if fmt not in output_formats:
            output_formats.append(fmt)
    if not output_formats:
        raise argparse.ArgumentTypeError(
            "output formats must be one or more of: pdf,png,jpg"
        )
    return output_formats


def parse_dpi(value):
    try:
        dpi = int(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError("dpi must be an integer >= 300.")
    if dpi < 300:
        raise argparse.ArgumentTypeError("dpi must be >= 300.")
    return dpi


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
    track_group = parser.add_mutually_exclusive_group(required=False)
    track_group.add_argument(
        "--gc",
        required=False,
        type=str,
        default=None,
        help="Precomputed GC-content track created by scripts/utils/create_gc_content.py. Mutually exclusive with --gene-content.",
    )
    parser.add_argument(
        "--kimura",
        required=False,
        type=str,
        default=None,
        help="RepeatMasker Kimura .divsum file.",
    )
    input_group = parser.add_mutually_exclusive_group(required=False)
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
    track_group.add_argument(
        "--gene-content",
        required=False,
        type=str,
        default=None,
        help="Precomputed gene-content track created by scripts/utils/create_gene_content.py. Mutually exclusive with --gc.",
    )
    parser.add_argument(
        "--figsize",
        required=False,
        type=parse_figsize,
        default=None,
        help=(
            "Figure size as W,H for Python general statistics bar plots only; "
            "does not affect karyoplots."
        ),
    )
    parser.add_argument(
        "--output-formats",
        required=False,
        type=parse_output_formats,
        default=["pdf"],
        help=(
            "Comma-separated output figure formats from pdf,png,jpg. "
            "Default: pdf."
        ),
    )
    parser.add_argument(
        "--dpi",
        required=False,
        type=parse_dpi,
        default=300,
        help="Raster output resolution in dots per inch. Must be >= 300. Default: 300.",
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
    args = parser.parse_args()
    if args.repeatmasker is None and args.edta is None:
        parser.error("Provide an input file using --repeatmasker or --edta.")
    if args.kimura is not None and args.edta is not None:
        parser.error("--kimura cannot be used with --edta.")
    if args.kimura is not None and args.repeatmasker is None:
        parser.error("--kimura requires --repeatmasker.")
    return args


def main():
    print("\nTEleVIZion is on!\n")

    args = parse_args()

    name = args.name
    genome = args.genome
    gc_content = args.gc

    if args.repeatmasker is not None:
        input_path = args.repeatmasker
        file_type = "RepeatMasker"
    elif args.edta is not None:
        input_path = args.edta
        file_type = "EDTA"

    kimura_file = args.kimura
    window_size = args.windowsize
    chrom_to_plot = args.chromtoplot if args.chromtoplot != "all" else "all"
    classes_order = args.classesorder.split(",") if args.classesorder is not None else None
    gene_content = args.gene_content
    per_chrom = args.perchromosome
    per_class = args.perclass
    output_formats = args.output_formats
    dpi = args.dpi

    if args.figsize is not None:
        width, height = args.figsize
    else:
        width, height = 10, 8

    palette = args.palette
    layout = args.layout
    zoom = args.zoom

    windows, chromosomes_to_plot, chrom_names = televizion_aggregation.build_windows(
        genome_file=genome,
        window_size=window_size,
        chrom_to_plot=chrom_to_plot,
    )


    print("\nInputs")
    print(f"- name: {name}")
    print(f"- genome: {genome}")
    print(f"- gc content: {gc_content}")
    print(f"- gene content: {gene_content}")
    print(f"- input ({file_type}): {input_path}")
    print(f"- kimura: {kimura_file}")

    print("\nOptions")
    print(f"- window size: {window_size}")
    print(f"- chromosomes to plot: {chrom_to_plot}")
    print(f"- classes order: {classes_order}")
    print(f"- per chromosome: {per_chrom}")
    print(f"- per class: {per_class}")
    print(f"- general statistics figure size: {width}, {height}")
    print(f"- output formats: {','.join(output_formats)}")
    print(f"- raster output dpi: {dpi}")
    print(f"- palette: {palette if palette is not None else 'default'}")
    print(f"- layout: {layout}")
    if zoom is not None:
        print(f"- karyotype zoom: {zoom}")
    else:
        print("- karyotype zoom: full selected chromosomes")
    print("\n")

    kimura_bed = None
    identity_bed = None

    kimura_div = None
    repeatmasker_identity = None
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
            (
                repeats,
                insertions,
                repeatmasker_identity,
            ) = televizion_io.parse_repeatmasker_annotations(
                repeatmasker_file=input_path,
                windows_dict=windows,
            )
    else:
        repeats, insertions = televizion_io.parse_edta_annotations(
            edta_file=input_path,
            windows_dict=windows,
        )

    agg_global, class_order, class_colors_hex, fam_colors = televizion_plotting_stats.build_color_maps(
        annotations_by_window=insertions,
        classes_order=classes_order,
        palette=palette,
    )

    televizion_plotting_stats.plot_repeat_family_bars(
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
        output_formats=output_formats,
        dpi=dpi,
    )

    reversed_classes, reversed_colors = televizion_aggregation.export_window_class_table(
        name=name,
        window_size=window_size,
        windows=windows,
        annotations_by_window=insertions,
        class_colors=class_colors_hex,
        class_order=class_order,
    )

    if kimura_file is not None and kimura_div is not None:
        kimura_bed = televizion_aggregation.export_window_kimura_table(
            name=name,
            window_size=window_size,
            kimura_bins=kimura_div,
            windows=windows,
            class_order=class_order,
        )

    if args.edta is not None:
        identity_bed = televizion_aggregation.export_window_identity_table(
            name=name,
            window_size=window_size,
            annotations_by_window=insertions,
            windows=windows,
            class_order=class_order,
        )
    elif repeatmasker_identity is not None:
        identity_bed = televizion_aggregation.export_window_identity_table(
            name=name,
            window_size=window_size,
            annotations_by_window=repeatmasker_identity,
            windows=windows,
            class_order=class_order,
        )

    televizion_plotting_stats.plot_karyotype_tracks(
        name=name,
        window_size=window_size,
        genome_file=genome,
        chromosomes=chromosomes_to_plot,
        gene_content=gene_content,
        gc_content=gc_content,
        classes=reversed_classes,
        colors=reversed_colors,
        plot_per_class=per_class,
        kimura_bed=kimura_bed,
        identity_bed=identity_bed,
        layout=layout,
        zoom=zoom,
        output_formats=output_formats,
        dpi=dpi,
    )


if __name__ == "__main__":
    main()
