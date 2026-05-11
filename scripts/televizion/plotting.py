"""
Plotting helpers for TEleVIZion.
"""

import os
import shlex
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_hex
import subprocess


# def build_color_maps(annotations_by_window, classes_order=None, palette=None):
#     print("# Running build_color_maps function\n")
#
#     records = []
#     for chrom_windows in annotations_by_window.values():
#         for classes in chrom_windows.values():
#             for rep_class, families in classes.items():
#                 for rep_family, metrics in families.items():
#                     records.append(
#                         {
#                             "rep_class": rep_class,
#                             "rep_family": rep_family,
#                             "count": metrics["count"],
#                             "length": metrics["length"],
#                         }
#                     )
#     df_global = pd.DataFrame(records)
#     agg_global = df_global.groupby(["rep_class", "rep_family"])[["length", "count"]].sum().reset_index()
#
#     if classes_order is None:
#         class_order = (
#             agg_global.groupby("rep_class")["length"].sum().sort_values(ascending=False).index.tolist()
#         )
#     else:
#         class_order = classes_order
#     n = max(len(class_order) - 1, 1)
#
#     fam_colors, class_colors, class_colors_hex = {}, {}, {}
#
#     if palette is None:
#         cmap = plt.get_cmap("turbo")
#         class_colors = {cls: cmap(i / n) for i, cls in enumerate(class_order)}
#         class_colors_hex = {cls: to_hex(rgba, keep_alpha=True) for cls, rgba in class_colors.items()}
#     else:
#         homemade_palette = {}
#         for line in open(palette).readlines():
#             desc, r, g, b, hex, cat = line[:-1].split("\t")
#             for category in cat.split(","):
#                 homemade_palette[category] = {
#                     "r": int(r),
#                     "g": int(g),
#                     "b": int(b),
#                     "hex": hex,
#                 }
#         for cls in class_order:
#             class_colors[cls] = [
#                 homemade_palette[cls]["r"] / 255,
#                 homemade_palette[cls]["g"] / 255,
#                 homemade_palette[cls]["b"] / 255,
#                 1,
#             ]
#             class_colors_hex[cls] = homemade_palette[cls]["hex"]
#
#     for cls in class_order:
#         fams = (
#             agg_global[agg_global["rep_class"] == cls]
#             .sort_values("length", ascending=False)["rep_family"]
#             .tolist()
#         )
#         alphas = np.linspace(1.0, 0.3, len(fams))
#         base = class_colors[cls]
#         fam_colors[cls] = {fam: (base[0], base[1], base[2], alpha) for fam, alpha in zip(fams, alphas)}
#
#     return (agg_global, class_order, class_colors_hex, fam_colors)


def build_class_color_maps(class_order, palette=None):
    """
    Build class-level color maps from a class order and optional palette.

    Returns:
        class_colors: RGBA color per class.
        class_colors_hex: hex color per class.
    """
    n = max(len(class_order) - 1, 1)
    class_colors = {}
    class_colors_hex = {}

    if palette is None:
        cmap = plt.get_cmap("turbo")
        class_colors = {cls: cmap(i / n) for i, cls in enumerate(class_order)}
        class_colors_hex = {
            cls: to_hex(rgba, keep_alpha=True) for cls, rgba in class_colors.items()
        }
    else:
        homemade_palette = {}
        with open(palette, "r") as handle:
            for line in handle:
                desc, r, g, b, hex, cat = line.rstrip("\n").split("\t")
                for category in cat.split(","):
                    homemade_palette[category] = {
                        "r": int(r),
                        "g": int(g),
                        "b": int(b),
                        "hex": hex,
                    }
        for cls in class_order:
            class_colors[cls] = [
                homemade_palette[cls]["r"] / 255,
                homemade_palette[cls]["g"] / 255,
                homemade_palette[cls]["b"] / 255,
                1,
            ]
            class_colors_hex[cls] = homemade_palette[cls]["hex"]

    return class_colors, class_colors_hex


def build_color_maps(annotations_by_window, classes_order=None, palette=None):
    """
    Build aggregate counts/lengths and color maps for repeat classes/families.

    Returns:
        agg_global: aggregated metrics per class/family.
        class_order: class ordering used in plots.
        class_colors_hex: hex color per class.
        fam_colors: RGBA colors per family.
    """
    print("# Running build_color_maps function\n")

    agg_map = {}
    for chrom_windows in annotations_by_window.values():
        for classes in chrom_windows.values():
            for rep_class, families in classes.items():
                for rep_family, metrics in families.items():
                    key = (rep_class, rep_family)
                    if key not in agg_map:
                        agg_map[key] = {"length": 0, "count": 0}
                    agg_map[key]["length"] += metrics["length"]
                    agg_map[key]["count"] += metrics["count"]

    records = [
        {
            "rep_class": rep_class,
            "rep_family": rep_family,
            "length": totals["length"],
            "count": totals["count"],
        }
        for (rep_class, rep_family), totals in agg_map.items()
    ]
    agg_global = pd.DataFrame(records)

    if classes_order is None:
        class_order = (
            agg_global.groupby("rep_class")["length"].sum().sort_values(ascending=False).index.tolist()
        )
    else:
        class_order = classes_order
    fam_colors = {}
    class_colors, class_colors_hex = build_class_color_maps(
        class_order=class_order,
        palette=palette,
    )

    for cls in class_order:
        fams = (
            agg_global[agg_global["rep_class"] == cls]
            .sort_values("length", ascending=False)["rep_family"]
            .tolist()
        )
        alphas = np.linspace(1.0, 0.3, len(fams))
        base = class_colors[cls]
        fam_colors[cls] = {fam: (base[0], base[1], base[2], alpha) for fam, alpha in zip(fams, alphas)}

    return (agg_global, class_order, class_colors_hex, fam_colors)


# def plot_repeat_family_bars(
#     name,
#     annotations_by_window,
#     chrom_names,
#     agg_global,
#     class_order,
#     class_colors_hex,
#     fam_colors,
#     per_chromosome=False,
#     width=10,
#     height=8,
# ):
#     print("# Running plot_repeat_family_bars function\n")
#
#     os.makedirs(f"analyses/{name}", exist_ok=True)
#
#     def _plot_stacked(agg_df, metric, ylabel, title, fname, width=10, height=8):
#         fig, ax = plt.subplots(figsize=(width, height))
#         x = np.arange(len(class_order))
#         xt = [cls.replace("_", "\n") for cls in class_order]
#         for i, cls in enumerate(class_order):
#             bottom = 0
#             total = 0
#             for fam, rgba in fam_colors[cls].items():
#                 mask = (agg_df["rep_class"] == cls) & (agg_df["rep_family"] == fam)
#                 if not mask.any():
#                     continue
#                 val = int(agg_df.loc[mask, metric].values[0])
#                 ax.bar(i, val, bottom=bottom, width=0.4, color=rgba, edgecolor="white", linewidth=0.3)
#                 bottom += val
#                 total += val
#             ax.text(i, total, f"{total:_}".replace("_", ","), ha="center", va="bottom")
#         handles = []
#         for cls in class_order:
#             for fam in reversed(list(fam_colors[cls].keys())):
#                 rgba = fam_colors[cls][fam]
#                 mask = (agg_df["rep_class"] == cls) & (agg_df["rep_family"] == fam)
#                 val = int(agg_df.loc[mask, metric].values[0]) if mask.any() else 0
#                 val_to_print = f"{val:_}".replace("_", ",")
#                 handles.append(mpatches.Patch(color=rgba, label=f"{cls} - {fam} ({val_to_print})"))
#         ax.set_xticks(x)
#         ax.set_xticklabels(xt, rotation=0)
#         ax.set_xlabel("Repeat class", weight="bold", labelpad=12)
#         ax.set_ylabel(ylabel, weight="bold", labelpad=12)
#         ax.set_title(title, weight="bold", pad=20)
#         ax.legend(handles=handles, title="Type", loc="upper left", bbox_to_anchor=(1.05, 1), ncol=1, fontsize="small")
#         plt.tight_layout()
#         fig.savefig(f"{fname}", bbox_inches="tight")
#         plt.close(fig)
#
#     def _plot_contiguous(agg_df, metric, ylabel, title, fname, width=width, height=height):
#         fig, ax = plt.subplots(figsize=(width, height))
#         families = []
#         values = []
#         colors = []
#         for cls in class_order:
#             for fam, rgba in fam_colors[cls].items():
#                 mask = (agg_df["rep_class"] == cls) & (agg_df["rep_family"] == fam)
#                 if not mask.any():
#                     continue
#                 val = int(agg_df.loc[mask, metric].values[0])
#                 families.append(f"{cls} - {fam}")
#                 values.append(val)
#                 colors.append(rgba)
#         y = np.arange(len(families))
#         ax.barh(y, values, color=colors, height=0.6)
#         for yi, val in zip(y, values):
#             ax.text(val + max(values) * 0.01, yi, f"{val:_}".replace("_", ","), va="center", ha="left")
#         ax.set_yticks(y)
#         ax.set_yticklabels(families)
#         ax.invert_yaxis()
#         ax.set_xlim(0, max(values) * 1.2)
#         ax.set_xlabel(ylabel, weight="bold", labelpad=12)
#         ax.set_ylabel("Repeat type", weight="bold", labelpad=12)
#         ax.set_title(title, weight="bold", pad=20)
#         legend_handles = [mpatches.Patch(color=class_colors_hex[cls], label=cls) for cls in class_order]
#         ax.legend(handles=legend_handles, title="Repeat class", loc="best", ncol=1, fontsize="small")
#         plt.tight_layout()
#         fig.savefig(f"{fname}", bbox_inches="tight")
#         plt.close(fig)
#
#     _plot_stacked(
#         agg_global,
#         "count",
#         "Insertion count",
#         f"Stacked Counts - Whole Genome - {name}",
#         f"analyses/{name}/{name}_whole_genome_stacked_counts_by_class.pdf",
#         width=width,
#         height=height,
#     )
#     _plot_stacked(
#         agg_global,
#         "length",
#         "Base pair span",
#         f"Stacked Lengths - Whole Genome - {name}",
#         f"analyses/{name}/{name}_whole_genome_stacked_lengths_by_class.pdf",
#         width=width,
#         height=height,
#     )
#     _plot_contiguous(
#         agg_global,
#         "count",
#         "Insertion count",
#         f"Counts by Type - Whole Genome - {name}",
#         f"analyses/{name}/{name}_whole_genome_contiguous_counts_by_class.pdf",
#         width=width,
#         height=height,
#     )
#     _plot_contiguous(
#         agg_global,
#         "length",
#         "Base pair span",
#         f"Lengths by Type - Whole Genome - {name}",
#         f"analyses/{name}/{name}_whole_genome_contiguous_lengths_by_class.pdf",
#         width=width,
#         height=height,
#     )
#
#     if per_chromosome:
#         for chrom, windows in annotations_by_window.items():
#             rec = []
#             for classes in windows.values():
#                 for rep_class, families in classes.items():
#                     for rep_family, metrics in families.items():
#                         rec.append(
#                             {
#                                 "rep_class": rep_class,
#                                 "rep_family": rep_family,
#                                 "count": metrics["count"],
#                                 "length": metrics["length"],
#                             }
#                         )
#             df_chr = pd.DataFrame(rec)
#             if df_chr.empty:
#                 continue
#             agg_chr = df_chr.groupby(["rep_class", "rep_family"])[["count", "length"]].sum().reset_index()
#             _plot_stacked(
#                 agg_chr,
#                 "count",
#                 "Insertion count",
#                 f"Stacked Counts in chromosome {chrom_names[chrom]} - {name}",
#                 f"analyses/{name}/{name}_{chrom_names[chrom]}_stacked_counts_by_class.pdf",
#                 width=width,
#                 height=height,
#             )
#             _plot_stacked(
#                 agg_chr,
#                 "length",
#                 "Base pair span",
#                 f"Stacked Lengths in chromosome {chrom_names[chrom]} - {name}",
#                 f"analyses/{name}/{name}_{chrom_names[chrom]}_stacked_lengths_by_class.pdf",
#                 width=width,
#                 height=height,
#             )
#             _plot_contiguous(
#                 agg_chr,
#                 "count",
#                 "Insertion count",
#                 f"Counts by Type in chromosome {chrom_names[chrom]} - {name}",
#                 f"analyses/{name}/{name}_{chrom_names[chrom]}_contiguous_counts_by_class.pdf",
#                 width=width,
#                 height=height,
#             )
#             _plot_contiguous(
#                 agg_chr,
#                 "length",
#                 "Base pair span",
#                 f"Lengths by Type in chromosome {chrom_names[chrom]} - {name}",
#                 f"analyses/{name}/{name}_{chrom_names[chrom]}_contiguous_lengths_by_class.pdf",
#                 width=width,
#                 height=height,
#             )


def plot_repeat_family_bars(
    name,
    annotations_by_window,
    chrom_names,
    agg_global,
    class_order,
    class_colors_hex,
    fam_colors,
    per_chromosome=False,
    width=10,
    height=8,
):
    """
    Plot stacked and contiguous repeat family bar charts for the genome.
    """
    print("# Running plot_repeat_family_bars function\n")

    os.makedirs(f"analyses/{name}", exist_ok=True)

    class_list = class_order
    fam_colors_map = fam_colors

    def _build_metric_lookup(agg_df, metric):
        lookup = {}
        for _, row in agg_df.iterrows():
            lookup[(row["rep_class"], row["rep_family"])] = int(row[metric])
        return lookup

    def _plot_stacked(agg_df, metric, ylabel, title, fname, width=10, height=8):
        value_lookup = _build_metric_lookup(agg_df, metric)
        fig, ax = plt.subplots(figsize=(width, height))
        x = np.arange(len(class_list))
        xt = [cls.replace("_", "\n") for cls in class_list]
        for i, cls in enumerate(class_list):
            bottom = 0
            total = 0
            for fam, rgba in fam_colors_map[cls].items():
                val = value_lookup.get((cls, fam), 0)
                if val == 0:
                    continue
                ax.bar(i, val, bottom=bottom, width=0.4, color=rgba, edgecolor="white", linewidth=0.3)
                bottom += val
                total += val
            ax.text(i, total, f"{total:_}".replace("_", ","), ha="center", va="bottom")
        handles = []
        for cls in class_list:
            for fam in reversed(list(fam_colors_map[cls].keys())):
                rgba = fam_colors_map[cls][fam]
                val = value_lookup.get((cls, fam), 0)
                val_to_print = f"{val:_}".replace("_", ",")
                handles.append(mpatches.Patch(color=rgba, label=f"{cls} - {fam} ({val_to_print})"))
        ax.set_xticks(x)
        ax.set_xticklabels(xt, rotation=0)
        ax.set_xlabel("Repeat class", weight="bold", labelpad=12)
        ax.set_ylabel(ylabel, weight="bold", labelpad=12)
        ax.set_title(title, weight="bold", pad=20)
        ax.legend(
            handles=handles,
            title="Type",
            loc="upper left",
            bbox_to_anchor=(1.05, 1),
            ncol=1,
            fontsize="small",
        )
        plt.tight_layout()
        fig.savefig(f"{fname}", bbox_inches="tight")
        plt.close(fig)

    def _plot_contiguous(agg_df, metric, ylabel, title, fname, width=width, height=height):
        value_lookup = _build_metric_lookup(agg_df, metric)
        fig, ax = plt.subplots(figsize=(width, height))
        families = []
        values = []
        colors = []
        for cls in class_list:
            for fam, rgba in fam_colors_map[cls].items():
                val = value_lookup.get((cls, fam), 0)
                if val == 0:
                    continue
                families.append(f"{cls} - {fam}")
                values.append(val)
                colors.append(rgba)
        y = np.arange(len(families))
        ax.barh(y, values, color=colors, height=0.6)
        max_val = max(values) if values else 0
        for yi, val in zip(y, values):
            ax.text(val + max_val * 0.01, yi, f"{val:_}".replace("_", ","), va="center", ha="left")
        ax.set_yticks(y)
        ax.set_yticklabels(families)
        ax.invert_yaxis()
        ax.set_xlim(0, max_val * 1.2 if max_val else 1)
        ax.set_xlabel(ylabel, weight="bold", labelpad=12)
        ax.set_ylabel("Repeat type", weight="bold", labelpad=12)
        ax.set_title(title, weight="bold", pad=20)
        legend_handles = [mpatches.Patch(color=class_colors_hex[cls], label=cls) for cls in class_list]
        ax.legend(handles=legend_handles, title="Repeat class", loc="best", ncol=1, fontsize="small")
        plt.tight_layout()
        fig.savefig(f"{fname}", bbox_inches="tight")
        plt.close(fig)

    _plot_stacked(
        agg_global,
        "count",
        "Insertion count",
        f"Stacked Counts - Whole Genome - {name}",
        f"analyses/{name}/{name}_whole_genome_stacked_counts_by_class.pdf",
        width=width,
        height=height,
    )
    _plot_stacked(
        agg_global,
        "length",
        "Base pair span",
        f"Stacked Lengths - Whole Genome - {name}",
        f"analyses/{name}/{name}_whole_genome_stacked_lengths_by_class.pdf",
        width=width,
        height=height,
    )
    _plot_contiguous(
        agg_global,
        "count",
        "Insertion count",
        f"Counts by Type - Whole Genome - {name}",
        f"analyses/{name}/{name}_whole_genome_contiguous_counts_by_class.pdf",
        width=width,
        height=height,
    )
    _plot_contiguous(
        agg_global,
        "length",
        "Base pair span",
        f"Lengths by Type - Whole Genome - {name}",
        f"analyses/{name}/{name}_whole_genome_contiguous_lengths_by_class.pdf",
        width=width,
        height=height,
    )

    if per_chromosome:
        for chrom, windows in annotations_by_window.items():
            rec = []
            for classes in windows.values():
                for rep_class, families in classes.items():
                    for rep_family, metrics in families.items():
                        rec.append(
                            {
                                "rep_class": rep_class,
                                "rep_family": rep_family,
                                "count": metrics["count"],
                                "length": metrics["length"],
                            }
                        )
            df_chr = pd.DataFrame(rec)
            if df_chr.empty:
                continue
            agg_chr = df_chr.groupby(["rep_class", "rep_family"])[["count", "length"]].sum().reset_index()
            _plot_stacked(
                agg_chr,
                "count",
                "Insertion count",
                f"Stacked Counts in chromosome {chrom_names[chrom]} - {name}",
                f"analyses/{name}/{name}_{chrom_names[chrom]}_stacked_counts_by_class.pdf",
                width=width,
                height=height,
            )
            _plot_stacked(
                agg_chr,
                "length",
                "Base pair span",
                f"Stacked Lengths in chromosome {chrom_names[chrom]} - {name}",
                f"analyses/{name}/{name}_{chrom_names[chrom]}_stacked_lengths_by_class.pdf",
                width=width,
                height=height,
            )
            _plot_contiguous(
                agg_chr,
                "count",
                "Insertion count",
                f"Counts by Type in chromosome {chrom_names[chrom]} - {name}",
                f"analyses/{name}/{name}_{chrom_names[chrom]}_contiguous_counts_by_class.pdf",
                width=width,
                height=height,
            )
            _plot_contiguous(
                agg_chr,
                "length",
                "Base pair span",
                f"Lengths by Type in chromosome {chrom_names[chrom]} - {name}",
                f"analyses/{name}/{name}_{chrom_names[chrom]}_contiguous_lengths_by_class.pdf",
                width=width,
                height=height,
            )


# def plot_karyotype_tracks(
#     name,
#     window_size,
#     genome_file,
#     chromosomes,
#     accessibility,
#     gc_content,
#     classes,
#     colors,
#     plot_per_class,
#     kimura_bed=None,
#     identity_bed=None,
# ):
#     print("# Running plot_karyotype_tracks function\n")
#
#     gc_arg = gc_content if gc_content is not None else "not_displayed"
#     accessibility_arg = accessibility if accessibility is not None else "not_displayed"
#     cmd = (
#         f"Rscript scripts/plot_chromosomes_kimura.R --name {name} "
#         f"--genome {genome_file} --chromosome {chromosomes} --accessibility {accessibility_arg} "
#         f"--gccontent {gc_arg} --input analyses/{name}/karyoplot_tables/{name}_{window_size}_repeat_classes.bed "
#         f"--classesorder {classes} --perclass {plot_per_class} --colorsorder {colors} "
#         f"--output analyses/{name}/{name}_{window_size} --kimura {kimura_bed} --identity {identity_bed}"
#     )
#     print(cmd)
#     subprocess.run(cmd.split(" "), check=True)


def plot_karyotype_tracks(
    name,
    window_size,
    genome_file,
    chromosomes,
    accessibility,
    gc_content,
    classes,
    colors,
    plot_per_class,
    kimura_bed=None,
    identity_bed=None,
    layout="horizontal",
    zoom=None,
    zoom_chromosome=None,
    zoom_start=None,
    zoom_end=None,
):
    """
    Run the R script to plot karyotype tracks for repeat annotations.
    """
    print("# Running plot_karyotype_tracks function\n")

    gc_arg = gc_content if gc_content is not None else "not_displayed"
    accessibility_arg = accessibility if accessibility is not None else "not_displayed"
    output_prefix = f"analyses/{name}/{name}_{window_size}"
    input_bed = f"analyses/{name}/karyoplot_tables/{name}_{window_size}_repeat_classes.bed"

    cmd = [
        "Rscript",
        "scripts/televizion/plot_chromosomes_kimura_optimized.R",
        "--name",
        str(name),
        "--genome",
        str(genome_file),
        "--chromosomes-order",
        str(chromosomes),
        "--accessibility",
        str(accessibility_arg),
        "--gc-content",
        str(gc_arg),
        "--classes-table",
        input_bed,
        "--classes-order",
        str(classes),
        "--per-class",
        str(plot_per_class),
        "--colors-order",
        str(colors),
        "--output",
        output_prefix,
        "--kimura-table",
        str(kimura_bed),
        "--identity-table",
        str(identity_bed),
        "--layout",
        str(layout),
    ]
    if zoom is not None:
        cmd.extend(["--zoom", str(zoom)])
    elif zoom_chromosome is not None:
        cmd.extend(
            [
                "--zoom-chromosome",
                str(zoom_chromosome),
                "--zoom-start",
                str(zoom_start),
                "--zoom-end",
                str(zoom_end),
            ]
        )
    print(shlex.join(cmd))
    subprocess.run(cmd, check=True)
