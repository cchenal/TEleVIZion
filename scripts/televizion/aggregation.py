"""
Windowing and aggregation utilities for TEleVIZion.
"""

import os
import pandas as pd
from intervaltree import IntervalTree


# def build_windows(genome_file, window_size, chrom_to_plot):
#     print("\n\n# Running build_windows function\n")
#
#     windows, chrom_names = {}, {}
#     for line in open(genome_file).readlines():
#         if not line.startswith("chr\tstart"):
#             fields = line[:-1].split("\t")
#             chrom, length, name = fields[0], int(fields[2]), fields[3]
#             chrom_names[chrom] = name
#             windows[chrom] = []
#             for start in range(0, length, window_size):
#                 end = min(start + window_size, length)
#                 window = f"{chrom}-{start + 1}-{end}"  # 1-based indexing
#                 windows[chrom].append(window)
#     if chrom_to_plot == "all":
#         tmp = []
#         for chrom in windows:
#             tmp.append(chrom)
#         chroms = ",".join(tmp)
#     else:
#         chroms = chrom_to_plot
#     return (windows, chroms, chrom_names)


def build_windows(genome_file, window_size, chrom_to_plot):
    """
    Build per-chromosome window labels from a genome metadata TSV.

    Returns:
        windows: dict of chromosome -> list of "chrom-start-end" window labels.
        chroms: comma-separated chromosome list for plotting.
        chrom_names: dict mapping chromosome ID to display name.
    """
    print("\n\n# Running build_windows function\n")

    windows, chrom_names = {}, {}
    with open(genome_file, "r") as handle:
        for line in handle:
            if not line or line.startswith("chr\tstart"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            chrom = fields[0]
            length = int(fields[2])
            name = fields[3]
            chrom_names[chrom] = name
            windows[chrom] = [
                f"{chrom}-{start + 1}-{min(start + window_size, length)}"
                for start in range(0, length, window_size)
            ]

    chroms = ",".join(windows.keys()) if chrom_to_plot == "all" else chrom_to_plot
    return (windows, chroms, chrom_names)


# def split_overlaps(window_annotations):
#     """
#     Takes a dictionary of interval entries and returns a dictionary of segments.
#
#     Keys: "start-end" string
#     Values: {
#         'count': number of overlapping intervals,
#         'entries': list of {id, rep_class, rep_family}
#     }
#     """
#     points = set()
#     for item in window_annotations.values():
#         points.add(item["start"])
#         points.add(item["end"])
#     sorted_points = sorted(points)
#
#     tree = IntervalTree()
#     for key, item in window_annotations.items():
#         tree.addi(item["start"], item["end"], key)
#
#     segments = {}
#
#     for i in range(len(sorted_points) - 1):
#         seg_start = sorted_points[i]
#         seg_end = sorted_points[i + 1]
#         if seg_start == seg_end:
#             continue
#
#         overlaps = tree[seg_start:seg_end]
#         key = f"{seg_start}-{seg_end}"
#         segment_data = {
#             "count": len(overlaps),
#             "entries": [],
#         }
#
#         for ov in overlaps:
#             info = window_annotations[ov.data]
#             if "match_identity" in info.keys():
#                 segment_data["entries"].append(
#                     {
#                         "id": ov.data,
#                         "rep_class": info["rep_class"],
#                         "rep_family": info["rep_family"],
#                         "rep_element": info["rep_element"],
#                         "match_identity": info["match_identity"],
#                     }
#                 )
#             else:
#                 segment_data["entries"].append(
#                     {
#                         "id": ov.data,
#                         "rep_class": info["rep_class"],
#                         "rep_family": info["rep_family"],
#                         "rep_element": info["rep_element"],
#                     }
#                 )
#         segments[key] = segment_data
#
#     return segments


def split_overlaps(window_annotations):
    """
    Takes a dictionary of interval entries and returns a dictionary of segments.

    Keys: "start-end" string
    Values: {
        'count': number of overlapping intervals,
        'entries': list of {id, rep_class, rep_family, rep_element, match_identity?}
    }
    """
    annotations = window_annotations
    points = set()
    for item in annotations.values():
        points.update((item["start"], item["end"]))
    sorted_points = sorted(points)

    tree = IntervalTree()
    for key, item in annotations.items():
        tree.addi(item["start"], item["end"], key)

    segments = {}
    for i in range(len(sorted_points) - 1):
        seg_start = sorted_points[i]
        seg_end = sorted_points[i + 1]
        if seg_start == seg_end:
            continue

        overlaps = tree[seg_start:seg_end]
        entries = []
        for ov in overlaps:
            info = annotations[ov.data]
            entry = {
                "id": ov.data,
                "rep_class": info["rep_class"],
                "rep_family": info["rep_family"],
                "rep_element": info["rep_element"],
            }
            if "divergence" in info:
                entry["divergence"] = info["divergence"]
            if "match_identity" in info:
                entry["match_identity"] = info["match_identity"]
            entries.append(entry)

        segments[f"{seg_start}-{seg_end}"] = {
            "count": len(overlaps),
            "entries": entries,
        }

    return segments


# def export_window_class_table(name, window_size, annotations_by_window, windows, class_colors, class_order):
#     print("# Running export_window_class_table function\n")
#
#     os.makedirs(f"analyses/{name}/karyoplot_tables", exist_ok=True)
#
#     output_bed = f"analyses/{name}/karyoplot_tables/{name}_{window_size}_repeat_classes.bed"
#
#     rows = []
#     for chrom, win_list in windows.items():
#         for win_label in win_list:
#             parts = win_label.split("-")
#             if len(parts) != 3:
#                 continue
#             _, start_s, end_s = parts
#             start = int(start_s) - 1
#             end = int(end_s)
#             window_len = end - start
#             window_karyoplotr = int(start + (window_len / 2))
#             row = {
#                 "chrom": chrom,
#                 "start": start,
#                 "end": end,
#                 "barycenter": window_karyoplotr,
#             }
#             for cls in class_order:
#                 fams = annotations_by_window.get(chrom, {}).get(win_label, {}).get(cls, {})
#                 cnt = sum(m["count"] for m in fams.values())
#                 lng = sum(m["length"] for m in fams.values())
#                 pct = (lng / window_len) if window_len > 0 else 0
#                 row[f"{cls}_count"] = cnt
#                 row[f"{cls}_length"] = lng
#                 row[f"{cls}_pct"] = round(pct, 3)
#             rows.append(row)
#
#     df = pd.DataFrame(rows)
#
#     count_cols = [f"{cls}_count" for cls in class_order]
#     stacked = df[count_cols].cumsum(axis=1)
#     stacked.columns = [f"{cls}_count_stacked" for cls in class_order]
#     df = pd.concat([df, stacked], axis=1)
#
#     pct_cols = [f"{cls}_pct" for cls in class_order]
#     stacked = df[pct_cols].cumsum(axis=1)
#     stacked.columns = [f"{cls}_pct_stacked" for cls in class_order]
#     df = pd.concat([df, stacked], axis=1)
#
#     cols = ["chrom", "start", "end", "barycenter"]
#     for cls in class_order:
#         cols.append(f"{cls}_count")
#         cols.append(f"{cls}_count_stacked")
#         cols.append(f"{cls}_length")
#         cols.append(f"{cls}_pct")
#         cols.append(f"{cls}_pct_stacked")
#
#     df = df[cols]
#     df.to_csv(output_bed, sep="\t", index=False)
#
#     reversed_classes = ",".join(list(reversed(class_order)))
#     tmp_colors = []
#     for cls in class_order:
#         tmp_colors.append(class_colors[cls])
#     reversed_colors = ",".join(list(reversed(tmp_colors)))
#     return (reversed_classes, reversed_colors)


def export_window_class_table(name, window_size, annotations_by_window, windows, class_colors, class_order):
    """
    Export per-window class summary metrics as a karyoplot table.

    Returns:
        reversed_classes: comma-separated class order for plotting (reversed).
        reversed_colors: comma-separated class colors matching reversed_classes.
    """
    print("# Running export_window_class_table function\n")

    os.makedirs(f"analyses/{name}/karyoplot_tables", exist_ok=True)

    output_bed = f"analyses/{name}/karyoplot_tables/{name}_{window_size}_repeat_classes.bed"

    rows = []
    annotations = annotations_by_window
    class_list = class_order
    for chrom, win_list in windows.items():
        chrom_windows = annotations.get(chrom, {})
        for win_label in win_list:
            parts = win_label.split("-")
            if len(parts) != 3:
                continue
            _, start_s, end_s = parts
            start = int(start_s) - 1
            end = int(end_s)
            window_len = end - start
            window_karyoplotr = int(start + (window_len / 2))
            row = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "barycenter": window_karyoplotr,
            }
            win_classes = chrom_windows.get(win_label, {})
            for cls in class_list:
                fams = win_classes.get(cls, {})
                cnt = 0
                lng = 0
                for metrics in fams.values():
                    cnt += metrics["count"]
                    lng += metrics["length"]
                pct = (lng / window_len) if window_len > 0 else 0
                row[f"{cls}_count"] = cnt
                row[f"{cls}_length"] = lng
                row[f"{cls}_pct"] = round(pct, 3)
            rows.append(row)

    df = pd.DataFrame(rows)

    count_cols = [f"{cls}_count" for cls in class_list]
    stacked = df[count_cols].cumsum(axis=1)
    stacked.columns = [f"{cls}_count_stacked" for cls in class_list]
    df = pd.concat([df, stacked], axis=1)

    pct_cols = [f"{cls}_pct" for cls in class_list]
    stacked = df[pct_cols].cumsum(axis=1)
    stacked.columns = [f"{cls}_pct_stacked" for cls in class_list]
    df = pd.concat([df, stacked], axis=1)

    cols = ["chrom", "start", "end", "barycenter"]
    for cls in class_list:
        cols.append(f"{cls}_count")
        cols.append(f"{cls}_count_stacked")
        cols.append(f"{cls}_length")
        cols.append(f"{cls}_pct")
        cols.append(f"{cls}_pct_stacked")

    df = df[cols]
    df.to_csv(output_bed, sep="\t", index=False)

    reversed_classes = ",".join(list(reversed(class_list)))
    reversed_colors = ",".join(list(reversed([class_colors[cls] for cls in class_list])))
    return (reversed_classes, reversed_colors)


# def export_window_kimura_table(name, window_size, kimura_bins, windows, class_order):
#     print("# Running export_window_kimura_table function\n")
#
#     os.makedirs(f"analyses/{name}/karyoplot_tables", exist_ok=True)
#     output_bed = f"analyses/{name}/karyoplot_tables/{name}_{window_size}_kimura.bed"
#
#     kimura_cats = ["0-10", "10-20", "20-30", "30-40", "40-70"]
#
#     rows = []
#     for chrom, win_list in windows.items():
#         for win_label in win_list:
#             parts = win_label.split("-")
#             if len(parts) != 3:
#                 continue
#             _, start_s, end_s = parts
#             start = int(start_s) - 1
#             end = int(end_s)
#             window_len = end - start
#             window_karyoplotr = int(start + (window_len / 2))
#             row = {
#                 "chrom": chrom,
#                 "start": start,
#                 "end": end,
#                 "barycenter": window_karyoplotr,
#             }
#             alls = []
#             for category in kimura_cats:
#                 alls.append(f"ALL_{category}_pct")
#                 row[f"ALL_{category}_count"] = 0
#                 row[f"ALL_{category}_pct"] = 0
#
#             for cls in class_order:
#                 cats = kimura_bins.get(chrom, {}).get(win_label, {}).get(cls, {})
#                 if len(cats) > 0:
#                     cnt_sum = sum(m for m in cats.values())
#                     for category in kimura_cats:
#                         row[f"{cls}_{category}_count"] = cats[category]
#                         row[f"ALL_{category}_count"] += cats[category]
#                         if cnt_sum != 0:
#                             row[f"{cls}_{category}_pct"] = round(cats[category] / cnt_sum, 4)
#                         else:
#                             row[f"{cls}_{category}_pct"] = 0
#                 else:
#                     for category in kimura_cats:
#                         row[f"{cls}_{category}_count"] = 0
#                         row[f"{cls}_{category}_pct"] = 0
#
#             cnt_sum = 0
#             for category in kimura_cats:
#                 cnt_sum += row[f"ALL_{category}_count"]
#             for category in kimura_cats:
#                 if cnt_sum != 0:
#                     row[f"ALL_{category}_pct"] = round(row[f"ALL_{category}_count"] / cnt_sum, 4)
#                 else:
#                     row[f"ALL_{category}_pct"] = 0
#
#             rows.append(row)
#
#     df = pd.DataFrame(rows)
#
#     stacked = df[alls].cumsum(axis=1)
#     stacked.columns = [f"ALL_{category}_pct_stacked" for category in kimura_cats]
#     df = pd.concat([df, stacked], axis=1)
#
#     for cls in class_order:
#         to_stack = []
#         for category in kimura_cats:
#             to_stack.append(f"{cls}_{category}_pct")
#         stacked = df[to_stack].cumsum(axis=1)
#         stacked.columns = [f"{cls}_{category}_pct_stacked" for category in kimura_cats]
#         df = pd.concat([df, stacked], axis=1)
#
#     cols = ["chrom", "start", "end", "barycenter"]
#     for category in kimura_cats:
#         cols.append(f"ALL_{category}_count")
#         cols.append(f"ALL_{category}_pct")
#         cols.append(f"ALL_{category}_pct_stacked")
#         df[f"ALL_{category}_pct_stacked"] = df[f"ALL_{category}_pct_stacked"].clip(upper=1)
#     for cls in class_order:
#         for category in kimura_cats:
#             cols.append(f"{cls}_{category}_count")
#             cols.append(f"{cls}_{category}_pct")
#             cols.append(f"{cls}_{category}_pct_stacked")
#             df[f"{cls}_{category}_pct_stacked"] = df[f"{cls}_{category}_pct_stacked"].clip(upper=1)
#
#     df = df[cols]
#     df.to_csv(output_bed, sep="\t", index=False)
#
#     return output_bed


def export_window_kimura_table(name, window_size, kimura_bins, windows, class_order):
    """
    Export per-window Kimura divergence summaries as a karyoplot table.

    Returns:
        output_bed: path to the generated BED-like TSV.
    """
    print("# Running export_window_kimura_table function\n")

    os.makedirs(f"analyses/{name}/karyoplot_tables", exist_ok=True)
    output_bed = f"analyses/{name}/karyoplot_tables/{name}_{window_size}_kimura.bed"

    kimura_cats = ["0-10", "10-20", "20-30", "30-40", "40-70"]
    all_pct_cols = [f"ALL_{category}_pct" for category in kimura_cats]

    rows = []
    class_list = class_order
    for chrom, win_list in windows.items():
        chrom_bins = kimura_bins.get(chrom, {})
        for win_label in win_list:
            parts = win_label.split("-")
            if len(parts) != 3:
                continue
            _, start_s, end_s = parts
            start = int(start_s) - 1
            end = int(end_s)
            window_len = end - start
            window_karyoplotr = int(start + (window_len / 2))
            row = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "barycenter": window_karyoplotr,
            }
            for category in kimura_cats:
                row[f"ALL_{category}_count"] = 0
                row[f"ALL_{category}_pct"] = 0

            win_bins = chrom_bins.get(win_label, {})
            for cls in class_list:
                cats = win_bins.get(cls, {})
                if cats:
                    cnt_sum = 0
                    for value in cats.values():
                        cnt_sum += value
                    for category in kimura_cats:
                        count = cats[category]
                        row[f"{cls}_{category}_count"] = count
                        row[f"ALL_{category}_count"] += count
                        row[f"{cls}_{category}_pct"] = round(count / cnt_sum, 4) if cnt_sum else 0
                else:
                    for category in kimura_cats:
                        row[f"{cls}_{category}_count"] = 0
                        row[f"{cls}_{category}_pct"] = 0

            total_count = 0
            for category in kimura_cats:
                total_count += row[f"ALL_{category}_count"]
            for category in kimura_cats:
                row[f"ALL_{category}_pct"] = (
                    round(row[f"ALL_{category}_count"] / total_count, 4) if total_count else 0
                )

            rows.append(row)

    df = pd.DataFrame(rows)

    stacked = df[all_pct_cols].cumsum(axis=1)
    stacked.columns = [f"ALL_{category}_pct_stacked" for category in kimura_cats]
    df = pd.concat([df, stacked], axis=1)

    for cls in class_list:
        to_stack = [f"{cls}_{category}_pct" for category in kimura_cats]
        stacked = df[to_stack].cumsum(axis=1)
        stacked.columns = [f"{cls}_{category}_pct_stacked" for category in kimura_cats]
        df = pd.concat([df, stacked], axis=1)

    cols = ["chrom", "start", "end", "barycenter"]
    for category in kimura_cats:
        cols.append(f"ALL_{category}_count")
        cols.append(f"ALL_{category}_pct")
        cols.append(f"ALL_{category}_pct_stacked")
        df[f"ALL_{category}_pct_stacked"] = df[f"ALL_{category}_pct_stacked"].clip(upper=1)
    for cls in class_list:
        for category in kimura_cats:
            cols.append(f"{cls}_{category}_count")
            cols.append(f"{cls}_{category}_pct")
            cols.append(f"{cls}_{category}_pct_stacked")
            df[f"{cls}_{category}_pct_stacked"] = df[f"{cls}_{category}_pct_stacked"].clip(upper=1)

    df = df[cols]
    df.to_csv(output_bed, sep="\t", index=False)

    return output_bed


def export_window_divergence_table(name, window_size, divergence_bins, windows, class_order):
    """
    Export per-window RepeatMasker divergence summaries as a karyoplot table.

    Returns:
        output_bed: path to the generated BED-like TSV.
    """
    print("# Running export_window_divergence_table function\n")

    os.makedirs(f"analyses/{name}/karyoplot_tables", exist_ok=True)
    output_bed = f"analyses/{name}/karyoplot_tables/{name}_{window_size}_divergence.bed"

    divergence_cats = ["0-10", "10-20", "20-30", "30-40", "40-70"]
    all_pct_cols = [f"ALL_{category}_pct" for category in divergence_cats]

    rows = []
    class_list = class_order
    for chrom, win_list in windows.items():
        chrom_bins = divergence_bins.get(chrom, {})
        for win_label in win_list:
            parts = win_label.split("-")
            if len(parts) != 3:
                continue
            _, start_s, end_s = parts
            start = int(start_s) - 1
            end = int(end_s)
            window_len = end - start
            window_karyoplotr = int(start + (window_len / 2))
            row = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "barycenter": window_karyoplotr,
            }
            for category in divergence_cats:
                row[f"ALL_{category}_count"] = 0
                row[f"ALL_{category}_pct"] = 0

            win_bins = chrom_bins.get(win_label, {})
            for cls in class_list:
                cats = win_bins.get(cls, {})
                if cats:
                    cnt_sum = 0
                    for value in cats.values():
                        cnt_sum += value
                    for category in divergence_cats:
                        count = cats.get(category, 0)
                        row[f"{cls}_{category}_count"] = count
                        row[f"ALL_{category}_count"] += count
                        row[f"{cls}_{category}_pct"] = round(count / cnt_sum, 4) if cnt_sum else 0
                else:
                    for category in divergence_cats:
                        row[f"{cls}_{category}_count"] = 0
                        row[f"{cls}_{category}_pct"] = 0

            total_count = 0
            for category in divergence_cats:
                total_count += row[f"ALL_{category}_count"]
            for category in divergence_cats:
                row[f"ALL_{category}_pct"] = (
                    round(row[f"ALL_{category}_count"] / total_count, 4) if total_count else 0
                )

            rows.append(row)

    df = pd.DataFrame(rows)

    stacked = df[all_pct_cols].cumsum(axis=1)
    stacked.columns = [f"ALL_{category}_pct_stacked" for category in divergence_cats]
    df = pd.concat([df, stacked], axis=1)

    for cls in class_list:
        to_stack = [f"{cls}_{category}_pct" for category in divergence_cats]
        stacked = df[to_stack].cumsum(axis=1)
        stacked.columns = [f"{cls}_{category}_pct_stacked" for category in divergence_cats]
        df = pd.concat([df, stacked], axis=1)

    cols = ["chrom", "start", "end", "barycenter"]
    for category in divergence_cats:
        cols.append(f"ALL_{category}_count")
        cols.append(f"ALL_{category}_pct")
        cols.append(f"ALL_{category}_pct_stacked")
        df[f"ALL_{category}_pct_stacked"] = df[f"ALL_{category}_pct_stacked"].clip(upper=1)
    for cls in class_list:
        for category in divergence_cats:
            cols.append(f"{cls}_{category}_count")
            cols.append(f"{cls}_{category}_pct")
            cols.append(f"{cls}_{category}_pct_stacked")
            df[f"{cls}_{category}_pct_stacked"] = df[f"{cls}_{category}_pct_stacked"].clip(upper=1)

    df = df[cols]
    df.to_csv(output_bed, sep="\t", index=False)

    return output_bed


# def export_window_identity_table(name, window_size, annotations_by_window, windows, class_order):
#     print("# Running export_window_identity_table function\n")
#
#     output_bed = f"analyses/{name}/karyoplot_tables/{name}_{window_size}_identity.bed"
#     os.makedirs(f"analyses/{name}/karyoplot_tables", exist_ok=True)
#
#     id_cats = ["1-0.9", "0.9-0.8", "0.8-0.7", "0.7-0.6", "0.6-0"]
#
#     identity_bins = {}
#
#     for chrom in annotations_by_window:
#         if chrom not in identity_bins:
#             identity_bins[chrom] = {}
#         for window in annotations_by_window[chrom]:
#             if window not in identity_bins[chrom]:
#                 identity_bins[chrom][window] = {}
#             for rep_class in annotations_by_window[chrom][window]:
#                 if rep_class not in identity_bins[chrom][window]:
#                     identity_bins[chrom][window][rep_class] = {
#                         "1-0.9": 0,
#                         "0.9-0.8": 0,
#                         "0.8-0.7": 0,
#                         "0.7-0.6": 0,
#                         "0.6-0": 0,
#                     }
#                 for rep_family in annotations_by_window[chrom][window][rep_class]:
#                     for i in annotations_by_window[chrom][window][rep_class][rep_family]["identity"]:
#                         if i >= 0.9:
#                             identity_bins[chrom][window][rep_class]["1-0.9"] += 1
#                         elif i >= 0.8:
#                             identity_bins[chrom][window][rep_class]["0.9-0.8"] += 1
#                         elif i >= 0.7:
#                             identity_bins[chrom][window][rep_class]["0.8-0.7"] += 1
#                         elif i >= 0.6:
#                             identity_bins[chrom][window][rep_class]["0.7-0.6"] += 1
#                         else:
#                             identity_bins[chrom][window][rep_class]["0.6-0"] += 1
#
#     rows = []
#     for chrom, win_list in windows.items():
#         for win_label in win_list:
#             parts = win_label.split("-")
#             if len(parts) != 3:
#                 continue
#             _, start_s, end_s = parts
#             start = int(start_s) - 1
#             end = int(end_s)
#             window_len = end - start
#             window_karyoplotr = int(start + (window_len / 2))
#             row = {
#                 "chrom": chrom,
#                 "start": start,
#                 "end": end,
#                 "barycenter": window_karyoplotr,
#             }
#             alls = []
#             for category in id_cats:
#                 alls.append(f"ALL_{category}_pct")
#                 row[f"ALL_{category}_count"] = 0
#                 row[f"ALL_{category}_pct"] = 0
#
#             for cls in class_order:
#                 cats = identity_bins.get(chrom, {}).get(win_label, {}).get(cls, {})
#                 if len(cats) > 0:
#                     cnt_sum = sum(m for m in cats.values())
#                     for category in id_cats:
#                         row[f"{cls}_{category}_count"] = cats[category]
#                         row[f"ALL_{category}_count"] += cats[category]
#                         if cnt_sum != 0:
#                             row[f"{cls}_{category}_pct"] = round(cats[category] / cnt_sum, 4)
#                         else:
#                             row[f"{cls}_{category}_pct"] = 0
#                 else:
#                     for category in id_cats:
#                         row[f"{cls}_{category}_count"] = 0
#                         row[f"{cls}_{category}_pct"] = 0
#
#             cnt_sum = 0
#             for category in id_cats:
#                 cnt_sum += row[f"ALL_{category}_count"]
#             for category in id_cats:
#                 if cnt_sum != 0:
#                     row[f"ALL_{category}_pct"] = round(row[f"ALL_{category}_count"] / cnt_sum, 4)
#                 else:
#                     row[f"ALL_{category}_pct"] = 0
#
#             rows.append(row)
#
#     df = pd.DataFrame(rows)
#
#     stacked = df[alls].cumsum(axis=1)
#     stacked.columns = [f"ALL_{category}_pct_stacked" for category in id_cats]
#     df = pd.concat([df, stacked], axis=1)
#
#     for cls in class_order:
#         to_stack = []
#         for category in id_cats:
#             to_stack.append(f"{cls}_{category}_pct")
#         stacked = df[to_stack].cumsum(axis=1)
#         stacked.columns = [f"{cls}_{category}_pct_stacked" for category in id_cats]
#         df = pd.concat([df, stacked], axis=1)
#
#     cols = ["chrom", "start", "end", "barycenter"]
#     for category in id_cats:
#         cols.append(f"ALL_{category}_count")
#         cols.append(f"ALL_{category}_pct")
#         cols.append(f"ALL_{category}_pct_stacked")
#         df[f"ALL_{category}_pct_stacked"] = df[f"ALL_{category}_pct_stacked"].clip(upper=1)
#     for cls in class_order:
#         for category in id_cats:
#             cols.append(f"{cls}_{category}_count")
#             cols.append(f"{cls}_{category}_pct")
#             cols.append(f"{cls}_{category}_pct_stacked")
#             df[f"{cls}_{category}_pct_stacked"] = df[f"{cls}_{category}_pct_stacked"].clip(upper=1)
#
#     df = df[cols]
#     df.to_csv(output_bed, sep="\t", index=False)
#
#     return output_bed


def export_window_identity_table(name, window_size, annotations_by_window, windows, class_order):
    """
    Export per-window identity summaries as a karyoplot table.

    Returns:
        output_bed: path to the generated BED-like TSV.
    """
    print("# Running export_window_identity_table function\n")

    output_bed = f"analyses/{name}/karyoplot_tables/{name}_{window_size}_identity.bed"
    os.makedirs(f"analyses/{name}/karyoplot_tables", exist_ok=True)

    id_cats = ["1-0.9", "0.9-0.8", "0.8-0.7", "0.7-0.6", "0.6-0"]
    all_pct_cols = [f"ALL_{category}_pct" for category in id_cats]

    identity_bins = None
    annotations = annotations_by_window

    for chrom, win_map in annotations.items():
        for win_label, class_map in win_map.items():
            for rep_class, rep_value in class_map.items():
                if isinstance(rep_value, dict) and rep_value:
                    if all(category in rep_value for category in id_cats):
                        identity_bins = annotations
                    else:
                        identity_bins = {}
                    break
            if identity_bins is not None:
                break
        if identity_bins is not None:
            break

    if identity_bins is None:
        identity_bins = {}
        for chrom in annotations:
            if chrom not in identity_bins:
                identity_bins[chrom] = {}
            for window in annotations[chrom]:
                if window not in identity_bins[chrom]:
                    identity_bins[chrom][window] = {}
                for rep_class in annotations[chrom][window]:
                    if rep_class not in identity_bins[chrom][window]:
                        identity_bins[chrom][window][rep_class] = {
                            "1-0.9": 0,
                            "0.9-0.8": 0,
                            "0.8-0.7": 0,
                            "0.7-0.6": 0,
                            "0.6-0": 0,
                        }
                    for rep_family in annotations[chrom][window][rep_class]:
                        identities = annotations[chrom][window][rep_class][rep_family]["identity"]
                        for ident in identities:
                            if ident >= 0.9:
                                identity_bins[chrom][window][rep_class]["1-0.9"] += 1
                            elif ident >= 0.8:
                                identity_bins[chrom][window][rep_class]["0.9-0.8"] += 1
                            elif ident >= 0.7:
                                identity_bins[chrom][window][rep_class]["0.8-0.7"] += 1
                            elif ident >= 0.6:
                                identity_bins[chrom][window][rep_class]["0.7-0.6"] += 1
                            else:
                                identity_bins[chrom][window][rep_class]["0.6-0"] += 1

    rows = []
    class_list = class_order
    for chrom, win_list in windows.items():
        chrom_bins = identity_bins.get(chrom, {})
        for win_label in win_list:
            parts = win_label.split("-")
            if len(parts) != 3:
                continue
            _, start_s, end_s = parts
            start = int(start_s) - 1
            end = int(end_s)
            window_len = end - start
            window_karyoplotr = int(start + (window_len / 2))
            row = {
                "chrom": chrom,
                "start": start,
                "end": end,
                "barycenter": window_karyoplotr,
            }
            for category in id_cats:
                row[f"ALL_{category}_count"] = 0
                row[f"ALL_{category}_pct"] = 0

            win_bins = chrom_bins.get(win_label, {})
            for cls in class_list:
                cats = win_bins.get(cls, {})
                if cats:
                    cnt_sum = 0
                    for value in cats.values():
                        cnt_sum += value
                    for category in id_cats:
                        count = cats[category]
                        row[f"{cls}_{category}_count"] = count
                        row[f"ALL_{category}_count"] += count
                        row[f"{cls}_{category}_pct"] = round(count / cnt_sum, 4) if cnt_sum else 0
                else:
                    for category in id_cats:
                        row[f"{cls}_{category}_count"] = 0
                        row[f"{cls}_{category}_pct"] = 0

            total_count = 0
            for category in id_cats:
                total_count += row[f"ALL_{category}_count"]
            for category in id_cats:
                row[f"ALL_{category}_pct"] = (
                    round(row[f"ALL_{category}_count"] / total_count, 4) if total_count else 0
                )

            rows.append(row)

    df = pd.DataFrame(rows)

    stacked = df[all_pct_cols].cumsum(axis=1)
    stacked.columns = [f"ALL_{category}_pct_stacked" for category in id_cats]
    df = pd.concat([df, stacked], axis=1)

    for cls in class_list:
        to_stack = [f"{cls}_{category}_pct" for category in id_cats]
        stacked = df[to_stack].cumsum(axis=1)
        stacked.columns = [f"{cls}_{category}_pct_stacked" for category in id_cats]
        df = pd.concat([df, stacked], axis=1)

    cols = ["chrom", "start", "end", "barycenter"]
    for category in id_cats:
        cols.append(f"ALL_{category}_count")
        cols.append(f"ALL_{category}_pct")
        cols.append(f"ALL_{category}_pct_stacked")
        df[f"ALL_{category}_pct_stacked"] = df[f"ALL_{category}_pct_stacked"].clip(upper=1)
    for cls in class_list:
        for category in id_cats:
            cols.append(f"{cls}_{category}_count")
            cols.append(f"{cls}_{category}_pct")
            cols.append(f"{cls}_{category}_pct_stacked")
            df[f"{cls}_{category}_pct_stacked"] = df[f"{cls}_{category}_pct_stacked"].clip(upper=1)

    df = df[cols]
    df.to_csv(output_bed, sep="\t", index=False)

    return output_bed
