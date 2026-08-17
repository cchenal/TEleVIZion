"""
I/O helpers: parse annotation files and read genome metadata.
"""

import os
import subprocess
import numpy as np

from .aggregation import split_overlaps

def compute_gc_content(fasta_file=None, gc_windows=10000):
    """
    Compute or reuse GC content windows from a FASTA file.

    Returns:
        gc_file: path to the GC windows TSV, or None if not computed.
    """
    print("# Running compute_gc_content function\n")

    if fasta_file is None:
        return None

    if not os.path.isfile(fasta_file):
        print(
            f"Warning! Provided .fasta file ({fasta_file}) doesn't exist -- GC content has not been calculated\n"
        )
        return None

    gc_file = f"{fasta_file[:-6]}_gc_windows_{gc_windows}.tsv"
    if not os.path.isfile(gc_file):
        cmd = [
            "python3",
            "scripts/utils/create_gc_content.py",
            fasta_file,
            "--window",
            str(gc_windows),
            "--out",
            gc_file,
        ]
        subprocess.run(cmd, check=True)
    else:
        print(f"GC file already exists! If you want to update it, please delete {gc_file}\n")

    return gc_file


# def parse_repeatmasker_kimura_bins(repeatmasker_kimura_file):
#     print("# Running parse_repeatmasker_kimura_bins function\n")
#
#     kimura_div = {}
#     for line in open(repeatmasker_kimura_file).readlines():
#         pseudo_fields = line[:-1].split()
#         fields = []
#         for f in pseudo_fields:
#             if len(f) > 0:
#                 fields.append(f)
#         if len(fields) > 0:
#             if (
#                 (not fields[0].isdigit())
#                 and (fields[0] not in ["Jukes/Cantor", "=================================================================", "File:", "Weighted", "Class", "-----", "Coverage", "Div"])
#             ):
#                 if fields[4] == "----":
#                     fields[4] = 0
#                 value = float(fields[4])
#                 if value < 10:
#                     cat = "0-10"
#                 elif value < 20:
#                     cat = "10-20"
#                 elif value < 30:
#                     cat = "20-30"
#                 elif value < 40:
#                     cat = "30-40"
#                 else:
#                     cat = "40-70"
#                 kimura_div[fields[1]] = {"pct_div": value, "category": cat}
#     return kimura_div


def parse_repeatmasker_kimura_bins(repeatmasker_kimura_file):
    """
    Parse RepeatMasker Kimura .divsum file into per-element bins.

    Returns:
        kimura_div: dict mapping element ID -> {pct_div, category}.
    """
    print("# Running parse_repeatmasker_kimura_bins function\n")

    kimura_div = {}
    headers = {
        "Jukes/Cantor",
        "=================================================================",
        "File:",
        "Weighted",
        "Class",
        "-----",
        "Coverage",
        "Div",
    }
    with open(repeatmasker_kimura_file, "r") as handle:
        for line in handle:
            fields = line.split()
            if not fields or fields[0].isdigit() or fields[0] in headers:
                continue
            if len(fields) <= 4:
                continue
            pct_raw = fields[4]
            if pct_raw == "----":
                pct_raw = "0"
            value = float(pct_raw)
            if value < 10:
                cat = "0-10"
            elif value < 20:
                cat = "10-20"
            elif value < 30:
                cat = "20-30"
            elif value < 40:
                cat = "30-40"
            else:
                cat = "40-70"
            kimura_div[fields[1]] = {"pct_div": value, "category": cat}
    return kimura_div


# def parse_repeatmasker_annotations(repeatmasker_file, windows_dict, kimura_dict=None):
#     print("# Running parse_repeatmasker_annotations function\n")
#
#     repeats = {}
#     for line in open(repeatmasker_file).readlines():
#         pseudo_fields = line[:-1].split()
#         fields = []
#         for f in pseudo_fields:
#             if len(f) > 0:
#                 fields.append(f)
#         if len(fields) > 1:
#             if fields[0][0].isnumeric():
#                 classification = fields[10].split("/")
#                 rep_class = classification[0]
#                 if len(classification) == 1:
#                     rep_family = "NA"
#                 else:
#                     if classification[1] == "":
#                         rep_family = "NA"
#                     else:
#                         rep_family = classification[1]
#                 if rep_class not in repeats:
#                     repeats[rep_class] = {}
#                 if rep_family not in repeats[rep_class]:
#                     repeats[rep_class][rep_family] = 0
#                 repeats[rep_class][rep_family] += 1
#
#     insertions = {}
#     for chrom in windows_dict:
#         insertions[chrom] = {}
#         for window in windows_dict[chrom]:
#             insertions[chrom][window] = {}
#             for rep_class in repeats:
#                 insertions[chrom][window][rep_class] = {}
#                 for rep_family in repeats[rep_class]:
#                     insertions[chrom][window][rep_class][rep_family] = {
#                         "count": 0,
#                         "length": 0,
#                         "divergence": [],
#                     }
#             insertions[chrom][window]["Ambiguous"] = {
#                 "NA": {"count": 0, "length": 0, "divergence": []}
#             }
#
#     if kimura_dict is not None:
#         kimura_div = {}
#         for chrom in windows_dict:
#             kimura_div[chrom] = {}
#             for window in windows_dict[chrom]:
#                 kimura_div[chrom][window] = {}
#                 for rep_class in repeats:
#                     kimura_div[chrom][window][rep_class] = {
#                         "0-10": 0,
#                         "10-20": 0,
#                         "20-30": 0,
#                         "30-40": 0,
#                         "40-70": 0,
#                     }
#
#     for chromosome in insertions:
#         tmp_insertions, counter = {}, 0
#         for line in open(repeatmasker_file).readlines():
#             pseudo_fields = line[:-1].split()
#             fields = []
#             for f in pseudo_fields:
#                 if len(f) > 0:
#                     fields.append(f)
#             if len(fields) > 1:
#                 if fields[0][0].isnumeric():
#                     chrom = fields[4]
#                     if chrom == chromosome:
#                         start = min(int(fields[5]), int(fields[6]))
#                         end = max(int(fields[5]), int(fields[6]))
#                         classification = fields[10].split("/")
#                         rep_class = classification[0]
#                         if len(classification) == 1:
#                             rep_family = "NA"
#                         else:
#                             if classification[1] == "":
#                                 rep_family = "NA"
#                             else:
#                                 rep_family = classification[1]
#                         tmp_insertions[str(counter)] = {
#                             "chrom": chrom,
#                             "start": int(start),
#                             "end": int(end),
#                             "rep_class": rep_class,
#                             "rep_family": rep_family,
#                             "rep_element": fields[9],
#                             "divergence": float(fields[1]),
#                         }
#                         counter += 1
#         for window in insertions[chromosome]:
#             coordinates = window.split("-")
#             window_insertions = {}
#             for tmp in tmp_insertions:
#                 if (
#                     tmp_insertions[tmp]["start"] < int(coordinates[2])
#                     and tmp_insertions[tmp]["end"] > int(coordinates[1])
#                 ):
#                     min_bp_to_consider = max(tmp_insertions[tmp]["start"], int(coordinates[1]))
#                     max_bp_to_consider = min(tmp_insertions[tmp]["end"], int(coordinates[2]))
#                     window_insertions[tmp] = {
#                         "chrom": tmp_insertions[tmp]["chrom"],
#                         "start": min_bp_to_consider,
#                         "end": max_bp_to_consider,
#                         "rep_class": tmp_insertions[tmp]["rep_class"],
#                         "rep_family": tmp_insertions[tmp]["rep_family"],
#                         "rep_element": tmp_insertions[tmp]["rep_element"],
#                         "divergence": tmp_insertions[tmp]["divergence"],
#                     }
#
#             segments = split_overlaps(window_annotations=window_insertions)
#             visited = []
#             for region, data in segments.items():
#                 r_start, r_end = region.split("-")
#                 len_bp_to_consider = int(r_end) - int(r_start)
#                 if data["count"] == 1:
#                     for entry in data["entries"]:
#                         if entry["id"] not in visited:
#                             insertions[chromosome][window][entry["rep_class"]][entry["rep_family"]][
#                                 "count"
#                             ] += 1
#                             visited.append(entry["id"])
#                             if (kimura_dict is not None) and (entry["rep_element"] in kimura_dict):
#                                 kimura_div[chromosome][window][entry["rep_class"]][
#                                     kimura_dict[entry["rep_element"]]["category"]
#                                 ] += 1
#                         insertions[chromosome][window][entry["rep_class"]][entry["rep_family"]][
#                             "length"
#                         ] += len_bp_to_consider
#                 elif data["count"] >= 2:
#                     rep_class_comp, rep_family_comp, divs = [], [], []
#                     for entry in data["entries"]:
#                         if entry["rep_class"] not in rep_class_comp:
#                             rep_class_comp.append(entry["rep_class"])
#                         if entry["rep_family"] not in rep_family_comp:
#                             rep_family_comp.append(entry["rep_family"])
#                         if entry["id"] not in visited:
#                             visited.append(entry["id"])
#                     if (len(rep_class_comp) == 1) and (len(rep_family_comp) == 1):
#                         insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]][
#                             "count"
#                         ] += 1
#                         if (kimura_dict is not None) and (entry["rep_element"] in kimura_dict):
#                             kimura_div[chromosome][window][entry["rep_class"]][
#                                 kimura_dict[entry["rep_element"]]["category"]
#                             ] += 1
#                         insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]][
#                             "length"
#                         ] += len_bp_to_consider
#                     else:
#                         insertions[chromosome][window]["Ambiguous"]["NA"]["count"] += 1
#                         if (kimura_dict is not None) and (entry["rep_element"] in kimura_dict):
#                             kimura_div[chromosome][window][entry["rep_class"]][
#                                 kimura_dict[entry["rep_element"]]["category"]
#                             ] += 1
#                         insertions[chromosome][window]["Ambiguous"]["NA"]["length"] += len_bp_to_consider
#
#     if kimura_dict is not None:
#         return (repeats, insertions, kimura_div)
#     return (repeats, insertions)

# def divergence_category(value):
#     """
#     Return the divergence category for a numeric value.

#     Categories:
#       0–10, 10–20, 20–30, 30–40, 40–70
#     """
#     if 0 <= value < 10:
#         return "0-10"
#     elif value < 20:
#         return "10-20"
#     elif value < 30:
#         return "20-30"
#     elif value < 40:
#         return "30-40"
#     else:
#         return "40-70"

def id_category(value):
    """
    Return the identity category for a numeric value.

    Categories:
      1-0.9, 0.9-0.8, 0.8-0.7, 0.7-0.6, 0.6-0
    """
    if value >= 0.9:
        return "1-0.9"
    elif value >= 0.8:
        return "0.9-0.8"
    elif value >= 0.7:
        return "0.8-0.7"
    elif value >= 0.6:
        return "0.7-0.6"
    else:
        return "0.6-0"

def parse_repeatmasker_annotations(repeatmasker_file, windows_dict, kimura_dict=None):
    """
    Parse RepeatMasker .out annotations into per-window insertion summaries.

    Returns:
        repeats: dict of repeat class -> family counts.
        insertions: per-window insertion metrics.
        kimura_div: optional per-window Kimura bins (if kimura_dict provided).
    """
    print("# Running parse_repeatmasker_annotations function\n")

    repeats = {}
    entries_by_chrom = {}
    counters = {}

    with open(repeatmasker_file, "r") as handle:
        for line in handle:
            fields = line.split()
            if len(fields) <= 1 or not fields[0][0].isnumeric():
                continue
            classification = fields[10].split("/")
            rep_class = classification[0]
            if len(classification) == 1 or classification[1] == "":
                rep_family = "NA"
            else:
                rep_family = classification[1]

            if rep_class not in repeats:
                repeats[rep_class] = {}
            if rep_family not in repeats[rep_class]:
                repeats[rep_class][rep_family] = 0
            repeats[rep_class][rep_family] += 1

            chrom = fields[4]
            start = min(int(fields[5]), int(fields[6]))
            end = max(int(fields[5]), int(fields[6]))
            counter = counters.get(chrom, 0)
            counters[chrom] = counter + 1
            entries_by_chrom.setdefault(chrom, []).append(
                {
                    "id": str(counter),
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "rep_class": rep_class,
                    "rep_family": rep_family,
                    "rep_element": fields[9],
                    "divergence": float(fields[1]),
                }
            )

    insertions = {}
    for chrom in windows_dict:
        insertions[chrom] = {}
        for window in windows_dict[chrom]:
            insertions[chrom][window] = {}
            for rep_class in repeats:
                insertions[chrom][window][rep_class] = {}
                for rep_family in repeats[rep_class]:
                    insertions[chrom][window][rep_class][rep_family] = {
                        "count": 0,
                        "length": 0
                    }
            insertions[chrom][window]["Ambiguous"] = {
                "NA": {"count": 0, "length": 0}
            }

    if kimura_dict is not None:
        kimura_div = {}
        for chrom in windows_dict:
            kimura_div[chrom] = {}
            for window in windows_dict[chrom]:
                kimura_div[chrom][window] = {}
                for rep_class in repeats:
                    kimura_div[chrom][window][rep_class] = {
                        "0-10": 0,
                        "10-20": 0,
                        "20-30": 0,
                        "30-40": 0,
                        "40-70": 0,
                    }
        missing_kimura_elements = set()
    else:
        # repeatmasker_div = {}
        # for chrom in windows_dict:
        #     repeatmasker_div[chrom] = {}
        #     for window in windows_dict[chrom]:
        #         repeatmasker_div[chrom][window] = {}
        #         for rep_class in repeats:
        #             repeatmasker_div[chrom][window][rep_class] = {
        #                 "0-10": 0,
        #                 "10-20": 0,
        #                 "20-30": 0,
        #                 "30-40": 0,
        #                 "40-70": 0,
        #             }
        repeatmasker_id = {}
        for chrom in windows_dict:
            repeatmasker_id[chrom] = {}
            for window in windows_dict[chrom]:
                repeatmasker_id[chrom][window] = {}
                for rep_class in repeats:
                    repeatmasker_id[chrom][window][rep_class] = {
                        "1-0.9": 0,
                        "0.9-0.8": 0,
                        "0.8-0.7": 0,
                        "0.7-0.6": 0,
                        "0.6-0": 0,
                    }
                repeatmasker_id[chrom][window]["Ambiguous"] = {
                    "1-0.9": 0,
                    "0.9-0.8": 0,
                    "0.8-0.7": 0,
                    "0.7-0.6": 0,
                    "0.6-0": 0,
                }

    def add_age_bin(chromosome, window, entry, target_class=None):
        rep_class = target_class or entry["rep_class"]
        if kimura_dict is not None:
            kimura_entry = kimura_dict.get(entry["rep_element"])
            if kimura_entry is None:
                missing_kimura_elements.add(entry["rep_element"])
                return
            kimura_div[chromosome][window][rep_class][
                kimura_entry["category"]
            ] += 1
            return

        repeatmasker_id[chromosome][window][rep_class][
            id_category((100 - float(entry["divergence"])) / 100)
        ] += 1

    for chromosome in insertions:
        chrom_entries = entries_by_chrom.get(chromosome, [])
        for window in insertions[chromosome]:
            coordinates = window.split("-")
            window_insertions = {}
            win_start = int(coordinates[1])
            win_end = int(coordinates[2])
            for entry in chrom_entries:
                if entry["start"] < win_end and entry["end"] > win_start:
                    min_bp_to_consider = max(entry["start"], win_start)
                    max_bp_to_consider = min(entry["end"], win_end)
                    window_insertions[entry["id"]] = {
                        "chrom": entry["chrom"],
                        "start": min_bp_to_consider,
                        "end": max_bp_to_consider,
                        "rep_class": entry["rep_class"],
                        "rep_family": entry["rep_family"],
                        "rep_element": entry["rep_element"],
                        "divergence": entry["divergence"],
                    }

            segments = split_overlaps(window_annotations=window_insertions)
            visited = []
            for region, data in segments.items():
                r_start, r_end = region.split("-")
                len_bp_to_consider = int(r_end) - int(r_start)
                if data["count"] == 1:
                    for entry in data["entries"]:
                        # print(entry, data["entries"])
                        if entry["id"] not in visited:
                            insertions[chromosome][window][entry["rep_class"]][entry["rep_family"]][
                                "count"
                            ] += 1
                            visited.append(entry["id"])
                            add_age_bin(chromosome, window, entry)
                        insertions[chromosome][window][entry["rep_class"]][entry["rep_family"]][
                            "length"
                        ] += len_bp_to_consider
                elif data["count"] >= 2:
                    rep_class_comp, rep_family_comp = [], []
                    for entry in data["entries"]:
                        if entry["rep_class"] not in rep_class_comp:
                            rep_class_comp.append(entry["rep_class"])
                        if entry["rep_family"] not in rep_family_comp:
                            rep_family_comp.append(entry["rep_family"])
                        if entry["id"] not in visited:
                            visited.append(entry["id"])
                    if (len(rep_class_comp) == 1) and (len(rep_family_comp) == 1):
                        insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]][
                            "count"
                        ] += 1
                        add_age_bin(chromosome, window, entry)
                        insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]][
                            "length"
                        ] += len_bp_to_consider
                    else:
                        insertions[chromosome][window]["Ambiguous"]["NA"]["count"] += 1
                        if kimura_dict is not None:
                            add_age_bin(chromosome, window, entry)
                        else:
                            add_age_bin(chromosome, window, entry, target_class="Ambiguous")
                        insertions[chromosome][window]["Ambiguous"]["NA"]["length"] += len_bp_to_consider

    if kimura_dict is not None:
        if missing_kimura_elements:
            print(
                "Warning! Kimura bins were missing for "
                f"{len(missing_kimura_elements)} repeat element(s); those insertions "
                "were excluded from the Kimura panel.\n"
            )
        return (repeats, insertions, kimura_div)
    else:
        return (repeats, insertions, repeatmasker_id)


# def parse_edta_annotations(edta_file, windows_dict):
#     print("# Running parse_edta_annotations function\n")
#
#     repeats = {}
#     for line in open(edta_file).readlines():
#         if line[0] != "#":
#             fields = line[:-1].split()
#             classification = fields[8].split("classification=")[1].split(";")[0].split("/")
#             rep_class = classification[0]
#             if len(classification) == 1:
#                 rep_family = "NA"
#             else:
#                 if classification[1] == "":
#                     rep_family = "NA"
#                 else:
#                     rep_family = classification[1]
#             if rep_class not in repeats:
#                 repeats[rep_class] = {}
#             if rep_family not in repeats[rep_class]:
#                 repeats[rep_class][rep_family] = 0
#             repeats[rep_class][rep_family] += 1
#
#     insertions = {}
#     for chrom in windows_dict:
#         insertions[chrom] = {}
#         for window in windows_dict[chrom]:
#             insertions[chrom][window] = {}
#             for rep_class in repeats:
#                 insertions[chrom][window][rep_class] = {}
#                 for rep_family in repeats[rep_class]:
#                     insertions[chrom][window][rep_class][rep_family] = {
#                         "count": 0,
#                         "length": 0,
#                         "identity": [],
#                     }
#             insertions[chrom][window]["Ambiguous"] = {
#                 "NA": {"count": 0, "length": 0, "identity": []}
#             }
#     for chromosome in insertions:
#         chrom_insertions, counter = {}, 0
#         for line in open(edta_file).readlines():
#             if line[0] != "#":
#                 fields = line[:-1].split("\t")
#                 chrom = fields[0]
#                 if chrom == chromosome:
#                     start = min(int(fields[3]), int(fields[4]))
#                     end = max(int(fields[3]), int(fields[4]))
#                     classification = fields[8].split("classification=")[1].split(";")[0].split("/")
#                     rep_class = classification[0]
#                     if len(classification) == 1:
#                         rep_family = "NA"
#                     else:
#                         if classification[1] == "":
#                             rep_family = "NA"
#                         else:
#                             rep_family = classification[1]
#                     rep_element = fields[8].split("Name=")[1].split(";")[0]
#                     identity = fields[8].split("identity=")[1].split(";")[0]
#                     if identity != "NA":
#                         identity = float(identity)
#                     chrom_insertions[str(counter)] = {
#                         "chrom": chrom,
#                         "start": int(start),
#                         "end": int(end),
#                         "rep_class": rep_class,
#                         "rep_family": rep_family,
#                         "rep_element": rep_element,
#                         "match_identity": identity,
#                     }
#                     counter += 1
#         for window in insertions[chromosome]:
#             coordinates = window.split("-")
#             window_insertions = {}
#             for tmp in chrom_insertions:
#                 if (
#                     chrom_insertions[tmp]["start"] < int(coordinates[2])
#                     and chrom_insertions[tmp]["end"] > int(coordinates[1])
#                 ):
#                     min_bp_to_consider = max(chrom_insertions[tmp]["start"], int(coordinates[1]))
#                     max_bp_to_consider = min(chrom_insertions[tmp]["end"], int(coordinates[2]))
#                     window_insertions[tmp] = {
#                         "chrom": chrom_insertions[tmp]["chrom"],
#                         "start": min_bp_to_consider,
#                         "end": max_bp_to_consider,
#                         "rep_class": chrom_insertions[tmp]["rep_class"],
#                         "rep_family": chrom_insertions[tmp]["rep_family"],
#                         "rep_element": chrom_insertions[tmp]["rep_element"],
#                         "match_identity": chrom_insertions[tmp]["match_identity"],
#                     }
#             segments = split_overlaps(window_annotations=window_insertions)
#             visited = []
#             for region, data in segments.items():
#                 r_start, r_end = region.split("-")
#                 len_bp_to_consider = int(r_end) - int(r_start)
#                 if data["count"] == 1:
#                     for entry in data["entries"]:
#                         if entry["id"] not in visited:
#                             insertions[chromosome][window][entry["rep_class"]][entry["rep_family"]][
#                                 "count"
#                             ] += 1
#                             visited.append(entry["id"])
#                             if entry["match_identity"] != "NA":
#                                 insertions[chromosome][window][entry["rep_class"]][entry["rep_family"]][
#                                     "identity"
#                                 ].append(entry["match_identity"])
#                         insertions[chromosome][window][entry["rep_class"]][entry["rep_family"]][
#                             "length"
#                         ] += len_bp_to_consider
#                 elif data["count"] >= 2:
#                     rep_class_comp, rep_family_comp, match_ids = [], [], []
#                     for entry in data["entries"]:
#                         if entry["rep_class"] not in rep_class_comp:
#                             rep_class_comp.append(entry["rep_class"])
#                         if entry["rep_family"] not in rep_family_comp:
#                             rep_family_comp.append(entry["rep_family"])
#                         if entry["id"] not in visited:
#                             visited.append(entry["id"])
#                         if entry["match_identity"] != "NA":
#                             match_ids.append(entry["match_identity"])
#                     if (len(rep_class_comp) == 1) and (len(rep_family_comp) == 1):
#                         insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]][
#                             "count"
#                         ] += 1
#                         insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]][
#                             "length"
#                         ] += len_bp_to_consider
#                         if entry["match_identity"] != "NA":
#                             insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]][
#                                 "identity"
#                             ].append(np.mean(match_ids))
#                     else:
#                         insertions[chromosome][window]["Ambiguous"]["NA"]["count"] += 1
#                         insertions[chromosome][window]["Ambiguous"]["NA"]["length"] += len_bp_to_consider
#                         if entry["match_identity"] != "NA":
#                             insertions[chromosome][window]["Ambiguous"]["NA"]["identity"].append(
#                                 np.mean(match_ids)
#                             )
#
#     return (repeats, insertions)


def parse_edta_annotations(edta_file, windows_dict):
    """
    Parse EDTA GFF3 annotations into per-window insertion summaries.

    Returns:
        repeats: dict of repeat class -> family counts.
        insertions: per-window insertion metrics.
    """
    print("# Running parse_edta_annotations function\n")

    def _parse_attributes(attr_field):
        attrs = {}
        for chunk in attr_field.split(";"):
            if "=" in chunk:
                key, value = chunk.split("=", 1)
                attrs[key] = value
        return attrs

    repeats = {}
    entries_by_chrom = {}
    counters = {}

    with open(edta_file, "r") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9:
                continue
            chrom = fields[0]
            start = min(int(fields[3]), int(fields[4]))
            end = max(int(fields[3]), int(fields[4]))
            attrs = _parse_attributes(fields[8])
            classification = attrs.get("classification", "").split("/")
            rep_class = classification[0] if classification else ""
            if len(classification) == 1 or classification[1] == "":
                rep_family = "NA"
            else:
                rep_family = classification[1]

            if rep_class not in repeats:
                repeats[rep_class] = {}
            if rep_family not in repeats[rep_class]:
                repeats[rep_class][rep_family] = 0
            repeats[rep_class][rep_family] += 1

            rep_element = attrs.get("Name", "")
            identity = attrs.get("identity", "NA")
            if identity != "NA":
                identity = float(identity)
            counter = counters.get(chrom, 0)
            counters[chrom] = counter + 1
            entries_by_chrom.setdefault(chrom, []).append(
                {
                    "id": str(counter),
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "rep_class": rep_class,
                    "rep_family": rep_family,
                    "rep_element": rep_element,
                    "match_identity": identity,
                }
            )

    insertions = {}
    for chrom in windows_dict:
        insertions[chrom] = {}
        for window in windows_dict[chrom]:
            insertions[chrom][window] = {}
            for rep_class in repeats:
                insertions[chrom][window][rep_class] = {}
                for rep_family in repeats[rep_class]:
                    insertions[chrom][window][rep_class][rep_family] = {
                        "count": 0,
                        "length": 0,
                        "identity": [],
                    }
            insertions[chrom][window]["Ambiguous"] = {
                "NA": {"count": 0, "length": 0, "identity": []}
            }

    for chromosome in insertions:
        chrom_entries = entries_by_chrom.get(chromosome, [])
        for window in insertions[chromosome]:
            coordinates = window.split("-")
            window_insertions = {}
            win_start = int(coordinates[1])
            win_end = int(coordinates[2])
            for entry in chrom_entries:
                if entry["start"] < win_end and entry["end"] > win_start:
                    min_bp_to_consider = max(entry["start"], win_start)
                    max_bp_to_consider = min(entry["end"], win_end)
                    window_insertions[entry["id"]] = {
                        "chrom": entry["chrom"],
                        "start": min_bp_to_consider,
                        "end": max_bp_to_consider,
                        "rep_class": entry["rep_class"],
                        "rep_family": entry["rep_family"],
                        "rep_element": entry["rep_element"],
                        "match_identity": entry["match_identity"],
                    }
            segments = split_overlaps(window_annotations=window_insertions)
            visited = []
            for region, data in segments.items():
                r_start, r_end = region.split("-")
                len_bp_to_consider = int(r_end) - int(r_start)
                if data["count"] == 1:
                    for entry in data["entries"]:
                        if entry["id"] not in visited:
                            insertions[chromosome][window][entry["rep_class"]][entry["rep_family"]][
                                "count"
                            ] += 1
                            visited.append(entry["id"])
                            if entry["match_identity"] != "NA":
                                insertions[chromosome][window][entry["rep_class"]][entry["rep_family"]][
                                    "identity"
                                ].append(entry["match_identity"])
                        insertions[chromosome][window][entry["rep_class"]][entry["rep_family"]][
                            "length"
                        ] += len_bp_to_consider
                elif data["count"] >= 2:
                    rep_class_comp, rep_family_comp, match_ids = [], [], []
                    for entry in data["entries"]:
                        if entry["rep_class"] not in rep_class_comp:
                            rep_class_comp.append(entry["rep_class"])
                        if entry["rep_family"] not in rep_family_comp:
                            rep_family_comp.append(entry["rep_family"])
                        if entry["id"] not in visited:
                            visited.append(entry["id"])
                        if entry["match_identity"] != "NA":
                            match_ids.append(entry["match_identity"])
                    if (len(rep_class_comp) == 1) and (len(rep_family_comp) == 1):
                        insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]][
                            "count"
                        ] += 1
                        insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]][
                            "length"
                        ] += len_bp_to_consider
                        if entry["match_identity"] != "NA":
                            insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]][
                                "identity"
                            ].append(np.mean(match_ids))
                    else:
                        insertions[chromosome][window]["Ambiguous"]["NA"]["count"] += 1
                        insertions[chromosome][window]["Ambiguous"]["NA"]["length"] += len_bp_to_consider
                        if entry["match_identity"] != "NA":
                            insertions[chromosome][window]["Ambiguous"]["NA"]["identity"].append(
                                np.mean(match_ids)
                            )

    return (repeats, insertions)
