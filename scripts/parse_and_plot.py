import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_hex
import pandas as pd
import subprocess
from intervaltree import IntervalTree

def define_windows(genome_file, window_size, chrom_to_plot):
    print("\n\n# Running define_windows function\n")

    windows, chrom_names = {}, {}
    for line in open(genome_file).readlines():
        if not line.startswith("chr\tstart"):
            fields = line[:-1].split("\t")
            chrom, length, name = fields[0], int(fields[2]), fields[3]
            chrom_names[chrom] = name
            windows[chrom] = []
            for start in range(0, length, window_size):
                end = min(start + window_size, length)
                window = f"{chrom}-{start + 1}-{end}"  # 1-based indexing
                windows[chrom].append(window)
    if chrom_to_plot == "all": 
        tmp = []
        for chrom in windows:
            tmp.append(chrom)
        chroms = ",".join(tmp)
    else:
        chroms = chrom_to_plot
    return(windows, chroms, chrom_names)

def parse_repeatmasker_kimura(repeatmasker_kimura_file):
    print("# Running parse_repeatmasker_kimura function\n")

    kimura_div = {}
    for line in open(repeatmasker_kimura_file).readlines():
        # Get rid of human friendly formatting 
        pseudo_fields = line[:-1].split()
        fields = []
        for f in pseudo_fields : 
            if len(f) > 0 :
                fields.append(f)
        if (len(fields) > 0):
            if (not fields[0].isdigit()) and (fields[0] not in ["Jukes/Cantor", "=================================================================", "File:", "Weighted", "Class", "-----", "Coverage", "Div"]):
                if fields[4] == "----":
                    fields[4] = 0
                value = float(fields[4])
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
                kimura_div[fields[1]] = {"pct_div":value, "category":cat}
    return(kimura_div)

def parse_repeatmasker_output(repeatmasker_file, windows_dict, kimura_dict = None):
    print("# Running parse_repeatmasker_output function\n")
    
    ## Retrieve all types of repeats in repeatmasker output file
    repeats = {}
    for line in open(repeatmasker_file).readlines():
        # Get rid of human friendly formatting 
        pseudo_fields = line[:-1].split()
        fields = []
        for f in pseudo_fields : 
            if len(f) > 0 :
                fields.append(f)
        if len(fields) > 1 : 
            if fields[0][0].isnumeric() : # Get rid of header
                # Parse the line
                classification = fields[10].split("/")
                rep_class = classification[0]
                if len(classification) == 1:
                    rep_family = "NA"
                else:
                    if classification[1] == "":
                        rep_family = "NA"
                    else:
                        rep_family = classification[1]
                # Update repeats dict
                if rep_class not in repeats:
                    repeats[rep_class] = {}
                if rep_family not in repeats[rep_class]:
                    repeats[rep_class][rep_family] = 0
                repeats[rep_class][rep_family] += 1
   
    ## Retrieve all insertions
    # Initialise insertions dict, by chromosomes and affiliated windows
    insertions = {}
    for chrom in windows_dict:
        insertions[chrom] = {}
        for window in windows_dict[chrom]:
            insertions[chrom][window] = {}
            for rep_class in repeats:
                insertions[chrom][window][rep_class] = {}
                for rep_family in repeats[rep_class]:
                    insertions[chrom][window][rep_class][rep_family] = {"count":0, "length":0}
            insertions[chrom][window]["Ambiguous"] = {"NA": {"count":0, "length":0}} 

    # Initialise kimura_div dict, by chromosomes and affiliated windows
    if kimura_dict != None:
        kimura_div = {}
        for chrom in windows_dict:
            kimura_div[chrom] = {}
            for window in windows_dict[chrom]:
                kimura_div[chrom][window] = {}
                for rep_class in repeats:
                    kimura_div[chrom][window][rep_class] = {"0-10":0, 
                                                            "10-20":0, 
                                                            "20-30":0, 
                                                            "30-40":0, 
                                                            "40-70":0} 

    # Fill insertions dict chromosome per chromosome to limit memory usage
    for chromosome in insertions:
        # print(chromosome)
        tmp_insertions, counter = {}, 0
        for line in open(repeatmasker_file).readlines():
            # Get rid of human friendly formatting 
            pseudo_fields = line[:-1].split()
            fields = []
            for f in pseudo_fields : 
                if len(f) > 0 :
                    fields.append(f)
            if len(fields) > 1 : 
                if fields[0][0].isnumeric() : # Get rid of header
                    # Parse the line
                    chrom = fields[4]
                    if chrom == chromosome:
                        start = min(int(fields[5]), int(fields[6]))
                        end = max(int(fields[5]), int(fields[6]))
                        classification = fields[10].split("/")
                        rep_class = classification[0]
                        if len(classification) == 1:
                            rep_family = "NA"                    
                        else:
                            if classification[1] == "":
                                rep_family = "NA"
                            else:
                                rep_family = classification[1]
                        # print(chrom, start, end, rep_class, rep_family)
                        tmp_insertions[str(counter)] = {
                            "chrom":chrom, 
                            "start":int(start), 
                            "end":int(end), 
                            "rep_class":rep_class, 
                            "rep_family":rep_family,
                            "rep_element":fields[9]}
                        counter += 1
        # Update insertions dict
        for window in insertions[chromosome]:
            coordinates = window.split("-")
            window_insertions = {}
            for tmp in tmp_insertions:
                if (tmp_insertions[tmp]["start"] < int(coordinates[2])) and (tmp_insertions[tmp]["end"] > int(coordinates[1])):
                    min_bp_to_consider = max(tmp_insertions[tmp]["start"], int(coordinates[1]))
                    max_bp_to_consider = min(tmp_insertions[tmp]["end"], int(coordinates[2]))
                    window_insertions[tmp] = {
                        "chrom":tmp_insertions[tmp]["chrom"], 
                        "start":min_bp_to_consider, 
                        "end":max_bp_to_consider, 
                        "rep_class":tmp_insertions[tmp]["rep_class"], 
                        "rep_family":tmp_insertions[tmp]["rep_family"],
                        "rep_element":tmp_insertions[tmp]["rep_element"]
                    }

            # Detect overlapping annotations to avoid length/percentages > 100%
            segments = detect_overlapping_annotations(window_insertions_dict = window_insertions)
            # for region, data in segments.items():
            #     print(f"Region {region}: {data['count']} sequence(s)")
            #     for entry in data["entries"]:
            #         print(f"  ID: {entry['id']}, Class: {entry['rep_class']}, Family: {entry['rep_family']}")
            visited = []
            for region, data in segments.items():
                r_start, r_end = region.split("-")
                len_bp_to_consider = int(r_end) - int(r_start)
                if data['count'] == 1:
                    for entry in data["entries"]:
                        if entry["id"] not in visited:
                            insertions[chromosome][window][entry['rep_class']][entry['rep_family']]["count"] += 1 
                            visited.append(entry["id"])
                            if (kimura_dict != None) and (entry['rep_element'] in kimura_dict):
                                kimura_div[chromosome][window][entry['rep_class']][kimura_dict[entry['rep_element']]['category']] += 1
                        insertions[chromosome][window][entry['rep_class']][entry['rep_family']]["length"] += len_bp_to_consider
                elif data['count'] >= 2:
                    rep_class_comp, rep_family_comp = [], []
                    for entry in data["entries"]:
                        if entry["rep_class"] not in rep_class_comp:
                            rep_class_comp.append(entry["rep_class"])
                        if entry["rep_family"] not in rep_family_comp:
                            rep_family_comp.append(entry["rep_family"])
                        if entry["id"] not in visited:
                            visited.append(entry["id"])
                    if (len(rep_class_comp) == 1) and (len(rep_family_comp) == 1):
                        # Same rep_class and rep_family
                        insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]]["count"] += 1 
                        if (kimura_dict != None) and (entry['rep_element'] in kimura_dict):
                            kimura_div[chromosome][window][entry['rep_class']][kimura_dict[entry['rep_element']]['category']] += 1
                        insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]]["length"] += len_bp_to_consider
                    else:
                        # Different rep_class and/or rep_family
                        insertions[chromosome][window]["Ambiguous"]["NA"]["count"] += 1 
                        if (kimura_dict != None) and (entry['rep_element'] in kimura_dict):
                            kimura_div[chromosome][window][entry['rep_class']][kimura_dict[entry['rep_element']]['category']] += 1
                        insertions[chromosome][window]["Ambiguous"]["NA"]["length"] += len_bp_to_consider

    if kimura_dict != None:
        # print(kimura_div)
        return(repeats, insertions, kimura_div)
    else:
        return(repeats, insertions)                    

def parse_edta_output(edta_file, windows_dict):
    print("# Running parse_edta_output function\n")
    
    ## Retrieve all types of repeats in EDTA output file
    repeats = {}
    for line in open(edta_file).readlines():
        if line[0] != "#":
            fields = line[:-1].split()
            classification = fields[8].split("classification=")[1].split(";")[0].split("/")
            rep_class = classification[0]
            if len(classification) == 1:
                rep_family = "NA"                    
            else:
                if classification[1] == "":
                    rep_family = "NA"
                else:
                    rep_family = classification[1]
            # Update repeats dict
            if rep_class not in repeats:
                repeats[rep_class] = {}
            if rep_family not in repeats[rep_class]:
                repeats[rep_class][rep_family] = 0
            repeats[rep_class][rep_family] += 1
   
    ## Retrieve all insertions
    # Initialize insertions dict, following the chromosomes and affiliated windows
    insertions = {}
    for chrom in windows_dict:
        insertions[chrom] = {}
        for window in windows_dict[chrom]:
            insertions[chrom][window] = {}
            for rep_class in repeats:
                insertions[chrom][window][rep_class] = {}
                for rep_family in repeats[rep_class]:
                    insertions[chrom][window][rep_class][rep_family] = {"count":0, "length":0}
            insertions[chrom][window]["Ambiguous"] = {"NA": {"count":0, "length":0}} # Can be refined later (DNA_LTR, DNA_DNA, ...).
    # Fill insertions dict chromosome per chromosome to limit memory usage
    for chromosome in insertions:
        # print(chromosome)
        chrom_insertions, counter = {}, 0
        for line in open(edta_file).readlines():
            if line[0] != "#":
                fields = line[:-1].split("\t")
                chrom = fields[0]
                if chrom == chromosome:
                    # print(fields)
                    start = min(int(fields[3]), int(fields[4]))
                    end = max(int(fields[3]), int(fields[4]))
                    classification = fields[8].split("classification=")[1].split(";")[0].split("/")
                    rep_class = classification[0]
                    if len(classification) == 1:
                        rep_family = "NA"                    
                    else:
                        if classification[1] == "":
                            rep_family = "NA"
                        else:
                            rep_family = classification[1]
                    chrom_insertions[str(counter)] = {
                        "chrom":chrom, 
                        "start":int(start), 
                        "end":int(end), 
                        "rep_class":rep_class, 
                        "rep_family":rep_family
                        }
                    counter += 1
        for window in insertions[chromosome]:
            coordinates = window.split("-")
            window_insertions = {}
            for tmp in chrom_insertions:
                if (chrom_insertions[tmp]["start"] < int(coordinates[2])) and (chrom_insertions[tmp]["end"] > int(coordinates[1])):
                    min_bp_to_consider = max(chrom_insertions[tmp]["start"], int(coordinates[1]))
                    max_bp_to_consider = min(chrom_insertions[tmp]["end"], int(coordinates[2]))
                    window_insertions[tmp] = {
                        "chrom":chrom_insertions[tmp]["chrom"], 
                        "start":min_bp_to_consider, 
                        "end":max_bp_to_consider, 
                        "rep_class":chrom_insertions[tmp]["rep_class"], 
                        "rep_family":chrom_insertions[tmp]["rep_family"]
                    }
            # Detect nested variations to avoid length/percentages > 100%
            segments = detect_overlapping_annotations(window_insertions_dict = window_insertions)
            # for region, data in segments.items():
            #     print(f"Region {region}: {data['count']} sequence(s)")
            #     for entry in data["entries"]:
            #         print(f"  ID: {entry['id']}, Class: {entry['rep_class']}, Family: {entry['rep_family']}")
            visited = []
            for region, data in segments.items():
                r_start, r_end = region.split("-")
                len_bp_to_consider = int(r_end) - int(r_start)
                if data['count'] == 1:
                    for entry in data["entries"]:
                        if entry["id"] not in visited:
                            insertions[chromosome][window][entry['rep_class']][entry['rep_family']]["count"] += 1 
                            visited.append(entry["id"])
                        insertions[chromosome][window][entry['rep_class']][entry['rep_family']]["length"] += len_bp_to_consider
                elif data['count'] >= 2:
                    rep_class_comp, rep_family_comp = [], []
                    for entry in data["entries"]:
                        if entry["rep_class"] not in rep_class_comp:
                            rep_class_comp.append(entry["rep_class"])
                        if entry["rep_family"] not in rep_family_comp:
                            rep_family_comp.append(entry["rep_family"])
                        if entry["id"] not in visited:
                            visited.append(entry["id"])
                    if (len(rep_class_comp) == 1) and (len(rep_family_comp) == 1):
                        # Same rep_class and rep_family
                        insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]]["count"] += 1 
                        insertions[chromosome][window][rep_class_comp[0]][rep_family_comp[0]]["length"] += len_bp_to_consider
                    else:
                        # Different rep_class and/or rep_family
                        insertions[chromosome][window]["Ambiguous"]["NA"]["count"] += 1 
                        insertions[chromosome][window]["Ambiguous"]["NA"]["length"] += len_bp_to_consider

    return(repeats, insertions)

def detect_overlapping_annotations(window_insertions_dict):
    """
    Takes a dictionary of interval entries and returns a dictionary of segments.
    
    Keys: "start-end" string
    Values: {
        'count': number of overlapping intervals,
        'entries': list of {id, rep_class, rep_family}
    }
    """
    points = set()
    for item in window_insertions_dict.values():
        points.add(item["start"])
        points.add(item["end"])
    sorted_points = sorted(points)

    tree = IntervalTree()
    for key, item in window_insertions_dict.items():
        tree.addi(item["start"], item["end"], key)

    segments = {}

    for i in range(len(sorted_points) - 1):
        seg_start = sorted_points[i]
        seg_end = sorted_points[i + 1]
        if seg_start == seg_end:
            continue

        overlaps = tree[seg_start:seg_end]
        key = f"{seg_start}-{seg_end}"
        segment_data = {
            "count": len(overlaps),
            "entries": []
        }

        for ov in overlaps:
            info = window_insertions_dict[ov.data]
            segment_data["entries"].append({
                "id": ov.data,
                "rep_class": info["rep_class"],
                "rep_family": info["rep_family"],
                "rep_element":info["rep_element"]
            })

        segments[key] = segment_data

    return(segments)

def manage_colors(insertions, classes_order, palette = None):
    """
    """
    print("# Running manage_colors function\n")

    # Global aggregation for color mapping
    records = []
    for chrom_windows in insertions.values():
        for classes in chrom_windows.values():
            for rep_class, families in classes.items():
                for rep_family, metrics in families.items():
                    records.append({
                        'rep_class': rep_class,
                        'rep_family': rep_family,
                        'count': metrics['count'],
                        'length': metrics['length']
                    })
    df_global = pd.DataFrame(records)
    agg_global = df_global.groupby(['rep_class','rep_family'])[['length','count']].sum().reset_index()

    # Determine class order
    if classes_order == None:
        class_order = agg_global.groupby('rep_class')['length'].sum().sort_values(ascending=False).index.tolist()
    else:
        class_order = classes_order
    n = max(len(class_order)-1,1)

    fam_colors, class_colors, class_colors_hex = {}, {}, {}

    if palette == None:
        cmap = plt.get_cmap('turbo')
        class_colors = {cls: cmap(i/n) for i,cls in enumerate(class_order)}
        # # Re‐cast to pure Python floats:
        # class_colors_py = {
        #     cls: tuple(float(x) for x in rgba)
        #     for cls, rgba in class_colors.items()
        # }
        # Convert to hex strings (including the alpha channel):
        class_colors_hex = {
            cls: to_hex(rgba, keep_alpha=True)
            for cls, rgba in class_colors.items()
        }

    else:
        # Make sure to know all classes present in your TE annotation file
        # scripts/palette.tsv should end by "\n"
        homemade_palette = {}
        for line in open(palette).readlines():
            desc, r, g, b, hex, cat = line[:-1].split("\t")
            for category in cat.split(","):
                homemade_palette[category] = {"r":int(r), "g":int(g), "b":int(b), "hex":hex}
        for cls in class_order:
            class_colors[cls] = [homemade_palette[cls]["r"]/255, homemade_palette[cls]["g"]/255, homemade_palette[cls]["b"]/255, 1]
            class_colors_hex[cls] = homemade_palette[cls]["hex"]

    # Family mapping
    for cls in class_order:
        fams = agg_global[agg_global['rep_class']==cls].sort_values('length',ascending=False)['rep_family'].tolist()
        alphas = np.linspace(1.0,0.3,len(fams))
        base = class_colors[cls]
        fam_colors[cls] = {fam:(base[0],base[1],base[2],alpha) for fam,alpha in zip(fams,alphas)}

    return(agg_global, class_order, class_colors_hex, fam_colors)

def plot_family_insertions(name, insertions, chrom_names, agg_global, class_order, class_colors_hex, fam_colors, per_chromosome=False, width=10, height=8):
    """
    """
    print("# Running plot_family_insertions function\n")

    os.makedirs('analyses', exist_ok=True)

    # Use Tahoma for all fonts
    # plt.rcParams['font.family'] = 'Tahoma'

    # # 1) Global aggregation for color mapping
    # records = []
    # for chrom_windows in insertions.values():
    #     for classes in chrom_windows.values():
    #         for rep_class, families in classes.items():
    #             for rep_family, metrics in families.items():
    #                 records.append({
    #                     'rep_class': rep_class,
    #                     'rep_family': rep_family,
    #                     'count': metrics['count'],
    #                     'length': metrics['length']
    #                 })
    # df_global = pd.DataFrame(records)
    # agg_global = df_global.groupby(['rep_class','rep_family'])[['count','length']].sum().reset_index()

    # # Determine class order
    # if classes_order == None:
    #     class_order = agg_global.groupby('rep_class')['count'].sum().sort_values(ascending=False).index.tolist()
    # else:
    #     class_order = classes_order

    # # Base colors per class
    # cmap = plt.get_cmap('turbo')
    # n = max(len(class_order)-1,1)
    # class_colors = {cls: cmap(i/n) for i,cls in enumerate(class_order)}
    # # # Re‐cast to pure Python floats:
    # # class_colors_py = {
    # #     cls: tuple(float(x) for x in rgba)
    # #     for cls, rgba in class_colors.items()
    # # }
    # # Convert to hex strings (including the alpha channel):
    # class_colors_hex = {
    #     cls: to_hex(rgba, keep_alpha=True)
    #     for cls, rgba in class_colors.items()
    # }

    # # Family alpha mapping
    # fam_colors = {}
    # for cls in class_order:
    #     fams = agg_global[agg_global['rep_class']==cls].sort_values('count',ascending=False)['rep_family'].tolist()
    #     alphas = np.linspace(1.0,0.3,len(fams))
    #     base = class_colors[cls]
    #     fam_colors[cls] = {fam:(base[0],base[1],base[2],alpha) for fam,alpha in zip(fams,alphas)}

    def _plot_stacked(agg_df, metric, ylabel, title, fname, width=10, height=8):
        fig, ax = plt.subplots(figsize=(width,height))
        x = np.arange(len(class_order))
        xt = [cls.replace('_','\n') for cls in class_order]
        for i,cls in enumerate(class_order):
            bottom=0
            total=0
            for fam,rgba in fam_colors[cls].items():
                mask=(agg_df['rep_class']==cls)&(agg_df['rep_family']==fam)
                if not mask.any(): continue
                val=int(agg_df.loc[mask,metric].values[0])
                ax.bar(i,val,bottom=bottom,width=0.4,color=rgba,edgecolor="white",linewidth=0.3)
                bottom+=val; total+=val
            ax.text(i,total,f"{total:_}".replace("_", ","),ha='center',va='bottom') 
        # legend
        handles=[]
        for cls in class_order:
            for fam in reversed(list(fam_colors[cls].keys())):
                rgba=fam_colors[cls][fam]
                mask=(agg_df['rep_class']==cls)&(agg_df['rep_family']==fam)
                val=int(agg_df.loc[mask,metric].values[0]) if mask.any() else 0
                val_to_print = f"{val:_}".replace("_", ",")
                handles.append(mpatches.Patch(color=rgba,label=f"{cls} - {fam} ({val_to_print})")) 
        ax.set_xticks(x); ax.set_xticklabels(xt,rotation=0)
        ax.set_xlabel('Repeat class',weight='bold',labelpad=12)
        ax.set_ylabel(ylabel,weight='bold',labelpad=12)
        ax.set_title(title,weight='bold',pad=20)
        ax.legend(handles=handles,title='Type',loc='upper left',bbox_to_anchor=(1.05,1),ncol=1,fontsize='small')
        plt.tight_layout()
        fig.savefig(f'{fname}',bbox_inches='tight')
        plt.close(fig)

    def _plot_contiguous(agg_df, metric, ylabel, title, fname, width=width, height=height):
        fig, ax = plt.subplots(figsize=(width,height))
        families=[]; values=[]; colors=[]
        for cls in class_order:
            for fam,rgba in fam_colors[cls].items():
                mask=(agg_df['rep_class']==cls)&(agg_df['rep_family']==fam)
                if not mask.any(): continue
                val=int(agg_df.loc[mask,metric].values[0])
                families.append(f"{cls} - {fam}")
                values.append(val)
                colors.append(rgba)
        y = np.arange(len(families))
        ax.barh(y, values, color=colors, height=0.6)
        # Annotate each bar with spacing and non-bold
        for yi, val in zip(y, values):
            ax.text(val + max(values)*0.01, yi, f"{val:_}".replace("_", ","), va='center', ha='left')
        ax.set_yticks(y); ax.set_yticklabels(families) # ,fontsize=8
        ax.invert_yaxis()
        # Extend x-axis limit to include labels
        ax.set_xlim(0, max(values)*1.2)
        ax.set_xlabel(ylabel,weight='bold',labelpad=12)
        ax.set_ylabel('Repeat type',weight='bold',labelpad=12)
        ax.set_title(title,weight='bold',pad=20)
        # legend of classes
        legend_handles = [mpatches.Patch(color=class_colors_hex[cls], label=cls)
                          for cls in class_order]
        ax.legend(handles=legend_handles, title='Repeat class',
                loc = "best", ncol=1, fontsize='small') # loc='upper left', bbox_to_anchor=(1.05,1)
        plt.tight_layout()
        fig.savefig(f'{fname}',bbox_inches='tight')
        plt.close(fig)

    # Generate plots if not already existing 
    _plot_stacked(agg_global,'count', 'Insertion count', f'Stacked Counts - Whole Genome - {name}', f'analyses/{name}_whole_genome_stacked_counts_by_class.pdf', width=width, height=height)
    _plot_stacked(agg_global,'length', 'Base pair span', f'Stacked Lengths - Whole Genome - {name}', f'analyses/{name}_whole_genome_stacked_lengths_by_class.pdf', width=width, height=height)
    _plot_contiguous(agg_global,'count', 'Insertion count', f'Counts by Type - Whole Genome - {name}', f'analyses/{name}_whole_genome_contiguous_counts_by_class.pdf', width=width, height=height)
    _plot_contiguous(agg_global,'length', 'Base pair span', f'Lengths by Type - Whole Genome - {name}', f'analyses/{name}_whole_genome_contiguous_lengths_by_class.pdf', width=width, height=height)

    if per_chromosome:
        for chrom, windows in insertions.items():
            rec=[] 
            for classes in windows.values():
                for rep_class,families in classes.items():
                    for rep_family,metrics in families.items():
                        rec.append({'rep_class':rep_class,'rep_family':rep_family,'count':metrics['count'],'length':metrics['length']})
            df_chr=pd.DataFrame(rec)
            if df_chr.empty: continue
            agg_chr=df_chr.groupby(['rep_class','rep_family'])[['count','length']].sum().reset_index()
            _plot_stacked(agg_chr,'count','Insertion count', f'Stacked Counts in chromosome {chrom_names[chrom]} - {name}', f'analyses/{name}_{chrom_names[chrom]}_stacked_counts_by_class.pdf', width=width, height=height)
            _plot_stacked(agg_chr,'length','Base pair span', f'Stacked Lengths in chromosome {chrom_names[chrom]} - {name}', f'analyses/{name}_{chrom_names[chrom]}_stacked_lengths_by_class.pdf', width=width, height=height)
            _plot_contiguous(agg_chr,'count','Insertion count', f'Counts by Type in chromosome {chrom_names[chrom]} - {name}', f'analyses/{name}_{chrom_names[chrom]}_contiguous_counts_by_class.pdf', width=width, height=height)
            _plot_contiguous(agg_chr,'length','Base pair span', f'Lengths by Type in chromosome {chrom_names[chrom]} - {name}', f'analyses/{name}_{chrom_names[chrom]}_contiguous_lengths_by_class.pdf', width=width, height=height)
    
    # return(class_order, class_colors_hex, fam_colors)

def export_window_class_bed(name, window_size, insertions, windows, class_colors, class_order):
    """
    """
    print("# Running export_window_class_bed function\n")

    os.makedirs('analyses/karyoplot_tables', exist_ok = True)

    output_bed = f'analyses/karyoplot_tables/{name}_{window_size}_repeat_classes.bed'

    # TO DO 
    # if (os.path.isfile(output_bed)) and (os.stat(output_bed).st_size != 0):
    #     print(f'{output_bed} already exists and isn’t empty. To force recalculation, please delete it first.')

    # Build rows
    rows = []
    for chrom, win_list in windows.items():
        for win_label in win_list:
            parts = win_label.split('-')
            if len(parts) != 3:
                continue
            _, start_s, end_s = parts
            start = int(start_s) - 1
            end = int(end_s)
            window_len = end - start
            window_karyoplotr = int(start + (window_len/2))
            row = {'chrom': chrom, 'start': start, 'end': end, 'barycenter': window_karyoplotr}
            for cls in class_order:
                fams = insertions.get(chrom, {}).get(win_label, {}).get(cls, {})
                cnt = sum(m['count'] for m in fams.values())
                lng = sum(m['length'] for m in fams.values())
                pct = (lng / window_len) if window_len > 0 else 0
                row[f'{cls}_count'] = cnt
                row[f'{cls}_length'] = lng
                row[f'{cls}_pct'] = round(pct, 3)
            rows.append(row)
    
    # Create DataFrame
    df = pd.DataFrame(rows)

    # Compute stacked counts
    count_cols = [f'{cls}_count' for cls in class_order]
    stacked = df[count_cols].cumsum(axis=1)
    stacked.columns = [f'{cls}_count_stacked' for cls in class_order]
    df = pd.concat([df, stacked], axis=1)

    # Compute stacked percentages
    pct_cols = [f'{cls}_pct' for cls in class_order]
    stacked = df[pct_cols].cumsum(axis=1)
    stacked.columns = [f'{cls}_pct_stacked' for cls in class_order]
    df = pd.concat([df, stacked], axis=1)

    # Build final column order: chrom, start, end, then per class:
    cols = ['chrom', 'start', 'end', 'barycenter']
    for cls in class_order:
        cols.append(f'{cls}_count')
        cols.append(f'{cls}_count_stacked')
        cols.append(f'{cls}_length')
        cols.append(f'{cls}_pct')
        cols.append(f'{cls}_pct_stacked')

    df = df[cols]
    df.to_csv(output_bed, sep='\t', index=False)

    reversed_classes = ",".join(list(reversed(class_order)))
    tmp_colors = []
    for cls in class_order:
        tmp_colors.append(class_colors[cls])
    reversed_colors = ",".join(list(reversed(tmp_colors)))
    return(reversed_classes, reversed_colors)

def export_window_kimura_bed(name, window_size, kimura_div, windows, class_order):
    """
    """
    print("# Running export_window_kimura_bed function\n")

    # os.makedirs('analyses/karyoplot_tables', exist_ok = True)

    output_bed = f'analyses/karyoplot_tables/{name}_{window_size}_kimura.bed'

    # TO DO 
    # if (os.path.isfile(output_bed)) and (os.stat(output_bed).st_size != 0):
    #     print(f'{output_bed} already exists and isn’t empty. To force recalculation, please delete it first.')

    kimura_cats = ["0-10", "10-20", "20-30", "30-40", "40-70"]

    # Build rows
    rows = []
    for chrom, win_list in windows.items():
        for win_label in win_list:
            parts = win_label.split('-')
            if len(parts) != 3:
                continue
            _, start_s, end_s = parts
            start = int(start_s) - 1
            end = int(end_s)
            window_len = end - start
            window_karyoplotr = int(start + (window_len/2))
            row = {'chrom': chrom, 'start': start, 'end': end, 'barycenter': window_karyoplotr}
            alls = []
            for category in kimura_cats:
                alls.append(f"ALL_{category}_pct")
                row[f"ALL_{category}_count"] = 0 
                row[f"ALL_{category}_pct"] = 0 

            ### Per class
            # not_to_plot = []
            for cls in class_order:
                cats = kimura_div.get(chrom, {}).get(win_label, {}).get(cls, {})
                if len(cats) > 0:
                    cnt_sum = sum(m for m in cats.values())
                    for category in kimura_cats:
                        row[f'{cls}_{category}_count'] = cats[category]
                        row[f"ALL_{category}_count"] += cats[category] # Prepare all classes
                        if cnt_sum != 0:
                            row[f'{cls}_{category}_pct'] = round(cats[category] / cnt_sum, 4)
                        else: 
                            row[f'{cls}_{category}_pct'] = 0
                else:
                    # not_to_plot.append(cls)
                    for category in kimura_cats:
                        row[f'{cls}_{category}_count'] = 0
                        row[f'{cls}_{category}_pct'] = 0

            ### All classes
            cnt_sum = 0
            for category in kimura_cats:
                cnt_sum += row[f"ALL_{category}_count"]
            for category in kimura_cats:
                if cnt_sum != 0:
                    row[f"ALL_{category}_pct"] = round(row[f"ALL_{category}_count"] / cnt_sum, 4)
                else: 
                    row[f'ALL_{category}_pct'] = 0
                
            rows.append(row)

    df = pd.DataFrame(rows)
    
    # Compute stacked percentages for ALL
    stacked = df[alls].cumsum(axis=1)
    stacked.columns = [f'ALL_{category}_pct_stacked' for category in kimura_cats]
    df = pd.concat([df, stacked], axis=1)

    # Compute stacked percentages per class
    # for cls in class_order:
    #     if cls not in not_to_plot:
    #         to_stack = []
    #         for category in kimura_cats:
    #             to_stack.append(f"{cls}_{category}_pct")
    #         stacked = df[to_stack].cumsum(axis=1)
    #         stacked.columns = [f'{cls}_{category}_pct_stacked' for category in kimura_cats]
    #         df = pd.concat([df, stacked], axis=1)
    for cls in class_order:
        to_stack = []
        for category in kimura_cats:
            to_stack.append(f"{cls}_{category}_pct")
        stacked = df[to_stack].cumsum(axis=1)
        stacked.columns = [f'{cls}_{category}_pct_stacked' for category in kimura_cats]
        df = pd.concat([df, stacked], axis=1)

    # Build final column order: chrom, start, end, then per class:
    cols = ['chrom', 'start', 'end', 'barycenter']
    for category in kimura_cats:
        cols.append(f'ALL_{category}_count')
        cols.append(f'ALL_{category}_pct')
        cols.append(f'ALL_{category}_pct_stacked')
        df[f'ALL_{category}_pct_stacked'] = df[f'ALL_{category}_pct_stacked'].clip(upper=1)
    # for cls in class_order:
    #     if cls not in not_to_plot:
    #         for category in kimura_cats:
    #             cols.append(f'{cls}_{category}_count')
    #             cols.append(f'{cls}_{category}_pct')
    #             cols.append(f'{cls}_{category}_pct_stacked')
    #             df[f'{cls}_{category}_pct_stacked'] = df[f'{cls}_{category}_pct_stacked'].clip(upper=1)
    for cls in class_order:
        for category in kimura_cats:
            cols.append(f'{cls}_{category}_count')
            cols.append(f'{cls}_{category}_pct')
            cols.append(f'{cls}_{category}_pct_stacked')
            df[f'{cls}_{category}_pct_stacked'] = df[f'{cls}_{category}_pct_stacked'].clip(upper=1)

    df = df[cols]
    df.to_csv(output_bed, sep='\t', index=False)

    return(output_bed)

def plot_karyoplots(name, window_size, genome_file, chromosomes, accessibility, classes, colors, plot_per_class, kimura_bed = None):
    """
    """
    print("# Running plot_karyoplots function\n")

    # Percentage and counts - Horizontal version
    # cmd = f'Rscript scripts/plot_chromosomes.R --name {name} --genome {genome_file} --chromosome {chromosomes} --accessibility {accessibility} --input analyses/karyoplot_tables/{name}_{window_size}_repeat_classes.bed --classesorder {classes} --perclass {plot_per_class} --colorsorder {colors} --output analyses/{name}_{window_size} --kimura {kimura_bed}'
    # print(cmd)
    # subprocess.run(cmd.split(" "), check=True)

    cmd = f'Rscript scripts/plot_chromosomes_kimura.R --name {name} --genome {genome_file} --chromosome {chromosomes} --accessibility {accessibility} --input analyses/karyoplot_tables/{name}_{window_size}_repeat_classes.bed --classesorder {classes} --perclass {plot_per_class} --colorsorder {colors} --output analyses/{name}_{window_size} --kimura {kimura_bed}'
    print(cmd)
    subprocess.run(cmd.split(" "), check=True)

    # Percentage and counts - Vertical version
    # cmd = "Rscript scripts/plot_chromosomes_vertical.R --genome " + genome + " --chromosome " + chromosomes_to_plot + " --accessibility " + accessibility + " --input analyses/karyoplotr_tables/rep_classes.bed --classesorder " + reversed_classes + " --colorsorder " + reversed_colors + " --output analyses/"


if __name__ == "__main__":
    
    ### PARSING ARGUMENTS 

    parser = argparse.ArgumentParser(
        description="Funny description to come"
    )

    parser.add_argument("--name", required=False, type=str, default="output")
    parser.add_argument("--genome", required=True, type=str)
    parser.add_argument("--repeatmasker", required=False, type=str, default=None) 
    parser.add_argument("--kimura", required=False, type=str, default=None) 
    parser.add_argument("--edta", required=False, type=str, default=None) 
    parser.add_argument("--windowsize", required=False, type=int, default=10000)
    parser.add_argument("--chromtoplot", required=False, type=str, default="all")
    parser.add_argument("--perchromosome", required=False, type=str, default=None)
    parser.add_argument("--classesorder", required=False, type=str, default=None)
    parser.add_argument("--perclass", required=False, type=str, default=None) 
    parser.add_argument("--accessibility", required=False, type=str, default=None)
    parser.add_argument("--figsize", required=False, type=str, default=None)
    parser.add_argument("--palette", required=False, type=str, default=None)

    args = parser.parse_args()

    print(f"\nTEleVIZion is on!\n")

    name = args.name
    print(f"Name: {name}")

    genome = args.genome
    print(f"Genome file: {genome}")

    if args.repeatmasker != None:
        file = args.repeatmasker 
        file_type = "RepeatMasker"
    elif args.edta != None:
        file = args.edta 
        file_type = "EDTA"
    else: 
        print("You should provide an input file using --repeatmasker or --edta")
        sys.exit()
    print(f"Input file ({file_type}): {file}")

    if args.kimura == None:
        kimura_file = None
    else:
        kimura_file = args.kimura 
    print(f"Additional Kimura features: {kimura_file}")

    window_size = args.windowsize
    print(f"Window size: {window_size}")

    if args.chromtoplot != "all":
        chrom_to_plot = args.chromtoplot
    else:
        chrom_to_plot = "all"
    print(f"Chromosomes to plot: {chrom_to_plot}")

    if args.classesorder != None:
        classes_order = args.classesorder.split(",")
    else:
        classes_order = None
    print(f"Classes order: {classes_order}")

    if args.accessibility != None:
        accessibility = args.accessibility 
    else:
        accessibility = "not_displayed"
    print(f"Accessibility: {accessibility}")

    if (args.perchromosome == None) or (args.perchromosome == "False"):
        per_chrom = False
    else:
        per_chrom = True
    print(f"Additional plots per chromosome: {per_chrom}")

    if (args.perclass == None) or (args.perclass == "False"):
        per_class = False
    else:
        per_class = True
    print(f"Additional plots per class: {per_class}")

    if args.figsize != None:
        tmp = args.figsize.split(",")
        width, height = int(tmp[0]), int(tmp[1])
    else:
        width, height = 10, 8
    print(f"Figure size for histograms: {width}, {height}")
    
    if args.palette == None:
        palette_type = "False, switch to default"
    else:
        palette_type = args.palette
    palette = args.palette
    print(f"Homemade palette: {palette_type}")


    ### CODE

    # Define windows 

    windows, chromosomes_to_plot, chrom_names = define_windows(
        genome_file = genome, 
        window_size = window_size,
        chrom_to_plot = chrom_to_plot
    )

    # Detect repeat types and actual insertions

    if args.repeatmasker != None:
        if kimura_file != None:
            kimura_dict = parse_repeatmasker_kimura(
                repeatmasker_kimura_file = kimura_file
            )
            repeats, insertions, kimura_div = parse_repeatmasker_output(
                repeatmasker_file = file,
                windows_dict = windows, 
                kimura_dict = kimura_dict
            )
        else:
            repeats, insertions = parse_repeatmasker_output(
                repeatmasker_file = file,
                windows_dict = windows
            )
        
    elif args.edta != None:
        repeats, insertions = parse_edta_output(
            edta_file = file,
            windows_dict = windows
        )

    # Plots (whole genome or list of chromosomes)

    agg_global, class_order, class_colors_hex, fam_colors = manage_colors(
        insertions = insertions, 
        classes_order = classes_order,
        palette = palette
        )
    
    plot_family_insertions(
        name = name,
        insertions = insertions, 
        per_chromosome = per_chrom,
        chrom_names = chrom_names,
        agg_global = agg_global,
        class_order = class_order,
        class_colors_hex = class_colors_hex,
        fam_colors = fam_colors,
        width = width,
        height = height
        )
           
    reversed_classes, reversed_colors = export_window_class_bed(
        name = name,
        window_size = window_size,
        windows = windows,
        insertions = insertions, 
        class_colors = class_colors_hex, 
        class_order = class_order
        ) 
    
    if kimura_file != None:
        kimura_bed = export_window_kimura_bed(
            name = name, 
            window_size = window_size, 
            kimura_div = kimura_div, 
            windows = windows, 
            class_order = class_order
            )
    else:
        kimura_bed = None
    
    plot_karyoplots(
        name = name, 
        window_size = window_size,
        genome_file = genome, 
        chromosomes = chromosomes_to_plot, 
        accessibility = accessibility, 
        classes = reversed_classes, 
        colors = reversed_colors,
        plot_per_class = per_class,
        kimura_bed = kimura_bed
    )


# Usage: TEleVIZion % python3 scripts/parse_and_plot.py --edta data/EDTA/GCA_943734845.1.sanitised.fa.mod.EDTA.TEanno.gff3 --genome data/EDTA/GCA_943734845.1_metadata.txt --windowsize 500000 --chromtoplot 2RL,3RL,X --classesorder unknown,Overlapping_annotations,LINE,MITE,DNA,LTR --name GCA_943734845.1
# python3 scripts/parse_and_plot.py --genome data/genomes/Ngousso/chromosomes_chr.tsv --accessibility data/genomes/Ngousso/accessibility.tsv --repeatmasker data/RepeatMasker/ngousso_chr.fasta.out --windowsize 500000 --chromtoplot Chr_2R,Chr_2L,Chr_3R,Chr_3L,Chr_X --name Ngousso --perchromosome False --classesorder Undetermined,Overlapping_annotations,Simple_repeat,Low_complexity,DNA,LINE,LTR,SINE