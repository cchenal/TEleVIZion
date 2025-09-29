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

def parse_repeatmasker_output(repeatmasker_file, windows_dict):
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
            insertions[chrom][window]["Overlapping_annotations"] = {"NA": {"count":0, "length":0}} 
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
                            "rep_family":rep_family}
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
                        "rep_family":tmp_insertions[tmp]["rep_family"]
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
                        insertions[chromosome][window]["Overlapping_annotations"]["NA"]["count"] += 1 
                        insertions[chromosome][window]["Overlapping_annotations"]["NA"]["length"] += len_bp_to_consider

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
            insertions[chrom][window]["Overlapping_annotations"] = {"NA": {"count":0, "length":0}} # Can be refined later (DNA_LTR, DNA_DNA, ...).
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
                        insertions[chromosome][window]["Overlapping_annotations"]["NA"]["count"] += 1 
                        insertions[chromosome][window]["Overlapping_annotations"]["NA"]["length"] += len_bp_to_consider

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
                "rep_family": info["rep_family"]
            })

        segments[key] = segment_data

    return(segments)

def plot_family_insertions(name, insertions, chrom_names_dict, per_chromosome=False, width=10, height=8):
    """
    """
    print("# Running plot_family_insertions function\n")

    os.makedirs('analyses', exist_ok=True)

    # Use Tahoma for all fonts
    # plt.rcParams['font.family'] = 'Tahoma'

    # 1) Global aggregation for color mapping
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
    agg_global = df_global.groupby(['rep_class','rep_family'])[['count','length']].sum().reset_index()

    # Determine class order
    if classes_order == None:
        class_order = agg_global.groupby('rep_class')['count'].sum().sort_values(ascending=False).index.tolist()
    else:
        class_order = classes_order

    # Base colors per class
    cmap = plt.get_cmap('turbo')
    n = max(len(class_order)-1,1)
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

    # Family alpha mapping
    fam_colors = {}
    for cls in class_order:
        fams = agg_global[agg_global['rep_class']==cls].sort_values('count',ascending=False)['rep_family'].tolist()
        alphas = np.linspace(1.0,0.3,len(fams))
        base = class_colors[cls]
        fam_colors[cls] = {fam:(base[0],base[1],base[2],alpha) for fam,alpha in zip(fams,alphas)}

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
            ax.text(i,total,f"{total:_}".replace("_", " "),ha='center',va='bottom') 
        # legend
        handles=[]
        for cls in class_order:
            for fam in reversed(list(fam_colors[cls].keys())):
                rgba=fam_colors[cls][fam]
                mask=(agg_df['rep_class']==cls)&(agg_df['rep_family']==fam)
                val=int(agg_df.loc[mask,metric].values[0]) if mask.any() else 0
                val_to_print = f"{val:_}".replace("_", " ")
                handles.append(mpatches.Patch(color=rgba,label=f"{fam} ({val_to_print})")) 
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
                families.append(fam)
                values.append(val)
                colors.append(rgba)
        y = np.arange(len(families))
        ax.barh(y, values, color=colors, height=0.6)
        # Annotate each bar with spacing and non-bold
        for yi, val in zip(y, values):
            ax.text(val + max(values)*0.01, yi, f"{val:_}".replace("_", " "), va='center', ha='left')
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
                loc = "lower right", ncol=1, fontsize='small') # loc='upper left', bbox_to_anchor=(1.05,1)
        plt.tight_layout()
        fig.savefig(f'{fname}',bbox_inches='tight')
        plt.close(fig)

    # Generate plots if not already existing 
    _plot_stacked(agg_global,'count', 'Insertion count', f'Stacked Counts - Whole Genome - {name}', f'analyses/{name}_whole_genome_stacked_counts_by_class.png', width=width, height=height)
    _plot_stacked(agg_global,'length', 'Base pair span', f'Stacked Lengths - Whole Genome - {name}', f'analyses/{name}_whole_genome_stacked_lengths_by_class.png', width=width, height=height)
    _plot_contiguous(agg_global,'count', 'Insertion count', f'Counts by Type - Whole Genome - {name}', f'analyses/{name}_whole_genome_contiguous_counts_by_class.png', width=width, height=height)
    _plot_contiguous(agg_global,'length', 'Base pair span', f'Lengths by Type - Whole Genome - {name}', f'analyses/{name}_whole_genome_contiguous_lengths_by_class.png', width=width, height=height)

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
            _plot_stacked(agg_chr,'count','Insertion count', f'Stacked Counts in chromosome {chrom_names[chrom]} - {name}', f'analyses/{name}_{chrom_names[chrom]}_stacked_counts_by_class.png', width=width, height=height)
            _plot_stacked(agg_chr,'length','Base pair span', f'Stacked Lengths in chromosome {chrom_names[chrom]} - {name}', f'analyses/{name}_{chrom_names[chrom]}_stacked_lengths_by_class.png', width=width, height=height)
            _plot_contiguous(agg_chr,'count','Insertion count', f'Counts by Type in chromosome {chrom_names[chrom]} - {name}', f'analyses/{name}_{chrom_names[chrom]}_contiguous_counts_by_class.png', width=width, height=height)
            _plot_contiguous(agg_chr,'length','Base pair span', f'Lengths by Type in chromosome {chrom_names[chrom]} - {name}', f'analyses/{name}_{chrom_names[chrom]}_contiguous_lengths_by_class.png', width=width, height=height)
    
    return(class_order, class_colors_hex, fam_colors)

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

def plot_karyoplots(name, window_size, genome_file, chromosomes, accessibility, classes, colors):
    """
    """
    print("# Running plot_karyoplots function\n")

    # Percentage and counts - Horizontal version
    cmd = f'Rscript scripts/plot_chromosomes.R --name {name} --genome {genome_file} --chromosome {chromosomes} --accessibility {accessibility} --input analyses/karyoplot_tables/{name}_{window_size}_repeat_classes.bed --classesorder {classes} --colorsorder {colors} --output analyses/{name}_{window_size}'
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
    parser.add_argument("--edta", required=False, type=str, default=None) 
    parser.add_argument("--windowsize", required=False, type=int, default=10000)
    parser.add_argument("--chromtoplot", required=False, type=str, default="all")
    parser.add_argument("--classesorder", required=False, type=str, default=None)
    parser.add_argument("--accessibility", required=False, type=str, default=None)
    parser.add_argument("--perchromosome", required=False, type=str, default=None)
    parser.add_argument("--figsize", required=False, type=str, default=None)

    args = parser.parse_args()

    name = args.name

    genome = args.genome

    if args.repeatmasker != None:
        file = args.repeatmasker 
    elif args.edta != None:
        file = args.edta 
    else: 
        print("You should provide an inut file using --repeatmasker or --edta")
        sys.exit()

    window_size = args.windowsize

    if args.chromtoplot != "all":
        chrom_to_plot = args.chromtoplot
    else:
        chrom_to_plot = "all"

    if args.classesorder != None:
        classes_order = args.classesorder.split(",")
    else:
        classes_order = None

    if args.accessibility != None:
        accessibility = args.accessibility 
    else:
        accessibility = "not_displayed"

    if (args.perchromosome == None) or (args.perchromosome == "False"):
        per_chrom = False
    else:
        per_chrom = True

    if args.figsize != None:
        tmp = args.figsize.split(",")
        width, height = int(tmp[0]), int(tmp[1])
    else:
        width, height = 10, 8


    ### CODE

    # Define windows 

    windows, chromosomes_to_plot, chrom_names = define_windows(
        genome_file = genome, 
        window_size = window_size,
        chrom_to_plot = chrom_to_plot
    )

    # Detect repeat types and actual insertions

    if args.repeatmasker != None:
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

    class_order, class_colors, fam_colors = plot_family_insertions(
        name = name,
        insertions = insertions, 
        per_chromosome = per_chrom,
        chrom_names_dict = chrom_names,
        width = width,
        height = height
        )
    
    reversed_classes, reversed_colors = export_window_class_bed(
        name = name,
        window_size = window_size,
        windows = windows,
        insertions = insertions, 
        class_colors = class_colors, 
        class_order = class_order
        ) 
    
    plot_karyoplots(
        name = name, 
        window_size = window_size,
        genome_file = genome, 
        chromosomes = chromosomes_to_plot, 
        accessibility = accessibility, 
        classes = reversed_classes, 
        colors = reversed_colors
    )


# Usage: TEleVIZion % python3 scripts/parse_and_plot.py --edta data/EDTA/GCA_943734845.1.sanitised.fa.mod.EDTA.TEanno.gff3 --genome data/EDTA/GCA_943734845.1_metadata.txt --windowsize 500000 --chromtoplot 2RL,3RL,X --classesorder unknown,Overlapping_annotations,LINE,MITE,DNA,LTR --name GCA_943734845.1
# python3 scripts/parse_and_plot.py --genome data/genomes/Ngousso/chromosomes_chr.tsv --accessibility data/genomes/Ngousso/accessibility.tsv --repeatmasker data/RepeatMasker/ngousso_chr.fasta.out --windowsize 500000 --chromtoplot Chr_2R,Chr_2L,Chr_3R,Chr_3L,Chr_X --name Ngousso --perchromosome False --classesorder Undetermined,Overlapping_annotations,Simple_repeat,Low_complexity,DNA,LINE,LTR,SINE