import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_hex
import pandas as pd
import subprocess

def define_windows(genome_file, window_size):
    windows = {}
    for line in open(genome_file).readlines():
        if not line.startswith("chr\tstart"):
            fields = line[:-1].split("\t")
            chrom, length = fields[0], int(fields[2])
            windows[chrom] = []
            for start in range(0, length, window_size):
                end = min(start + window_size, length)
                window = f"{chrom}-{start + 1}-{end}"  # 1-based indexing
                windows[chrom].append(window)
    return(windows)

def parse_repeatmasker_output(repeatmasker_file, windows_dict):
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
            for tmp in tmp_insertions:
                if (tmp_insertions[tmp]["start"] < int(coordinates[2])) and (tmp_insertions[tmp]["end"] > int(coordinates[1])):
                    min_bp_to_consider = max(tmp_insertions[tmp]["start"], int(coordinates[1]))
                    max_bp_to_consider = min(tmp_insertions[tmp]["end"], int(coordinates[2]))
                    len_bp_to_consider = max_bp_to_consider - min_bp_to_consider
                    insertions[chromosome][window][tmp_insertions[tmp]["rep_class"]][tmp_insertions[tmp]["rep_family"]]["count"] += 1 
                    insertions[chromosome][window][tmp_insertions[tmp]["rep_class"]][tmp_insertions[tmp]["rep_family"]]["length"] += len_bp_to_consider
    return(repeats, insertions)                    

def parse_edta_output(edta_file, windows_dict):
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
    # Fill insertions dict chromosome per chromosome to limit memory usage
    for chromosome in insertions:
        # print(chromosome)
        tmp_insertions, counter = {}, 0
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
                    tmp_insertions[str(counter)] = {
                        "chrom":chrom, 
                        "start":int(start), 
                        "end":int(end), 
                        "rep_class":rep_class, 
                        "rep_family":rep_family}
                    # print(tmp_insertions[str(counter)])
                    counter += 1
        # Update insertions dict
        for window in insertions[chromosome]:
            coordinates = window.split("-")
            for tmp in tmp_insertions:
                if (tmp_insertions[tmp]["start"] < int(coordinates[2])) and (tmp_insertions[tmp]["end"] > int(coordinates[1])):
                    min_bp_to_consider = max(tmp_insertions[tmp]["start"], int(coordinates[1]))
                    max_bp_to_consider = min(tmp_insertions[tmp]["end"], int(coordinates[2]))
                    len_bp_to_consider = max_bp_to_consider - min_bp_to_consider
                    insertions[chromosome][window][tmp_insertions[tmp]["rep_class"]][tmp_insertions[tmp]["rep_family"]]["count"] += 1 
                    insertions[chromosome][window][tmp_insertions[tmp]["rep_class"]][tmp_insertions[tmp]["rep_family"]]["length"] += len_bp_to_consider
    return(repeats, insertions)

def plot_family_insertions(insertions, per_chromosome=False):
    """
    Creates both stacked and contiguous non-stacked bar charts of TE insertions by class/family.
    - per_chromosome=False: global genome-wide charts
    - per_chromosome=True: per-chromosome charts
    Saves PNGs to analyses/ and returns fam_colors mapping.
    """
    # os.makedirs('analyses', exist_ok=True)

    # Use Tahoma for all fonts
    plt.rcParams['font.family'] = 'Tahoma'

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
    class_order = agg_global.groupby('rep_class')['count'].sum().sort_values(ascending=False).index.tolist()

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

    def _plot_stacked(agg_df, metric, ylabel, title, fname):
        fig, ax = plt.subplots(figsize=(9,8))
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
            ax.text(i,total,str(total),ha='center',va='bottom')
        # legend
        handles=[]
        for cls in class_order:
            for fam in reversed(list(fam_colors[cls].keys())):
                rgba=fam_colors[cls][fam]
                mask=(agg_df['rep_class']==cls)&(agg_df['rep_family']==fam)
                val=int(agg_df.loc[mask,metric].values[0]) if mask.any() else 0
                handles.append(mpatches.Patch(color=rgba,label=f"{fam} ({val})"))
        ax.set_xticks(x); ax.set_xticklabels(xt,rotation=0)
        ax.set_xlabel('Repeat class',weight='bold',labelpad=12)
        ax.set_ylabel(ylabel,weight='bold',labelpad=12)
        ax.set_title(title,weight='bold',pad=20)
        ax.legend(handles=handles,title='Families',loc='upper left',bbox_to_anchor=(1.05,1),ncol=1,fontsize='small')
        plt.tight_layout()
        fig.savefig(f'analyses/{fname}',bbox_inches='tight')
        plt.close(fig)

    def _plot_contiguous(agg_df, metric, ylabel, title, fname):
        fig, ax = plt.subplots(figsize=(9,8))
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
            ax.text(val + max(values)*0.01, yi, str(val), va='center', ha='left')
        ax.set_yticks(y); ax.set_yticklabels(families,fontsize=8)
        ax.invert_yaxis()
        # Extend x-axis limit to include labels
        ax.set_xlim(0, max(values)*1.1)
        ax.set_xlabel(ylabel,weight='bold',labelpad=12)
        ax.set_ylabel('Repeat family',weight='bold',labelpad=12)
        ax.set_title(title,weight='bold',pad=20)
        plt.tight_layout()
        fig.savefig(f'analyses/{fname}',bbox_inches='tight')
        plt.close(fig)

    # Generate plots
    _plot_stacked(agg_global,'count','Insertion count','Stacked Counts - Whole Genome','stacked_counts_by_class.png')
    _plot_stacked(agg_global,'length','Base pair span','Stacked Lengths - Whole Genome','stacked_lengths_by_class.png')
    _plot_contiguous(agg_global,'count','Insertion count','Counts by Family - Whole Genome','contiguous_counts_by_class.png')
    _plot_contiguous(agg_global,'length','Base pair span','Lengths by Family - Whole Genome','contiguous_lengths_by_class.png')

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
            _plot_stacked(agg_chr,'count','Insertion count',f'Stacked Counts in chromosome {chrom}',f'{chrom}_stacked_counts_by_class.png')
            _plot_stacked(agg_chr,'length','Base pair span',f'Stacked Lengths in chromosome {chrom}',f'{chrom}_stacked_lengths_by_class.png')
            _plot_contiguous(agg_chr,'count','Insertion count',f'Counts by Family in chromosome {chrom}',f'{chrom}_contiguous_counts_by_class.png')
            _plot_contiguous(agg_chr,'length','Base pair span',f'Lengths by Family in chromosome {chrom}',f'{chrom}_contiguous_lengths_by_class.png')
    
    return(class_colors_hex, fam_colors)

def export_window_class_bed(insertions, windows, output_bed, class_colors, class_order=None):
    """
    Export a window-based summary BED file with per-class counts and lengths,
    optionally percentage of window span, and cumulative stacked counts.

    Args:
      insertions: dict mapping insertions[chrom][window_label][rep_class][rep_family] = {'count','length'}
      windows: dict mapping chrom -> list of window_label "chr-start-end"
      output_bed: path to write the summary
      class_order: optional list defining the output order of classes
                   (default = all detected, sorted)
    """
    # os.makedirs(os.path.dirname(output_bed) or '.', exist_ok=True)

    # Discover all classes present
    all_classes = sorted({
        rep_class
        for chrom_data in insertions.values()
        for window_dict in chrom_data.values()
        for rep_class in window_dict.keys()
    })

    # Determine final class list and order
    if class_order:
        # keep only classes that actually exist
        classes = [c for c in class_order if c in all_classes]
    else:
        # default: include all, sorted
        classes = all_classes

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
            for cls in classes:
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
    count_cols = [f'{cls}_count' for cls in classes]
    stacked = df[count_cols].cumsum(axis=1)
    stacked.columns = [f'{cls}_count_stacked' for cls in classes]
    df = pd.concat([df, stacked], axis=1)

    # Compute stacked percentages
    pct_cols = [f'{cls}_pct' for cls in classes]
    stacked = df[pct_cols].cumsum(axis=1)
    stacked.columns = [f'{cls}_pct_stacked' for cls in classes]
    df = pd.concat([df, stacked], axis=1)

    # Build final column order: chrom, start, end, then per class:
    cols = ['chrom', 'start', 'end', 'barycenter']
    for cls in classes:
        cols.append(f'{cls}_count')
        cols.append(f'{cls}_count_stacked')
        cols.append(f'{cls}_length')
        cols.append(f'{cls}_pct')
        cols.append(f'{cls}_pct_stacked')

    df = df[cols]
    df.to_csv(output_bed, sep='\t', index=False)

    reversed_classes = ",".join(list(reversed(classes)))
    tmp_colors = []
    for cls in classes:
        tmp_colors.append(class_colors[cls])
    reversed_colors = ",".join(list(reversed(tmp_colors)))
    return(reversed_classes, reversed_colors)



if __name__ == "__main__":
    
    ### PARSING ARGUMENTS 

    parser = argparse.ArgumentParser(
        description="Funny description to come"
    )

    parser.add_argument("--genome", required=True, type=str)
    parser.add_argument("--repeatmasker", required=False, type=str, default=None) 
    parser.add_argument("--edta", required=False, type=str, default=None) 
    parser.add_argument("--windowsize", required=False, type=int, default=10000)
    parser.add_argument("--classesorder", required=False, type=str, default=None)

    args = parser.parse_args()

    genome = args.genome
    if args.repeatmasker != None:
        file = args.repeatmasker 
    elif args.edta != None:
        file = args.edta 
    else: 
        print("You shhould provide an inut file using --repeatmasker or --edta")
        sys.exit()
    window_size = args.windowsize
    if args.classesorder != None:
        classes_order = args.classesorder.split(",")
    else:
        classes_order = None

    ### CODE

    # Define windows 

    windows = define_windows(
        genome_file = genome, 
        window_size = window_size
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

    class_colors, fam_colors = plot_family_insertions(insertions, per_chromosome=False)
    # reversed_classes, reversed_colors = export_window_class_bed(insertions, windows, "analyses/karyoplotr_tables/rep_classes.bed", class_colors, class_order=classes_order) 
    # cmd = "Rscript scripts/plot_chromosomes.R --genome data/genomes/Ngousso/chromosomes_chr.tsv --chromosome Chr_2R,Chr_2L,Chr_3R,Chr_3L,Chr_X --accessibility data/genomes/Ngousso/accessibility.tsv --input analyses/karyoplotr_tables/rep_classes.bed --classesorder " + reversed_classes + " --colorsorder " + reversed_colors + " --output analyses/"
    # subprocess.run(cmd.split(" "), check=True)

# Usage: python3 scripts/parse_and_plot.py --repeatmasker data/RepeatMasker/ngousso_chr.fasta.out --genome data/genomes/Ngousso/chromosomes_chr.tsv --windowsize 500000 --classesorder Simple_repeat,Low_complexity,DNA,LINE,LTR,SINE,Undetermined