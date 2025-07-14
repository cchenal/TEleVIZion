import os
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
                ax.bar(i,val,bottom=bottom,width=0.4,color=rgba)
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
            _plot_stacked(agg_chr,'count','Insertion count',f'Stacked Counts in chr {chrom}',f'{chrom}_stacked_counts_by_class.png')
            _plot_stacked(agg_chr,'length','Total bases impacted',f'Stacked Lengths in chr {chrom}',f'{chrom}_stacked_lengths_by_class.png')
            _plot_contiguous(agg_chr,'count','Insertion count',f'Counts by Family in chr {chrom}',f'{chrom}_contiguous_counts_by_class.png')
            _plot_contiguous(agg_chr,'length','Total bases impacted',f'Lengths by Family in chr {chrom}',f'{chrom}_contiguous_lengths_by_class.png')
    else:
        _plot_stacked(agg_global,'count','Insertion count','Stacked Counts','stacked_counts_by_class.png')
        _plot_stacked(agg_global,'length','Total bases impacted','Stacked Lengths','stacked_lengths_by_class.png')
        _plot_contiguous(agg_global,'count','Insertion count','Counts by Family','contiguous_counts_by_class.png')
        _plot_contiguous(agg_global,'length','Total bases impacted','Lengths by Family','contiguous_lengths_by_class.png')

    return fam_colors

# def export_insertions_for_karyoploter_fam(insertions, fam_colors, output_bed):
#     """
#     Flatten per-window, per-class insertion lengths and export to a BED-like file
#     enriched with colors for karyoploteR.
    
#     Assumes `insertions[chrom][window_label][rep_class][rep_family] = {'count', 'length'}` 
#     where window_label is "chr-start-end".
    
#     Args:
#       insertions: nested dict of insertions by chromosome and window_label.
#       fam_colors: dict mapping rep_class -> rep_family -> RGBA tuple.
#       output_bed: path to write the tab-delimited file
      
#     Generates a file with columns:
#       chrom  start   end   rep_class   rep_family   length    color
#     """
#     os.makedirs(os.path.dirname(output_bed) or '.', exist_ok=True)
#     rows = []
#     for chrom, win_dict in insertions.items():
#         for win_label, class_dict in win_dict.items():
#             # win_label expected "chr-start-end"
#             parts = win_label.split('-')
#             if len(parts) != 3:
#                 continue
#             _, start_str, end_str = parts
#             start = int(start_str) - 1  # convert to 0-based
#             end = int(end_str)
#             for rep_class, families in class_dict.items():
#                 for rep_family, metrics in families.items():
#                     length = metrics.get('length', 0)
#                     rgba = fam_colors.get(rep_class, {}).get(rep_family, (0,0,0,1))
#                     color = to_hex(rgba, keep_alpha=True)
#                     rows.append({
#                         'chrom':      chrom,
#                         'start':      start,
#                         'end':        end,
#                         'rep_class':  rep_class,
#                         'rep_family': rep_family,
#                         'length':     length,
#                         'color':      color
#                     })
#     df = pd.DataFrame(rows)
#     df.to_csv(output_bed, sep='\t', index=False)

def export_class_insertions_for_karyoploter(insertions, fam_colors, output_bed):
    """
    Flatten per-window, per-class total insertion lengths and export to a BED-like file
    enriched with class colors (using the first family's color).
    
    Args:
      insertions: nested dict of insertions by chromosome and window_label.
                  insertions[chrom][window_label][rep_class][rep_family] = {'count', 'length'}
      fam_colors: dict mapping rep_class -> rep_family -> RGBA tuple
      output_bed: path to write the tab-delimited file
      
    Generates a file with columns:
      chrom  start   end   rep_class   length    color
    where length is total across all families and color is the hex RGBA of the first family.
    """
    # Ensure output directory exists
    os.makedirs(os.path.dirname(output_bed) or '.', exist_ok=True)
    
    rows = []
    for chrom, win_dict in insertions.items():
        for win_label, class_dict in win_dict.items():
            parts = win_label.split('-')
            if len(parts) != 3:
                continue
            _, start_str, end_str = parts
            start = int(start_str) - 1  # convert to 0-based for BED
            end = int(end_str)
            
            for rep_class, families in class_dict.items():
                # Sum lengths across families
                total_length = sum(metrics.get('length', 0) for metrics in families.values())
                # Choose color of the first family in fam_colors[rep_class]
                fam_list = list(fam_colors.get(rep_class, {}).keys())
                if fam_list:
                    first_family = fam_list[0]
                    rgba = fam_colors[rep_class][first_family]
                else:
                    rgba = (0, 0, 0, 1)
                color = to_hex(rgba, keep_alpha=True)
                
                rows.append({
                    'chrom':     chrom,
                    'start':     start,
                    'end':       end,
                    'rep_class': rep_class,
                    'length':    total_length,
                    'color':     color
                })
    
    df = pd.DataFrame(rows)
    df.to_csv(output_bed, sep='\t', index=False)





if __name__ == "__main__":
    
    ### PARSING ARGUMENTS 

    parser = argparse.ArgumentParser(
        description="Funny description to come"
    )
    parser.add_argument("--genome", required=True, type=str)
    parser.add_argument("--repeatmasker", required=True, type=str) # To be changed when several tools available
    parser.add_argument("--windowsize", required=False, type=int, default=10000)
    args = parser.parse_args()

    genome = args.genome
    file = args.repeatmasker # to be modified when several tools will be available
    window_size = args.windowsize

    ### CODE

    # Define windows 

    windows = define_windows(
        genome_file = genome, 
        window_size = window_size
    )

    # Detect repeats and repeat types

    repeats, insertions = parse_repeatmasker_output(
        repeatmasker_file = file,
        windows_dict = windows
    )

    # Plots (whole genome or list of chromosomes)

    colors = plot_family_insertions(insertions, per_chromosome=False)
    # export_class_insertions_for_karyoploter(insertions, colors, "analyses/insertions_by_window_class.bed")
    # cmd = "Rscript scripts/plot_chromosomes.R --genome data/genomes/Ngousso/chromosomes_chr.tsv --chromosome Chr_2R,Chr_2L,Chr_3R,Chr_3L,Chr_X --accessibility data/genomes/Ngousso/accessibility.tsv --out analyses/karyoplotr.png".split(" ")
    # subprocess.run(cmd, check=True)
