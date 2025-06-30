import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.colors import to_hex
import pandas as pd

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

def plot_total_repeat_content(repeats_dict, insertions_dict):
    # Get all unique rep_family for stacking and consistent order
    all_rep_family = set()
    for rep_family in repeats_dict.values():
        all_rep_family.update(rep_family.keys())
    all_rep_family = sorted(all_rep_family)

    # Initialize insertions_dict for stacking
    rep_class_labels = list(repeats_dict.keys())
    count_matrix = np.zeros((len(all_rep_family), len(rep_class_labels)))
    length_matrix = np.zeros((len(all_rep_family), len(rep_class_labels)))

    # Fill in the matrices
    for j, rep_class in enumerate(rep_class_labels):
        for i, rep_family in enumerate(all_rep_family):
            if rep_family in repeats_dict[rep_class]:
                for chrom in insertions_dict:
                    for window in insertions_dict[chrom]:
                        count_matrix[i, j] += insertions_dict[chrom][window][rep_class][rep_family]['count']
                        length_matrix[i, j] += insertions_dict[chrom][window][rep_class][rep_family]['length']

    # Plotting
    fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 10), sharex=True,  gridspec_kw={'wspace': 0.1}, layout="constrained") # , gridspec_kw={'wspace': 0.4, 'hspace': 0.4}
    fig.suptitle('Repeat content - whole genome', fontsize=16, weight="bold")

    # Stack plot for counts
    bottom_counts = np.zeros(len(rep_class_labels))
    for i, rep_family in enumerate(all_rep_family):
        axes[0].bar(rep_class_labels, count_matrix[i], bottom=bottom_counts, label=rep_family, width=0.4)
        bottom_counts += count_matrix[i]
    axes[0].set_title('Stacked counts per class by family')
    axes[0].set_ylabel('Count')
    # leg = axes[0].legend() # bbox_to_anchor=(1.05, -1), loc='bottom left'
    # axes[0].tick_params(axis='x', rotation=45)

    # Stack plot for lengths
    bottom_lengths = np.zeros(len(rep_class_labels))
    for i, rep_family in enumerate(all_rep_family):
        axes[1].bar(rep_class_labels, length_matrix[i], bottom=bottom_lengths, label=rep_family, width=0.4)
        bottom_lengths += length_matrix[i]
    axes[1].set_title('Stacked lengths per class by family')
    axes[1].set_ylabel('Length')
    leg = axes[1].legend(bbox_to_anchor=(1.09, 0), loc='lower left', fontsize="small")
    axes[1].tick_params(axis='x', rotation=90)

    # leg.set_in_layout(True)
    
    # plt.tight_layout()
    plt.savefig("analyses/test_total.png") #  bbox_extra_artists=[leg]

    # # Re-defining the updated plotting function after kernel reset
    # import matplotlib.pyplot as plt
    # import numpy as np
    # from matplotlib.colors import to_rgba

    
    # # Get all unique rep_family for stacking and consistent order
    # all_rep_family = set()
    # for rep_family in repeats_dict.values():
    #     all_rep_family.update(rep_family.keys())
    # all_rep_family = sorted(all_rep_family)

    # # Initialize insertions_dict for stacking
    # rep_class_labels = list(repeats_dict.keys())
    # count_matrix = np.zeros((len(all_rep_family), len(rep_class_labels)))
    # length_matrix = np.zeros((len(all_rep_family), len(rep_class_labels)))

    # # Fill in the matrices
    # for j, rep_class in enumerate(rep_class_labels):
    #     for i, rep_family in enumerate(all_rep_family):
    #         if rep_family in repeats_dict[rep_class]:
    #             for chrom in insertions_dict:
    #                 for window in insertions_dict[chrom]:
    #                     count_matrix[i, j] += insertions_dict[chrom][window][rep_class][rep_family]['count']
    #                     length_matrix[i, j] += insertions_dict[chrom][window][rep_class][rep_family]['length']

    # # Generate colors for families, 'NA' is grey
    # base_cmap = plt.get_cmap('tab20', len(all_rep_family))
    # family_colors = {}
    # for i, fam in enumerate(all_rep_family):
    #     if fam == 'NA':
    #         family_colors[fam] = (0.6, 0.6, 0.6, 1.0)  # grey RGBA
    #     else:
    #         family_colors[fam] = base_cmap(i)

    # def adjust_brightness(color, factor):
    #     r, g, b, a = to_rgba(color)
    #     return (min(1, r * factor), min(1, g * factor), min(1, b * factor), a)

    # # Assign colors with brightness adjustment for each class
    # bar_colors_counts = np.zeros((len(all_rep_family), len(rep_class_labels), 4))
    # bar_colors_lengths = np.zeros((len(all_rep_family), len(rep_class_labels), 4))
    # for j, rep_class in enumerate(rep_class_labels):
    #     brightness = 0.6 + 0.4 * (j / len(rep_class_labels))
    #     for i, fam in enumerate(all_rep_family):
    #         base_color = family_colors[fam]
    #         if fam == 'NA':
    #             color = base_color
    #         else:
    #             color = adjust_brightness(base_color, brightness)
    #         bar_colors_counts[i, j] = color
    #         bar_colors_lengths[i, j] = color

    # # Plotting
    # fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(6, 10), sharex=True,
    #                         gridspec_kw={'wspace': 0.1}, layout="constrained")
    # fig.suptitle('Repeat content - whole genome', fontsize=16, weight="bold")

    # # Stack plot for counts
    # bottom_counts = np.zeros(len(rep_class_labels))
    # for i, rep_family in enumerate(all_rep_family):
    #     axes[0].bar(rep_class_labels, count_matrix[i], bottom=bottom_counts,
    #                 color=bar_colors_counts[i], label=rep_family, width=0.4)
    #     bottom_counts += count_matrix[i]
    # axes[0].set_title('Stacked counts per class by family')
    # axes[0].set_ylabel('Count')

    # # Stack plot for lengths
    # bottom_lengths = np.zeros(len(rep_class_labels))
    # for i, rep_family in enumerate(all_rep_family):
    #     axes[1].bar(rep_class_labels, length_matrix[i], bottom=bottom_lengths,
    #                 color=bar_colors_lengths[i], label=rep_family, width=0.4)
    #     bottom_lengths += length_matrix[i]
    # axes[1].set_title('Stacked lengths per class by family')
    # axes[1].set_ylabel('Length')
    # axes[1].tick_params(axis='x', rotation=90)

    # # Shared legend
    # handles, labels = axes[1].get_legend_handles_labels()
    # fig.legend(handles, labels, bbox_to_anchor=(1.09, 0.5), loc='center left', fontsize="small", title="Families")

    # plt.savefig("analyses/test_total.png")
    # plt.close()

def plot_family_insertions_by_class(insertions):
    """
    Generate two color‐coded bar charts for transposable element insertions,
    grouping repeat families contiguously by repeat class:
      1. Total insertion counts per repeat family
      2. Total bases impacted per repeat family
    Bars are colored by repeat class with a legend.
    """
    # Flatten nested dict into records
    records = []
    for chrom_windows in insertions.values():
        for classes in chrom_windows.values():
            for rep_class, families in classes.items():
                for rep_family, metrics in families.items():
                    records.append({
                        'rep_class':  rep_class,
                        'rep_family': rep_family,
                        'count':      metrics['count'],
                        'length':     metrics['length']
                    })

    # Create DataFrame and aggregate per class-family
    df = pd.DataFrame(records)
    agg = (
        df.groupby(['rep_class', 'rep_family'])[['count', 'length']]
          .sum()
          .reset_index()
    )

    # Determine class order by total count desc
    class_order = (
        agg.groupby('rep_class')['count']
           .sum()
           .sort_values(ascending=False)
           .index
           .tolist()
    )

    # Build overall family order and corresponding values
    order = []
    counts = []
    lengths = []
    bar_classes = []
    for cls in class_order:
        fams = (
            agg[agg['rep_class'] == cls]
            .sort_values('count', ascending=False)
            ['rep_family']
            .tolist()
        )
        for fam in fams:
            order.append(fam)
            bar_classes.append(cls)
            row = agg[(agg['rep_class'] == cls) & (agg['rep_family'] == fam)]
            counts.append(int(row['count']))
            lengths.append(int(row['length']))

    # Assign colors to classes
    cmap = plt.get_cmap('tab10')
    class_colors = {cls: cmap(i % 10) for i, cls in enumerate(class_order)}

    # Create legend handles
    legend_handles = [mpatches.Patch(color=class_colors[cls], label=cls) 
                      for cls in class_order]

    # Plot total counts
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(range(len(order)), counts, color=[class_colors[c] for c in bar_classes], edgecolor='black')
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=8)
    ax.set_xlabel('Repeat family (grouped by class)')
    ax.set_ylabel('Insertion count')
    ax.set_title('Total Insertions per Repeat Family by Class')
    ax.legend(handles=legend_handles, title='Repeat class', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig("analyses/test_counts.png")

    # Plot total lengths
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.bar(range(len(order)), lengths, color=[class_colors[c] for c in bar_classes], edgecolor='black')
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels(order, rotation=90, fontsize=8)
    ax.set_xlabel('Repeat family (grouped by class)')
    ax.set_ylabel('Total bases impacted')
    ax.set_title('Total Base Impact per Repeat Family by Class')
    ax.legend(handles=legend_handles, title='Repeat class', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()
    plt.savefig("analyses/test_length.png")

plt.rcParams['font.family'] = 'Tahoma'

def plot_family_insertions_stacked_by_class(insertions):
    
    # Flatten nested dict into records
    records = []
    for chrom_windows in insertions.values():
        for classes in chrom_windows.values():
            for rep_class, families in classes.items():
                for rep_family, metrics in families.items():
                    records.append({
                        'rep_class':  rep_class,
                        'rep_family': rep_family,
                        'count':      metrics['count'],
                        'length':     metrics['length']
                    })

    # Aggregate
    df = pd.DataFrame(records)
    agg = df.groupby(['rep_class', 'rep_family'])[['count', 'length']].sum().reset_index()

    # Determine class order
    class_order = (
        agg.groupby('rep_class')['count']
           .sum()
           .sort_values(ascending=False)
           .index
           .tolist()
    )

    # Base colors from rainbow
    cmap = plt.get_cmap('turbo')
    class_colors = {
        cls: cmap(i / max(len(class_order)-1, 1))
        for i, cls in enumerate(class_order)
    }

    # Replace underscores in class labels with newline for x-axis
    x_labels = [cls.replace('_', '\n') for cls in class_order]
    x = np.arange(len(class_order))

    def _plot(metric, fname, ylabel, title):
        fig, ax = plt.subplots(figsize=(9, 8))
        # Plot stacks
        for i, cls in enumerate(class_order):
            fams = agg[agg['rep_class'] == cls]['rep_family'].tolist()
            alphas = np.linspace(1.0, 0.3, len(fams))
            bottom = 0
            total = 0
            for fam, alpha in zip(fams, alphas):
                val = agg.loc[
                    (agg['rep_class'] == cls) & (agg['rep_family'] == fam),
                    metric
                ].values[0]
                ax.bar(i, val, bottom=bottom, width=0.4,
                       color=class_colors[cls], alpha=alpha)
                bottom += val
                total += val
            ax.text(i, total, str(total), ha='center', va='bottom')

        # Build legend with family values
        legend_handles = []
        for cls in class_order:
            fams = agg[agg['rep_class'] == cls]['rep_family'].tolist()
            alphas = np.linspace(1.0, 0.3, len(fams))
            for fam, alpha in zip(reversed(fams), reversed(alphas)):
                val = agg.loc[
                    (agg['rep_class'] == cls) & (agg['rep_family'] == fam),
                    metric
                ].values[0]
                label = f"{fam} ({val})"
                legend_handles.append(
                    mpatches.Patch(color=class_colors[cls], alpha=alpha, label=label)
                )

        ax.set_xticks(x)
        ax.set_xticklabels(x_labels, rotation=0)
        ax.set_xlabel('Repeat class', weight='bold', labelpad=12)
        ax.set_ylabel(ylabel, weight='bold', labelpad=12)
        ax.set_title(title, weight='bold', pad=20)
        ax.legend(handles=legend_handles, title='Families',
                  loc='upper left', bbox_to_anchor=(1.05, 1),
                  ncol=1, fontsize='small')
        plt.tight_layout()
        fig.savefig(f'analyses/{fname}', bbox_inches='tight')

    # Generate plots
    _plot('count', 'stacked_counts_by_class.png',
          'Insertion count', 'Stacked Insertion Counts per Class by Family')
    _plot('length', 'stacked_lengths_by_class.png',
          'Total bases impacted', 'Stacked Base Impact per Class by Family')

def plot_family_insertions_stacked_by_class_per_chromosome(insertions):

    # Common colormap
    cmap = plt.get_cmap('turbo')

    for chrom, windows in insertions.items():
        # Flatten for this chromosome
        records = []
        for classes in windows.values():
            for rep_class, families in classes.items():
                for rep_family, metrics in families.items():
                    records.append({
                        'rep_class':  rep_class,
                        'rep_family': rep_family,
                        'count':      metrics['count'],
                        'length':     metrics['length']
                    })
        df = pd.DataFrame(records)
        if df.empty:
            continue

        # Aggregate and determine order
        agg = df.groupby(['rep_class', 'rep_family'])[['count', 'length']].sum().reset_index()
        class_order = (
            agg.groupby('rep_class')['count']
               .sum()
               .sort_values(ascending=False)
               .index
               .tolist()
        )
        # Assign class colors
        class_colors = {
            cls: cmap(i / max(len(class_order)-1, 1))
            for i, cls in enumerate(class_order)
        }
        # Prepare x-axis labels
        x_labels = [cls.replace('_', '\n') for cls in class_order]
        x = np.arange(len(class_order))

        def _plot(metric, fname, ylabel, title):
            fig, ax = plt.subplots(figsize=(9, 8))
            # Plot stacks
            for i, cls in enumerate(class_order):
                fams = agg[agg['rep_class'] == cls]['rep_family'].tolist()
                alphas = np.linspace(1.0, 0.3, len(fams))
                bottom = 0
                total = 0
                for fam, alpha in zip(fams, alphas):
                    val = agg.loc[
                        (agg['rep_class'] == cls) & (agg['rep_family'] == fam),
                        metric
                    ].values[0]
                    ax.bar(i, val, bottom=bottom, width=0.4,
                           color=class_colors[cls], alpha=alpha)
                    bottom += val
                    total += val
                ax.text(i, total, str(total), ha='center', va='bottom')

            # Build legend
            legend_handles = []
            for cls in class_order:
                fams = agg[agg['rep_class'] == cls]['rep_family'].tolist()
                alphas = np.linspace(1.0, 0.3, len(fams))
                for fam, alpha in zip(reversed(fams), reversed(alphas)):
                    val = agg.loc[
                        (agg['rep_class'] == cls) & (agg['rep_family'] == fam),
                        metric
                    ].values[0]
                    label = f"{fam} ({val})"
                    legend_handles.append(
                        mpatches.Patch(color=class_colors[cls], alpha=alpha, label=label)
                    )

            ax.set_xticks(x)
            ax.set_xticklabels(x_labels, rotation=0)
            ax.set_xlabel('Repeat class', weight='bold', labelpad=12)
            ax.set_ylabel(ylabel, weight='bold', labelpad=12)
            ax.set_title(title, weight='bold', pad=20)
            ax.legend(handles=legend_handles, title='Families',
                      loc='upper left', bbox_to_anchor=(1.05, 1),
                      ncol=1, fontsize='small')
            plt.tight_layout()
            fig.savefig(f'analyses/{fname}', bbox_inches='tight')
            plt.close(fig)

        # Create per-chromosome plots
        _plot('count',
              f'{chrom}_stacked_counts.png',
              'Insertion count',
              f'Stacked Insertion Counts per Class by Family in chromosome {chrom}')
        _plot('length',
              f'{chrom}_stacked_lengths.png',
              'Total bases impacted',
              f'Stacked Base Impact per Class by Family in chromosome {chrom}')






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

    windows = define_windows(
        genome_file = genome, 
        window_size = window_size
    )

    repeats, insertions = parse_repeatmasker_output(
        repeatmasker_file = file,
        windows_dict = windows
    )


    # print(repeats)
    # plot_total_repeat_content(repeats, insertions)
    # plot_family_insertions_by_class(insertions)
    # plot_family_insertions_stacked_by_class(insertions)
    plot_family_insertions_stacked_by_class_per_chromosome(insertions)