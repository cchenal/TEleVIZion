import argparse

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
    file = args.repeatmasker
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

    print(insertions["Chr_X"]['Chr_X-1-500000'])