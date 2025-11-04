vars = {}

for line in open("/Users/cc54/Documents/Lab/Analyses/TEleVIZion/data/genomes/Ngousso/NRR_threshold_3pairs.tsv").readlines():
    fields = line[:-1].split('\t')
    if fields[0] != "variation":
        vars[fields[0]] = {"chr":fields[2], "coord":int(fields[3])}

with open("/Users/cc54/Documents/Lab/Analyses/TEleVIZion/data/genomes/Ngousso/HMchr_NRR_threshold_3pairs.fasta.out", "w") as out:
    for line in open("/Users/cc54/Documents/Lab/Analyses/TEleVIZion/data/genomes/Ngousso/HM_NRR_threshold_3pairs.fasta.out").readlines():
        # print(line)
        pseudo_fields = line[:-1].split()
        fields = []
        counter = 0
        for f in pseudo_fields:
            if f != "" :
                fields.append(f)
                counter += 1
        if counter > 0:
            if fields[0].isnumeric() :
                var = fields[4]
                # len = int(fields[6]) - int(fields[5]) + 1
                fields[4] = vars[var]["chr"]
                fields[5] = str(vars[var]["coord"])
                # fields[6] = str(vars[var]["coord"] + len)
                fields[6] = str(vars[var]["coord"] + 1)
                # print(fields)
                out.write("\t".join(fields) + "\n")
            else:
                out.write(line)
        else:
            out.write(line)
        