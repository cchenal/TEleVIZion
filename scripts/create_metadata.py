import os
import subprocess

# chrom_of_interest = ["2R", "2L", "3R", "3L", "X"]

# for f in os.listdir("data/EDTA/"):
#     if f.endswith("_scaff_lengths.tsv"):
#         name = f.split("_scaff")[0]
#         if os.path.isfile("data/EDTA/" + name + ".sanitised.fa.mod.EDTA.TEanno.gff3"):
#             print(name)
#             tmp = {}
#             for line in open("data/EDTA/" + f).readlines():
#                 if line[0] != "#":
#                     chrom, length = line[:-1].split("\t")
#                     if chrom in chrom_of_interest:
#                         tmp[chrom] = length
#             if len(tmp) == 5 :
#                 with open("data/EDTA/" + name + "_partial_metadata.txt", "w") as out:
#                     out.write("chr\tstart\tend\tname\tgieStain\n")
#                     for chr in tmp:
#                         out.write("\t".join([chr, str(1), tmp[chr], chr, "gneg"]) + "\n")
#                 cmd = f'python3 scripts/parse_and_plot.py --edta data/EDTA/{name}.sanitised.fa.mod.EDTA.TEanno.gff3 --genome data/EDTA/{name}_partial_metadata.txt --windowsize 500000 --chromtoplot 2R,2L,3R,3L,X --name {name}'
#                 subprocess.run(cmd.split(" "), check=True)
#                 #--classesorder unknown,LINE,MITE,DNA,LTR

chrom_of_interest = ["chr2R", "chr2L", "chr3R", "chr3L", "chrX"]

for f in os.listdir("data/EDTA/"):
    if f.endswith("_scaff_lengths.tsv"):
        name = f.split("_scaff")[0]
        if os.path.isfile("data/EDTA/" + name + ".sanitised.fa.mod.EDTA.TEanno.gff3"):
            print(name)
            tmp = {}
            for line in open("data/EDTA/" + f).readlines():
                if line[0] != "#":
                    chrom, length = line[:-1].split("\t")
                    if chrom in chrom_of_interest:
                        tmp[chrom] = length
            if len(tmp) == 5 :
                with open("data/EDTA/" + name + "_partial_metadata.txt", "w") as out:
                    out.write("chr\tstart\tend\tname\tgieStain\n")
                    for chr in tmp:
                        out.write("\t".join([chr, str(1), tmp[chr], chr, "gneg"]) + "\n")
                cmd = f'python3 scripts/parse_and_plot.py --edta data/EDTA/{name}.sanitised.fa.mod.EDTA.TEanno.gff3 --genome data/EDTA/{name}_partial_metadata.txt --windowsize 500000 --chromtoplot chr2R,chr2L,chr3R,chr3L,chrX --name {name}'
                subprocess.run(cmd.split(" "), check=True)
                #--classesorder unknown,LINE,MITE,DNA,LTR


chrom_of_interest = ["chrX", "chr2", "chr3"]

for f in os.listdir("data/EDTA/"):
    if f.endswith("_scaff_lengths.tsv"):
        name = f.split("_scaff")[0]
        if os.path.isfile("data/EDTA/" + name + ".sanitised.fa.mod.EDTA.TEanno.gff3"):
            print(name)
            tmp = {}
            for line in open("data/EDTA/" + f).readlines():
                if line[0] != "#":
                    chrom, length = line[:-1].split("\t")
                    if chrom in chrom_of_interest:
                        tmp[chrom] = length
            if len(tmp) == 3 :
                with open("data/EDTA/" + name + "_partial_metadata.txt", "w") as out:
                    out.write("chr\tstart\tend\tname\tgieStain\n")
                    for chr in tmp:
                        out.write("\t".join([chr, str(1), tmp[chr], chr, "gneg"]) + "\n")
                cmd = f'python3 scripts/parse_and_plot.py --edta data/EDTA/{name}.sanitised.fa.mod.EDTA.TEanno.gff3 --genome data/EDTA/{name}_partial_metadata.txt --windowsize 500000 --chromtoplot chr2,chr3,chrX --name {name}'
                subprocess.run(cmd.split(" "), check=True)


chrom_of_interest = ["2RL", "3RL", "X"]

for f in os.listdir("data/EDTA/"):
    if f.endswith("_scaff_lengths.tsv"):
        name = f.split("_scaff")[0]
        if os.path.isfile("data/EDTA/" + name + ".sanitised.fa.mod.EDTA.TEanno.gff3"):
            print(name)
            tmp = {}
            for line in open("data/EDTA/" + f).readlines():
                if line[0] != "#":
                    chrom, length = line[:-1].split("\t")
                    if chrom in chrom_of_interest:
                        tmp[chrom] = length
            if len(tmp) == 3 :
                with open("data/EDTA/" + name + "_partial_metadata.txt", "w") as out:
                    out.write("chr\tstart\tend\tname\tgieStain\n")
                    for chr in tmp:
                        out.write("\t".join([chr, str(1), tmp[chr], chr, "gneg"]) + "\n")
                cmd = f'python3 scripts/parse_and_plot.py --edta data/EDTA/{name}.sanitised.fa.mod.EDTA.TEanno.gff3 --genome data/EDTA/{name}_partial_metadata.txt --windowsize 500000 --chromtoplot 2RL,3RL,X --name {name}'
                subprocess.run(cmd.split(" "), check=True)


