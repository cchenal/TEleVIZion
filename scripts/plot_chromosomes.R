cat("\n### IMPORTING LIBRARIES ###\n\n")

library(karyoploteR)
library(optparse)
 
cat("\n### PARSING ARGUMENTS ###\n\n")

# option_list = list(
#   make_option(c("-f", "--file"), type="character", default=NULL, help="dataset file name", metavar="character"), 
#   make_option(c("-o", "--out"), type="character", default="out.txt", help="output file name [default= %default]", metavar="character")
# ); 

option_list = list(
  make_option(c("-g", "--genome"), type="character", default=NULL, help="genome file name, tabulated: [chr, start, end, name, gieStain]", metavar="character"), 
  make_option(c("-c", "--chromosomes"), type="character", default=NULL, help="chromosome order, use comma", metavar="character"), 
  make_option(c("-a", "--accessibility"), type="character", default=NULL, help="genome accessibility file name, tabulated: [chr, start, end, name, itemRgb]", metavar="character"), 
  make_option(c("-o", "--out"), type="character", default=NULL, help="output file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$genome) | is.null(opt$out)){
    print_help(opt_parser)
}

if (is.null(opt$genome)){
  stop("Missing argument: -g <input file>", call.=FALSE)
}

if (is.null(opt$out)){
  stop("Missing argument: -o <output file>", call.=FALSE)
}


cat("\n### CODE ###\n\n")

### Retrieve chromosomes informations

genome <- toGRanges(opt$genome) # NB: opt$genome should have a "\n" as last character

if (is.character(opt$chromosomes)){
    chromosome_order <- unlist(strsplit(x = opt$chromosomes, split = ","))
} else {
    chromosome_order <- genome@seqinfo@seqnames
}

if (is.character(opt$accessibility)){
    # acc_table <- read.table(file = opt$accessibility)
    # accessibility <- toGRanges(
    #     seqnames = acc_table$chromosome,
    #     ranges = IRanges(start = acc_table$start, end = acc_table$end),
    #     name = acc_table$name,
    #     itemRgb = acc_table$itemRgb)
    accessibility <- toGRanges(opt$accessibility)
}

### Plotting 

png(opt$out, width = 29.7, height = 8, units = "cm", res = 300)

# Graphical parameters 
pp <- getDefaultPlotParams(plot.type=4)
pp$data1inmargin <- 25 
pp$leftmargin <- 0.085
pp$rightmargin <- 0.095
pp$ideogramlateralmargin <- 0.01
pp$data2inmargin <- 10

# Plotting the outline of each chromosome
kp <- plotKaryotype(
    genome = genome, 
    cytobands = genome, 
    chromosomes = chromosome_order, 
    plot.type = 4, 
    plot.params = pp
)

# Painting loci according to their accessibility
if (is.character(opt$accessibility)){
    kpAddLabels(
        karyoplot = kp, 
        labels = "Accessibility", 
        data.panel = "ideogram", 
        cex = 0.8, 
        label.margin = 0.01)
    kpPlotRegions(
        karyoplot = kp, 
        data = accessibility, 
        col = accessibility$itemRgb, 
        data.panel = "ideogram")
}

# Adding the scale along each chromosome 
kpAddBaseNumbers(
    karyoplot = kp, 
    tick.dist = 10000000, 
    tick.len = 3, 
    tick.col = "black" #, cex = 2
)

# kpAddCytobandLabels(
#     karyoplot = kp
# )

# Add the plot's scale
# kpAxis(
#     karyoplot = kp, 
#     data.panel = 1, 
#     ymin = 0, 
#     ymax = 100, 
#     r0 = 0, 
#     r1 = 1, 
#     side = "right", 
#     numticks = 5, 
#     cex = 0.8
# )

# dev.off()