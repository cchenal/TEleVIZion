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
  make_option(c("-a", "--accessibility"), type="character", default="not_displayed", help="genome accessibility file name, tabulated: [chr, start, end, name, itemRgb]", metavar="character"), 
  make_option(c("-i", "--input"), type="character", default=NULL, help="input file name, tabulated: [chr, start, end, <areas>]", metavar="character"), 
  make_option(c("-k", "--classesorder"), type="character", default=NULL, help="classes order, should be reversed for display purposes, use comma", metavar="character"), 
  make_option(c("-l", "--colorsorder"), type="character", default=NULL, help="colors order, should be reversed for display purposes, use comma", metavar="character"), 
  make_option(c("-o", "--output"), type="character", default=NULL, help="output prefix file name [default= %default]", metavar="character")
); 
 
opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

if (is.null(opt$genome) | is.null(opt$input) | is.null(opt$output)){
    print_help(opt_parser)
}

if (is.null(opt$genome)){
  stop("Missing argument: -g <input file>", call.=FALSE)
}

if (is.null(opt$input)){
  stop("Missing argument: -i <input file>", call.=FALSE)
}

if (is.null(opt$classesorder)){
  stop("Missing argument: -k <list,of,classes>", call.=FALSE)
} else {
  classes_order <- unlist(strsplit(x = opt$classesorder, split = ","))
}

if (is.null(opt$colorsorder)){
  stop("Missing argument: -l <list,of,colors>", call.=FALSE)
} else {
  colors_order <- unlist(strsplit(x = opt$colorsorder, split = ","))
}

if (is.null(opt$output)){
  stop("Missing argument: -o <output prefix>", call.=FALSE)
}


cat("\n### CODE ###\n\n")

### Get the data 

data <- read.csv(opt$input, sep="\t")


### Retrieve chromosomes informations

genome <- toGRanges(opt$genome) # NB: opt$genome should have a "\n" as last character

if (is.character(opt$chromosomes)){
    chromosome_order <- unlist(strsplit(x = opt$chromosomes, split = ","))
} else {
    chromosome_order <- genome@seqinfo@seqnames
}

# if (is.character(opt$accessibility)){
#     accessibility <- toGRanges(opt$accessibility)
# }


if (opt$accessibility != "not_displayed"){
  accessibility <- toGRanges(opt$accessibility)
}

### Percentage of repeats along the genome

file <- paste0(opt$output, "karyoplot_stacked_percentage_by_class.png")
png(file, width = 29.7, height = 8, units = "cm", res = 300)

# Graphical parameters 
pp <- getDefaultPlotParams(plot.type=4)
pp$data1inmargin <- 12 
pp$leftmargin <- 0.085
pp$rightmargin <- 0.195
pp$ideogramlateralmargin <- 0.01
pp$data2inmargin <- 10

# Plotting the outline of each chromosome
kp <- plotKaryotype(
    genome = genome, 
    # cytobands = genome, 
    chromosomes = chromosome_order, 
    plot.type = 4, 
    plot.params = pp,
    cex = 0.8
)

# Painting loci according to their accessibility
if (opt$accessibility != "not_displayed"){
    kpAddLabels(
        karyoplot = kp, 
        labels = "Accessibility", 
        data.panel = "ideogram", 
        cex = 0.8, 
        label.margin = 0.01,
        family="Tahoma")
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
    tick.col = "black",
    family="Tahoma", 
    cex = 0.5
)

# Background
kpDataBackground(kp, data.panel = 1, col="grey96")


# Area
kpAddLabels(kp, labels="Percentage\nof repeated\ncontent along\nthe genome", data.panel=1, r0=0, r1=1, cex = 0.8, family="Tahoma")
columns <- paste0(classes_order, "_pct_stacked")
for (i in 1:length(columns)){
    kpArea(kp, chr=data$chrom, x=data$barycenter, y=data[[columns[i]]], col=colors_order[i], border="NA", r0=0, r1=1)
}    

# Add the plot's scale
kpAxis(
    karyoplot = kp, 
    data.panel = 1, 
    # ymin = 0, 
    # ymax = 100, 
    r0 = 0, 
    r1 = 1, 
    side = "right", 
    numticks = 5,
    labels = c("0%", "25%", "50%", "75%", "100%"), 
    cex = 0.8
)

# Legend
twidth  <- max(strwidth(classes_order, units="user"))

legend(
  x       = 0.87,
  y       = 1.25,
  legend  = classes_order,
  fill    = colors_order,
  border  = "grey5",
  bty     = "o",
  box.lwd = 0.3,
  title   = "Repeat class",
  cex     = 0.8,
  text.width = twidth,
  xpd = TRUE
)

# if (is.character(opt$accessibility)){
if (opt$accessibility != "not_displayed"){
    legend(
    x       = 0.87,
    y       = 0,
    legend  = c("Low", "High"),
    fill    = c("grey5", "grey95"),
    border  = "grey5", #NA
    bty     = "o",
    box.lwd = 0.3,
    title   = "Accessibility",
    cex     = 0.8,
    text.width = twidth,
    xpd = TRUE
    )
}

dev.off()




### Number of insertions along the genome 

file <- paste0(opt$output, "karyoplot_stacked_counts_by_class.png")
png(file, width = 29.7, height = 8, units = "cm", res = 300)

# Graphical parameters 
pp <- getDefaultPlotParams(plot.type=4)
pp$data1inmargin <- 12 
pp$leftmargin <- 0.085
pp$rightmargin <- 0.195
pp$ideogramlateralmargin <- 0.01
pp$data2inmargin <- 10

# Plotting the outline of each chromosome
kp <- plotKaryotype(
    genome = genome, 
    # cytobands = genome, 
    chromosomes = chromosome_order, 
    plot.type = 4, 
    plot.params = pp,
    cex = 0.8
)

# Painting loci according to their accessibility
if (opt$accessibility != "not_displayed"){
    kpAddLabels(
        karyoplot = kp, 
        labels = "Accessibility", 
        data.panel = "ideogram", 
        cex = 0.8, 
        label.margin = 0.01,
        family="Tahoma")
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
    tick.col = "black",
    family="Tahoma", 
    cex = 0.5
)

# Background
kpDataBackground(kp, data.panel = 1, col="grey96")


# Area
kpAddLabels(kp, labels="Number of\ninsertions\nalong the\ngenome", data.panel=1, r0=0, r1=1, cex = 0.8, family="Tahoma")
columns <- paste0(classes_order, "_count_stacked")
overall_max <- max(data[, columns])
for (i in 1:length(columns)){
    kpArea(kp, chr=data$chrom, x=data$barycenter, y=data[[columns[i]]]/overall_max, col=colors_order[i], border="NA", r0=0, r1=1)
}    

# Add the plot's scale
kpAxis(
    karyoplot = kp, 
    data.panel = 1, 
    ymin = 0, 
    ymax = overall_max, 
    r0 = 0, 
    r1 = 1, 
    side = "right", 
    numticks = 5,
    cex = 0.8
)

# Legend
twidth  <- max(strwidth(classes_order, units="user"))

legend(
  x       = 0.87,
  y       = 1.25,
  legend  = classes_order,
  fill    = colors_order,
  border  = "grey5",
  bty     = "o",
  box.lwd = 0.3,
  title   = "Repeat class",
  cex     = 0.8,
  text.width = twidth,
  xpd = TRUE
)

if (opt$accessibility != "not_displayed"){
    legend(
    x       = 0.87,
    y       = 0,
    legend  = c("Low", "High"),
    fill    = c("grey5", "grey95"),
    border  = "grey5", #NA
    bty     = "o",
    box.lwd = 0.3,
    title   = "Accessibility",
    cex     = 0.8,
    text.width = twidth,
    xpd = TRUE
    )
}
