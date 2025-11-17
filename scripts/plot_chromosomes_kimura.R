# cat("\n### IMPORTING LIBRARIES ###\n\n")

# library(karyoploteR)
# library(optparse)
 
# cat("\n### PARSING ARGUMENTS ###\n\n")

# # option_list = list(
# #   make_option(c("-f", "--file"), type="character", default=NULL, help="dataset file name", metavar="character"), 
# #   make_option(c("-o", "--out"), type="character", default="out.txt", help="output file name [default= %default]", metavar="character")
# # ); 

# option_list = list(
#   make_option(c("-n", "--name"), type="character", default=NULL, help="prefix to use for outputs", metavar="character"), 
#   make_option(c("-g", "--genome"), type="character", default=NULL, help="genome file name, tabulated: [chr, start, end, name, gieStain]", metavar="character"), 
#   make_option(c("-c", "--chromosomes"), type="character", default=NULL, help="chromosome order, use comma", metavar="character"), 
#   make_option(c("-a", "--accessibility"), type="character", default="not_displayed", help="genome accessibility file name, tabulated: [chr, start, end, name, itemRgb]", metavar="character"), 
#   make_option(c("-i", "--input"), type="character", default=NULL, help="input file name, tabulated: [chr, start, end, <areas>]", metavar="character"), 
#   make_option(c("-k", "--classesorder"), type="character", default=NULL, help="classes order, should be reversed for display purposes, use comma", metavar="character"), 
#   make_option(c("-p", "--perclass"), type="character", default=NULL, help="plot additional karyoplots", metavar="character"), 
#   make_option(c("-l", "--colorsorder"), type="character", default=NULL, help="colors order, should be reversed for display purposes, use comma", metavar="character"), 
#   make_option(c("-o", "--output"), type="character", default=NULL, help="output prefix file name [default= %default]", metavar="character"),
#   make_option(c("-m", "--kimura"), type="character", default=NULL, help="Additional panel [RepeatMasker only}", metavar="character")
# ); 
 
# opt_parser = OptionParser(option_list=option_list);
# opt = parse_args(opt_parser);

# if (is.null(opt$genome) | is.null(opt$input) | is.null(opt$output)){
#     print_help(opt_parser)
# }

# if (is.null(opt$genome)){
#   stop("Missing argument: -g <input file>", call.=FALSE)
# }

# if (is.null(opt$input)){
#   stop("Missing argument: -i <input file>", call.=FALSE)
# }

# if (is.null(opt$classesorder)){
#   stop("Missing argument: -k <list,of,classes>", call.=FALSE)
# } else {
#   classes_order <- unlist(strsplit(x = opt$classesorder, split = ","))
# }

# if (is.null(opt$colorsorder)){
#   stop("Missing argument: -l <list,of,colors>", call.=FALSE)
# } else {
#   colors_order <- unlist(strsplit(x = opt$colorsorder, split = ","))
# }

# if (is.null(opt$output)){
#   stop("Missing argument: -o <output prefix>", call.=FALSE)
# }

# if (! is.null(opt$perclass)){
#   if (is.character(opt$perclass) && opt$perclass %in% c("False", "None")){
#     opt$perclass <- NULL
#   }
# }

# if (! is.null(opt$kimura)){
#   if (is.character(opt$kimura) && opt$kimura %in% c("False", "None")){
#     opt$perclass <- NULL
#   }
# }

# cat("\n### CODE ###\n\n")

# ### Get the data 

# data <- read.csv(opt$input, sep="\t")
# data_kimura <- read.csv(opt$kimura, sep="\t")
# name <- opt$name

# ### Retrieve chromosomes informations

# genome <- toGRanges(opt$genome) # NB: opt$genome should have a "\n" as last character

# if (is.character(opt$chromosomes)){
#     chromosome_order <- unlist(strsplit(x = opt$chromosomes, split = ","))
# } else {
#     chromosome_order <- genome@seqinfo@seqnames
# }
# names_vec <- mcols(genome)$name
# df_names <- data.frame(chr = as.character(seqnames(genome)), name = names_vec)
# ordered_names <- unlist(lapply(chromosome_order, function(chr){df_names$name[df_names$chr == chr]}))

# if (opt$accessibility != "not_displayed"){
#   accessibility <- toGRanges(opt$accessibility)
# }

# ### Kimura parameters

# kimura_cats <- c("40-70", "30-40", "20-30", "10-20", "0-10")
# kimura_colors <- c("#2c699a", "#0db39e", "#83e377", "#efea5a", "#f29e4c")



# ### Percentage of repeats along the genome


# ##### All classes

# file <- paste0(opt$output, "_karyoplot_stacked_percentage_by_class.pdf")
# # png(file, width = 29.7, height = 10, units = "cm", res = 300)

# # Graphical parameters 
# pp <- getDefaultPlotParams(plot.type=3)
# pp$data1inmargin <- 14 
# # pp$leftmargin <- 0.085
# pp$topmargin <- 100
# pp$leftmargin <- 0.02
# pp$rightmargin <- 0.24
# pp$ideogramlateralmargin <- 0.01
# pp$data2inmargin <- 20
# pp$data2height <- pp$data1height / 3

# # Plotting the outline of each chromosome
# kp <- plotKaryotype(
#     genome = genome, 
#     # cytobands = genome, 
#     chromosomes = chromosome_order, 
#     plot.type = 3, 
#     plot.params = pp,
#     labels.plotter = NULL
# )
# kpAddChromosomeNames(kp, chr.names = ordered_names, cex = 0.8)

# # Painting loci according to their accessibility
# if (opt$accessibility != "not_displayed"){
#     # kpAddLabels(
#     #     karyoplot = kp, 
#     #     labels = "Accessibility", 
#     #     data.panel = "ideogram", 
#     #     cex = 0.8, 
#     #     label.margin = 0.01) # family="Tahoma")
#     kpPlotRegions(
#         karyoplot = kp, 
#         data = accessibility, 
#         col = accessibility$itemRgb, 
#         data.panel = "ideogram")
# }

# # Adding the scale along each chromosome 
# kpAddBaseNumbers(
#     karyoplot = kp, 
#     tick.dist = 10000000, 
#     tick.len = 3, 
#     tick.col = "black", # family="Tahoma" 
#     cex = 0.5
# )

# # Background
# kpDataBackground(kp, data.panel = 1, col="grey96")
# kpDataBackground(kp, data.panel = 2, col="grey96")


# # Area
# ## TE
# columns <- paste0(classes_order, "_pct_stacked")
# for (i in 1:length(columns)){
#     kpArea(kp, chr=data$chrom, x=data$barycenter, y=data[[columns[i]]], col=colors_order[i], border="NA", r0=0, r1=1, data.panel = 1)
# }    

# kpAxis(
#     karyoplot = kp, 
#     data.panel = 1, 
#     r0 = 0, 
#     r1 = 1, 
#     side = "right", 
#     numticks = 5,
#     labels = c("0%", "25%", "50%", "75%", "100%"), 
#     cex = 0.7
# )

# ## K2p
# columns <- paste0("ALL_", gsub("-", ".", kimura_cats), "_pct_stacked")
# for (i in 1:length(columns)){
#     kpArea(kp, chr=data$chrom, x=data$barycenter, y=data_kimura[[columns[i]]], col=kimura_colors[i], border="NA", r0=1, r1=0, data.panel = 2)
# }    

# kpAxis(
#     karyoplot = kp, 
#     data.panel = 2, 
#     r0 = 1, 
#     r1 = 0, 
#     side = "right", 
#     numticks = 3,
#     labels = c("0%", "50%", "100%"), 
#     cex = 0.7
# )

# # Legend
# twidth  <- max(strwidth(classes_order, units="user"))

# legend(
#   x       = 0.82,
#   y       = 1,
#   legend  = classes_order,
#   fill    = colors_order,
#   border  = "grey5",
#   bty     = "o",
#   box.lwd = 0.3,
#   title   = "Repeat class",
#   cex     = 0.8,
#   text.width = twidth,
#   xpd = TRUE
# )

# legend(
#   x       = 0.82,
#   y       = 0.1,
#   legend  = kimura_cats,
#   fill    = kimura_colors,
#   border  = "grey5",
#   bty     = "o",
#   box.lwd = 0.3,
#   title   = "K2p",
#   cex     = 0.8,
#   text.width = twidth,
#   xpd = TRUE
# )

# # # if (is.character(opt$accessibility)){
# # if (opt$accessibility != "not_displayed"){
# #     legend(
# #     x       = 0.82,
# #     y       = -0.25,
# #     legend  = c("Low", "High"),
# #     fill    = c("grey5", "grey95"),
# #     border  = "grey5", #NA
# #     bty     = "o",
# #     box.lwd = 0.3,
# #     title   = "Accessibility",
# #     cex     = 0.8,
# #     text.width = twidth,
# #     xpd = TRUE
# #     )
# # }

# # Title
# title(paste0("Percentage of repeated content along chromosomes - ", name), cex.main = 0.8, line = 2.5)

# dev.off()




# ##### Per class

# if (!is.null(opt$perclass)){
#   for (cls in classes_order) {

#     file <- paste0(opt$output, "_karyoplot_percentage_", cls, ".pdf")
#     # png(file, width = 29.7, height = 10, units = "cm", res = 300)

#     # Graphical parameters 
#     pp <- getDefaultPlotParams(plot.type=3)
#     pp$data1inmargin <- 14 
#     # pp$leftmargin <- 0.085
#     pp$topmargin <- 100
#     pp$leftmargin <- 0.02
#     pp$rightmargin <- 0.24
#     pp$ideogramlateralmargin <- 0.01
#     pp$data2inmargin <- 20
#     pp$data2height <- pp$data1height / 3

#     # Plotting the outline of each chromosome
#     kp <- plotKaryotype(
#         genome = genome, 
#         # cytobands = genome, 
#         chromosomes = chromosome_order, 
#         plot.type = 3, 
#         plot.params = pp,
#         labels.plotter = NULL
#     )
#     kpAddChromosomeNames(kp, chr.names = ordered_names, cex = 0.8)

#     # Painting loci according to their accessibility
#     if (opt$accessibility != "not_displayed"){
#         # kpAddLabels(
#         #     karyoplot = kp, 
#         #     labels = "Accessibility", 
#         #     data.panel = "ideogram", 
#         #     cex = 0.8, 
#         #     label.margin = 0.01) # family="Tahoma")
#         kpPlotRegions(
#             karyoplot = kp, 
#             data = accessibility, 
#             col = accessibility$itemRgb, 
#             data.panel = "ideogram")
#     }

#     # Adding the scale along each chromosome 
#     kpAddBaseNumbers(
#         karyoplot = kp, 
#         tick.dist = 10000000, 
#         tick.len = 3, 
#         tick.col = "black", # family="Tahoma" 
#         cex = 0.5
#     )

#     # Background
#     kpDataBackground(kp, data.panel = 1, col="grey96")
#     kpDataBackground(kp, data.panel = 2, col="grey96")

#     # Area
#     ## TE
#     columns <- paste0(classes_order, "_pct")
#     for (i in 1:length(columns)){
#       if (columns[i] == paste0(cls, "_pct")) {
#         kpArea(kp, chr=data$chrom, x=data$barycenter, y=data[[columns[i]]], col=colors_order[i], border="NA", r0=0, r1=1)
#       }
#     }

#     kpAxis(
#         karyoplot = kp, 
#         data.panel = 1, 
#         r0 = 0, 
#         r1 = 1, 
#         side = "right", 
#         numticks = 5,
#         labels = c("0%", "25%", "50%", "75%", "100%"), 
#         cex = 0.7
#     )

#     ## K2p
#     columns <- paste0(cls, "_", gsub("-", ".", kimura_cats), "_pct_stacked")
#     for (i in 1:length(columns)){
#         kpArea(kp, chr=data$chrom, x=data$barycenter, y=data_kimura[[columns[i]]], col=kimura_colors[i], border="NA", r0=1, r1=0, data.panel = 2)
#     }    

#     kpAxis(
#         karyoplot = kp, 
#         data.panel = 2, 
#         r0 = 1, 
#         r1 = 0, 
#         side = "right", 
#         numticks = 3,
#         labels = c("0%", "50%", "100%"), 
#         cex = 0.7
#     )

#     # Legend
#     twidth  <- max(strwidth(classes_order, units="user"))

#     legend(
#       x       = 0.82,
#       y       = 1,
#       legend  = classes_order,
#       fill    = colors_order,
#       border  = "grey5",
#       bty     = "o",
#       box.lwd = 0.3,
#       title   = "Repeat class",
#       cex     = 0.8,
#       text.width = twidth,
#       xpd = TRUE
#     )

#     legend(
#       x       = 0.82,
#       y       = 0.1,
#       legend  = kimura_cats,
#       fill    = kimura_colors,
#       border  = "grey5",
#       bty     = "o",
#       box.lwd = 0.3,
#       title   = "K2p",
#       cex     = 0.8,
#       text.width = twidth,
#       xpd = TRUE
#     )


#     # # if (is.character(opt$accessibility)){
#     # if (opt$accessibility != "not_displayed"){
#     #     legend(
#     #     x       = 0.82,
#     #     y       = -0.25,
#     #     legend  = c("Low", "High"),
#     #     fill    = c("grey5", "grey95"),
#     #     border  = "grey5", #NA
#     #     bty     = "o",
#     #     box.lwd = 0.3,
#     #     title   = "Accessibility",
#     #     cex     = 0.8,
#     #     text.width = twidth,
#     #     xpd = TRUE
#     #     )
#     # }

#     # Title
#     title(paste0("Percentage of ", cls, " content along chromosomes - ", name), cex.main = 0.8, line = 2.5)

#     dev.off()
#   }
# }





# ### Number of insertions along the genome 

# ##### All classes

# file <- paste0(opt$output, "_karyoplot_stacked_counts_by_class.pdf")
# # png(file, width = 29.7, height = 10, units = "cm", res = 300)

# # Graphical parameters 
# pp <- getDefaultPlotParams(plot.type=3)
# pp$data1inmargin <- 14 
# # pp$leftmargin <- 0.085
# pp$topmargin <- 100
# pp$leftmargin <- 0.02
# pp$rightmargin <- 0.24
# pp$ideogramlateralmargin <- 0.01
# pp$data2inmargin <- 20
# pp$data2height <- pp$data1height / 3

# # Plotting the outline of each chromosome
# kp <- plotKaryotype(
#     genome = genome, 
#     # cytobands = genome, 
#     chromosomes = chromosome_order, 
#     plot.type = 3, 
#     plot.params = pp,
#     labels.plotter = NULL
# )
# kpAddChromosomeNames(kp, chr.names = ordered_names, cex = 0.8)

# # Painting loci according to their accessibility
# if (opt$accessibility != "not_displayed"){
#     # kpAddLabels(
#     #     karyoplot = kp, 
#     #     labels = "Accessibility", 
#     #     data.panel = "ideogram", 
#     #     cex = 0.8, 
#     #     label.margin = 0.01) # family="Tahoma")
#     kpPlotRegions(
#         karyoplot = kp, 
#         data = accessibility, 
#         col = accessibility$itemRgb, 
#         data.panel = "ideogram")
# }

# # Adding the scale along each chromosome 
# kpAddBaseNumbers(
#     karyoplot = kp, 
#     tick.dist = 10000000, 
#     tick.len = 3, 
#     tick.col = "black", # family="Tahoma" 
#     cex = 0.5
# )

# # Background
# kpDataBackground(kp, data.panel = 1, col="grey96")
# kpDataBackground(kp, data.panel = 2, col="grey96")

# # Area
# # TE
# columns <- paste0(classes_order, "_count_stacked")
# overall_max <- max(data[, columns])
# for (i in 1:length(columns)){
#     kpArea(kp, chr=data$chrom, x=data$barycenter, y=data[[columns[i]]]/overall_max, col=colors_order[i], border="NA", r0=0, r1=1, data.panel = 1)
# }    

# kpAxis(
#     karyoplot = kp, 
#     data.panel = 1, 
#     ymin = 0, 
#     ymax = overall_max, 
#     r0 = 0, 
#     r1 = 1, 
#     side = "right", 
#     numticks = 5,
#     cex = 0.7
# )

# ## K2p
# columns <- paste0("ALL_", gsub("-", ".", kimura_cats), "_pct_stacked")
# for (i in 1:length(columns)){
#     kpArea(kp, chr=data$chrom, x=data$barycenter, y=data_kimura[[columns[i]]], col=kimura_colors[i], border="NA", r0=1, r1=0, data.panel = 2)
# }    

# kpAxis(
#     karyoplot = kp, 
#     data.panel = 2, 
#     r0 = 1, 
#     r1 = 0, 
#     side = "right", 
#     numticks = 3,
#     labels = c("0%", "50%", "100%"), 
#     cex = 0.7
# )

# # Legend
# twidth  <- max(strwidth(classes_order, units="user"))

# legend(
#   x       = 0.82,
#   y       = 1,
#   legend  = classes_order,
#   fill    = colors_order,
#   border  = "grey5",
#   bty     = "o",
#   box.lwd = 0.3,
#   title   = "Repeat class",
#   cex     = 0.8,
#   text.width = twidth,
#   xpd = TRUE
# )

# legend(
#   x       = 0.82,
#   y       = 0.1,
#   legend  = kimura_cats,
#   fill    = kimura_colors,
#   border  = "grey5",
#   bty     = "o",
#   box.lwd = 0.3,
#   title   = "K2p",
#   cex     = 0.8,
#   text.width = twidth,
#   xpd = TRUE
# )

# # if (opt$accessibility != "not_displayed"){
# #     legend(
# #     x       = 0.82,
# #     y       = -0.25,
# #     legend  = c("Low", "High"),
# #     fill    = c("grey5", "grey95"),
# #     border  = "grey5", #NA
# #     bty     = "o",
# #     box.lwd = 0.3,
# #     title   = "Accessibility",
# #     cex     = 0.8,
# #     text.width = twidth,
# #     xpd = TRUE
# #     )
# # }

# # Title
# title(paste0("Number of insertions along chromosomes - ", name), cex.main = 0.8, line = 2.5)


# ##### Per class 

# if (!is.null(opt$perclass)){
#   for (cls in classes_order) {
#     file <- paste0(opt$output, "_karyoplot_counts_", cls, ".pdf")
#     # png(file, width = 29.7, height = 10, units = "cm", res = 300)

#     # Graphical parameters 
#     pp <- getDefaultPlotParams(plot.type=3)
#     pp$data1inmargin <- 14 
#     # pp$leftmargin <- 0.085
#     pp$topmargin <- 100
#     pp$leftmargin <- 0.02
#     pp$rightmargin <- 0.24
#     pp$ideogramlateralmargin <- 0.01
#     pp$data2inmargin <- 20
#     pp$data2height <- pp$data1height / 3

#     # Plotting the outline of each chromosome
#     kp <- plotKaryotype(
#         genome = genome, 
#         # cytobands = genome, 
#         chromosomes = chromosome_order, 
#         plot.type = 3, 
#         plot.params = pp,
#         labels.plotter = NULL
#     )
#     kpAddChromosomeNames(kp, chr.names = ordered_names, cex = 0.8)

#     # Painting loci according to their accessibility
#     if (opt$accessibility != "not_displayed"){
#         # kpAddLabels(
#         #     karyoplot = kp, 
#         #     labels = "Accessibility", 
#         #     data.panel = "ideogram", 
#         #     cex = 0.8, 
#         #     label.margin = 0.01) # family="Tahoma")
#         kpPlotRegions(
#             karyoplot = kp, 
#             data = accessibility, 
#             col = accessibility$itemRgb, 
#             data.panel = "ideogram")
#     }

#     # Adding the scale along each chromosome 
#     kpAddBaseNumbers(
#         karyoplot = kp, 
#         tick.dist = 10000000, 
#         tick.len = 3, 
#         tick.col = "black", # family="Tahoma" 
#         cex = 0.5
#     )

#     # Background
#     kpDataBackground(kp, data.panel = 1, col="grey96")
#     kpDataBackground(kp, data.panel = 2, col="grey96")

#     # Area
#     ## TE
#     columns <- paste0(classes_order, "_count")
#     overall_max <- max(data[, columns])
#     for (i in 1:length(columns)){
#       if (columns[i] == paste0(cls, "_count")) {
#         kpArea(kp, chr=data$chrom, x=data$barycenter, y=data[[columns[i]]]/overall_max, col=colors_order[i], border="NA", r0=0, r1=1)
#       }
#     }    

#     kpAxis(
#         karyoplot = kp, 
#         data.panel = 1, 
#         ymin = 0, 
#         ymax = overall_max, 
#         r0 = 0, 
#         r1 = 1, 
#         side = "right", 
#         numticks = 5,
#         cex = 0.7
#     )

#     ## K2p
#     columns <- paste0(cls, "_", gsub("-", ".", kimura_cats), "_pct_stacked")
#     for (i in 1:length(columns)){
#         kpArea(kp, chr=data$chrom, x=data$barycenter, y=data_kimura[[columns[i]]], col=kimura_colors[i], border="NA", r0=1, r1=0, data.panel = 2)
#     }    

#     kpAxis(
#         karyoplot = kp, 
#         data.panel = 2, 
#         r0 = 1, 
#         r1 = 0, 
#         side = "right", 
#         numticks = 3,
#         labels = c("0%", "50%", "100%"), 
#         cex = 0.7
#     )

#     # Legend
#     twidth  <- max(strwidth(classes_order, units="user"))

#     legend(
#       x       = 0.82,
#       y       = 1,
#       legend  = classes_order,
#       fill    = colors_order,
#       border  = "grey5",
#       bty     = "o",
#       box.lwd = 0.3,
#       title   = "Repeat class",
#       cex     = 0.8,
#       text.width = twidth,
#       xpd = TRUE
#     )

#     legend(
#       x       = 0.82,
#       y       = 0.1,
#       legend  = kimura_cats,
#       fill    = kimura_colors,
#       border  = "grey5",
#       bty     = "o",
#       box.lwd = 0.3,
#       title   = "K2p",
#       cex     = 0.8,
#       text.width = twidth,
#       xpd = TRUE
#     )

#     # if (opt$accessibility != "not_displayed"){
#     #     legend(
#     #     x       = 0.82,
#     #     y       = -0.25,
#     #     legend  = c("Low", "High"),
#     #     fill    = c("grey5", "grey95"),
#     #     border  = "grey5", #NA
#     #     bty     = "o",
#     #     box.lwd = 0.3,
#     #     title   = "Accessibility",
#     #     cex     = 0.8,
#     #     text.width = twidth,
#     #     xpd = TRUE
#     #     )
#     # }

#     # Title
#     title(paste0("Number of insertions of ", cls, " along chromosomes - ", name), cex.main = 0.8, line = 2.5)

#   }
# }

cat("\n### IMPORTING LIBRARIES ###\n\n")
library(karyoploteR)
library(optparse)

cat("\n### PARSING ARGUMENTS ###\n\n")

option_list <- list(
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="prefix to use for outputs", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="genome file name, tabulated: [chr, start, end, name, gieStain]", metavar="character"),
  make_option(c("-c", "--chromosomes"), type="character", default=NULL,
              help="chromosome order, use comma", metavar="character"),
  make_option(c("-a", "--accessibility"), type="character", default="not_displayed",
              help="genome accessibility file name, tabulated: [chr, start, end, name, itemRgb]", metavar="character"),
  make_option(c("-i", "--input"), type="character", default=NULL,
              help="input file name, tabulated: [chr, start, end, <areas>]", metavar="character"),
  make_option(c("-k", "--classesorder"), type="character", default=NULL,
              help="classes order, should be reversed for display purposes, use comma", metavar="character"),
  make_option(c("-p", "--perclass"), type="character", default=NULL,
              help="plot additional karyoplots", metavar="character"),
  make_option(c("-l", "--colorsorder"), type="character", default=NULL,
              help="colors order, should be reversed for display purposes, use comma", metavar="character"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="output prefix file name [default= %default]", metavar="character"),
  make_option(c("-m", "--kimura"), type="character", default=NULL,
              help="Additional panel [RepeatMasker only]", metavar="character"),
  make_option(c("-e", "--edtaidentity"), type="character", default=NULL,
              help="Additional panel [EDTA only]", metavar="character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# ---- Required args
if (is.null(opt$genome) || is.null(opt$input) || is.null(opt$output)) {
  print_help(opt_parser)
}
if (is.null(opt$genome))      stop("Missing argument: -g <input file>", call. = FALSE)
if (is.null(opt$input))       stop("Missing argument: -i <input file>", call. = FALSE)
if (is.null(opt$classesorder)) stop("Missing argument: -k <list,of,classes>", call. = FALSE)
if (is.null(opt$colorsorder))  stop("Missing argument: -l <list,of,colors>", call. = FALSE)
if (is.null(opt$output))       stop("Missing argument: -o <output prefix>", call. = FALSE)

classes_order <- unlist(strsplit(opt$classesorder, ","))
colors_order  <- unlist(strsplit(opt$colorsorder, ","))

# ---- Normalize flags that come as strings
if (!is.null(opt$perclass)) {
  if (is.character(opt$perclass) && opt$perclass %in% c("False", "None")) {
    opt$perclass <- NULL
  }
}
if (!is.null(opt$kimura)) {
  if (is.character(opt$kimura) && opt$kimura %in% c("False", "None")) {
    opt$kimura <- NULL   # (fixed: previously nulled perclass by mistake)
  }
}
if (!is.null(opt$edtaidentity)) {
  if (is.character(opt$edtaidentity) && opt$edtaidentity %in% c("False", "None")) {
    opt$edtaidentity <- NULL   # (fixed: previously nulled perclass by mistake)
  }
}

cat("\n### CODE ###\n\n")

# ---- Data
data        <- read.csv(opt$input,  sep = "\t")
data_kimura <- if (!is.null(opt$kimura)) read.csv(opt$kimura, sep = "\t") else NULL
data_identity <- if (!is.null(opt$edtaidentity)) read.csv(opt$edtaidentity, sep = "\t") else NULL
name        <- opt$name

# ---- Genome / chromosomes
genome <- toGRanges(opt$genome)  # requires trailing newline in file

chromosome_order <- if (is.character(opt$chromosomes)) {
  unlist(strsplit(opt$chromosomes, ","))
} else {
  genome@seqinfo@seqnames
}

names_vec   <- mcols(genome)$name
df_names    <- data.frame(chr = as.character(seqnames(genome)), name = names_vec)
ordered_names <- unlist(lapply(chromosome_order, function(chr) df_names$name[df_names$chr == chr]))

accessibility <- NULL
if (!identical(opt$accessibility, "not_displayed")) {
  accessibility <- toGRanges(opt$accessibility)
}

# ---- Kimura
kimura_cats   <- c("40-70", "30-40", "20-30", "10-20", "0-10")
kimura_colors <- c("#2c699a", "#0db39e", "#83e377", "#efea5a", "#f29e4c")

# ---- Identity
identity_cats   <- c("0.5-0", "0.75-0.5", "0.9-0.75", "1-0.9")
identity_colors <- c("#0db39e", "#83e377", "#efea5a", "#f29e4c") # "#2c699a", 

# ---------------------------
# Helpers to avoid repetition
# ---------------------------

make_pp <- function() {
  pp <- getDefaultPlotParams(plot.type = 3)
  pp$data1inmargin <- 14
  pp$topmargin <- 100
  pp$leftmargin <- 0.04
  pp$rightmargin <- 0.24
  pp$ideogramlateralmargin <- 0.01
  pp$data2inmargin <- 20
  pp$data2height <- pp$data1height / 3
  pp
}

plot_karyo_base <- function() {
  pp <- make_pp()
  kp <- plotKaryotype(
    genome       = genome,
    chromosomes  = chromosome_order,
    plot.type    = 3,
    plot.params  = pp,
    labels.plotter = NULL,
    cytobands = genome
  )
  kpAddChromosomeNames(kp, chr.names = ordered_names, cex = 0.8)

  if (!is.null(accessibility)) {
    kpPlotRegions(karyoplot = kp, data = accessibility,
                  col = accessibility$itemRgb, data.panel = "ideogram")
  }

  kpAddBaseNumbers(karyoplot = kp, tick.dist = 10000000, tick.len = 3,
                   tick.col = "black", cex = 0.5)

  # kpDataBackground(kp, data.panel = 1, col = "grey98")  
  # kpDataBackground(kp, data.panel = 2, col = "grey98")  # , col = "grey96"

  kpAddLabels(kp, labels = "TE", data.panel = 1, cex=0.8)
  if (!is.null(data_kimura)) {
  kpAddLabels(kp, labels = "K2p", data.panel = 2, cex=0.8)
  }
  if (!is.null(data_identity)) {
    kpAddLabels(kp, labels = "Identity", data.panel = 2, cex=0.8)
  }
  

  kp
}

add_axis_pct <- function(kp, panel = 1, r0 = 0, r1 = 1, labels = c("0%", "25%", "50%", "75%", "100%")) {
  kpAxis(karyoplot = kp, data.panel = panel, r0 = r0, r1 = r1,
         side = "right", numticks = length(labels), labels = labels, cex = 0.7)
}

add_axis_counts <- function(kp, overall_max) {
  kpAxis(karyoplot = kp, data.panel = 1, ymin = 0, ymax = overall_max,
         r0 = 0, r1 = 1, side = "right", numticks = 5, cex = 0.7)
}

draw_k2p_panel <- function(kp, dat) {
  cols <- paste0("ALL_", gsub("-", ".", kimura_cats), "_pct_stacked")
  for (i in seq_along(cols)) {
    kpArea(kp, chr = data$chrom, x = data$barycenter, y = dat[[cols[i]]],
           col = kimura_colors[i], border = "NA", r0 = 1, r1 = 0, data.panel = 2)
  }
  add_axis_pct(kp, panel = 2, r0 = 1, r1 = 0, labels = c("0%", "50%", "100%"))
}

draw_identity_panel <- function(kp, dat) {
  cols <- paste0("ALL_", gsub("-", ".", identity_cats), "_pct_stacked")
  for (i in seq_along(cols)) {
    kpArea(kp, chr = data$chrom, x = data$barycenter, y = dat[[cols[i]]],
           col = identity_colors[i], border = "NA", r0 = 1, r1 = 0, data.panel = 2)
  }
  add_axis_pct(kp, panel = 2, r0 = 1, r1 = 0, labels = c("0%", "50%", "100%"))
}

add_legends <- function() {
  twidth <- max(strwidth(classes_order, units = "user"))
  legend(x = 0.82, y = 1, legend = classes_order, fill = colors_order,
         border = "grey5", bty = "o", box.lwd = 0.3, title = "Repeat class",
         cex = 0.8, text.width = twidth, xpd = TRUE)
  if (!is.null(data_kimura)) {
    legend(x = 0.82, y = 0.1, legend = kimura_cats, fill = kimura_colors,
         border = "grey5", bty = "o", box.lwd = 0.3, title = "K2p",
         cex = 0.8, text.width = twidth, xpd = TRUE)
  }
  if (!is.null(data_identity)) {
  legend(x = 0.82, y = 0.1, legend = identity_cats, fill = identity_colors,
         border = "grey5", bty = "o", box.lwd = 0.3, title = "Identity",
         cex = 0.8, text.width = twidth, xpd = TRUE)
  }
}

# =========================================
# Percentage of repeats along the genome
# =========================================

# ---- All classes (stacked % by class + K2p)
file <- paste0(opt$output, "_karyoplot_stacked_percentage_by_class.pdf")
# png(file, width = 29.7, height = 10, units = "cm", res = 300)
pdf(file, width = 11.417, height = 3.937)
on.exit(dev.off(), add = TRUE)

kp <- plot_karyo_base()

# TE stacked percentages
columns <- paste0(classes_order, "_pct_stacked")
for (i in seq_along(columns)) {
  kpArea(kp, chr = data$chrom, x = data$barycenter, y = data[[columns[i]]],
         col = colors_order[i], border = "NA", r0 = 0, r1 = 1, data.panel = 1)
}
add_axis_pct(kp, panel = 1)

# K2p stacked percentages
if (!is.null(data_kimura)) {
  draw_k2p_panel(kp, data_kimura)
}

if (!is.null(data_identity)) {
  draw_identity_panel(kp, data_identity)
}

add_legends()
title(paste0("Percentage of repeated content along chromosomes - ", name),
      cex.main = 0.8, line = 2.5)
dev.off()

# ---- Per class (% + K2p), if requested
if (!is.null(opt$perclass)) {
  for (cls in classes_order) {
    file <- paste0(opt$output, "_karyoplot_percentage_", cls, ".pdf")
    # png(file, width = 29.7, height = 10, units = "cm", res = 300)
    pdf(file, width = 11.417, height = 3.937)
    on.exit(dev.off(), add = TRUE)

    kp <- plot_karyo_base()

    # Only the selected class percentage
    cols_all <- paste0(classes_order, "_pct")
    for (i in seq_along(cols_all)) {
      if (identical(cols_all[i], paste0(cls, "_pct"))) {
        kpArea(kp, chr = data$chrom, x = data$barycenter, y = data[[cols_all[i]]],
               col = colors_order[i], border = "NA", r0 = 0, r1 = 1)
      }
    }
    add_axis_pct(kp, panel = 1)

    # K2p stacked for that class
    if (!is.null(data_kimura)) {
      columns <- paste0(cls, "_", gsub("-", ".", kimura_cats), "_pct_stacked")
      for (i in seq_along(columns)) {
        kpArea(kp, chr = data$chrom, x = data$barycenter, y = data_kimura[[columns[i]]],
               col = kimura_colors[i], border = "NA", r0 = 1, r1 = 0, data.panel = 2)
      }
      add_axis_pct(kp, panel = 2, r0 = 1, r1 = 0, labels = c("0%", "50%", "100%"))
    }

    if (!is.null(data_identity)) {
      columns <- paste0(cls, "_", gsub("-", ".", identity_cats), "_pct_stacked")
      for (i in seq_along(columns)) {
        kpArea(kp, chr = data$chrom, x = data$barycenter, y = data_identity[[columns[i]]],
               col = identity_colors[i], border = "NA", r0 = 1, r1 = 0, data.panel = 2)
      }
      add_axis_pct(kp, panel = 2, r0 = 1, r1 = 0, labels = c("0%", "50%", "100%"))
    }

    add_legends()
    title(paste0("Percentage of ", cls, " content along chromosomes - ", name),
          cex.main = 0.8, line = 2.5)
    dev.off()
  }
}

# =========================================
# Number of insertions along the genome
# =========================================

# ---- All classes (stacked counts + K2p)
file <- paste0(opt$output, "_karyoplot_stacked_counts_by_class.pdf")
# png(file, width = 29.7, height = 10, units = "cm", res = 300)
pdf(file, width = 11.417, height = 3.937)
on.exit(dev.off(), add = TRUE)

kp <- plot_karyo_base()

columns <- paste0(classes_order, "_count_stacked")
overall_max <- max(data[, columns])
for (i in seq_along(columns)) {
  kpArea(kp, chr = data$chrom, x = data$barycenter, y = data[[columns[i]]]/overall_max,
         col = colors_order[i], border = "NA", r0 = 0, r1 = 1, data.panel = 1)
}
add_axis_counts(kp, overall_max)

if (!is.null(data_kimura)) {
  draw_k2p_panel(kp, data_kimura)
}

if (!is.null(data_identity)) {
  draw_identity_panel(kp, data_identity)
}

add_legends()
title(paste0("Number of insertions along chromosomes - ", name),
      cex.main = 0.8, line = 2.5)
dev.off()

# ---- Per class (counts + K2p), if requested
if (!is.null(opt$perclass)) {
  for (cls in classes_order) {
    file <- paste0(opt$output, "_karyoplot_counts_", cls, ".pdf")
    # png(file, width = 29.7, height = 10, units = "cm", res = 300)
    pdf(file, width = 11.417, height = 3.937)
    on.exit(dev.off(), add = TRUE)

    kp <- plot_karyo_base()

    cols_all <- paste0(classes_order, "_count")
    overall_max <- max(data[, cols_all])
    for (i in seq_along(cols_all)) {
      if (identical(cols_all[i], paste0(cls, "_count"))) {
        kpArea(kp, chr = data$chrom, x = data$barycenter, y = data[[cols_all[i]]]/overall_max,
               col = colors_order[i], border = "NA", r0 = 0, r1 = 1)
      }
    }
    add_axis_counts(kp, overall_max)

    if (!is.null(data_kimura)) {
      columns <- paste0(cls, "_", gsub("-", ".", kimura_cats), "_pct_stacked")
      for (i in seq_along(columns)) {
        kpArea(kp, chr = data$chrom, x = data$barycenter, y = data_kimura[[columns[i]]],
               col = kimura_colors[i], border = "NA", r0 = 1, r1 = 0, data.panel = 2)
      }
      add_axis_pct(kp, panel = 2, r0 = 1, r1 = 0, labels = c("0%", "50%", "100%"))
    }

    add_legends()
    title(paste0("Number of insertions of ", cls, " along chromosomes - ", name),
          cex.main = 0.8, line = 2.5)
    dev.off()
  }
}
