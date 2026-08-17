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
  make_option(c("-b", "--gccontent"), type="character", default="not_displayed",
              help="genome gc content file name, tabulated: [chr, start, end, name, itemRgb, gc_content]", metavar="character"),
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
  make_option(c("-e", "--identity"), type="character", default=NULL,
              help="Additional panel [EDTA and RepeatMasker without K2p]", metavar="character")
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
if (!is.null(opt$gccontent)) {
  if (is.character(opt$gccontent) && opt$gccontent %in% c("not_displayed")) {
    gccontent <- NULL   
  } else {
    gccontent <- toGRanges(opt$gccontent)
  }
}

accessibility <- NULL
if (!identical(opt$accessibility, "not_displayed")) {
  accessibility <- toGRanges(opt$accessibility)
}
if (!is.null(opt$kimura)) {
  if (is.character(opt$kimura) && opt$kimura %in% c("False", "None")) {
    opt$kimura <- NULL   
  }
}
if (!is.null(opt$identity)) {
  if (is.character(opt$identity) && opt$identity %in% c("False", "None")) {
    opt$identity <- NULL  
  }
}

cat("\n### CODE ###\n\n")

# ---- Data
data        <- read.csv(opt$input,  sep = "\t")
data_kimura <- if (!is.null(opt$kimura)) read.csv(opt$kimura, sep = "\t") else NULL
data_identity <- if (!is.null(opt$identity)) read.csv(opt$identity, sep = "\t") else NULL
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

# ---- Kimura
kimura_cats   <- c("40-70", "30-40", "20-30", "10-20", "0-10")
kimura_colors <- c("#2c699a", "#0db39e", "#83e377", "#efea5a", "#f29e4c")

# ---- Identity
identity_cats   <- c("0.6-0", "0.7-0.6", "0.8-0.7", "0.9-0.8", "1-0.9")
identity_colors <- c("#2c699a", "#0db39e", "#83e377", "#efea5a", "#f29e4c") 

# ---------------------------
# Helpers to avoid repetition
# ---------------------------

make_pp <- function() {
  pp <- getDefaultPlotParams(plot.type = 3)
  pp$data1inmargin <- 14
  pp$topmargin <- 100
  pp$leftmargin <- 0.06
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
  if (!is.null(gccontent)) {
    kpPlotRegions(karyoplot = kp, data = gccontent,
                  col = gccontent$itemRgb, data.panel = "ideogram")
  }

  kpAddBaseNumbers(karyoplot = kp, tick.dist = 10000000, tick.len = 3,
                   tick.col = "black", cex = 0.5)

  # kpDataBackground(kp, data.panel = 1, col = "grey98")  
  # kpDataBackground(kp, data.panel = 2, col = "grey98")  # , col = "grey96"

  kpAddLabels(kp, labels = "Repeats", data.panel = 1, cex=0.8)
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

# add_gc_colorbar <- function(left=0.88, right=0.98, bottom=0.15, top=0.90,
#                             at=seq(0.2, 0.8, by=0.1),
#                             title="GC content (%)") {
#   op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)

#   par(new = TRUE, fig = c(left, right, bottom, top), mar = c(2, 1, 2, 4))

#   vals <- seq(0.2, 0.8, length.out = 256)
#   image(x = 1,
#         y = vals,
#         z = t(as.matrix(vals)),
#         col = rev(hcl.colors(256, "Spectral")),
#         axes = FALSE, xlab = "", ylab = "")

#   axis(4, at = at, labels = sprintf("%f", at*100), las = 1)
#   mtext(title, side = 4, line = 2.5, cex = 0.9)
#   box()
# }

add_gc_colorbar_horizontal <- function(
  left = 0.8, right = 0.91, bottom = 0.3, top = 0.48,
  cmap_min = 0.2, cmap_max = 0.8,
  n = 256,
  title = "",
  label_pos = seq(0.2, 0.8, by = 0.1)
) {
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)

  # panel for the colorbar
  par(new = TRUE, fig = c(left, right, bottom, top), mar = c(2, 1, 1, 1), mgp = c(1.8, 0.2, 0))

  # gradient domain + colors
  x <- seq(cmap_min, cmap_max, length.out = n)   # horizontal axis
  y <- c(0, 1)                                   # two-row strip
  cols <- rev(hcl.colors(n, "Spectral"))

  # z must be nrow = length(x), ncol = length(y)
  z <- cbind(x, x)   # n x 2 matrix (duplicates the gradient)

  image(
    x = x, y = y, z = z,
    col = cols,
    axes = FALSE,
    xlab = title, ylab = "",
    useRaster = TRUE
  )

  axis(1, at = label_pos, labels = sprintf("%.0f", label_pos * 100), cex.axis = 0.8, tck = -0.5)
  box()
}


add_legends <- function() {
  twidth <- max(strwidth(classes_order, units = "user"))
  legend(x = 0.82, y = 1.2, legend = classes_order, fill = colors_order,
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
  if (!is.null(gccontent)){
    tmp <- c(0, 1)
    legend(x = 0.82, y = 0.45, legend = tmp, fill = "white",
      border = "white", bty = "o", box.lwd = 0.3, title = "GC content (%)", title.col = "grey5", text.col = "white",
      cex = 0.8, text.width = twidth, xpd = TRUE)
    add_gc_colorbar_horizontal()
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
