cat("\n### IMPORTING LIBRARIES ###\n\n")

library(karyoploteR)
library(optparse)

cat("\n### PARSING ARGUMENTS ###\n\n")

option_list <- list(
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="Output prefix for generated files.", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="Genome table: [chr, start, end, name, gieStain].", metavar="character"),
  make_option(c("-c", "--chromosomes-order", "--chromosomes"), type="character", default=NULL,
              help="Chromosome order as comma-separated list.", metavar="character", dest="chromosomes"),
  make_option(c("--accessibility"), type="character", default="not_displayed",
              help="Accessibility table: [chr, start, end, name, itemRgb] or 'not_displayed'.", metavar="character"),
  make_option(c("--gc-content", "--gccontent"), type="character", default="not_displayed",
              help="GC content table: [chr, start, end, name, itemRgb, gc_content] or 'not_displayed'.", metavar="character", dest="gccontent"),
  make_option(c("-i", "--classes-table", "--input"), type="character", default=NULL,
              help="Per-window class table: [chr, start, end, <areas>].", metavar="character", dest="input"),
  make_option(c("-k", "--classes-order", "--classesorder"), type="character", default=NULL,
              help="Class order (comma-separated, reversed for display).", metavar="character", dest="classesorder"),
  make_option(c("--per-class", "--perclass"), type="character", default=NULL,
              help="Generate per-class plots (use 'True').", metavar="character", dest="perclass"),
  make_option(c("-l", "--colors-order", "--colorsorder"), type="character", default=NULL,
              help="Class color order (comma-separated, reversed for display).", metavar="character", dest="colorsorder"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output prefix path.", metavar="character"),
  make_option(c("--kimura-table", "--kimura"), type="character", default=NULL,
              help="Kimura table (RepeatMasker only).", metavar="character", dest="kimura"),
  make_option(c("--identity-table", "--identity"), type="character", default=NULL,
              help="Identity table (EDTA or RepeatMasker without K2p).", metavar="character", dest="identity")
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
data <- read.csv(opt$input,  sep = "\t")
data_kimura <- if (!is.null(opt$kimura)) read.csv(opt$kimura, sep = "\t") else NULL
data_identity <- if (!is.null(opt$identity)) read.csv(opt$identity, sep = "\t") else NULL
name <- opt$name

# ---- Genome / chromosomes
genome <- toGRanges(opt$genome)

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
kimura_cols_all <- paste0("ALL_", gsub("-", ".", kimura_cats), "_pct_stacked")

# ---- Identity
identity_cats   <- c("0.6-0", "0.7-0.6", "0.8-0.7", "0.9-0.8", "1-0.9")
identity_colors <- c("#2c699a", "#0db39e", "#83e377", "#efea5a", "#f29e4c")
identity_cols_all <- paste0("ALL_", gsub("-", ".", identity_cats), "_pct_stacked")

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

draw_k2p_panel <- function(kp, dat, cols) {
  for (i in seq_along(cols)) {
    kpArea(kp, chr = data$chrom, x = data$barycenter, y = dat[[cols[i]]],
           col = kimura_colors[i], border = "NA", r0 = 1, r1 = 0, data.panel = 2)
  }
  add_axis_pct(kp, panel = 2, r0 = 1, r1 = 0, labels = c("0%", "50%", "100%"))
}

draw_identity_panel <- function(kp, dat, cols) {
  for (i in seq_along(cols)) {
    kpArea(kp, chr = data$chrom, x = data$barycenter, y = dat[[cols[i]]],
           col = identity_colors[i], border = "NA", r0 = 1, r1 = 0, data.panel = 2)
  }
  add_axis_pct(kp, panel = 2, r0 = 1, r1 = 0, labels = c("0%", "50%", "100%"))
}

add_gc_colorbar_horizontal <- function(
  left = 0.8, right = 0.91, bottom = 0.3, top = 0.48,
  cmap_min = 0.2, cmap_max = 0.8,
  n = 256,
  title = "",
  label_pos = seq(0.2, 0.8, by = 0.1)
) {
  op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
  par(new = TRUE, fig = c(left, right, bottom, top), mar = c(2, 1, 1, 1), mgp = c(1.8, 0.2, 0))
  x <- seq(cmap_min, cmap_max, length.out = n)
  y <- c(0, 1)
  cols <- rev(hcl.colors(n, "Spectral"))
  z <- cbind(x, x)
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
pdf(file, width = 11.417, height = 3.937)

kp <- plot_karyo_base()

columns_pct <- paste0(classes_order, "_pct_stacked")
for (i in seq_along(columns_pct)) {
  kpArea(kp, chr = data$chrom, x = data$barycenter, y = data[[columns_pct[i]]],
         col = colors_order[i], border = "NA", r0 = 0, r1 = 1, data.panel = 1)
}
add_axis_pct(kp, panel = 1)

if (!is.null(data_kimura)) {
  draw_k2p_panel(kp, data_kimura, kimura_cols_all)
}

if (!is.null(data_identity)) {
  draw_identity_panel(kp, data_identity, identity_cols_all)
}

add_legends()
title(paste0("Percentage of repeated content along chromosomes - ", name),
      cex.main = 0.8, line = 2.5)

dev.off()

# ---- Per class (% + K2p), if requested
if (!is.null(opt$perclass)) {
  cols_all_pct <- paste0(classes_order, "_pct")
  for (cls in classes_order) {
    file <- paste0(opt$output, "_karyoplot_percentage_", cls, ".pdf")
    pdf(file, width = 11.417, height = 3.937)

    kp <- plot_karyo_base()

    for (i in seq_along(cols_all_pct)) {
      if (identical(cols_all_pct[i], paste0(cls, "_pct"))) {
        kpArea(kp, chr = data$chrom, x = data$barycenter, y = data[[cols_all_pct[i]]],
               col = colors_order[i], border = "NA", r0 = 0, r1 = 1)
      }
    }
    add_axis_pct(kp, panel = 1)

    if (!is.null(data_kimura)) {
      columns <- paste0(cls, "_", gsub("-", ".", kimura_cats), "_pct_stacked")
      draw_k2p_panel(kp, data_kimura, columns)
    }

    if (!is.null(data_identity)) {
      columns <- paste0(cls, "_", gsub("-", ".", identity_cats), "_pct_stacked")
      draw_identity_panel(kp, data_identity, columns)
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
pdf(file, width = 11.417, height = 3.937)

kp <- plot_karyo_base()

columns_cnt <- paste0(classes_order, "_count_stacked")
overall_max <- max(data[, columns_cnt])
for (i in seq_along(columns_cnt)) {
  kpArea(kp, chr = data$chrom, x = data$barycenter, y = data[[columns_cnt[i]]]/overall_max,
         col = colors_order[i], border = "NA", r0 = 0, r1 = 1, data.panel = 1)
}
add_axis_counts(kp, overall_max)

if (!is.null(data_kimura)) {
  draw_k2p_panel(kp, data_kimura, kimura_cols_all)
}

if (!is.null(data_identity)) {
  draw_identity_panel(kp, data_identity, identity_cols_all)
}

add_legends()
title(paste0("Number of insertions along chromosomes - ", name),
      cex.main = 0.8, line = 2.5)

dev.off()

# ---- Per class (counts + K2p), if requested
if (!is.null(opt$perclass)) {
  cols_all_cnt <- paste0(classes_order, "_count")
  overall_max <- max(data[, cols_all_cnt])
  for (cls in classes_order) {
    file <- paste0(opt$output, "_karyoplot_counts_", cls, ".pdf")
    pdf(file, width = 11.417, height = 3.937)

    kp <- plot_karyo_base()

    for (i in seq_along(cols_all_cnt)) {
      if (identical(cols_all_cnt[i], paste0(cls, "_count"))) {
        kpArea(kp, chr = data$chrom, x = data$barycenter, y = data[[cols_all_cnt[i]]]/overall_max,
               col = colors_order[i], border = "NA", r0 = 0, r1 = 1)
      }
    }
    add_axis_counts(kp, overall_max)

    if (!is.null(data_kimura)) {
      columns <- paste0(cls, "_", gsub("-", ".", kimura_cats), "_pct_stacked")
      draw_k2p_panel(kp, data_kimura, columns)
    }

    add_legends()
    title(paste0("Number of insertions of ", cls, " along chromosomes - ", name),
          cex.main = 0.8, line = 2.5)
    dev.off()
  }
}
