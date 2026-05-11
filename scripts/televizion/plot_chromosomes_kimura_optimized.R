cat("\n### IMPORTING LIBRARIES ###\n\n")

library(karyoploteR)
library(optparse)

cat("\n### PARSING ARGUMENTS ###\n\n")

option_list <- list(
  make_option(c("-n", "--name"), type="character", default=NULL,
              help="Output prefix for generated files.", metavar="character"),
  make_option(c("-g", "--genome"), type="character", default=NULL,
              help="Genome table: [chr, start, end, name, gieStain].", metavar="character"),
  make_option(c("-c", "--chromosomes-order"), type="character", default=NULL,
              help="Chromosome order as comma-separated list.", metavar="character", dest="chromosomes"),
  make_option(c("-a", "--accessibility"), type="character", default="not_displayed",
              help="Accessibility table: [chr, start, end, name, itemRgb] or 'not_displayed'.", metavar="character"),
  make_option(c("-b", "--gc-content"), type="character", default="not_displayed",
              help="GC content table: [chr, start, end, name, itemRgb, gc_content] or 'not_displayed'.", metavar="character", dest="gccontent"),
  make_option(c("-i", "--classes-table"), type="character", default=NULL,
              help="Per-window class table: [chr, start, end, <areas>].", metavar="character", dest="input"),
  make_option(c("-k", "--classes-order"), type="character", default=NULL,
              help="Class order (comma-separated, reversed for display).", metavar="character", dest="classesorder"),
  make_option(c("-p", "--per-class"), type="character", default=NULL,
              help="Generate per-class plots (use 'True').", metavar="character", dest="perclass"),
  make_option(c("-l", "--colors-order"), type="character", default=NULL,
              help="Class color order (comma-separated, reversed for display).", metavar="character", dest="colorsorder"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output prefix path.", metavar="character"),
  make_option(c("--output-formats"), type="character", default="pdf",
              help="Comma-separated output figure formats: pdf,png,jpg. Default: pdf.", metavar="character", dest="output_formats"),
  make_option(c("--dpi"), type="integer", default=300,
              help="Raster output resolution for png/jpg. Must be >= 300. Default: 300.", metavar="integer"),
  make_option(c("-m", "--kimura-table"), type="character", default=NULL,
              help="Kimura table (RepeatMasker only).", metavar="character", dest="kimura"),
  make_option(c("-e", "--identity-table"), type="character", default=NULL,
              help="Identity table (EDTA or RepeatMasker without K2p).", metavar="character", dest="identity"),
  make_option(c("--layout"), type="character", default="horizontal",
              help="Layout mode: horizontal (plot.type=3) or vertical (plot.type=2).", metavar="character"),
  make_option(c("--zoom"), type="character", default=NULL,
              help="Optional zoom region as chr:start-end. Default: plot all selected chromosomes.", metavar="character")
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

parse_zoom_region <- function(region) {
  region <- trimws(region)
  match <- regexec("^([^:]+):([0-9,]+)-([0-9,]+)$", region)
  pieces <- regmatches(region, match)[[1]]
  if (length(pieces) != 4) {
    stop("Invalid --zoom value. Use chr:start-end, for example chr1:1000000-2000000.", call. = FALSE)
  }
  list(
    chromosome = pieces[2],
    start = as.numeric(gsub(",", "", pieces[3])),
    end = as.numeric(gsub(",", "", pieces[4]))
  )
}

build_zoom_request <- function(opt) {
  if (!is.null(opt$zoom)) {
    return(parse_zoom_region(opt$zoom))
  }
  NULL
}

validate_zoom_request <- function(zoom_request, genome) {
  if (is.null(zoom_request)) {
    return(NULL)
  }
  if (is.na(zoom_request$start) || is.na(zoom_request$end)) {
    stop("Zoom start and end must be numeric.", call. = FALSE)
  }
  zoom_request$start <- as.numeric(zoom_request$start)
  zoom_request$end <- as.numeric(zoom_request$end)
  if (zoom_request$start != floor(zoom_request$start) || zoom_request$end != floor(zoom_request$end)) {
    stop("Zoom start and end must be whole-number coordinates.", call. = FALSE)
  }
  if (zoom_request$start < 1 || zoom_request$end < zoom_request$start) {
    stop("Zoom coordinates must satisfy 1 <= start <= end.", call. = FALSE)
  }

  genome_chromosomes <- as.character(seqnames(genome))
  if (!zoom_request$chromosome %in% genome_chromosomes) {
    stop(paste0("Zoom chromosome '", zoom_request$chromosome, "' is not present in the genome table."), call. = FALSE)
  }

  chrom_idx <- which(genome_chromosomes == zoom_request$chromosome)[1]
  chrom_start <- start(genome)[chrom_idx]
  chrom_end <- end(genome)[chrom_idx]
  if (zoom_request$start < chrom_start || zoom_request$end > chrom_end) {
    stop(
      paste0(
        "Zoom coordinates for ", zoom_request$chromosome, " must be within ",
        chrom_start, "-", chrom_end, "."
      ),
      call. = FALSE
    )
  }
  zoom_request
}

subset_table_to_zoom <- function(dat, zoom_region) {
  if (is.null(dat) || is.null(zoom_region)) {
    return(dat)
  }
  keep <- dat$chrom == zoom_region$chromosome &
    dat$end >= zoom_region$start &
    dat$start <= zoom_region$end
  dat[keep, , drop = FALSE]
}

subset_granges_to_zoom <- function(gr, zoom_region) {
  if (is.null(gr) || is.null(zoom_region)) {
    return(gr)
  }
  keep <- as.character(seqnames(gr)) == zoom_region$chromosome &
    end(gr) >= zoom_region$start &
    start(gr) <= zoom_region$end
  gr[keep]
}

format_coord <- function(x) {
  format(x, scientific = FALSE, trim = TRUE)
}

sanitize_filename <- function(x) {
  gsub("[^A-Za-z0-9._-]+", "_", x)
}

parse_output_formats <- function(value) {
  formats <- tolower(trimws(unlist(strsplit(value, ","))))
  formats <- formats[formats != ""]
  allowed_formats <- c("pdf", "png", "jpg")
  invalid_formats <- setdiff(formats, allowed_formats)
  if (length(formats) == 0 || length(invalid_formats) > 0) {
    stop("Invalid --output-formats value. Use one or more of: pdf,png,jpg.", call. = FALSE)
  }
  unique(formats)
}

output_formats <- parse_output_formats(opt$output_formats)
output_dpi <- as.integer(opt$dpi)
if (is.na(output_dpi) || output_dpi < 300) {
  stop("--dpi must be an integer >= 300.", call. = FALSE)
}

# ---- Data
data <- read.csv(opt$input,  sep = "\t")
data_kimura <- if (!is.null(opt$kimura)) read.csv(opt$kimura, sep = "\t") else NULL
data_identity <- if (!is.null(opt$identity)) read.csv(opt$identity, sep = "\t") else NULL
name <- opt$name

# ---- Genome / chromosomes
genome <- toGRanges(opt$genome)
zoom_region <- validate_zoom_request(build_zoom_request(opt), genome)

chromosome_order <- if (is.character(opt$chromosomes)) {
  unlist(strsplit(opt$chromosomes, ","))
} else {
  genome@seqinfo@seqnames
}
if (!is.null(zoom_region)) {
  chromosome_order <- zoom_region$chromosome
  data <- subset_table_to_zoom(data, zoom_region)
  data_kimura <- subset_table_to_zoom(data_kimura, zoom_region)
  data_identity <- subset_table_to_zoom(data_identity, zoom_region)
  accessibility <- subset_granges_to_zoom(accessibility, zoom_region)
  gccontent <- subset_granges_to_zoom(gccontent, zoom_region)
  if (nrow(data) == 0) {
    stop("The zoom region contains no plotted windows. Use a wider interval or a smaller --windowsize.", call. = FALSE)
  }
}

zoom_granges <- if (!is.null(zoom_region)) {
  GRanges(
    seqnames = zoom_region$chromosome,
    ranges = IRanges(start = zoom_region$start, end = zoom_region$end)
  )
} else {
  NULL
}
zoom_label <- if (!is.null(zoom_region)) {
  paste0(zoom_region$chromosome, ":", format_coord(zoom_region$start), "-", format_coord(zoom_region$end))
} else {
  NULL
}
zoom_suffix <- if (!is.null(zoom_region)) {
  paste0(
    "_zoom_",
    sanitize_filename(zoom_region$chromosome),
    "_",
    format_coord(zoom_region$start),
    "_",
    format_coord(zoom_region$end)
  )
} else {
  ""
}

names_vec   <- mcols(genome)$name
df_names    <- data.frame(chr = as.character(seqnames(genome)), name = names_vec)
ordered_names <- unlist(lapply(chromosome_order, function(chr) df_names$name[df_names$chr == chr]))

layout_mode <- tolower(opt$layout)
if (!layout_mode %in% c("horizontal", "vertical")) {
  stop("Invalid --layout value. Use 'horizontal' or 'vertical'.", call. = FALSE)
}
plot_type <- if (layout_mode == "vertical") 2 else 3
n_chr <- length(chromosome_order)
plot_width <- 11.417
plot_height <- if (layout_mode == "vertical") 1.2 * n_chr + 2 else 3.937

plot_file <- function(stem, output_format) {
  paste0(opt$output, zoom_suffix, stem, ".", output_format)
}

open_plot_device <- function(file, output_format) {
  if (output_format == "pdf") {
    pdf(file, width = plot_width, height = plot_height)
  } else if (output_format == "png") {
    png(file, width = plot_width, height = plot_height, units = "in", res = output_dpi)
  } else if (output_format == "jpg") {
    jpeg(file, width = plot_width, height = plot_height, units = "in", res = output_dpi, quality = 95)
  }
}

write_plot <- function(stem, plotter) {
  for (output_format in output_formats) {
    file <- plot_file(stem, output_format)
    open_plot_device(file, output_format)
    plotter()
    dev.off()
  }
}

plot_title <- function(text) {
  if (is.null(zoom_label)) {
    paste0(text, " - ", name)
  } else {
    paste0(text, " (", zoom_label, ") - ", name)
  }
}

safe_max <- function(x) {
  value <- suppressWarnings(max(x, na.rm = TRUE))
  if (!is.finite(value) || value <= 0) {
    return(1)
  }
  value
}

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
  pp <- getDefaultPlotParams(plot.type = plot_type)
  pp$data1inmargin <- if (layout_mode == "vertical") 20 else 14
  pp$topmargin <- 100
  pp$leftmargin <- 0.06
  pp$rightmargin <- if (layout_mode == "vertical") 0.28 else 0.24
  pp$ideogramlateralmargin <- 0.01
  pp$data2inmargin <- if (layout_mode == "vertical") 32 else 20
  pp$data2height <- pp$data1height / 3
  if (layout_mode == "vertical") {
    pp$ideogramheight <- pp$data2height / 2
    pp$data2outmargin <- 40
  }
  pp
}

plot_karyo_base <- function() {
  pp <- make_pp()
  plot_args <- list(
    genome       = genome,
    chromosomes  = chromosome_order,
    plot.type    = plot_type,
    plot.params  = pp,
    labels.plotter = NULL,
    cytobands = genome
  )
  if (!is.null(zoom_granges)) {
    plot_args$zoom <- zoom_granges
  }
  kp <- do.call(plotKaryotype, plot_args)
  kpAddChromosomeNames(kp, chr.names = ordered_names, cex = 0.8)

  if (!is.null(accessibility)) {
    kpPlotRegions(karyoplot = kp, data = accessibility,
                  col = accessibility$itemRgb, data.panel = "ideogram")
  }
  if (!is.null(gccontent)) {
    kpPlotRegions(karyoplot = kp, data = gccontent,
                  col = gccontent$itemRgb, data.panel = "ideogram")
  }

  tick_dist <- if (is.null(zoom_region)) 10000000 else nice_step(zoom_region$end - zoom_region$start + 1)
  kpAddBaseNumbers(karyoplot = kp, tick.dist = tick_dist, tick.len = 4,
                   tick.col = "black", cex = 0.5)

  if (layout_mode == "vertical") {
    kpAddLabels(kp, labels = "Repeats", data.panel = 1, cex=0.6, side="left", srt=90, pos=3) # kpAddLabels(kp, labels = "Repeats", data.panel = 1, cex=0.5, side="right", srt=90, pos=3, label.margin = 0.05)
  } else {
    kpAddLabels(kp, labels = "Repeats", data.panel = 2, cex=0.8)
  }
  if (!is.null(data_kimura)) {
    if (layout_mode == "vertical") {
      kpAddLabels(kp, labels = "K2p", data.panel = 2, cex=0.6, side="left", srt=90, pos=3) # kpAddLabels(kp, labels = "K2p", data.panel = 2, cex=0.5, side="right", srt=90, pos=3, label.margin = 0.05)
    } else {
      kpAddLabels(kp, labels = "K2p", data.panel = 2, cex=0.8)
    }
  }
  if (!is.null(data_identity)) {
    if (layout_mode == "vertical") {
      kpAddLabels(kp, labels = "Identity", data.panel = 2, cex=0.6, side="left", srt=90, pos=3) # kpAddLabels(kp, labels = "Identity", data.panel = 2, cex=0.5, side="right", srt=90, pos=3, label.margin = 0.05)
    } else {
      kpAddLabels(kp, labels = "Identity", data.panel = 2, cex=0.8)
    }
  }

  kp
}

add_axis_pct <- function(kp, panel = 1, r0 = 0, r1 = 1, labels = c("0%", "25%", "50%", "75%", "100%")) {
  kpAxis(karyoplot = kp, data.panel = panel, r0 = r0, r1 = r1,
         side = "right", numticks = length(labels), labels = labels, cex = 0.5)
}

nice_step <- function(x, nticks = 5) {
  if (x <= 0) return(1)
  raw <- x / (nticks - 1)
  mag <- 10 ^ floor(log10(raw))
  nice <- c(1, 2, 5, 10)
  step <- nice[which.min(abs(raw / mag - nice))] * mag
  step
}

add_axis_counts <- function(kp, overall_max, nticks = 5) {
  step <- nice_step(overall_max, nticks)
  ticks <- seq(0, floor(overall_max / step) * step, by = step)
  kpAxis(karyoplot = kp, data.panel = 1, ymin = 0, ymax = overall_max,
         r0 = 0, r1 = 1, side = "right", numticks = length(ticks), tick.pos = ticks, labels = ticks, cex = 0.5)
}

draw_k2p_panel <- function(kp, dat, cols) {
  for (i in seq_along(cols)) {
    kpArea(kp, chr = dat$chrom, x = dat$barycenter, y = dat[[cols[i]]],
           col = kimura_colors[i], border = "NA", r0 = 1, r1 = 0, data.panel = 2)
  }
  add_axis_pct(kp, panel = 2, r0 = 1, r1 = 0, labels = c("0%", "50%", "100%"))
}

draw_identity_panel <- function(kp, dat, cols) {
  for (i in seq_along(cols)) {
    kpArea(kp, chr = dat$chrom, x = dat$barycenter, y = dat[[cols[i]]],
           col = identity_colors[i], border = "NA", r0 = 1, r1 = 0, data.panel = 2)
  }
  add_axis_pct(kp, panel = 2, r0 = 1, r1 = 0, labels = c("0%", "50%", "100%"))
}

add_gc_colorbar_horizontal <- function(
  legs_width = 1, top = 1, bottom = 1,
  cmap_min = 0.2, cmap_max = 0.8,
  n = 256,
  title = "",
  label_pos = seq(0.2, 0.8, by = 0.1)
) {
  if (layout_mode == "vertical") {
    left <- x_legs - 0.02
    right <- x_legs + (0.78 * legs_width)
    # print(left)
    # print(right)
    # print(top)
    # print(bottom)
  } else {
    left <- 0.8
    right <- 0.91
    bottom <- 0.3
    top <- 0.48
  }
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

# Legends 
x_legs <- 0.82
if (layout_mode == "vertical") {
  y_rep <- 1
} else {
  y_rep <- 1.2
}

add_legends <- function() {
  twidth <- max(strwidth(classes_order, units = "user"))
  rep_leg <- legend(x = x_legs, y = y_rep, legend = classes_order, fill = colors_order,
          border = "grey5", bty = "o", box.lwd = 0.3, title = "Repeat class",
          cex = 0.8, text.width = twidth, xpd = TRUE) 
  y_div <- rep_leg$rect$top - (1.05 * rep_leg$rect$h)
  legs_width <- rep_leg$rect$w
  if (!is.null(data_kimura)) {
    div_leg <- legend(x = x_legs, y = y_div, legend = kimura_cats, fill = kimura_colors,
           border = "grey5", bty = "o", box.lwd = 0.3, title = "K2p",
           cex = 0.6, text.width = twidth, xpd = TRUE)
  }
  if (!is.null(data_identity)) {
    div_leg <- legend(x = x_legs, y = y_div, legend = identity_cats, fill = identity_colors,
           border = "grey5", bty = "o", box.lwd = 0.3, title = "Identity",
           cex = 0.8, text.width = twidth, xpd = TRUE)
  }
  if (exists("div_leg", inherits = FALSE)) {
    y_gc <- div_leg$rect$top - (1.1 * div_leg$rect$h)
  } else {
    y_gc <- y_div
  }

  if (!is.null(gccontent)){
    tmp <- c(0, 1)
    gc_leg <- legend(x = 0.82, y = y_gc, legend = tmp, fill = "white",
           border = "white", bty = "o", box.lwd = 0.3, title = "GC content (%)", title.col = "grey5", text.col = "white",
           cex = 0.8, text.width = twidth, xpd = TRUE)
    add_gc_colorbar_horizontal(legs_width = legs_width, top = gc_leg$rect$top - 0.024, bottom = gc_leg$rect$top - gc_leg$rect$h - 0.024)
  }
}

# =========================================
# Percentage of repeats along the genome
# =========================================

# ---- All classes (stacked % by class + K2p)
write_plot("_karyoplot_stacked_percentage_by_class", function() {
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
  title(plot_title("Percentage of repeated content along chromosomes"),
        cex.main = 0.8, line = 2.5)
})

# ---- Per class (% + K2p), if requested
if (!is.null(opt$perclass)) {
  cols_all_pct <- paste0(classes_order, "_pct")
  for (cls in classes_order) {
    write_plot(paste0("_karyoplot_percentage_", cls), function() {
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
      title(plot_title(paste0("Percentage of ", cls, " content along chromosomes")),
            cex.main = 0.8, line = 2.5)
    })
  }
}

# =========================================
# Number of insertions along the genome
# =========================================

# ---- All classes (stacked counts + K2p)
write_plot("_karyoplot_stacked_counts_by_class", function() {
  kp <- plot_karyo_base()

  columns_cnt <- paste0(classes_order, "_count_stacked")
  overall_max <- safe_max(data[, columns_cnt])
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
  title(plot_title("Number of insertions along chromosomes"),
        cex.main = 0.8, line = 2.5)
})

# ---- Per class (counts + K2p), if requested
if (!is.null(opt$perclass)) {
  cols_all_cnt <- paste0(classes_order, "_count")
  overall_max <- safe_max(data[, cols_all_cnt])
  for (cls in classes_order) {
    write_plot(paste0("_karyoplot_counts_", cls), function() {
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
      title(plot_title(paste0("Number of insertions of ", cls, " along chromosomes")),
            cex.main = 0.8, line = 2.5)
    })
  }
}
