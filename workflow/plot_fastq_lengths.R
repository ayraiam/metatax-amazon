# make the script accept infile, outdir, and label
args    <- commandArgs(trailingOnly = TRUE)
infile  <- if (length(args) >= 1) args[1] else "results/lengths/all_lengths.tsv"
outdir  <- if (length(args) >= 2) args[2] else dirname(infile)
label   <- if (length(args) >= 3) args[3] else ""

suffix  <- if (nzchar(label)) paste0("_", label) else ""
png_out <- file.path(outdir, paste0("read_length_boxplots", suffix, ".png"))
pdf_out <- file.path(outdir, paste0("read_length_boxplots", suffix, ".pdf"))

make_boxplot <- function(){
  
  suppressPackageStartupMessages({
    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package 'data.table' not found. Install it in the libsQC env.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' not found. Install it in the libsQC env.")
    }
  })
  
  ### >>> CHANGED: use arg-provided infile/outdir
  if (!file.exists(infile)) {
    stop(sprintf("Input not found: %s (did you run plot_fastq_length_boxplots to build it?)", infile))
  }
  
  dt <- read.table(infile, header = TRUE, sep = "\t")
  dt$length <- as.integer(dt$length)
  
  # If samples are like "Group/Sample", split; otherwise single group is fine
  if (any(grepl("/", dt$sample, fixed = TRUE))) {
    parts <- strsplit(as.character(dt$sample), "/", fixed = TRUE)
    dt$group  <- vapply(parts, `[`, character(1), 1)
    dt$sample <- vapply(parts, function(x) if (length(x) > 1) x[2] else x[1], character(1))
    dt$group  <- factor(dt$group)
  } else {
    dt$group <- factor("All")
  }
  
  meds <- aggregate(length ~ sample, dt, function(x) stats::median(x, na.rm = TRUE))
  ord  <- meds[order(meds$length), "sample"]
  dt$sample <- factor(dt$sample, levels = ord)
  
  y_cap <- unname(quantile(dt$length, 0.99, na.rm = TRUE))
  
  p <- ggplot2::ggplot(dt, ggplot2::aes(sample, length, group = sample)) +
    ggplot2::geom_violin(fill = "grey60", color = NA, alpha = 0.4, scale = "width", trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.25, outlier.shape = NA, linewidth = 0.3, fill = "white") +
    ggplot2::coord_cartesian(ylim = c(0, y_cap)) +
    ggplot2::scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    ggplot2::labs(title = "Read-length distribution per FASTQ (violin + box, outliers hidden)",
                  x = "FASTQ file (base name)", y = "Read length (bp)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1, size = 5),
                   panel.grid.minor = ggplot2::element_blank())
  
  # With per-group runs, there will typically be one facet only; harmless to keep
  if (length(levels(dt$group)) > 1) {
    p <- p + ggplot2::facet_wrap(~ group, scales = "free_x", ncol = 1)
  }
  
  # cap width to avoid the 50-inch error (and allow override)
  n <- length(levels(dt$sample))
  w <- min(40, max(8, n * 0.25))  # 8–40 inches
  ggplot2::ggsave(png_out, p, width = w, height = 6, dpi = 200, limitsize = FALSE)
  ggplot2::ggsave(pdf_out,  p, width = w, height = 6, limitsize = FALSE)
  
  message(sprintf("Wrote:\n  %s\n  %s", png_out, pdf_out))
}

make_boxplot()
