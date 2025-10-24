make_boxplot <- function(){
  
  # workflow/plot_fastq_lengths.R
  suppressPackageStartupMessages({
    if (!requireNamespace("data.table", quietly = TRUE)) {
      stop("Package 'data.table' not found. Install it in the libsQC env.")
    }
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
      stop("Package 'ggplot2' not found. Install it in the libsQC env.")
    }
  })
  
  infile  <- "results/lengths/all_lengths.tsv"
  png_out <- "results/lengths/read_length_boxplots.png"
  pdf_out <- "results/lengths/read_length_boxplots.pdf"
  
  if (!file.exists(infile)) {
    stop(sprintf("Input not found: %s (did you run plot_fastq_length_boxplots to build it?)", infile))
  }
  
  # Expect two columns: sample, length
  dt <- read.table(infile, header = TRUE, sep = "\t")
  dt$length <- as.integer(dt$length)
  
  # reorder by median (optional)
  meds <- aggregate(length ~ sample, dt, function(x) stats::median(x, na.rm = TRUE))
  ord  <- meds[order(meds$length), "sample"]
  dt$sample <- factor(dt$sample, levels = ord)
  
  # choose a zoom cap (tweak 0.99 → 0.995 or 0.98 as you like)
  y_cap <- unname(quantile(dt$length, 0.99, na.rm = TRUE))
  
  p <- ggplot2::ggplot(dt, ggplot2::aes(sample, length, group = sample)) +
    ggplot2::geom_violin(fill = "grey60", color = NA, alpha = 0.4, scale = "width", trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.25, outlier.shape = NA, linewidth = 0.3, fill = "white") +  # hide outlier dots
    ggplot2::coord_cartesian(ylim = c(0, y_cap)) +                                               # zoom to central range
    ggplot2::scale_y_continuous(labels = function(x) format(x, big.mark = ",", scientific = FALSE)) +
    ggplot2::labs(title = "Read-length distribution per FASTQ (violin + box, outliers hidden)",
         x = "FASTQ file (base name)", y = "Read length (bp)") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
          panel.grid.minor = ggplot2::element_blank())
  
  n <- length(levels(dt$sample))
  w <- max(12, n * 0.25)
  ggplot2::ggsave("results/lengths/read_length_boxplots.png", p, width = w, height = 6, dpi = 200)
  ggplot2::ggsave("results/lengths/read_length_boxplots.pdf",  p, width = w, height = 6)
  
  message(sprintf("Wrote:\n  %s\n  %s", png_out, pdf_out))
}

# Calling the function
make_boxplot()
