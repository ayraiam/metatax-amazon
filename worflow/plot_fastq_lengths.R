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
  
  # Order samples by median length (base R)
  meds <- aggregate(length ~ sample, data = dt, FUN = function(x) stats::median(x, na.rm = TRUE))
  ord  <- meds[order(meds$length), "sample"]
  dt$sample <- factor(dt$sample, levels = ord)
  
  p <- ggplot2::ggplot(dt, ggplot2::aes(x = sample, y = length)) +
    ggplot2::geom_violin(fill = "grey60", color = NA, alpha = 0.4, scale = "width", trim = TRUE) +
    ggplot2::geom_boxplot(width = 0.25, outlier.alpha = 0.2, linewidth = 0.2, fill = "white") +
    ggplot2::labs(
      title = "Read-length distribution per FASTQ (one box per file)",
      x = "FASTQ file (base name)",
      y = "Read length (bp)"
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1),
      panel.grid.minor = ggplot2::element_blank()
    )
  
  # Width scales with number of boxes so ~70 boxes stay readable
  n <- length(levels(dt$sample))
  w <- max(12, n * 0.25)
  
  ggplot2::ggsave(png_out, p, width = w, height = 6, dpi = 200)
  ggplot2::ggsave(pdf_out, p, width = w, height = 6)
  
  message(sprintf("Wrote:\n  %s\n  %s", png_out, pdf_out))
}

# Calling the function
make_boxplot()
