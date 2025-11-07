suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript plot_emu_genus_stacks.R <abundance_combined.tsv> <outdir>")
}
abund_file <- args[1]
outdir <- args[2]
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

dt <- fread(abund_file, sep = "\t", header = TRUE, na.strings = c("", "NA"))
setnames(dt, tolower(names(dt)))

# Expect columns at least: file, rank, taxon, abundance
req <- c("file","rank","taxon","abundance")
if (!all(req %in% names(dt))) {
  stop(sprintf("Missing required columns. Saw: %s\nNeed: %s",
               paste(names(dt), collapse=", "),
               paste(req, collapse=", ")))
}

# <<< CHANGED: early exit if no rows after read
if (nrow(dt) == 0) {
  message("No rows in abundance table; nothing to plot. Exiting quietly.")
  quit(save="no", status=0)
}

# Coerce abundance numeric
suppressWarnings(dt[, abundance := as.numeric(abundance)])

# ---- Genus extraction (unchanged) --------------------------------------
get_genus <- function(taxon, rank){
  if (is.na(taxon) || taxon == "") return(NA_character_)
  if (!is.na(rank) && tolower(rank) == "genus") return(taxon)
  x <- unlist(strsplit(taxon, ";", fixed = TRUE))
  x <- x[nchar(x) > 0]
  last <- if (length(x)) x[length(x)] else taxon
  last <- sub("^[a-z]__","", last, perl = TRUE)
  pieces <- strsplit(last, "\\s+")[[1]]
  if (length(pieces) >= 1) pieces[1] else last
}

dt[, genus := mapply(get_genus, taxon, rank)]
dt <- dt[!is.na(genus) & !is.na(abundance)]

# guard again if trimming removes all rows
if (nrow(dt) == 0) {
  message("No genus rows after filtering; nothing to plot. Exiting quietly.")
  quit(save="no", status=0)
}

# ---- Normalization per file (defensive) ---------------------------------
# robust total computation and NA protection
dt[, total := sum(abundance[is.finite(abundance)], na.rm=TRUE), by = file]
dt[, rel := {
  tot <- unique(total)
  if (!length(tot) || is.na(tot) || !is.finite(tot) || tot == 0) NA_real_
  else if (tot > 90 && tot < 110) abundance/100 else abundance/tot
}, by = file]
dt <- dt[is.finite(rel) & !is.na(rel) & rel >= 0]

# ---- Site facet from file name ------------------------------------------
# handle both basenames and paths
dt[, file_base := basename(file)]
dt[, file_base := sub("\\.(fastq|fq)(\\.gz)?$", "", file_base, ignore.case = TRUE)]
dt[, site := sub("(_I_|_II_|_III_).*", "", file_base, perl = TRUE)]
dt[, file_factor := factor(file_base, levels = unique(file_base))]

# ---- Keep top N genera per file, collapse others ------------------------
N <- 20
dt[, rank_in_file := frank(-rel, ties.method = "min"), by = file_base]
top_dt    <- dt[rank_in_file <= N, .(file_base, site, file_factor, genus, rel)]
other_dt  <- dt[rank_in_file >  N, .(rel = sum(rel)), by = .(file_base, site, file_factor)]
if (nrow(other_dt)) other_dt[, genus := "Other"]

plot_dt <- rbindlist(list(top_dt, other_dt), use.names = TRUE, fill = TRUE)

# Put Other at end
genus_levels <- unique(plot_dt[genus != "Other"][order(-rel), genus])
plot_dt[, genus := factor(genus, levels = c(genus_levels, "Other"))]

# ---- Facet in chunks of up to 20 sites ---------------------------------
sites  <- unique(plot_dt$site)
if (!length(sites)) {
  message("No sites to facet; nothing to plot.")
  quit(save="no", status=0)
}
chunks <- split(sites, ceiling(seq_along(sites) / 20))

plot_chunk <- function(sv, idx){
  d <- plot_dt[site %in% sv]
  if (!nrow(d)) return(invisible())
  p <- ggplot(d, aes(x = file_factor, y = rel, fill = genus)) +
    geom_col(width = 0.9) +
    facet_wrap(~ site, scales = "free_x") +
    labs(x = NULL, y = "Relative abundance (genus, %)", fill = "Genus") +
    scale_y_continuous(labels = function(z) 100*z,
                       expand = expansion(mult = c(0, 0.02))) +
    theme_bw(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      panel.grid.major.x = element_blank(),
      strip.background = element_rect(fill = "grey90"),
      legend.position = "bottom"
    )
  ggsave(file.path(outdir, sprintf("emu_genus_stacks_%02d.png", idx)), p, width = 16, height = 9, dpi = 200)
  ggsave(file.path(outdir, sprintf("emu_genus_stacks_%02d.pdf", idx)), p, width = 16, height = 9)
}

for (i in seq_along(chunks)) plot_chunk(chunks[[i]], i)
