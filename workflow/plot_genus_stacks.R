suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(RColorBrewer)
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

# --- Standardize columns -------------------------------------------------
has_genus <- "genus" %in% names(dt)
if (!has_genus && !"taxon" %in% names(dt)) {
  for (alt in c("name","taxonomy","clade_name","lineage","species")) {
    if (alt %in% names(dt)) { setnames(dt, alt, "taxon"); break }
  }
}
if (!"file" %in% names(dt) || (!has_genus && !"taxon" %in% names(dt))) {
  stop(sprintf(
    "Missing required columns. Saw: %s\nNeed: file AND (genus OR taxon/lineage alias).",
    paste(names(dt), collapse=", ")
  ))
}

if (!"abundance" %in% names(dt)) dt[, abundance := NA_real_]
suppressWarnings(dt[, abundance := as.numeric(abundance)])
if (("count" %in% names(dt)) && (!any(is.finite(dt$abundance)) || sum(dt$abundance, na.rm = TRUE) == 0)) {
  suppressWarnings(dt[, abundance := as.numeric(count)])
}

# ---- Genus choice: prefer provided column, fallback to parse ------------
if (has_genus) {
  dt[, genus := as.character(genus)]
} else {
  if (!"rank" %in% names(dt)) dt[, rank := NA_character_]
  get_genus <- function(taxon, rank){
    if (is.na(taxon) || taxon == "") return(NA_character_)
    if (!missing(rank) && !is.null(rank) && !is.na(rank) && tolower(rank) == "genus") return(taxon)
    x <- unlist(strsplit(taxon, ";", fixed = TRUE))
    x <- x[nchar(x) > 0]
    last <- if (length(x)) x[length(x)] else taxon
    last <- sub("^[a-z]__","", last, perl = TRUE)
    pieces <- strsplit(last, "\\s+")[[1]]
    if (length(pieces) >= 1) pieces[1] else last
  }
  dt[, genus := mapply(get_genus, taxon, rank)]
}

# Map missing/blank genus to "No genus"; keep only finite abundances
dt[is.na(genus) | !nzchar(genus), genus := "No genus"]
dt <- dt[is.finite(abundance)]

# ---- Collapse to one row per file × genus -------------------------------
dt <- dt[, .(abundance = sum(as.numeric(abundance), na.rm = TRUE)), by = .(file, genus)]

# ---- Normalization per file (abundance already 0–1) ---------------------
dt[, total := sum(abundance, na.rm = TRUE), by = file]
dt[, rel := fifelse(total > 0.99 & total < 1.01,
                    abundance,
                    fifelse(total > 0, abundance/total, NA_real_)),
   by = file]
dt <- dt[is.finite(rel) & !is.na(rel) & rel >= 0]

if (nrow(dt) == 0L) {
  message("No genus rows after filtering; nothing to plot. Exiting quietly.")
  quit(save = "no", status = 0)
}

# ---- Facets and file_base ----------------------------------------------
file_base <- basename(dt$file)
file_base <- sub("\\.(fastq|fq)(\\.gz)?$", "", file_base, ignore.case = TRUE)
dt[, file_base := file_base]

# Site = prefix before _I_/_II_/_III_
dt[, site := sub("(_I_|_II_|_III_).*", "", file_base, perl = TRUE)]

# Replicate symbol extracted as I / II / III
dt[, rep := fifelse(grepl("(_I_|_II_|_III_)", file_base),
                    sub(".*_(I{1,3})_.*", "\\1", file_base),
                    NA_character_)]

# Build replicate labels within each site (Ia, Ib for duplicates)
lab_map_dt <- unique(dt[, .(file_base, site, rep)])
lab_map_dt[, dupN := .N, by = .(site, rep)]
lab_map_dt[dupN == 1, rep_label := rep]
lab_map_dt[dupN >  1, rep_label := paste0(rep, letters[seq_len(.N)]), by = .(site, rep)]

# Order within site by replicate then letter
lab_map_dt[, rep_rank := match(rep, c("I","II","III"))]
setorder(lab_map_dt, site, rep_rank, rep_label)

# Map: file_base -> x-axis label
lab_map <- setNames(lab_map_dt$rep_label, lab_map_dt$file_base)

# x factor levels (global order)
order_levels <- lab_map_dt[, file_base]
dt[, file_factor := factor(file_base, levels = order_levels)]

# ---- Strict Top-N + 'Other' (deterministic tie-break by genus) ----------
N <- 20
setorder(dt, file_base, -rel, genus)   
dt[, idx := seq_len(.N), by = file_base]

top_dt   <- dt[idx <= N, .(file_base, site, file_factor, genus, rel)]
other_dt <- dt[idx >  N, .(rel = sum(rel)), by = .(file_base, site, file_factor)]
if (nrow(other_dt)) other_dt[, genus := "Other"]

plot_dt <- rbindlist(list(top_dt, other_dt), use.names = TRUE, fill = TRUE)

# Genus levels: keep "Other" last
genus_levels <- unique(plot_dt[genus != "Other"][order(-rel), genus])
plot_dt[, genus := factor(genus, levels = c(genus_levels, "Other"))]

# ---- Color palette selection --------------------------------------------
make_diverging_palette <- function(n_colors, has_other = FALSE) {
  palette_names <- c("Spectral","RdYlGn","RdYlBu","RdGy","RdBu","PuOr","PRGn","PiYG","BrBG")
  info <- RColorBrewer::brewer.pal.info
  cols <- character(0)
  needed <- n_colors
  for (p in palette_names) {
    if (!(p %in% rownames(info))) next
    maxn <- info[p, "maxcolors"]
    take <- min(maxn, max(0, needed))
    if (take > 0) {
      cols <- c(cols, RColorBrewer::brewer.pal(take, p)[seq_len(take)])
      needed <- n_colors - length(cols)
      if (needed <= 0) break
    }
  }
  if (length(cols) < n_colors) cols <- rep(cols, length.out = n_colors)
  cols[seq_len(n_colors)]
}

lev <- levels(plot_dt$genus)
has_other <- length(lev) > 0 && tail(lev, 1) == "Other"
k <- length(lev) - ifelse(has_other, 1L, 0L)
base_cols <- if (k > 0) make_diverging_palette(k) else character(0)
fill_vals <- if (has_other) c(base_cols, "grey80") else base_cols
names(fill_vals) <- lev

# ---- Plotting ------------------------------------------------------------
sites  <- unique(plot_dt$site)
chunks <- split(sites, ceiling(seq_along(sites) / 20))

plot_chunk <- function(sites_vec, idx){
  d <- plot_dt[site %in% sites_vec]
  if (!nrow(d)) return(invisible())
  p <- ggplot(d, aes(x = file_factor, y = rel, fill = genus)) +
    geom_col(width = 0.7) +  
    facet_wrap(~ site, scales = "free_x") +
    labs(x = NULL, y = "Relative abundance (genus, %)", fill = "Genus") +
    scale_y_continuous(labels = function(z) 100*z,
                       expand = expansion(mult = c(0, 0.02))) +
    scale_x_discrete(labels = lab_map) +
    scale_fill_manual(values = fill_vals, drop = FALSE) +
    theme_classic(base_size = 10) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      strip.background = element_rect(fill = "white"),
      strip.text = element_text(size = 12, face = "bold"),
      legend.position = "right",                      
      legend.text = element_text(size = 8),           
      legend.title = element_text(size = 9)           
    )
  ggsave(file.path(outdir, sprintf("emu_genus_stacks_%02d.png", idx)), p, width = 16, height = 9, dpi = 200)
  ggsave(file.path(outdir, sprintf("emu_genus_stacks_%02d.pdf", idx)), p, width = 16, height = 9)
}

for (i in seq_along(chunks)) plot_chunk(chunks[[i]], i)
