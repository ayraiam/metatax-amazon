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

# --- Standardize columns -------------------------------------------------
# We expect: file, taxon (or alias), and either abundance or count
if (!"taxon" %in% names(dt)) {                                    
  for (alt in c("name","taxonomy","clade_name","lineage","species")) {  
    if (alt %in% names(dt)) { setnames(dt, alt, "taxon"); break }       
  }                                                                     
}                                                                       
if (!all(c("file","taxon") %in% names(dt))) {                           
  stop(sprintf("Missing required columns. Saw: %s\nNeed at least: file and taxon (or aliases)",  
               paste(names(dt), collapse=", ")))                        
}                                                                       

# Ensure an abundance column exists (fallback to count when needed)       
if (!"abundance" %in% names(dt)) dt[, abundance := NA_real_]            
suppressWarnings(dt[, abundance := as.numeric(abundance)])              
if (("count" %in% names(dt)) && (!any(is.finite(dt$abundance)) || sum(dt$abundance, na.rm = TRUE) == 0)) {  
  suppressWarnings(dt[, abundance := as.numeric(count)])                
}                                                                       

# ---- Genus extraction ----------------------------------------------------
# If rank == genus use 'taxon' as-is. Otherwise try best-effort parse.
get_genus <- function(taxon, rank){
  if (is.na(taxon) || taxon == "") return(NA_character_)
  if (!is.null(rank) && !is.na(rank) && tolower(rank) == "genus") return(taxon)
  
  # Split hierarchical strings like "k__;p__;...;g__Escherichia;s__coli" OR "Bacteria;Proteobacteria;...;Escherichia coli"
  x <- unlist(strsplit(taxon, ";", fixed = TRUE))
  x <- x[nchar(x) > 0]
  last <- if (length(x)) x[length(x)] else taxon
  
  # Strip common prefixes like "g__"
  last <- sub("^[a-z]__","", last, perl = TRUE)
  
  # If it looks like "Escherichia coli", take first token
  pieces <- strsplit(last, "\\s+")[[1]]
  if (length(pieces) >= 1) pieces[1] else last
}

# rank may be missing entirely; handle safely
if (!"rank" %in% names(dt)) dt[, rank := NA_character_]               
dt[, genus := mapply(get_genus, taxon, rank)]
dt <- dt[!is.na(genus) & !is.na(abundance)]

# ---- Normalization per file ---------------------------------------------
# Convert to proportions per file (0-1). Heuristics:
# - If total ~100, treat as percent.
# - Else divide by sum (counts).
dt[, total := sum(abundance, na.rm=TRUE), by = file]
dt[, rel := {
  if (!is.finite(total) || total == 0) NA_real_
  else if (total > 90 && total < 110) abundance/100
  else abundance/total
}, by = file]
dt <- dt[!is.na(rel)]

if (nrow(dt) == 0L) {
  message("No genus rows after filtering; nothing to plot. Exiting quietly.")
  quit(save = "no", status = 0)
}

# ---- Site facet from file name ------------------------------------------
file_base <- basename(dt$file)
# drop extensions if present
file_base <- sub("\\.(fastq|fq)(\\.gz)?$", "", file_base, ignore.case = TRUE)
dt[, file_base := file_base]

# Everything before _I_/_II_/_III_
dt[, site := sub("(_I_|_II_|_III_).*", "", file_base, perl = TRUE)]

# Preserve original file order on x-axis
dt[, file_factor := factor(file_base, levels = unique(file_base))]

# ---- Keep top N genera per file, collapse others into a single "Other" ---
N <- 20
# Rank within file by relative abundance
dt[, rank_in_file := frank(-rel, ties.method = "min"), by = file_base]

top_dt    <- dt[rank_in_file <= N, .(file_base, site, file_factor, genus, rel)]
other_dt  <- dt[rank_in_file >  N, .(rel = sum(rel)), by = .(file_base, site, file_factor)]
if (nrow(other_dt)) other_dt[, genus := "Other"]

plot_dt <- rbindlist(list(
  top_dt,
  other_dt
), use.names = TRUE, fill = TRUE)

# Ensure genus is a factor with "Other" at end
genus_levels <- unique(plot_dt[genus != "Other"][order(-rel), genus])
plot_dt[, genus := factor(genus, levels = c(genus_levels, "Other"))]

# ---- Facet in chunks of up to 20 sites ----------------------------------
sites  <- unique(plot_dt$site)
chunks <- split(sites, ceiling(seq_along(sites) / 20))

plot_chunk <- function(sites_vec, idx){
  d <- plot_dt[site %in% sites_vec]
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
  
  ggsave(file.path(outdir, sprintf("emu_genus_stacks_%02d.png", idx)), p,
         width = 16, height = 9, dpi = 200)
  ggsave(file.path(outdir, sprintf("emu_genus_stacks_%02d.pdf", idx)), p,
         width = 16, height = 9)
}

for (i in seq_along(chunks)) plot_chunk(chunks[[i]], i)
