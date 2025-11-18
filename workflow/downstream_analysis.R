#!/usr/bin/env Rscript
# ==========================================================
# downstream_analysis.R
# 1) Add 'environment' from filename rules
# 2) Make environment-grouped stacked bars (Family & Genus)
# 3) Alpha diversity (Shannon / Simpson) + pairwise Wilcoxon
# 4) Beta diversity (Bray–Curtis PCoA, Genus)
# ==========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(ggplot2)
  library(vegan)
  library(tidyr)
  library(dplyr)
  library(cowplot)
  library(RColorBrewer)
  library(ggbeeswarm)
})

# ---------- Input arguments ----------
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript downstream_analysis.R <input.tsv> <output_dir> [prefix]")
infile  <- args[1]
outdir  <- args[2]
prefix  <- ifelse(length(args) >= 3, args[3], "downstream")
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

# ==========================================================
# 0) If infile does NOT exist, merge batch tables_b*/abundance_combined.tsv
# ==========================================================
if (!file.exists(infile)) {
  message(">>> Input file not found: ", infile)
  message(">>> Trying to merge batch abundance_combined.tsv files from results/tables_b*/ ...")
  
  # Try to infer the 'results' dir from the expected infile path:
  # e.g. .../results/tables/abundance_combined.tsv -> .../results
  results_dir <- dirname(dirname(infile))
  if (!dir.exists(results_dir)) {
    # fallback: use current working directory + 'results'
    results_dir <- file.path(getwd(), "results")
  }
  message(">>> Searching for batch tables under: ", results_dir)
  
  # 1) find all abundance_combined.tsv files recursively
  all_abund <- list.files(
    path       = results_dir,
    pattern    = "^abundance_combined\\.tsv$",
    recursive  = TRUE,
    full.names = TRUE
  )
  # 2) keep only those inside a tables_b* directory
  batch_files <- all_abund[grepl("/tables_b", all_abund)]
  
  if (length(batch_files) == 0L) {
    stop("Could not find infile or any batch results/tables_b*/abundance_combined.tsv")
  }
  
  message(">>> Found ", length(batch_files), " batch abundance tables. Merging...")
  dt_list <- lapply(
    batch_files,
    function(f) data.table::fread(f, sep = "\t", header = TRUE, na.strings = c("", "NA"))
  )
  merged_dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  data.table::setnames(merged_dt, tolower(names(merged_dt)))
  
  # Ensure target directory exists (e.g. results/tables/)
  dir.create(dirname(infile), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(merged_dt, file = infile, sep = "\t")
  message(">>> Merged table written to: ", infile)
}

# ---------- Helper: add environment column ----------
add_environment <- function(dt) {
  stopifnot("file" %in% names(dt))
  
  # Normalize possible typos like "LO1" -> "L01" or "lO2" -> "L02"
  dt[, file_norm := toupper(file)]
  dt[, file_norm := gsub("LO([0-9])", "L0\\1", file_norm)]
  dt[, file_norm := gsub("IO([0-9])", "L0\\1", file_norm)]
  dt[, file_norm := gsub("O([0-9])", "0\\1", file_norm)]  # generic fix
  
  dt[, environment := NA_character_]
  
  rules <- c(
    "^CAMP"     = "Campina",
    "L02_500"   = "Floresta", "L02_1500" = "Floresta", "L02_2500" = "Floresta",
    "L02_2900"  = "Igarape",  "L02_3500" = "Floresta", "L02_4500" = "Floresta",
    "L02_4350"  = "Igarape",  "NS2_550"  = "Igarape",  "L01_4500" = "Floresta",
    "L01_3500"  = "Floresta", "L01_3050" = "Igarape",  "L01_2500" = "Floresta",
    "L01_1500"  = "Floresta", "TRAV_0"   = "Igarape",  "L01_500"  = "Floresta",
    "NS1_50"    = "Igarape",
    "^PENEIRA"  = "Peneira"
  )
  
  # Match using normalized filenames
  for (pat in names(rules)) {
    dt[str_detect(file_norm, pat), environment := rules[[pat]]]
  }
  
  dt[, file_norm := NULL]
  dt
}

# ---------- Read + cleaning ----------
read_and_shape <- function(path){
  dt <- fread(path, sep = "\t", header = TRUE, na.strings = c("", "NA"))
  setnames(dt, tolower(names(dt)))
  if (!"file" %in% names(dt)) stop("Input must contain 'file' column")
  dt <- add_environment(dt)
  
  # numeric abundance
  if (!"abundance" %in% names(dt)) dt[, abundance := NA_real_]
  suppressWarnings(dt[, abundance := as.numeric(abundance)])
  if ("estimated_counts" %in% names(dt) && sum(dt$abundance, na.rm = TRUE) == 0)
    dt[, abundance := as.numeric(estimated_counts)]
  dt <- dt[is.finite(abundance)]
  
  if (!"genus" %in% names(dt))  dt[, genus  := NA_character_]
  if (!"family" %in% names(dt)) dt[, family := NA_character_]
  
  # derive genus if missing
  if (all(is.na(dt$genus)) && "name" %in% names(dt))
    dt[, genus := sub("\\s.*$", "", name)]
  dt[is.na(genus) | genus == "", genus := "No genus"]
  
  drop_pat <- "(?i)^(unmapped|mapped_unclassified)$"
  dt <- dt[!str_detect(coalesce(genus, ""), drop_pat)]
  list(raw = dt)
}

# ---------- Matrix builder ----------
build_matrix <- function(dt, tax_col){
  stopifnot(tax_col %in% names(dt))
  long <- dt[, .(file, environment, taxon = get(tax_col), value = abundance)]
  long[taxon == "" | is.na(taxon), taxon := paste0("No ", tax_col)]
  long <- long[, .(value = sum(value, na.rm = TRUE)), by = .(file, environment, taxon)]
  wide <- long |>
    select(file, taxon, value) |>
    pivot_wider(names_from = taxon, values_from = value, values_fill = 0) |>
    as.data.table()
  meta <- unique(long[, .(file, environment)])
  setkey(meta, file); setkey(wide, file)
  mat <- as.data.frame(wide); rownames(mat) <- mat$file; mat$file <- NULL
  list(mat = mat, meta = meta, long = long)
}

# ---------- Helpers for plotting ----------
collapse_topN_by_env <- function(long_df, N = 20){
  tmp <- copy(long_df)
  tmp[, total := sum(value), by = file]
  tmp[, rel := fifelse(total > 0, value / total, 0)]
  ranks <- tmp[, .(mean_rel = mean(rel, na.rm = TRUE)),
               by = .(environment, taxon)] |>
    arrange(environment, desc(mean_rel), taxon) |>
    group_by(environment) |>
    mutate(rank = row_number()) |>
    as.data.table()
  tmp <- merge(tmp, ranks[, .(environment, taxon, rank)],
               by = c("environment", "taxon"), all.x = TRUE)
  tmp[, group := ifelse(rank <= N, taxon, "Other")]
  tmp[, .(rel = sum(rel)), by = .(environment, file, group)] |>
    setnames("group", "taxon")
}

palette_with_other_first <- function(levels_vec){
  lev <- unique(levels_vec)
  ordered <- c("Other", setdiff(lev, "Other"))
  
  # Build a long discrete palette: Set1 -> Set2 -> Set3
  get_pal <- function(name){
    maxn <- RColorBrewer::brewer.pal.info[name, "maxcolors"]
    RColorBrewer::brewer.pal(maxn, name)
  }
  pool <- c(get_pal("Set1"), get_pal("Set2"), get_pal("Set3"))
  
  k <- max(0, length(ordered) - 1)  # non-"Other" colors needed
  cols_non_other <- if (k == 0) character(0) else rep(pool, length.out = k)
  
  vals <- c("grey80", cols_non_other)  # "Other" first and grey
  names(vals) <- ordered
  list(levels = ordered, values = vals)
}

legend_label_wrap <- function(x) {
  # replace one or more spaces with a line break
  vapply(x, function(s) if (is.na(s)) s else gsub("\\s+", "\n", s), character(1))
}

make_env_stacks <- function(dt_raw, rank_col, out_png, out_pdf, N = 20, title_rank = "Genus"){
  if (!rank_col %in% names(dt_raw)) return(invisible(NULL))
  if (all(is.na(dt_raw[[rank_col]])) || all(dt_raw[[rank_col]] == "")) return(invisible(NULL))
  
  shaped <- build_matrix(dt_raw, tax_col = rank_col)
  long  <- shaped$long
  collapsed <- collapse_topN_by_env(long, N = N)
  
  # order samples within each environment by file
  collapsed[, file_fac := factor(file, levels = unique(file)), by = environment]
  
  # 'Other' first and legend title set to rank
  lvl_info <- palette_with_other_first(collapsed$taxon)
  collapsed[, taxon := factor(taxon, levels = lvl_info$levels)]
  
  # Determine if facet labels should be shown (hide for Genus)
  facet_text_size <- if (tolower(title_rank) == "genus") 0 else 12
  facet_strip <- if (tolower(title_rank) == "genus") element_blank() else element_text(size = facet_text_size)
  
  p <- ggplot(collapsed, aes(x = file_fac, y = rel, fill = taxon)) +
    geom_col(width = 0.98) +  # bars nearly touching
    facet_grid(~ environment, scales = "free_x", space = "free_x") +
    scale_y_continuous(labels = function(z) 100 * z,
                       expand = expansion(mult = c(0, 0.02))) +
    scale_fill_manual(
      values = lvl_info$values,
      breaks = lvl_info$levels,
      labels = legend_label_wrap(lvl_info$levels),
      name   = title_rank,
      drop   = FALSE
    ) +
    guides(fill = guide_legend(title = title_rank)) +
    labs(x = NULL, y = "Relative abundance (%)") +
    theme_classic(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.spacing.x = unit(0.02, "lines"),
      strip.text = facet_strip,
      strip.background = element_rect(fill = "white", colour = NA),
      axis.title = element_text(face = "plain"),
      axis.text  = element_text(face = "plain"),
      legend.title = element_text(face = "plain"),
      legend.text  = element_text(face = "plain"),
      panel.border = element_blank(),
      legend.position = "right"
    )
  
  ggsave(out_png, p, width = 14, height = 5, dpi = 300)
  ggsave(out_pdf, p, width = 14, height = 5)
  p
}

# ==========================================================
# 0) Read + environment + cleaning
# ==========================================================
read_obj <- read_and_shape(infile)
dt_raw   <- read_obj$raw

fwrite(
  dt_raw,
  file = file.path(outdir, paste0(prefix, "_with_environment.tsv")),
  sep  = "\t"
)

# ==========================================================
# 1) STACKED BARS
# ==========================================================
p_family <- make_env_stacks(
  dt_raw, "family",
  file.path(outdir, paste0(prefix, "_stacks_family.png")),
  file.path(outdir, paste0(prefix, "_stacks_family.pdf")),
  N = 20, title_rank = "Family"
)
p_genus <- make_env_stacks(
  dt_raw, "genus",
  file.path(outdir, paste0(prefix, "_stacks_genus.png")),
  file.path(outdir, paste0(prefix, "_stacks_genus.pdf")),
  N = 20, title_rank = "Genus"
)

if (!is.null(p_family) && !is.null(p_genus)) {
  combined <- plot_grid(
    p_family,  # keep legend
    p_genus,   # keep legend
    ncol = 1, rel_heights = c(1, 1), align = "v"
  )
  ggsave(file.path(outdir, paste0(prefix, "_stacks_family_genus_grid.png")),
         combined, width = 22, height = 10, dpi = 300)
  ggsave(file.path(outdir, paste0(prefix, "_stacks_family_genus_grid.pdf")),
         combined, width = 22, height = 10)
}

# ==========================================================
# 2) ALPHA DIVERSITY (Shannon & Simpson + pairwise Wilcoxon)
# ==========================================================
mx_rel  <- build_matrix(dt_raw, "genus")
mat_rel <- mx_rel$mat
meta    <- mx_rel$meta

# Normalize rows to relative abundance (safe even if already proportions)
row_sums <- rowSums(mat_rel, na.rm = TRUE); row_sums[row_sums == 0] <- 1
rel <- sweep(mat_rel, 1, row_sums, "/")

# Metrics
shannon <- vegan::diversity(rel, index = "shannon")
simpson <- vegan::diversity(rel, index = "simpson")  # 1 - D in vegan

# Save table
alpha_df <- data.table::data.table(
  file    = rownames(mat_rel),
  Shannon = as.numeric(shannon),
  Simpson = as.numeric(simpson)
)
alpha_df <- merge(alpha_df, meta, by = "file", all.x = TRUE)

# Make environment a factor for stable x positions
alpha_df[, environment := factor(environment)]
env_levels <- levels(alpha_df$environment)

data.table::fwrite(
  alpha_df,
  file = file.path(outdir, paste0(prefix, "_alpha_diversity.tsv")),
  sep  = "\t"
)

# ---------- Pairwise Wilcoxon tests (Mann–Whitney) ----------
pairwise_wilcox <- function(df, value_col, metric_name) {
  df <- df[!is.na(environment) & !is.na(get(value_col))]
  envs <- levels(df$environment)
  envs <- envs[envs %in% unique(df$environment)]
  if (length(envs) < 2) {
    return(data.table(env1 = character(0), env2 = character(0),
                      metric = character(0), p_value = numeric(0)))
  }
  
  cmb <- t(combn(envs, 2))
  res_list <- apply(cmb, 1, function(pair_env) {
    e1 <- pair_env[1]; e2 <- pair_env[2]
    x <- df[environment == e1][[value_col]]
    y <- df[environment == e2][[value_col]]
    p <- tryCatch(
      wilcox.test(x, y)$p.value,
      error = function(e) NA_real_
    )
    data.table(env1 = e1, env2 = e2, metric = metric_name, p_value = p)
  })
  rbindlist(res_list)
}

pw_sh <- pairwise_wilcox(alpha_df, "Shannon", "Shannon")
pw_sp <- pairwise_wilcox(alpha_df, "Simpson", "Simpson")

pairwise_all <- rbindlist(list(pw_sh, pw_sp))
data.table::fwrite(
  pairwise_all,
  file = file.path(outdir, paste0(prefix, "_alpha_pairwise_wilcox.tsv")),
  sep  = "\t"
)

# ---------- Helpers for plotting bars + stars ----------
p_to_star <- function(p) {
  if (!is.finite(p) || is.na(p) || p > 0.05) "" else
    if (p <= 0.001) "***" else if (p <= 0.01) "**" else "*"
}

# Significant comparisons only (p <= 0.05)
sig_sh <- pw_sh[is.finite(p_value) & p_value <= 0.05]
sig_sp <- pw_sp[is.finite(p_value) & p_value <= 0.05]

# Add x positions for each env (factor index)
env_to_x <- function(e) which(env_levels == e)

if (nrow(sig_sh) > 0) {
  sig_sh[, x1 := env_to_x(env1)]
  sig_sh[, x2 := env_to_x(env2)]
  sig_sh[, star := vapply(p_value, p_to_star, character(1))]
}

if (nrow(sig_sp) > 0) {
  sig_sp[, x1 := env_to_x(env1)]
  sig_sp[, x2 := env_to_x(env2)]
  sig_sp[, star := vapply(p_value, p_to_star, character(1))]
}

# y-positions for bars
sh_max <- max(alpha_df$Shannon, na.rm = TRUE)
sh_range <- diff(range(alpha_df$Shannon, na.rm = TRUE))
if (!is.finite(sh_range) || sh_range == 0) sh_range <- 1
base_sh <- sh_max + 0.08 * sh_range
step_sh <- 0.06 * sh_range

sp_max <- max(alpha_df$Simpson, na.rm = TRUE)
sp_range <- diff(range(alpha_df$Simpson, na.rm = TRUE))
if (!is.finite(sp_range) || sp_range == 0) sp_range <- 0.1
base_sp <- sp_max + 0.08 * sp_range
step_sp <- 0.06 * sp_range

if (nrow(sig_sh) > 0) {
  sig_sh[, y := base_sh + step_sh * (seq_len(.N) - 1L)]
}
if (nrow(sig_sp) > 0) {
  sig_sp[, y := base_sp + step_sp * (seq_len(.N) - 1L)]
}

# Common theme + colors
theme_base <- theme_classic(base_size = 12) + theme(panel.grid = element_blank())

env_colors <- c(
  "Campina"  = "#FFCC00",
  "Floresta" = "#99CC33",
  "Igarape"  = "#3399FF",
  "Peneira"  = "#FF9900"
)

# ==========================================================
# Plots with bars + asterisks
# ==========================================================
p_sh <- ggplot(alpha_df, aes(x = environment, y = Shannon, fill = environment, color = environment)) +
  geom_violin(alpha = 0.25, linewidth = 0, position = position_dodge(width = 0.75), show.legend = FALSE) +
  geom_quasirandom(shape = 21, size = 1, dodge.width = 0.75,
                   alpha = 0.5, color = "black", show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.9,
               color = "black", fill = "white", show.legend = FALSE) +
  # bars + stars for significant comparisons
  { if (nrow(sig_sh) > 0) 
    geom_segment(data = sig_sh,
                 aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.4) 
    else NULL } +
  { if (nrow(sig_sh) > 0)
    geom_text(data = sig_sh,
              aes(x = (x1 + x2)/2, y = y + step_sh * 0.2, label = star),
              inherit.aes = FALSE, vjust = 0, size = 3)
    else NULL } +
  expand_limits(y = if (nrow(sig_sh) > 0) max(sig_sh$y + step_sh * 0.4) else sh_max) +
  labs(x = "\nEnvironment", y = "Shannon (H')\n", title = "") +
  scale_fill_manual(values = env_colors, drop = FALSE) +
  scale_color_manual(values = env_colors, drop = FALSE) +
  theme_base
ggsave(file.path(outdir, paste0(prefix, "_alpha_shannon_env.png")),
       p_sh, width = 3, height = 5, dpi = 300)

p_sp <- ggplot(alpha_df, aes(x = environment, y = Simpson, fill = environment, color = environment)) +
  geom_violin(alpha = 0.25, linewidth = 0, position = position_dodge(width = 0.75), show.legend = FALSE) +
  geom_quasirandom(shape = 21, size = 1, dodge.width = 0.75,
                   alpha = 0.5, color = "black", show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.9,
               color = "black", fill = "white", show.legend = FALSE) +
  # bars + stars for significant comparisons
  { if (nrow(sig_sp) > 0) 
    geom_segment(data = sig_sp,
                 aes(x = x1, xend = x2, y = y, yend = y),
                 inherit.aes = FALSE, linewidth = 0.4) 
    else NULL } +
  { if (nrow(sig_sp) > 0)
    geom_text(data = sig_sp,
              aes(x = (x1 + x2)/2, y = y + step_sp * 0.2, label = star),
              inherit.aes = FALSE, vjust = 0, size = 3)
    else NULL } +
  expand_limits(y = if (nrow(sig_sp) > 0) max(sig_sp$y + step_sp * 0.4) else sp_max) +
  labs(x = "\nEnvironment", y = "Simpson (1 - D)\n", title = "") +
  scale_fill_manual(values = env_colors, drop = FALSE) +
  scale_color_manual(values = env_colors, drop = FALSE) +
  theme_base
ggsave(file.path(outdir, paste0(prefix, "_alpha_simpson_env.png")),
       p_sp, width = 3, height = 5, dpi = 300)

# ==========================================================
# 3) BETA DIVERSITY – Bray–Curtis PCoA (Genus)
# ==========================================================
bray <- vegan::vegdist(rel, method = "bray")

pcoa <- cmdscale(bray, k = 2, eig = TRUE)
eig  <- pcoa$eig
eig[eig < 0] <- 0  # guard against tiny negative numerical eigenvalues
var_expl <- eig / sum(eig)
pc1_lab <- sprintf("PC1 (%.1f%%)", 100 * var_expl[1])
pc2_lab <- sprintf("PC2 (%.1f%%)", 100 * var_expl[2])

pcoa_df <- data.table::data.table(
  file = rownames(rel),
  PC1  = pcoa$points[, 1],
  PC2  = pcoa$points[, 2]
)
pcoa_df <- merge(pcoa_df, meta, by = "file", all.x = TRUE)
data.table::fwrite(
  pcoa_df,
  file = file.path(outdir, paste0(prefix, "_pcoa_braycurtis.tsv")),
  sep  = "\t"
)

p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = environment)) +
  geom_point(size = 2.5, alpha = 0.9) +
  scale_color_manual(values = env_colors, na.value = "grey70", drop = FALSE) +
  labs(
    x = pc1_lab,
    y = pc2_lab,
    title = ""  # Bray–Curtis PCoA (Genus)
  ) +
  theme_base
ggsave(file.path(outdir, paste0(prefix, "_pcoa_braycurtis_env.png")),
       p_pcoa, width = 5, height = 4, dpi = 300)

message(">>> Done. Outputs in: ", outdir)
