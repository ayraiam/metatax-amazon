#!/usr/bin/env Rscript
# ==========================================================
# downstream_analysis.R
# 1) Add 'environment' from filename rules
# 2) Make environment-grouped stacked bars (Family & Genus)
# 3) Alpha diversity (Shannon / Simpson)
# 4) Beta diversity (Bray–Curtis PCoA, Genus)
# 5) Beta-diversity stats: PERMANOVA + betadisper
# 6) NEW: Per-code genus CLR concordance scatter plots (Peneira vs Floresta)
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

### downstream mode + numeric source (passed from runall.sh)
MODE <- toupper(Sys.getenv("MODE", "16S"))
USE_COUNTS <- as.integer(Sys.getenv("USE_COUNTS", "1"))
message(">>> MODE = ", MODE)
message(">>> USE_COUNTS = ", USE_COUNTS)

# ==========================================================
# 0) If infile does NOT exist, merge batch tables_b*/abundance_combined.tsv
# ==========================================================
if (!file.exists(infile)) {
  message(">>> Input file not found: ", infile)
  message(">>> Trying to merge batch abundance_combined.tsv files from results/tables_b*/ ...")
  
  results_dir <- dirname(dirname(infile))
  if (!dir.exists(results_dir)) {
    results_dir <- file.path(getwd(), "results")
  }
  message(">>> Searching for batch tables under: ", results_dir)
  
  all_abund <- list.files(
    path       = results_dir,
    pattern    = "^abundance_combined\\.tsv$",
    recursive  = TRUE,
    full.names = TRUE
  )
  
  tables_dir_name <- basename(dirname(infile))
  pattern_tables  <- paste0("/", tables_dir_name, "_b")
  message(">>> Using batch pattern filter: ", pattern_tables)
  
  batch_files <- all_abund[grepl(pattern_tables, all_abund)]
  
  if (length(batch_files) == 0L) {
    stop(
      "Could not find infile or any batch ", tables_dir_name,
      "_b*/abundance_combined.tsv under: ", results_dir
    )
  }
  
  message(">>> Found ", length(batch_files), " batch abundance tables. Merging...")
  dt_list <- lapply(
    batch_files,
    function(f) data.table::fread(f, sep = "\t", header = TRUE, na.strings = c("", "NA"))
  )
  merged_dt <- data.table::rbindlist(dt_list, use.names = TRUE, fill = TRUE)
  data.table::setnames(merged_dt, tolower(names(merged_dt)))
  
  dir.create(dirname(infile), showWarnings = FALSE, recursive = TRUE)
  data.table::fwrite(merged_dt, file = infile, sep = "\t")
  message(">>> Merged table written to: ", infile)
}

# ---------- Helper: add environment column ----------
add_environment <- function(dt) {
  stopifnot("file" %in% names(dt))
  
  dt[, file_norm := toupper(file)]
  dt[, file_norm := gsub("LO([0-9])", "L0\\1", file_norm)]
  dt[, file_norm := gsub("IO([0-9])", "L0\\1", file_norm)]
  dt[, file_norm := gsub("O([0-9])", "0\\1", file_norm)]
  
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
  
  for (pat in names(rules)) {
    dt[str_detect(file_norm, pat), environment := rules[[pat]]]
  }
  
  dt[, file_norm := NULL]
  dt
}

### extract numeric code from filename (0500 -> 500)
add_code <- function(dt) {
  stopifnot("file" %in% names(dt))
  f <- toupper(dt$file)
  
  code_chr <- rep(NA_character_, length(f))
  
  # PENEIRA_0500_...
  idx_p <- str_detect(f, "^PENEIRA_")
  if (any(idx_p)) {
    code_chr[idx_p] <- str_match(f[idx_p], "^PENEIRA_([0-9]+)_")[, 2]
  }
  
  # L01_500_... or L02_1500_...
  idx_l <- str_detect(f, "^L0[12]_")
  if (any(idx_l)) {
    code_chr[idx_l] <- str_match(f[idx_l], "^L0[12]_([0-9]+)_")[, 2]
  }
  
  dt[, code := suppressWarnings(as.integer(code_chr))]
  dt
}

# ---------- Read + cleaning ----------
read_and_shape <- function(path){
  dt <- fread(path, sep = "\t", header = TRUE, na.strings = c("", "NA"))
  setnames(dt, tolower(names(dt)))
  
  ### normalize possible "estimated counts" column name
  if ("estimated counts" %in% names(dt)) setnames(dt, "estimated counts", "estimated_counts")
  
  if (!"file" %in% names(dt)) stop("Input must contain 'file' column")
  dt <- add_environment(dt)
  dt <- add_code(dt)  
  
  # numeric abundance + counts
  if (!"abundance" %in% names(dt)) dt[, abundance := NA_real_]
  suppressWarnings(dt[, abundance := as.numeric(abundance)])
  
  if (!"estimated_counts" %in% names(dt)) dt[, estimated_counts := NA_real_]
  suppressWarnings(dt[, estimated_counts := as.numeric(estimated_counts)])
  
  ### unify numeric column used downstream
  dt[, value := if (USE_COUNTS == 1) estimated_counts else abundance]
  
  dt <- dt[is.finite(value)]
  
  if (!"genus" %in% names(dt))  dt[, genus  := NA_character_]
  if (!"family" %in% names(dt)) dt[, family := NA_character_]
  
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
  ### use unified 'value' instead of 'abundance'
  long <- dt[, .(file, environment, taxon = get(tax_col), value = value)]
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
  
  get_pal <- function(name){
    maxn <- RColorBrewer::brewer.pal.info[name, "maxcolors"]
    RColorBrewer::brewer.pal(maxn, name)
  }
  pool <- c(get_pal("Set1"), get_pal("Set2"), get_pal("Set3"))
  
  k <- max(0, length(ordered) - 1)
  cols_non_other <- if (k == 0) character(0) else rep(pool, length.out = k)
  
  vals <- c("grey80", cols_non_other)
  names(vals) <- ordered
  list(levels = ordered, values = vals)
}

legend_label_wrap <- function(x) {
  vapply(x, function(s) if (is.na(s)) s else gsub("\\s+", "\n", s), character(1))
}

make_env_stacks <- function(dt_raw, rank_col, out_png, out_pdf, N = 20, title_rank = "Genus"){
  if (!rank_col %in% names(dt_raw)) return(invisible(NULL))
  if (all(is.na(dt_raw[[rank_col]])) || all(dt_raw[[rank_col]] == "")) return(invisible(NULL))
  
  shaped <- build_matrix(dt_raw, tax_col = rank_col)
  long  <- shaped$long
  collapsed <- collapse_topN_by_env(long, N = N)
  
  collapsed[, file_fac := factor(file, levels = unique(file)), by = environment]
  
  lvl_info <- palette_with_other_first(collapsed$taxon)
  collapsed[, taxon := factor(taxon, levels = lvl_info$levels)]
  
  facet_text_size <- if (tolower(title_rank) == "genus") 0 else 12
  facet_strip <- if (tolower(title_rank) == "genus") element_blank() else element_text(size = facet_text_size)
  
  p <- ggplot(collapsed, aes(x = file_fac, y = rel, fill = taxon)) +
    geom_col(width = 0.98) +
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
    p_family,
    p_genus,
    ncol = 1, rel_heights = c(1, 1), align = "v"
  )
  ggsave(file.path(outdir, paste0(prefix, "_stacks_family_genus_grid.png")),
         combined, width = 22, height = 10, dpi = 300)
  ggsave(file.path(outdir, paste0(prefix, "_stacks_family_genus_grid.pdf")),
         combined, width = 22, height = 10)
}

# ==========================================================
# 2) ALPHA DIVERSITY (Shannon & Simpson only)
# ==========================================================
mx_rel  <- build_matrix(dt_raw, "genus")
mat_rel <- mx_rel$mat
meta    <- mx_rel$meta

row_sums <- rowSums(mat_rel, na.rm = TRUE); row_sums[row_sums == 0] <- 1
rel <- sweep(mat_rel, 1, row_sums, "/")

shannon <- vegan::diversity(rel, index = "shannon")
simpson <- vegan::diversity(rel, index = "simpson")

alpha_df <- data.table::data.table(
  file    = rownames(mat_rel),
  Shannon = as.numeric(shannon),
  Simpson = as.numeric(simpson)
)
alpha_df <- merge(alpha_df, meta, by = "file", all.x = TRUE)

alpha_df[, environment := factor(
  environment,
  levels = c("Campina", "Floresta", "Igarape", "Peneira")
)]

data.table::fwrite(
  alpha_df,
  file = file.path(outdir, paste0(prefix, "_alpha_diversity.tsv")),
  sep  = "\t"
)

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

pval_to_stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "ns")
         )
  )
}

build_sig_df <- function(sig_pw, df, value_col) {
  if (nrow(sig_pw) == 0) return(NULL)
  
  env_levels <- levels(df$environment)
  rng <- range(df[[value_col]], na.rm = TRUE)
  y_min <- rng[1]; y_max <- rng[2]
  y_step <- 0.05 * (y_max - y_min)
  if (is.na(y_step) || y_step == 0) y_step <- 0.1
  
  sig_pw <- copy(sig_pw)[order(p_value)]
  sig_pw[, idx := seq_len(.N)]
  sig_pw[, `:=`(
    x    = match(env1, env_levels),
    xend = match(env2, env_levels),
    y    = y_max + idx * y_step,
    label = pval_to_stars(p_value)
  )]
  sig_pw
}

pw_shannon <- pairwise_wilcox(alpha_df, "Shannon", "Shannon")
pw_simpson <- pairwise_wilcox(alpha_df, "Simpson", "Simpson")
pw_all <- rbind(pw_shannon, pw_simpson)

fwrite(
  pw_all,
  file = file.path(outdir, paste0(prefix, "_alpha_pairwise_wilcox.tsv")),
  sep  = "\t"
)

sig_shannon <- pw_shannon[!is.na(p_value) & p_value < 0.05]
sig_simpson <- pw_simpson[!is.na(p_value) & p_value < 0.05]

sig_sh_df <- build_sig_df(sig_shannon, alpha_df, "Shannon")
sig_sp_df <- build_sig_df(sig_simpson, alpha_df, "Simpson")

theme_base <- theme_classic(base_size = 12) + theme(panel.grid = element_blank())

env_colors <- c(
  "Campina"  = "#FFCC00",
  "Floresta" = "#99CC33",
  "Igarape"  = "#3399FF",
  "Peneira"  = "#FF9900"
)

p_sh <- ggplot(alpha_df, aes(x = environment, y = Shannon, fill = environment, color = environment)) +
  geom_violin(alpha = 0.25, linewidth = 0, position = position_dodge(width = 0.75), show.legend = FALSE) +
  geom_quasirandom(shape = 21, size = 1, dodge.width = 0.75,
                   alpha = 0.5, color = "black", show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.9, color = "black", fill = "white", show.legend = FALSE) +
  labs(x = "\nEnvironment", y = "Shannon (H')\n", title = "") +
  scale_fill_manual(values = env_colors) +
  scale_color_manual(values = env_colors) +
  theme_base

if (!is.null(sig_sh_df) && nrow(sig_sh_df) > 0) {
  p_sh <- p_sh +
    geom_segment(
      data = sig_sh_df,
      aes(x = x, xend = xend, y = y, yend = y),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = sig_sh_df,
      aes(x = (x + xend) / 2, y = y, label = label),
      vjust = -0.3, size = 3, inherit.aes = FALSE
    )
}

ggsave(file.path(outdir, paste0(prefix, "_alpha_shannon_env.png")),
       p_sh, width = 3, height = 5, dpi = 300)

p_sp <- ggplot(alpha_df, aes(x = environment, y = Simpson, fill = environment, color = environment)) +
  geom_violin(alpha = 0.25, linewidth = 0, position = position_dodge(width = 0.75), show.legend = FALSE) +
  geom_quasirandom(shape = 21, size = 1, dodge.width = 0.75,
                   alpha = 0.5, color = "black", show.legend = FALSE) +
  geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.9, color = "black", fill = "white", show.legend = FALSE) +
  labs(x = "\nEnvironment", y = "Simpson (1 - D)\n", title = "") +
  scale_fill_manual(values = env_colors) +
  scale_color_manual(values = env_colors) +
  theme_base

if (!is.null(sig_sp_df) && nrow(sig_sp_df) > 0) {
  p_sp <- p_sp +
    geom_segment(
      data = sig_sp_df,
      aes(x = x, xend = xend, y = y, yend = y),
      inherit.aes = FALSE
    ) +
    geom_text(
      data = sig_sp_df,
      aes(x = (x + xend) / 2, y = y, label = label),
      vjust = -0.3, size = 3, inherit.aes = FALSE
    )
}

ggsave(file.path(outdir, paste0(prefix, "_alpha_simpson_env.png")),
       p_sp, width = 3, height = 5, dpi = 300)

# ==========================================================
# 3) BETA DIVERSITY – Bray–Curtis PCoA (Genus)
# ==========================================================
bray <- vegan::vegdist(rel, method = "bray")

pcoa <- cmdscale(bray, k = 2, eig = TRUE)
eig  <- pcoa$eig
eig[eig < 0] <- 0
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
  scale_color_manual(values = env_colors, na.value = "grey70") +
  labs(
    x = pc1_lab,
    y = pc2_lab,
    title = ""
  ) +
  theme_base
ggsave(file.path(outdir, paste0(prefix, "_pcoa_braycurtis_env.png")),
       p_pcoa, width = 5, height = 4, dpi = 300)

# ==========================================================
# 4) BETA STATS – PERMANOVA + betadisper
# ==========================================================
set.seed(2025)
meta$environment <- factor(meta$environment)

perm <- vegan::adonis2(bray ~ environment,
                       data = meta,
                       permutations = 999)
perm_df <- as.data.frame(perm)
data.table::fwrite(
  as.data.table(perm_df, keep.rownames = "term"),
  file = file.path(outdir, paste0(prefix, "_beta_permanova.tsv")),
  sep  = "\t"
)

bd <- vegan::betadisper(bray, meta$environment)

bd_anova   <- as.data.frame(anova(bd))
bd_perm    <- vegan::permutest(bd, permutations = 999)
bd_perm_df <- as.data.frame(bd_perm$tab)

data.table::fwrite(
  as.data.table(bd_anova, keep.rownames = "term"),
  file = file.path(outdir, paste0(prefix, "_beta_betadisper_anova.tsv")),
  sep  = "\t"
)
data.table::fwrite(
  as.data.table(bd_perm_df, keep.rownames = "term"),
  file = file.path(outdir, paste0(prefix, "_beta_betadisper_permutest.tsv")),
  sep  = "\t"
)

# ==========================================================
# 5) PER-CODE GENUS CLR CONCORDANCE (Peneira vs Floresta)
#     - One scatter plot per code
#     - x = CLR(Floresta), y = CLR(Peneira)
#     - dots = genera
# ==========================================================
dt_pf <- dt_raw[environment %in% c("Peneira", "Floresta") & !is.na(code)]

target_codes <- c(500, 1500, 2500, 3500, 4500)
dt_pf <- dt_pf[code %in% target_codes]

mx_g <- build_matrix(dt_pf, "genus")
mat_g <- mx_g$mat
meta_g <- as.data.table(mx_g$meta)

# attach code to meta (file -> code)
meta_g <- merge(meta_g, unique(dt_pf[, .(file, code)]), by = "file", all.x = TRUE)

# CLR: relative -> +pseudocount -> log -> center
row_sums_g <- rowSums(mat_g, na.rm = TRUE); row_sums_g[row_sums_g == 0] <- 1
rel_g <- sweep(mat_g, 1, row_sums_g, "/")
pseudocount <- 1e-6
rel_g <- rel_g + pseudocount
log_g <- log(rel_g)
clr_g <- sweep(log_g, 1, rowMeans(log_g), "-")

corr_dir <- file.path(outdir, paste0(prefix, "_code_concordance_", MODE))
dir.create(corr_dir, showWarnings = FALSE, recursive = TRUE)

for (cc in target_codes) {
  files_p <- meta_g[environment == "Peneira" & code == cc, file]
  files_f <- meta_g[environment == "Floresta" & code == cc, file]
  
  if (length(files_p) == 0 || length(files_f) == 0) next
  
  vP <- colMeans(clr_g[files_p, , drop = FALSE], na.rm = TRUE)
  vF <- colMeans(clr_g[files_f, , drop = FALSE], na.rm = TRUE)
  
  df <- data.table(
    genus = names(vP),
    CLR_Floresta = as.numeric(vF[names(vP)]),
    CLR_Peneira  = as.numeric(vP)
  )
  df <- df[is.finite(CLR_Floresta) & is.finite(CLR_Peneira)]
  
  ct <- suppressWarnings(cor.test(df$CLR_Floresta, df$CLR_Peneira, method = "pearson"))
  r <- unname(ct$estimate)
  r2 <- r^2
  pval <- ct$p.value
  
  df[, resid_abs := abs(CLR_Peneira - CLR_Floresta)]
  top_lab <- df[order(-resid_abs)][1:min(5, .N)]
  
  p <- ggplot(df, aes(x = CLR_Floresta, y = CLR_Peneira)) +
    geom_point(alpha = 0.75, size = 1.7) +
    geom_smooth(method = "lm", se = FALSE, linewidth = 0.6, linetype = "dashed") +
    geom_text(data = top_lab, aes(label = genus), size = 3, vjust = -0.6) +
    labs(
      title = paste0("Code ", cc, " (", MODE, ")"),
      subtitle = sprintf("Pearson R² = %.2f, p = %.2g", r2, pval),
      x = "CLR abundance (Floresta)",
      y = "CLR abundance (Peneira)"
    ) +
    theme_classic(base_size = 12)
  
  ggsave(file.path(corr_dir, paste0(prefix, "_code_", cc, "_scatter.png")),
         p, width = 5.2, height = 4.2, dpi = 300)
  ggsave(file.path(corr_dir, paste0(prefix, "_code_", cc, "_scatter.pdf")),
         p, width = 5.2, height = 4.2)
  
  fwrite(df[, .(genus, CLR_Floresta, CLR_Peneira)],
         file = file.path(corr_dir, paste0(prefix, "_code_", cc, "_clr.tsv")),
         sep = "\t")
}

message(">>> Done. Outputs in: ", outdir)
message(">>> Concordance scatter plots in: ", corr_dir)
