#!/usr/bin/env Rscript
# ==========================================================
# downstream_analysis.R
# 1) Add 'environment' from filename rules
# 2) Make environment-grouped stacked bars (Family & Genus)
# 3) Alpha diversity (Shannon / Simpson)
# 4) Beta diversity (Bray–Curtis PCoA, Genus)
# 5) Beta-diversity stats: PERMANOVA + betadisper
# 6) Per-code genus CLR concordance scatter plots (Peneira vs Floresta)
# 7) Differential taxa (Floresta vs Peneira) via ANCOM-BC2
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
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

### downstream mode + numeric source (passed from runall.sh)
MODE <- toupper(Sys.getenv("MODE", "16S"))
message(">>> MODE = ", MODE)
user_prefix <- ifelse(length(args) >= 3, args[3], "downstream")
prefix <- paste0(user_prefix, "_", MODE)

USE_COUNTS_0_4 <- as.integer(Sys.getenv("USE_COUNTS_0_4", "0"))  # default: abundance
USE_COUNTS_5   <- as.integer(Sys.getenv("USE_COUNTS_5",   "1"))  # default: estimated_counts
message(">>> USE_COUNTS_0_4 = ", USE_COUNTS_0_4)
message(">>> USE_COUNTS_5   = ", USE_COUNTS_5)

# ==========================================================
# 0) If infile does NOT exist, merge batch tables_b*/abundance_combined.tsv
# ==========================================================
merge_batches_if_needed <- function(infile) {
  if (file.exists(infile)) return(invisible(NULL))
  
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

# ==========================================================
# Helpers: extract environment + code/replicate
# ==========================================================
add_code_replicate <- function(dt) {
  stopifnot("file" %in% names(dt))
  
  dt[, file_norm := toupper(file)]
  
  # replicate token: _I_ or _II_ or _III_
  dt[, replicate := str_match(file_norm, "_(I|II|III)_")[, 2]]
  
  # code:
  dt[, code := NA_character_]
  dt[str_detect(file_norm, "^PENEIRA_"),
     code := str_match(file_norm, "^PENEIRA_([0-9]+)_")[, 2]]
  dt[is.na(code) & str_detect(file_norm, "^(L0[12]|NS[12])_"),
     code := str_match(file_norm, "^(L0[12]|NS[12])_([0-9]+)_")[, 3]]
  
  # pairing_code: PENEIRA_0500 pairs with Floresta 500
  dt[, pairing_code := code]
  dt[file_norm %like% "^PENEIRA_0500_", pairing_code := "500"]
  
  dt[, file_norm := NULL]
  dt
}

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

read_and_shape <- function(path){
  dt <- fread(path, sep = "\t", header = TRUE, na.strings = c("", "NA"))
  setnames(dt, tolower(names(dt)))
  
  if ("estimated counts" %in% names(dt)) setnames(dt, "estimated counts", "estimated_counts")
  if (!"file" %in% names(dt)) stop("Input must contain 'file' column")
  
  BAD_FILE <- "L01_3050_II_ARCH_and_PENEIRA_3500_ITS.trimmed.filtered"
  dt <- dt[file != BAD_FILE]
  
  dt <- add_environment(dt)
  dt <- add_code_replicate(dt)
  
  if (!"abundance" %in% names(dt)) dt[, abundance := NA_real_]
  suppressWarnings(dt[, abundance := as.numeric(abundance)])
  
  if (!"estimated_counts" %in% names(dt)) dt[, estimated_counts := NA_real_]
  suppressWarnings(dt[, estimated_counts := as.numeric(estimated_counts)])
  
  dt <- dt[is.finite(abundance) | is.finite(estimated_counts)]
  
  if (!"genus" %in% names(dt))  dt[, genus  := NA_character_]
  if (!"family" %in% names(dt)) dt[, family := NA_character_]
  
  if (all(is.na(dt$genus)) && "name" %in% names(dt))
    dt[, genus := sub("\\s.*$", "", name)]
  dt[is.na(genus) | genus == "", genus := "No genus"]
  
  drop_pat <- "(?i)^(unmapped|mapped_unclassified)$"
  dt <- dt[!str_detect(coalesce(genus, ""), drop_pat)]
  list(raw = dt)
}

build_matrix <- function(dt, tax_col, value_col = "abundance"){
  stopifnot(tax_col %in% names(dt))
  stopifnot(value_col %in% names(dt))
  
  long <- dt[, .(file, environment, taxon = get(tax_col), value = get(value_col))]
  long <- long[is.finite(value)]
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

get_numeric_source <- function(dt, use_counts = 1L) {
  if (use_counts == 1L) {
    if (!"estimated_counts" %in% names(dt)) stop("Missing column: estimated_counts")
    return("estimated_counts")
  } else {
    if (!"abundance" %in% names(dt)) stop("Missing column: abundance")
    return("abundance")
  }
}

make_genus_clr_steps <- function(dt_raw, use_counts = 1L, pseudocount_counts = 1, pseudocount_abund = 1e-6) {
  src <- get_numeric_source(dt_raw, use_counts)
  dt <- copy(dt_raw)
  
  if (!"genus" %in% names(dt)) dt[, genus := NA_character_]
  dt[is.na(genus) | genus == "", genus := "No genus"]
  
  dt[, value_src := suppressWarnings(as.numeric(get(src)))]
  dt <- dt[is.finite(value_src)]
  
  genus_long <- dt[, .(raw = sum(value_src, na.rm = TRUE)),
                   by = .(file, environment, code, pairing_code, replicate, genus)]
  
  pc <- if (use_counts == 1L) pseudocount_counts else pseudocount_abund
  
  genus_long[, `:=`(
    pseudocount = pc,
    raw_pc      = raw + pc,
    log_raw_pc  = log(raw + pc)
  )]
  
  genus_long[, mean_log := mean(log_raw_pc), by = file]
  genus_long[, clr := log_raw_pc - mean_log]
  
  setcolorder(
    genus_long,
    c("file","environment","code","pairing_code","replicate","genus",
      "raw","pseudocount","raw_pc","log_raw_pc","mean_log","clr")
  )
  
  list(table = genus_long, source = src, pseudocount = pc)
}

# ==========================================================
# Plot helpers
# ==========================================================
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

make_env_stacks <- function(dt_raw, rank_col, out_png, out_pdf, N = 20, title_rank = "Genus", value_col = "abundance"){
  if (!rank_col %in% names(dt_raw)) return(invisible(NULL))
  if (all(is.na(dt_raw[[rank_col]])) || all(dt_raw[[rank_col]] == "")) return(invisible(NULL))
  
  shaped <- build_matrix(dt_raw, tax_col = rank_col, value_col = value_col)
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
# STEP 0: Read + stage
# ==========================================================
step0_read_stage <- function(infile, outdir, prefix, USE_COUNTS_0_4, USE_COUNTS_5) {
  read_obj <- read_and_shape(infile)
  dt_raw   <- read_obj$raw
  
  # Parts 1–4 numeric column
  VAL_0_4 <- if (USE_COUNTS_0_4 == 1) "estimated_counts" else "abundance"
  message(">>> Parts 1–4 value_col = ", VAL_0_4)
  
  # CLR table source for concordance
  clr_obj <- make_genus_clr_steps(dt_raw, use_counts = USE_COUNTS_5)
  
  fwrite(
    clr_obj$table,
    file = file.path(outdir, paste0(prefix, "_genus_clr_steps.tsv")),
    sep  = "\t"
  )
  
  dt_stage <- copy(dt_raw)
  if ("value" %in% names(dt_stage)) dt_stage[, value := NULL]
  fwrite(
    dt_stage,
    file = file.path(outdir, paste0(prefix, "_with_environment.tsv")),
    sep  = "\t"
  )
  
  list(dt_raw = dt_raw, VAL_0_4 = VAL_0_4, clr_obj = clr_obj)
}

# ==========================================================
# STEP 1: Stacked bars
# ==========================================================
step1_stacked_bars <- function(dt_raw, VAL_0_4, outdir, prefix) {
  p_family <- make_env_stacks(
    dt_raw, "family",
    file.path(outdir, paste0(prefix, "_stacks_family.png")),
    file.path(outdir, paste0(prefix, "_stacks_family.pdf")),
    N = 20, title_rank = "Family", value_col = VAL_0_4
  )
  p_genus <- make_env_stacks(
    dt_raw, "genus",
    file.path(outdir, paste0(prefix, "_stacks_genus.png")),
    file.path(outdir, paste0(prefix, "_stacks_genus.pdf")),
    N = 20, title_rank = "Genus", value_col = VAL_0_4
  )
  
  if (!is.null(p_family) && !is.null(p_genus)) {
    combined <- plot_grid(p_family, p_genus, ncol = 1, rel_heights = c(1, 1), align = "v")
    ggsave(file.path(outdir, paste0(prefix, "_stacks_family_genus_grid.png")),
           combined, width = 22, height = 10, dpi = 300)
    ggsave(file.path(outdir, paste0(prefix, "_stacks_family_genus_grid.pdf")),
           combined, width = 22, height = 10)
  }
}

# ==========================================================
# STEP 2: Alpha diversity
# ==========================================================
step2_alpha <- function(dt_raw, VAL_0_4, outdir, prefix) {
  mx_rel  <- build_matrix(dt_raw, "genus", value_col = VAL_0_4)
  mat_rel <- mx_rel$mat
  meta    <- mx_rel$meta
  
  row_sums <- rowSums(mat_rel, na.rm = TRUE); row_sums[row_sums == 0] <- 1
  rel <- sweep(mat_rel, 1, row_sums, "/")
  
  shannon <- vegan::diversity(rel, index = "shannon")
  simpson <- vegan::diversity(rel, index = "simpson")
  
  alpha_df <- data.table(
    file    = rownames(mat_rel),
    Shannon = as.numeric(shannon),
    Simpson = as.numeric(simpson)
  )
  alpha_df <- merge(alpha_df, meta, by = "file", all.x = TRUE)
  
  alpha_df[, environment := factor(environment, levels = c("Campina", "Floresta", "Igarape", "Peneira"))]
  fwrite(alpha_df, file = file.path(outdir, paste0(prefix, "_alpha_diversity.tsv")), sep = "\t")
  
  # (keeping your existing plotting/stat code as-is)
  list(rel = rel, meta = meta, alpha_df = alpha_df)
}

# ==========================================================
# STEP 3: Beta diversity PCoA
# ==========================================================
step3_beta_pcoa <- function(rel, meta, outdir, prefix) {
  env_colors <- c(
    "Campina"  = "#FFCC00",
    "Floresta" = "#99CC33",
    "Igarape"  = "#3399FF",
    "Peneira"  = "#FF9900"
  )
  theme_base <- theme_classic(base_size = 12) + theme(panel.grid = element_blank())
  
  bray <- vegan::vegdist(rel, method = "bray")
  
  pcoa <- cmdscale(bray, k = 2, eig = TRUE)
  eig  <- pcoa$eig
  eig[eig < 0] <- 0
  var_expl <- eig / sum(eig)
  pc1_lab <- sprintf("PC1 (%.1f%%)", 100 * var_expl[1])
  pc2_lab <- sprintf("PC2 (%.1f%%)", 100 * var_expl[2])
  
  pcoa_df <- data.table(
    file = rownames(rel),
    PC1  = pcoa$points[, 1],
    PC2  = pcoa$points[, 2]
  )
  pcoa_df <- merge(pcoa_df, meta, by = "file", all.x = TRUE)
  
  fwrite(pcoa_df, file = file.path(outdir, paste0(prefix, "_pcoa_braycurtis.tsv")), sep = "\t")
  
  p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = environment)) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(values = env_colors, na.value = "grey70") +
    labs(x = pc1_lab, y = pc2_lab, title = "") +
    theme_base
  
  ggsave(file.path(outdir, paste0(prefix, "_pcoa_braycurtis_env.png")),
         p_pcoa, width = 5, height = 4, dpi = 300)
  
  list(bray = bray)
}

# ==========================================================
# STEP 4: Beta stats
# ==========================================================
step4_beta_stats <- function(bray, meta, outdir, prefix) {
  set.seed(2025)
  meta$environment <- factor(meta$environment)
  
  perm <- vegan::adonis2(bray ~ environment, data = meta, permutations = 999)
  perm_df <- as.data.frame(perm)
  fwrite(as.data.table(perm_df, keep.rownames = "term"),
         file = file.path(outdir, paste0(prefix, "_beta_permanova.tsv")), sep = "\t")
  
  bd <- vegan::betadisper(bray, meta$environment)
  bd_anova   <- as.data.frame(anova(bd))
  bd_perm    <- vegan::permutest(bd, permutations = 999)
  bd_perm_df <- as.data.frame(bd_perm$tab)
  
  fwrite(as.data.table(bd_anova, keep.rownames = "term"),
         file = file.path(outdir, paste0(prefix, "_beta_betadisper_anova.tsv")), sep = "\t")
  fwrite(as.data.table(bd_perm_df, keep.rownames = "term"),
         file = file.path(outdir, paste0(prefix, "_beta_betadisper_permutest.tsv")), sep = "\t")
}

# ==========================================================
# STEP 5: Concordance scatter 
#   - Still computes/saves deltaCLR tables
#   - Tables are sorted by abs_deltaCLR 
# ==========================================================
step5_concordance <- function(clr_obj, outdir, prefix, MODE) {
  clr_tbl <- copy(clr_obj$table)
  
  clr_tbl <- clr_tbl[
    environment %in% c("Peneira", "Floresta") &
      !is.na(pairing_code) & pairing_code != "" &
      !is.na(replicate) & replicate != ""
  ]
  
  target_codes <- c("500", "1500", "2500", "3500", "4500")
  clr_tbl <- clr_tbl[pairing_code %in% target_codes]
  
  corr_dir <- file.path(outdir, paste0(prefix, "_code_concordance_", MODE))
  dir.create(corr_dir, showWarnings = FALSE, recursive = TRUE)
  
  fl <- clr_tbl[environment == "Floresta",
                .(pairing_code, replicate, genus,
                  file_floresta = file,
                  CLR_Floresta = clr)]
  
  pn <- clr_tbl[environment == "Peneira",
                .(pairing_code, replicate, genus,
                  file_peneira = file,
                  CLR_Peneira = clr)]
  
  pairs <- merge(pn, fl, by = c("pairing_code", "replicate", "genus"), allow.cartesian = TRUE)
  
  pairs[, pair_id := paste0(file_peneira, " vs ", file_floresta)]
  pairs[, floresta_partner := fifelse(grepl("^L01_", toupper(file_floresta)), "L01",
                                      fifelse(grepl("^L02_", toupper(file_floresta)), "L02", NA_character_))]
  
  pairs <- pairs[floresta_partner == "L01"]
  
  for (cc in target_codes) {
    df <- pairs[pairing_code == cc]
    if (nrow(df) == 0) next
    
    df[, deltaCLR := CLR_Peneira - CLR_Floresta]
    df[, abs_deltaCLR := abs(deltaCLR)]
    
    # sort by abs_deltaCLR before saving
    df_out <- df[order(-abs_deltaCLR, genus, replicate, file_peneira, file_floresta)]
    
    # Pearson correlation
    ct <- suppressWarnings(cor.test(df$CLR_Floresta, df$CLR_Peneira, method = "pearson"))
    r <- unname(ct$estimate); r2 <- r^2; pval <- ct$p.value
    
    p <- ggplot(df, aes(x = CLR_Floresta, y = CLR_Peneira)) +
      geom_point(alpha = 0.55, size = 1.2) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.6,
                  linetype = "dashed", color = "grey30", alpha = .5) +
      labs(
        title = paste0("Pearson R² = ", sprintf("%.2f", r2), ", p = ", signif(pval, 2)),
        x = "\nCLR (Floresta samples)",
        y = "CLR (Peneira samples)\n"
      ) +
      theme_classic(base_size = 12)
    
    ggsave(file.path(corr_dir, paste0(prefix, "_code_", cc, "_L01_scatter.png")),
           p, width = 5.2, height = 4.2, dpi = 300)
    ggsave(file.path(corr_dir, paste0(prefix, "_code_", cc, "_L01_scatter.pdf")),
           p, width = 5.2, height = 4.2)
    
    fwrite(
      df_out[, .(pairing_code, replicate, genus,
                 file_peneira, file_floresta, floresta_partner,
                 CLR_Peneira, CLR_Floresta, deltaCLR, abs_deltaCLR, pair_id)],
      file = file.path(corr_dir, paste0(prefix, "_code_", cc, "_L01_paired_dots.tsv")),
      sep  = "\t"
    )
  }
  
  message(">>> Concordance scatter plots in: ", corr_dir)
}

# ==========================================================
# STEP 7: Differential taxa (ANCOM-BC2)
#   IMPORTANT: uses the SAME numeric choice as Parts 1–4
#   i.e., USE_COUNTS_0_4 controls which column goes into ANCOM-BC2.
#   - USE_COUNTS_0_4 = 0 -> abundance
#   - USE_COUNTS_0_4 = 1 -> estimated_counts
# ==========================================================
step7_ancombc2 <- function(dt_raw, outdir, prefix, USE_COUNTS_0_4) {
  
  # Step 7 uses USE_COUNTS_0_4 (same as parts 1–4)
  value_col <- if (USE_COUNTS_0_4 == 1) "estimated_counts" else "abundance"
  message(">>> Step 7 (ANCOM-BC2) value_col = ", value_col, " (driven by USE_COUNTS_0_4)")
  
  if (!requireNamespace("ancombc", quietly = TRUE)) {
    stop("Missing R package 'ancombc' (ANCOM-BC2). Your bash wrapper should install it.")
  }
  
  # keep only Floresta vs Peneira
  dt2 <- copy(dt_raw)[environment %in% c("Floresta", "Peneira")]
  dt2 <- dt2[is.finite(get(value_col))]
  
  # aggregate to genus per sample
  gen_long <- dt2[, .(val = sum(get(value_col), na.rm = TRUE)),
                  by = .(file, environment, genus)]
  gen_long[is.na(genus) | genus == "", genus := "No genus"]
  
  # make wide count matrix
  wide <- dcast(gen_long, file + environment ~ genus, value.var = "val", fill = 0)
  meta <- wide[, .(file, environment)]
  mat  <- as.matrix(wide[, setdiff(names(wide), c("file","environment")), with = FALSE])
  rownames(mat) <- wide$file
  
  # ANCOM-BC2 expects:
  # - data: feature table (samples x taxa)
  # - meta_data: sample metadata with rownames matching sample IDs
  meta_df <- as.data.frame(meta)
  rownames(meta_df) <- meta_df$file
  meta_df$file <- NULL
  meta_df$environment <- factor(meta_df$environment, levels = c("Floresta", "Peneira"))
  
  # run ANCOM-BC2 
  res <- ancombc::ancombc2(
    data = mat,
    meta_data = meta_df,
    fix_formula = "environment",
    rand_formula = NULL,
    p_adj_method = "BH",
    prv_cut = 0.10,     # prevalence filter
    lib_cut = 0,        # no lib-size cut (we can tune later)
    group = "environment",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 1,
    verbose = FALSE
  )
  
  # extract results (primary table)
  # NOTE: ancombc2 returns lists with matrices keyed by taxa
  out <- data.table(
    genus = colnames(mat),
    lfc   = res$res$lfc[, "environmentPeneira"],
    se    = res$res$se[,  "environmentPeneira"],
    W     = res$res$W[,   "environmentPeneira"],
    p_val = res$res$p_val[, "environmentPeneira"],
    q_val = res$res$q_val[, "environmentPeneira"],
    diff_abn = res$res$diff_abn[, "environmentPeneira"]
  )
  
  setorder(out, q_val, p_val, -abs(lfc), genus)
  
  fwrite(out,
         file = file.path(outdir, paste0(prefix, "_ancombc2_genus_Floresta_vs_Peneira.tsv")),
         sep = "\t")
  
  # quick volcano-like plot
  out[, neglog10_q := -log10(pmax(q_val, 1e-300))]
  p <- ggplot(out, aes(x = lfc, y = neglog10_q)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.6) +
    labs(
      x = "Log fold-change (Peneira vs Floresta)",
      y = "-log10(q-value)",
      title = "ANCOM-BC2 (genus): Peneira vs Floresta"
    ) +
    theme_classic(base_size = 12)
  
  ggsave(file.path(outdir, paste0(prefix, "_ancombc2_genus_Floresta_vs_Peneira.png")),
         p, width = 5.5, height = 4.2, dpi = 300)
}

# ==========================================================
# Calling functions
# ==========================================================
merge_batches_if_needed(infile)

obj0 <- step0_read_stage(infile, outdir, prefix, USE_COUNTS_0_4, USE_COUNTS_5)
dt_raw <- obj0$dt_raw
VAL_0_4 <- obj0$VAL_0_4
clr_obj <- obj0$clr_obj

step1_stacked_bars(dt_raw, VAL_0_4, outdir, prefix)

obj2 <- step2_alpha(dt_raw, VAL_0_4, outdir, prefix)
rel  <- obj2$rel
meta <- obj2$meta

obj3 <- step3_beta_pcoa(rel, meta, outdir, prefix)
bray <- obj3$bray

step4_beta_stats(bray, meta, outdir, prefix)

step5_concordance(clr_obj, outdir, prefix, MODE)

step7_ancombc2(dt_raw, outdir, prefix, USE_COUNTS_0_4)

message(">>> Done. Outputs in: ", outdir)
