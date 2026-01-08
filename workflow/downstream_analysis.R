#!/usr/bin/env Rscript
# ==========================================================
# downstream_analysis.R
# 1) Add 'environment' from filename rules
# 2) Make environment-grouped stacked bars (Family & Genus)
# 3) Alpha diversity (Shannon / Simpson)
# 4) Beta diversity (Bray–Curtis PCoA, Genus)
# 5) Beta-diversity stats: PERMANOVA + betadisper
# 6) Per-code genus CLR concordance scatter plots (Peneira vs Floresta)
# 7) Differential taxa (Floresta vs Peneira) via ANCOM-BC2 (L01 only)
# 8) Differential taxa per code (Floresta CODE vs Peneira CODE) via ANCOM-BC2 (L01 only)
# 9) Heatmap (pheatmap) of ANCOM-BC2 significant genera (L01 Floresta vs Peneira)
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

#global output roots (always under results/)
RESULTS_ROOT <- dirname(outdir)  # assumes outdir is inside results/...
TABLES_DIR   <- file.path(RESULTS_ROOT, "tables")
PLOTS_DIR    <- file.path(RESULTS_ROOT, "plots")
dir.create(TABLES_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR,  showWarnings = FALSE, recursive = TRUE)
message(">>> TABLES_DIR = ", TABLES_DIR)
message(">>> PLOTS_DIR  = ", PLOTS_DIR)

### downstream mode + numeric source (passed from runall.sh)
MODE <- toupper(Sys.getenv("MODE", "16S"))
message(">>> MODE = ", MODE)
user_prefix <- ifelse(length(args) >= 3, args[3], "downstream")
prefix <- paste0(user_prefix, "_", MODE)

USE_COUNTS_0_4 <- as.integer(Sys.getenv("USE_COUNTS_0_4", "0"))  # default: abundance
USE_COUNTS_5   <- as.integer(Sys.getenv("USE_COUNTS_5",   "1"))  # default: estimated_counts

#Dedicated flag for ANCOM-BC2 input
USE_COUNTS_ANCOM <- as.integer(Sys.getenv("USE_COUNTS_ANCOM", "1"))  # default: estimated_counts

message(">>> USE_COUNTS_0_4 = ", USE_COUNTS_0_4)
message(">>> USE_COUNTS_5   = ", USE_COUNTS_5)
message(">>> USE_COUNTS_ANCOM = ", USE_COUNTS_ANCOM)  

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

#normalize sample IDs BEFORE any downstream logic
### - forces LO1/IO1/O1 etc -> L01
### - forces LO2/IO2/O2 etc -> L02
### - also uppercases and trims spaces
normalize_file_id <- function(x) {
  x <- toupper(trimws(x))
  x <- gsub("LO([0-9])", "L0\\1", x)  # LO1 -> L01
  x <- gsub("IO([0-9])", "L0\\1", x)  # IO1 -> L01
  x <- gsub("L0I", "L01", x)
  x <- gsub("L0O", "L00", x)
  x <- gsub("O([0-9])", "0\\1", x)
  x <- gsub("\\bL1_", "L01_", x)
  x <- gsub("\\bL2_", "L02_", x)
  x
}

add_code_replicate <- function(dt) {
  stopifnot("file" %in% names(dt))
  
  dt[, file_norm := toupper(file)]  # file is already normalized in read_and_shape()
  
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
  dt[file_norm %like% "^PENEIRA_0500_", pairing_code := "500"]  # (already correct)
  
  dt[, file_norm := NULL]
  dt
}

add_environment <- function(dt) {
  stopifnot("file" %in% names(dt))
  
  dt[, file_norm := toupper(file)]  # file is already normalized in read_and_shape()
  
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
  
  # keep original file, then normalize file IDs BEFORE anything else
  dt[, file_raw := file]
  dt[, file := normalize_file_id(file)]
  
  # ==========================================================
  # Quick debug + robust BAD_FILE exclusion (Option A)
  # ==========================================================
  BAD_FILE <- "L01_3050_II_ARCH_and_PENEIRA_3500_ITS.trimmed.filtered"
  BAD_FILE_NORM <- normalize_file_id(BAD_FILE)
  
  message(">>> BAD_FILE (raw)  = ", BAD_FILE)
  message(">>> BAD_FILE (norm) = ", BAD_FILE_NORM)
  message(">>> Any exact match in file?  ", any(dt$file == BAD_FILE))
  message(">>> Any substring match in file? ", any(str_detect(dt$file, fixed(BAD_FILE_NORM))))
  message(">>> Any substring match in file_raw? ", any(str_detect(toupper(dt$file_raw), fixed(toupper(BAD_FILE)))))
  
  message(">>> Examples containing '3050_II' (from file):")
  ex <- unique(dt$file[str_detect(dt$file, "3050_II")])
  if (length(ex) == 0) {
    message(">>>   (none found)")
  } else {
    print(head(ex, 10))
  }
  
  n_before <- nrow(dt)
  dt <- dt[
    !str_detect(file, fixed(BAD_FILE_NORM)) &
      !str_detect(toupper(file_raw), fixed(toupper(BAD_FILE)))
  ]
  n_after <- nrow(dt)
  message(">>> BAD_FILE exclusion removed rows: ", n_before - n_after)
  # ==========================================================
  
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
step0_read_stage <- function(infile, outdir, prefix, USE_COUNTS_0_4, USE_COUNTS_5, TABLES_DIR) {
  read_obj <- read_and_shape(infile)
  dt_raw   <- read_obj$raw
  
  VAL_0_4 <- if (USE_COUNTS_0_4 == 1) "estimated_counts" else "abundance"
  message(">>> Parts 1–4 value_col = ", VAL_0_4)
  
  clr_obj <- make_genus_clr_steps(dt_raw, use_counts = USE_COUNTS_5)
  
  fwrite(
    clr_obj$table,
    file = file.path(TABLES_DIR, paste0(prefix, "_genus_clr_steps.tsv")),
    sep  = "\t"
  )
  
  dt_stage <- copy(dt_raw)
  if ("value" %in% names(dt_stage)) dt_stage[, value := NULL]
  
  fwrite(
    dt_stage,
    file = file.path(TABLES_DIR, paste0(prefix, "_with_environment.tsv")),
    sep  = "\t"
  )
  
  list(dt_raw = dt_raw, VAL_0_4 = VAL_0_4, clr_obj = clr_obj)
}

# ==========================================================
# STEP 1: Stacked bars
# ==========================================================
step1_stacked_bars <- function(dt_raw, VAL_0_4, outdir, prefix, PLOTS_DIR) {
  p_family <- make_env_stacks(
    dt_raw, "family",
    file.path(PLOTS_DIR, paste0(prefix, "_stacks_family.png")),
    file.path(PLOTS_DIR, paste0(prefix, "_stacks_family.pdf")),
    N = 20, title_rank = "Family", value_col = VAL_0_4
  )
  p_genus <- make_env_stacks(
    dt_raw, "genus",
    file.path(PLOTS_DIR, paste0(prefix, "_stacks_genus.png")),
    file.path(PLOTS_DIR, paste0(prefix, "_stacks_genus.pdf")),
    N = 20, title_rank = "Genus", value_col = VAL_0_4
  )
  
  if (!is.null(p_family) && !is.null(p_genus)) {
    combined <- plot_grid(p_family, p_genus, ncol = 1, rel_heights = c(1, 1), align = "v")
    ggsave(file.path(PLOTS_DIR, paste0(prefix, "_stacks_family_genus_grid.png")),
           combined, width = 22, height = 10, dpi = 300)
    ggsave(file.path(PLOTS_DIR, paste0(prefix, "_stacks_family_genus_grid.pdf")),
           combined, width = 22, height = 10)
  }
}

# ==========================================================
# STEP 2: Alpha diversity (Shannon & Simpson only)
# ==========================================================
step2_alpha <- function(dt_raw, VAL_0_4, outdir, prefix, TABLES_DIR, PLOTS_DIR) {
  mx_rel  <- build_matrix(dt_raw, "genus", value_col = VAL_0_4)
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
  
  alpha_df[, environment := factor(environment, levels = c("Campina", "Floresta", "Igarape", "Peneira"))]
  
  data.table::fwrite(alpha_df, file = file.path(TABLES_DIR, paste0(prefix, "_alpha_diversity.tsv")), sep  = "\t")
  
  # ---- pairwise tests + plots (unchanged) ----
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
      p <- tryCatch(wilcox.test(x, y)$p.value, error = function(e) NA_real_)
      data.table(env1 = e1, env2 = e2, metric = metric_name, p_value = p)
    })
    rbindlist(res_list)
  }
  
  pval_to_stars <- function(p) {
    ifelse(p < 0.001, "***",
           ifelse(p < 0.01, "**",
                  ifelse(p < 0.05, "*", "ns")))
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
  
  fwrite(pw_all, file = file.path(TABLES_DIR, paste0(prefix, "_alpha_pairwise_wilcox.tsv")), sep = "\t")
  
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
    geom_quasirandom(shape = 21, size = 1, dodge.width = 0.75, alpha = 0.5, color = "black", show.legend = FALSE) +
    geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.9, color = "black", fill = "white", show.legend = FALSE) +
    labs(x = "\nEnvironment", y = "Shannon (H')\n", title = "") +
    scale_fill_manual(values = env_colors) +
    scale_color_manual(values = env_colors) +
    theme_base
  
  if (!is.null(sig_sh_df) && nrow(sig_sh_df) > 0) {
    p_sh <- p_sh +
      geom_segment(data = sig_sh_df, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE) +
      geom_text(data = sig_sh_df, aes(x = (x + xend) / 2, y = y, label = label),
                vjust = -0.3, size = 3, inherit.aes = FALSE)
  }
  
  ggsave(file.path(PLOTS_DIR, paste0(prefix, "_alpha_shannon_env.png")), p_sh, width = 3, height = 5, dpi = 300)
  
  p_sp <- ggplot(alpha_df, aes(x = environment, y = Simpson, fill = environment, color = environment)) +
    geom_violin(alpha = 0.25, linewidth = 0, position = position_dodge(width = 0.75), show.legend = FALSE) +
    geom_quasirandom(shape = 21, size = 1, dodge.width = 0.75, alpha = 0.5, color = "black", show.legend = FALSE) +
    geom_boxplot(outlier.shape = NA, width = 0.3, alpha = 0.9, color = "black", fill = "white", show.legend = FALSE) +
    labs(x = "\nEnvironment", y = "Simpson (1 - D)\n", title = "") +
    scale_fill_manual(values = env_colors) +
    scale_color_manual(values = env_colors) +
    theme_base
  
  if (!is.null(sig_sp_df) && nrow(sig_sp_df) > 0) {
    p_sp <- p_sp +
      geom_segment(data = sig_sp_df, aes(x = x, xend = xend, y = y, yend = y), inherit.aes = FALSE) +
      geom_text(data = sig_sp_df, aes(x = (x + xend) / 2, y = y, label = label),
                vjust = -0.3, size = 3, inherit.aes = FALSE)
  }
  
  ggsave(file.path(PLOTS_DIR, paste0(prefix, "_alpha_simpson_env.png")), p_sp, width = 3, height = 5, dpi = 300)
  
  list(rel = rel, meta = meta, alpha_df = alpha_df)
}

# ==========================================================
# STEP 3: Beta diversity PCoA
# ==========================================================
step3_beta_pcoa <- function(rel, meta, outdir, prefix, TABLES_DIR, PLOTS_DIR) {
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
  
  fwrite(pcoa_df, file = file.path(TABLES_DIR, paste0(prefix, "_pcoa_braycurtis.tsv")), sep = "\t")
  
  p_pcoa <- ggplot(pcoa_df, aes(x = PC1, y = PC2, color = environment)) +
    geom_point(size = 2.5, alpha = 0.9) +
    scale_color_manual(values = env_colors, na.value = "grey70") +
    labs(x = pc1_lab, y = pc2_lab, title = "") +
    theme_base
  
  ggsave(file.path(PLOTS_DIR, paste0(prefix, "_pcoa_braycurtis_env.png")), p_pcoa, width = 5, height = 4, dpi = 300)
  
  list(bray = bray)
}

# ==========================================================
# STEP 4: Beta stats
# ==========================================================
step4_beta_stats <- function(bray, meta, outdir, prefix, TABLES_DIR) {
  set.seed(2025)
  meta$environment <- factor(meta$environment)
  
  perm <- vegan::adonis2(bray ~ environment, data = meta, permutations = 999)
  perm_df <- as.data.frame(perm)
  
  fwrite(as.data.table(perm_df, keep.rownames = "term"),
         file = file.path(TABLES_DIR, paste0(prefix, "_beta_permanova.tsv")), sep = "\t")
  
  bd <- vegan::betadisper(bray, meta$environment)
  bd_anova   <- as.data.frame(anova(bd))
  bd_perm    <- vegan::permutest(bd, permutations = 999)
  bd_perm_df <- as.data.frame(bd_perm$tab)
  
  fwrite(as.data.table(bd_anova, keep.rownames = "term"),
         file = file.path(TABLES_DIR, paste0(prefix, "_beta_betadisper_anova.tsv")), sep = "\t")
  
  fwrite(as.data.table(bd_perm_df, keep.rownames = "term"),
         file = file.path(TABLES_DIR, paste0(prefix, "_beta_betadisper_permutest.tsv")), sep = "\t")
}

# ==========================================================
# STEP 5: Concordance scatter
# ==========================================================
step5_concordance <- function(clr_obj, outdir, prefix, MODE, TABLES_DIR, PLOTS_DIR) {
  clr_tbl <- copy(clr_obj$table)
  
  clr_tbl <- clr_tbl[
    environment %in% c("Peneira", "Floresta") &
      !is.na(pairing_code) & pairing_code != "" &
      !is.na(replicate) & replicate != ""
  ]
  
  target_codes <- c("500", "1500", "2500", "3500", "4500")
  clr_tbl <- clr_tbl[pairing_code %in% target_codes]
  
  corr_dir <- file.path(PLOTS_DIR, paste0(prefix, "_code_concordance_", MODE))
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
    
    df_out <- df[order(-abs_deltaCLR, genus, replicate, file_peneira, file_floresta)]
    
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
      file = file.path(TABLES_DIR, paste0(prefix, "_code_", cc, "_L01_paired_dots.tsv")),
      sep  = "\t"
    )
  }
  
  message(">>> Concordance scatter plots in: ", corr_dir)
}

# ==========================================================
# STEP 6: pooled (all codes) genus CLR concordance scatter
#   - each dot = (pairing_code x replicate x genus)
#   - x: Floresta CLR (L01 only)
#   - y: Peneira  CLR (PENEIRA_* only)
#   - pools all codes together in one scatter
# ==========================================================
step6c_allcodes_clr_concordance <- function(
    clr_obj, outdir, prefix, MODE,
    TABLES_DIR, PLOTS_DIR,
    enforce_L01_floresta = TRUE,
    enforce_peneira_only = TRUE,
    label_top_n = 0,              # set >0 if you want to label most discordant points
    color_by = c("pairing_code", "replicate", "none"),
    point_alpha = 0.55,
    point_size  = 1.1
) {
  color_by <- match.arg(color_by)
  
  if (is.null(clr_obj) || is.null(clr_obj$table)) {
    stop("step6c_allcodes_clr_concordance: clr_obj$table missing.")
  }
  
  dt <- data.table::as.data.table(clr_obj$table)
  
  req <- c("file","environment","genus","clr","pairing_code","replicate")
  miss <- setdiff(req, names(dt))
  if (length(miss) > 0) stop("step6c_allcodes_clr_concordance: missing columns: ", paste(miss, collapse = ", "))
  
  # keep only target groups + finite clr
  dt <- dt[environment %in% c("Floresta","Peneira")]
  dt <- dt[is.finite(clr)]
  dt[is.na(genus) | genus == "", genus := "No genus"]
  dt[, pairing_code := as.character(pairing_code)]
  dt[, replicate := as.character(replicate)]
  
  # enforce the same sample inclusion logic you used elsewhere
  if (enforce_L01_floresta) {
    dt <- dt[
      (environment == "Floresta" & grepl("^L01_", file)) |
        (environment == "Peneira")
    ]
  }
  if (enforce_peneira_only) {
    dt <- dt[
      (environment == "Floresta") |
        (environment == "Peneira" & grepl("^PENEIRA_", file))
    ]
  }
  
  # split and pair
  flo <- dt[environment == "Floresta",
            .(pairing_code, replicate, genus, clr_flo = clr, file_flo = file)]
  pen <- dt[environment == "Peneira",
            .(pairing_code, replicate, genus, clr_pen = clr, file_pen = file)]
  
  paired <- merge(
    flo, pen,
    by = c("pairing_code","replicate","genus"),
    all = FALSE
  )
  
  if (nrow(paired) < 50) {
    message(">>> Step 6C: too few paired points (n=", nrow(paired), "); skipping.")
    return(invisible(NULL))
  }
  
  # correlations (pooled)
  pear  <- suppressWarnings(stats::cor(paired$clr_flo, paired$clr_pen, method = "pearson"))
  spear <- suppressWarnings(stats::cor(paired$clr_flo, paired$clr_pen, method = "spearman"))
  
  # optional labeling: biggest absolute deviation from y=x
  paired[, delta := clr_pen - clr_flo]
  paired[, abs_delta := abs(delta)]
  data.table::setorder(paired, -abs_delta)
  
  # plotting mapping
  if (color_by == "pairing_code") {
    aes_col <- ggplot2::aes(color = pairing_code)
  } else if (color_by == "replicate") {
    aes_col <- ggplot2::aes(color = replicate)
  } else {
    aes_col <- ggplot2::aes()
  }
  
  p <- ggplot2::ggplot(paired, ggplot2::aes(x = clr_flo, y = clr_pen)) +
    ggplot2::geom_point(aes_col, alpha = point_alpha, size = point_size) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed", alpha = 0.6) +
    ggplot2::geom_vline(xintercept = 0, linetype = "dotted", alpha = 0.35) +
    ggplot2::geom_hline(yintercept = 0, linetype = "dotted", alpha = 0.35) +
    ggplot2::labs(
      title = "Genus CLR concordance (pooled across all codes): Floresta vs Peneira",
      subtitle = paste0(
        "Each dot = (code × replicate × genus). ",
        "Pearson r=", signif(pear, 3), " | Spearman ρ=", signif(spear, 3),
        " | MODE=", MODE, " | source=", clr_obj$source, " | pc=", clr_obj$pseudocount
      ),
      x = "Floresta CLR (per sample, per genus)",
      y = "Peneira CLR (paired sample, per genus)",
      color = if (color_by == "none") NULL else color_by
    ) +
    ggplot2::theme_classic(base_size = 12) +
    ggplot2::theme(legend.position = if (color_by == "none") "none" else "right")
  
  # optional point labels
  if (label_top_n > 0) {
    lab <- paired[1:min(label_top_n, .N)]
    lab[, label := paste0(genus, " | code=", pairing_code, " | rep=", replicate)]
    if (requireNamespace("ggrepel", quietly = TRUE)) {
      p <- p + ggrepel::geom_text_repel(
        data = lab,
        ggplot2::aes(label = label),
        size = 2.8,
        max.overlaps = Inf,
        box.padding = 0.25,
        point.padding = 0.2
      )
    } else {
      p <- p + ggplot2::geom_text(
        data = lab,
        ggplot2::aes(label = label),
        size = 2.6,
        vjust = -0.6
      )
    }
  }
  
  # save
  out_base <- file.path(PLOTS_DIR, paste0(prefix, "_allcodes_CLR_concordance"))
  ggplot2::ggsave(paste0(out_base, ".png"), p, width = 6.6, height = 6.0, dpi = 300)
  ggplot2::ggsave(paste0(out_base, ".pdf"), p, width = 6.6, height = 6.0)
  
  # save table (useful for debugging)
  data.table::fwrite(
    paired,
    file = file.path(TABLES_DIR, paste0(prefix, "_allcodes_CLR_concordance_points.tsv")),
    sep = "\t"
  )
  
  message(">>> Step 6C: wrote pooled all-codes CLR concordance plot to: ", out_base, ".{png,pdf}")
  invisible(list(points = paired, plot = p, pearson = pear, spearman = spear))
}

# ==========================================================
# STEP 7: Differential taxa (ANCOM-BC2) via PHYLOSEQ
#   #L01-only Floresta, keep Peneira; exclude L02 entirely
# ==========================================================
#Signature now uses USE_COUNTS_ANCOM
step7_ancombc2 <- function(dt_raw, outdir, prefix, USE_COUNTS_ANCOM, TABLES_DIR, PLOTS_DIR) {
  value_col <- if (USE_COUNTS_ANCOM == 1) "estimated_counts" else "abundance"
  message(">>> Step 7 (ANCOM-BC2) value_col = ", value_col, " (driven by USE_COUNTS_ANCOM)")
  
  suppressPackageStartupMessages({
    if (!requireNamespace("ANCOMBC", quietly = TRUE) && !requireNamespace("ancombc", quietly = TRUE)) {
      stop("Missing R package ANCOMBC/ancombc.")
    }
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      stop("Missing R package 'phyloseq'.")
    }
    library(phyloseq)
  })
  
  pkg <- if (requireNamespace("ANCOMBC", quietly = TRUE)) "ANCOMBC" else "ancombc"
  message(">>> Using ANCOM-BC2 package namespace: ", pkg)
  
  # keep only Floresta vs Peneira + numeric finite
  dt2 <- data.table::copy(dt_raw)[environment %in% c("Floresta", "Peneira")]
  dt2 <- dt2[is.finite(get(value_col))]
  
  #enforce L01-only for Floresta; keep all Peneira; drop all L02
  dt2 <- dt2[
    (environment == "Floresta" & grepl("^L01_", file)) |
      (environment == "Peneira"  & grepl("^PENEIRA_", file))
  ]
  message(">>> Step 7: after L01-only Floresta filter: samples=", length(unique(dt2$file)))
  
  # aggregate to genus per sample
  gen_long <- dt2[, .(val = sum(get(value_col), na.rm = TRUE)),
                  by = .(file, environment, genus)]
  gen_long[is.na(genus) | genus == "", genus := "No genus"]
  
  # wide sample x genus
  wide <- data.table::dcast(gen_long, file + environment ~ genus, value.var = "val", fill = 0)
  
  # ----- metadata -----
  meta_df <- as.data.frame(wide[, .(file, environment)])
  colnames(meta_df)[colnames(meta_df) == "environment"] <- "env_group"
  
  meta_df$env_group <- factor(meta_df$env_group, levels = c("Floresta", "Peneira"))
  
  rownames(meta_df) <- meta_df$file
  meta_df$file <- NULL
  
  # ----- OTU table (taxa x samples) -----
  mat <- as.matrix(wide[, setdiff(names(wide), c("file", "environment")), with = FALSE])
  rownames(mat) <- wide$file
  otu <- t(mat)  # taxa x samples
  storage.mode(otu) <- "numeric"
  
  common <- intersect(colnames(otu), rownames(meta_df))
  if (length(common) < 2) stop(">>> Step 7 ERROR: too few common samples between otu and meta.")
  otu <- otu[, common, drop = FALSE]
  meta_df <- meta_df[common, , drop = FALSE]
  
  ps <- phyloseq(
    otu_table(otu, taxa_are_rows = TRUE),
    sample_data(meta_df)
  )
  
  message(">>> Step 7: phyloseq built. taxa=", ntaxa(ps), " samples=", nsamples(ps))
  message(">>> Step 7: sample vars in phyloseq = ", paste(colnames(as(sample_data(ps), "data.frame")), collapse = ", "))
  
  ancombc2_fun <- get("ancombc2", envir = asNamespace(pkg))
  
  res <- ancombc2_fun(
    data = ps,
    fix_formula = "env_group",
    rand_formula = NULL,
    p_adj_method = "BH",
    prv_cut = 0.10,
    lib_cut = 0,
    group = "env_group",
    struc_zero = TRUE,
    neg_lb = TRUE,
    alpha = 0.05,
    n_cl = 1,
    verbose = TRUE
  )
  
  save_prefix <- paste0(prefix, "_ancombc2_raw")
  saveRDS(res, file = file.path(TABLES_DIR, paste0(save_prefix, ".rds")))
  message(">>> Saved ANCOMBC2 res object to: ", file.path(TABLES_DIR, paste0(save_prefix, ".rds")))
  
  res_tbl <- if (!is.null(res$res)) res$res else res
  
  if (!is.null(res_tbl)) {
    res_dt <- data.table::as.data.table(res_tbl)
    
    if (!"taxon" %in% names(res_dt)) {
      if (nrow(res_dt) > 0 && !is.null(rownames(res_tbl))) {
        res_dt[, taxon := rownames(res_tbl)]
      } else {
        setnames(res_dt, 1, "taxon")
      }
    }
    
    fwrite(res_dt, file = file.path(TABLES_DIR, paste0(prefix, "_ancombc2_full_table.tsv")), sep = "\t")
    
    lvl <- levels(meta_df$env_group)
    baseline_level <- lvl[1]
    other_level    <- lvl[2]
    message(">>> Step 7: baseline = ", baseline_level, " ; coefficient level = ", other_level)
    
    esc_other <- gsub("([^A-Za-z0-9])", "\\\\\\1", other_level)
    pat_lfc   <- paste0("^lfc_.*", esc_other, "$")
    pat_se    <- paste0("^se_.*",  esc_other, "$")
    pat_W     <- paste0("^W_.*",   esc_other, "$")
    pat_p     <- paste0("^p_.*",   esc_other, "$")
    pat_q     <- paste0("^q_.*",   esc_other, "$")
    pat_diff  <- paste0("^diff_.*", esc_other, "$")
    pat_dr    <- paste0("^diff_robust_.*", esc_other, "$")
    pat_pss   <- paste0("^passed_ss_.*", esc_other, "$")
    
    pick1 <- function(pat, cols) {
      hit <- grep(pat, cols, value = TRUE)
      if (length(hit) == 0) return(NA_character_)
      hit[1]
    }
    
    cols <- names(res_dt)
    col_lfc  <- pick1(pat_lfc, cols)
    col_se   <- pick1(pat_se, cols)
    col_W    <- pick1(pat_W, cols)
    col_p    <- pick1(pat_p, cols)
    col_q    <- pick1(pat_q, cols)
    col_diff <- pick1(pat_diff, cols)
    col_dr   <- pick1(pat_dr, cols)
    col_pss  <- pick1(pat_pss, cols)
    
    message(">>> Step 7: picked columns:")
    message("    lfc  = ", col_lfc)
    message("    se   = ", col_se)
    message("    W    = ", col_W)
    message("    p    = ", col_p)
    message("    q    = ", col_q)
    message("    diff = ", col_diff)
    message("    diff_robust = ", col_dr)
    message("    passed_ss   = ", col_pss)
    
    if (is.na(col_lfc) || is.na(col_q)) {
      stop("ANCOM-BC2 table does not contain expected lfc/q columns for the non-baseline level. See *_ancombc2_full_table.tsv in TABLES_DIR.")
    }
    
    out <- data.table(
      genus = res_dt$taxon,
      lfc   = as.numeric(res_dt[[col_lfc]]),
      se    = if (!is.na(col_se))   as.numeric(res_dt[[col_se]])   else NA_real_,
      W     = if (!is.na(col_W))    as.numeric(res_dt[[col_W]])    else NA_real_,
      p_val = if (!is.na(col_p))    as.numeric(res_dt[[col_p]])    else NA_real_,
      q_val = as.numeric(res_dt[[col_q]]),
      diff_abn     = if (!is.na(col_diff)) as.integer(res_dt[[col_diff]]) else NA_integer_,
      diff_robust  = if (!is.na(col_dr))   as.integer(res_dt[[col_dr]])   else NA_integer_,
      passed_ss    = if (!is.na(col_pss))  as.integer(res_dt[[col_pss]])  else NA_integer_
    )
    
    out[, lfc_direction := ifelse(lfc > 0,
                                  paste0("Higher in ", other_level),
                                  paste0("Higher in ", baseline_level))]
    
    out[, abs_lfc := abs(lfc)]
    data.table::setorder(out, q_val, p_val, -abs_lfc, genus)
    
    out_tsv <- file.path(TABLES_DIR, paste0(prefix, "_ancombc2_genus_", baseline_level, "_vs_", other_level, ".tsv"))
    fwrite(out, file = out_tsv, sep = "\t")
    message(">>> Step 7: wrote ANCOM-BC2 summary table: ", out_tsv)
    
    out[, neglog10_q := -log10(pmax(q_val, 1e-300))]
    out[, significance := ifelse(q_val < 0.05, "q < 0.05", "ns")]
    
    p <- ggplot(out, aes(x = lfc, y = neglog10_q, color = significance)) +
      geom_point(alpha = 0.6, size = 1.2) +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.6) +
      geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
      scale_color_manual(values = c("q < 0.05" = "red", "ns" = "grey60")) +
      labs(
        x = paste0("Log fold-change (", other_level, " vs ", baseline_level, ")"),
        y = "-log10(q-value)",
        title = paste0("ANCOM-BC2 (genus): ", other_level, " vs ", baseline_level),
        color = "Significance"
      ) +
      theme_classic(base_size = 12) +
      theme(legend.position = "bottom")
    
    out_png <- file.path(PLOTS_DIR, paste0(prefix, "_ancombc2_genus_", baseline_level, "_vs_", other_level, ".png"))
    ggsave(out_png, p, width = 6, height = 5, dpi = 300)
    
    message(">>> Step 7: ANCOM-BC2 completed successfully. Found ",
            sum(out$q_val < 0.05, na.rm = TRUE), " significant genera at q < 0.05.")
    return(invisible(out))
  }
  
  stop("Step 7: Could not locate a single ANCOM-BC2 results table in the returned object. (You should still have the saved .rds in TABLES_DIR to inspect.)")
}

# ==========================================================
# STEP 8: Differential taxa per pairing_code (ANCOM-BC2)
#   #L01-only Floresta, keep Peneira; drop L02
# ==========================================================
#Signature now uses USE_COUNTS_ANCOM
step8_ancombc2_by_code <- function(dt_raw, outdir, prefix, USE_COUNTS_ANCOM, TABLES_DIR, PLOTS_DIR,
                                   target_codes = c("500", "1500", "2500", "3500", "4500")) {
  
  if (!"pairing_code" %in% names(dt_raw)) {
    stop(">>> Step 8 ERROR: dt_raw has no 'pairing_code' column. (Check add_code_replicate())")
  }
  
  dt_raw[, pairing_code := as.character(pairing_code)]
  
  message(">>> Step 8: Running ANCOM-BC2 per code: ", paste(target_codes, collapse = ", "))
  
  for (cc in target_codes) {
    
    dt_cc <- data.table::copy(dt_raw)[environment %in% c("Floresta", "Peneira") & pairing_code == cc]
    
    #enforce L01-only for Floresta; keep all Peneira; drop all L02
    dt_cc <- dt_cc[
      (environment == "Floresta" & grepl("^L01_", file)) |
        (environment == "Peneira"  & grepl("^PENEIRA_", file))
    ]
    
    n_samp <- length(unique(dt_cc$file))
    n_env  <- length(unique(dt_cc$environment))
    if (n_samp < 4 || n_env < 2) {
      message(">>> Step 8: code=", cc, " skipped (samples=", n_samp, ", env_groups=", n_env, ")")
      next
    }
    
    prefix_cc <- paste0(prefix, "_code_", cc)
    
    message(">>> Step 8: code=", cc, " -> running ANCOM-BC2 with prefix: ", prefix_cc)
    
    step7_ancombc2(
      dt_raw = dt_cc,
      outdir = outdir,
      prefix = prefix_cc,
      USE_COUNTS_ANCOM = USE_COUNTS_ANCOM,   
      TABLES_DIR = TABLES_DIR,
      PLOTS_DIR = PLOTS_DIR
    )
  }
  
  message(">>> Step 8: done.")
}

# ==========================================================
# STEP 9: Heatmap (pheatmap) of ANCOM-BC2 significant genera
# ==========================================================

step9_heatmap_ancom_sig <- function(dt_raw, ancom_out, outdir, prefix, USE_COUNTS_ANCOM, TABLES_DIR, PLOTS_DIR,
                                    q_cut = 0.05, top_n = 40, pseudocount = 1) {
  
  # Install/load pheatmap
  if (!requireNamespace("pheatmap", quietly = TRUE)) {
    message(">>> Installing CRAN package 'pheatmap'...")
    install.packages("pheatmap", repos = "https://cloud.r-project.org")
  }
  suppressPackageStartupMessages(library(pheatmap))
  
  if (is.null(ancom_out) || nrow(ancom_out) == 0) {
    message(">>> Step 9: ancom_out is empty; skipping heatmap.")
    return(invisible(NULL))
  }
  
  value_col <- if (USE_COUNTS_ANCOM == 1) "estimated_counts" else "abundance"
  message(">>> Step 9 (heatmap) value_col = ", value_col)
  
  sig <- data.table::as.data.table(ancom_out)
  sig <- sig[is.finite(q_val) & q_val < q_cut]
  if (nrow(sig) == 0) {
    message(">>> Step 9: No genera pass q < ", q_cut, "; skipping heatmap.")
    return(invisible(NULL))
  }
  
  data.table::setorder(sig, q_val, -abs_lfc, genus)
  if (!is.null(top_n) && nrow(sig) > top_n) sig <- sig[1:top_n]
  
  keep_genera <- unique(sig$genus)
  
  # Filter to the same comparison set used in Step 7
  dt2 <- data.table::copy(dt_raw)[environment %in% c("Floresta", "Peneira")]
  dt2 <- dt2[is.finite(get(value_col))]
  dt2 <- dt2[
    (environment == "Floresta" & grepl("^L01_", file)) |
      (environment == "Peneira"  & grepl("^PENEIRA_", file))
  ]
  
  gen_long <- dt2[, .(val = sum(get(value_col), na.rm = TRUE)),
                  by = .(file, environment, code, pairing_code, replicate, genus)]
  gen_long[is.na(genus) | genus == "", genus := "No genus"]
  gen_long <- gen_long[genus %in% keep_genera]
  
  if (nrow(gen_long) == 0) {
    message(">>> Step 9: No data for selected genera after filtering; skipping heatmap.")
    return(invisible(NULL))
  }
  
  wide <- data.table::dcast(gen_long, file + environment + code + pairing_code + replicate ~ genus,
                            value.var = "val", fill = 0)
  
  meta <- as.data.frame(wide[, .(file, environment, code, pairing_code, replicate)])
  rownames(meta) <- meta$file
  
  mat_samp_x_gen <- as.matrix(wide[, setdiff(names(wide), c("file","environment","code","pairing_code","replicate")), with = FALSE])
  rownames(mat_samp_x_gen) <- wide$file
  
  # CLR per sample: log(count + pc) - mean(log(count + pc))
  logm <- log(mat_samp_x_gen + pseudocount)
  clr  <- logm - rowMeans(logm)
  
  # pheatmap expects rows=features, cols=samples (typical)
  hm <- t(clr)  # genera x samples
  
  # Order samples: Floresta then Peneira; within by pairing_code/code/replicate if available
  meta$environment <- factor(meta$environment, levels = c("Floresta", "Peneira"))
  meta$pairing_code <- as.character(meta$pairing_code)
  meta$code <- as.character(meta$code)
  meta$replicate <- as.character(meta$replicate)
  
  ord <- order(meta$environment, meta$pairing_code, meta$code, meta$replicate, rownames(meta))
  meta_ord <- meta[ord, , drop = FALSE]
  hm <- hm[, rownames(meta_ord), drop = FALSE]
  
  ann_col <- meta_ord[, c("environment","pairing_code","replicate"), drop = FALSE]
  colnames(ann_col) <- c("Environment", "Code", "Replicate")
  
  out_png <- file.path(PLOTS_DIR, paste0(prefix, "_heatmap_ancom_sig_clr.png"))
  out_pdf <- file.path(PLOTS_DIR, paste0(prefix, "_heatmap_ancom_sig_clr.pdf"))
  
  # Save PNG
  png(out_png, width = 1600, height = 1200, res = 150)
  pheatmap::pheatmap(
    hm,
    annotation_col = ann_col,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    fontsize_row = 8,
    fontsize_col = 6,
    main = paste0("ANCOM-BC2 significant genera (q<", q_cut, "), CLR heatmap (", value_col, ")")
  )
  dev.off()
  
  # Save PDF
  pdf(out_pdf, width = 12, height = 9)
  pheatmap::pheatmap(
    hm,
    annotation_col = ann_col,
    cluster_rows = TRUE,
    cluster_cols = FALSE,
    scale = "row",
    fontsize_row = 8,
    fontsize_col = 6,
    main = paste0("ANCOM-BC2 significant genera (q<", q_cut, "), CLR heatmap (", value_col, ")")
  )
  dev.off()
  
  # Also export the underlying matrix (for reproducibility)
  fwrite(
    data.table::as.data.table(hm, keep.rownames = "genus"),
    file = file.path(TABLES_DIR, paste0(prefix, "_heatmap_matrix_clr.tsv")),
    sep = "\t"
  )
  fwrite(
    data.table::as.data.table(meta_ord, keep.rownames = "file"),
    file = file.path(TABLES_DIR, paste0(prefix, "_heatmap_sample_annotation.tsv")),
    sep = "\t"
  )
  
  message(">>> Step 9: Heatmap written to: ", out_png, " and ", out_pdf)
  invisible(list(matrix = hm, annotation = meta_ord, sig = sig))
}

### OpenTree helpers (used in Step 10)

filter_otl_matches <- function(mdt) {
  # expects columns: search_string, ott_id; optionally flags
  stopifnot(all(c("search_string","ott_id") %in% names(mdt)))
  
  mdt <- mdt[!is.na(ott_id)]
  
  # Drop known-bad flags if present
  if ("flags" %in% names(mdt)) {
    bad_pat <- "(pruned_ott_id|deprecated_ott_id|forwarded_ott_id|blocked|taxon_incomplete)"
    mdt <- mdt[
      is.na(flags) | flags == "" | !grepl(bad_pat, flags, ignore.case = TRUE)
    ]
  }
  
  # Keep first/best hit per search_string (you already do this)
  mdt <- mdt[, .SD[1], by = search_string]
  
  mdt
}

safe_tol_induced_subtree <- function(ott_ids, label_format = "name") {
  ott_ids <- unique(ott_ids)
  
  repeat {
    if (length(ott_ids) < 2) {
      stop("Too few OTT IDs to build induced subtree after filtering.")
    }
    
    tr <- tryCatch(
      rotl::tol_induced_subtree(ott_ids = ott_ids, label_format = label_format),
      error = function(e) e
    )
    
    if (!inherits(tr, "error")) return(tr)
    
    msg <- conditionMessage(tr)
    
    # Try to extract an ott id from the error message and drop it, then retry
    hit <- regmatches(msg, regexec("ott[0-9]+", msg))[[1]]
    if (length(hit) >= 1) {
      bad <- hit[1]
      bad_int <- suppressWarnings(as.integer(gsub("^ott", "", bad)))
      message(">>> Step 10: OpenTree error for ", bad, " ; dropping and retrying.")
      ott_ids <- ott_ids[ott_ids != bad_int]
      next
    }
    
    # Unknown error: stop
    stop(tr)
  }
}
### -----------------------------------------------------------------------

# ==========================================================
# STEP 10: Real phylogeny tree (Open Tree of Life) + tip icons
#   - builds an induced subtree from OTL using genus names
#   - tip icons:
#       * genera selected from ANCOM-BC2 "res" summary: color by LFC (continuous)
#       * genera selected from ANCOM-BC2 "zero_ind": binary presence/absence category
# ==========================================================
step10_otl_phylogeny_lfc_zeroind <- function(
    res7_summary,              # data.table from step7_ancombc2() return (genus, lfc, q_val, abs_lfc, ...)
    ancombc2_raw_rds,          # path to prefix_ancombc2_raw.rds saved by step7_ancombc2
    outdir, prefix,
    PLOTS_DIR, TABLES_DIR,
    q_cut = 0.05,
    top_n_each_dir = 10,
    zero_n_each_dir = 10,
    seed = 1,
    use_only_qsig_for_top = TRUE,   # if TRUE, top lists are q<q_cut; if none pass, falls back to smallest q
    min_unique_tips = 8
) {
  #No runtime installs on HPC; require packages instead
  if (!requireNamespace("rotl", quietly = TRUE)) {
    stop("Missing R package 'rotl'. Install it in the downstream env before running Step 10.")
  }
  if (!requireNamespace("ape", quietly = TRUE)) {
    stop("Missing R package 'ape'. Install it in the downstream env before running Step 10.")
  }
  suppressPackageStartupMessages({
    library(rotl)
    library(ape)
  })
  
  if (is.null(res7_summary) || nrow(res7_summary) == 0) {
    message(">>> Step 10: res7_summary is empty; skipping phylogeny plot.")
    return(invisible(NULL))
  }
  if (!file.exists(ancombc2_raw_rds)) {
    stop(">>> Step 10: ancombc2_raw_rds not found: ", ancombc2_raw_rds)
  }
  
  # ---------- 1) pick genera from "res" summary ----------
  res_dt <- data.table::as.data.table(res7_summary)
  res_dt <- res_dt[is.finite(lfc) & is.finite(q_val)]
  if (!"abs_lfc" %in% names(res_dt)) res_dt[, abs_lfc := abs(lfc)]
  
  # Prefer q-significant for top lists, but fall back if none
  if (use_only_qsig_for_top) {
    sig_only <- res_dt[q_val < q_cut]
    if (nrow(sig_only) > 0) {
      base_for_top <- sig_only
    } else {
      base_for_top <- res_dt  # fallback: smallest q
    }
  } else {
    base_for_top <- res_dt
  }
  
  top_peneira <- base_for_top[lfc > 0][order(q_val, -abs(lfc))][1:min(top_n_each_dir, .N)]
  top_floresta <- base_for_top[lfc < 0][order(q_val, -abs(lfc))][1:min(top_n_each_dir, .N)]
  
  top_res <- data.table::rbindlist(list(top_peneira, top_floresta), use.names = TRUE, fill = TRUE)
  top_res <- unique(top_res, by = "genus")
  top_res[, source := "res"]
  top_res[, icon_type := "lfc"]  # continuous
  
  # ---------- 2) pull "zero_ind" from raw ANCOM-BC2 object ----------
  raw_obj <- readRDS(ancombc2_raw_rds)
  
  # Try multiple likely locations (depends on ANCOMBC/ancombc version)
  zero_tbl <- NULL
  if (is.list(raw_obj)) {
    if (!is.null(raw_obj$zero_ind)) zero_tbl <- raw_obj$zero_ind
    if (is.null(zero_tbl) && !is.null(raw_obj$res) && is.list(raw_obj$res) && !is.null(raw_obj$res$zero_ind)) zero_tbl <- raw_obj$res$zero_ind
  }
  
  if (is.null(zero_tbl)) {
    message(">>> Step 10: could not find 'zero_ind' inside the ANCOM-BC2 raw object. Skipping zero_ind sampling.")
    zero_dt <- data.table::data.table()
  } else {
    zero_dt <- data.table::as.data.table(zero_tbl)
    
    # Standardize the taxon column name to "genus" if possible
    if (!"genus" %in% names(zero_dt)) {
      if ("taxon" %in% names(zero_dt)) {
        data.table::setnames(zero_dt, "taxon", "genus")
      } else if (nrow(zero_dt) > 0 && !is.null(rownames(zero_tbl))) {
        zero_dt[, genus := rownames(zero_tbl)]
      } else {
        # last resort: assume first col is taxon
        data.table::setnames(zero_dt, 1, "genus")
      }
    }
  }
  
  # We need two groups from zero_ind:
  #   - present in Peneira and NOT in Floresta
  #   - present in Floresta and NOT in Peneira
  #
  # Because zero_ind schemas vary, we implement a robust "best-effort":
  # - look for columns that include both group names and look binary (0/1 or TRUE/FALSE)
  # - otherwise, skip with message
  zero_pick <- data.table::data.table()
  if (nrow(zero_dt) > 0) {
    cn <- names(zero_dt)
    
    # Find candidate columns
    # Examples seen in wild: "zero_ind_Floresta", "zero_ind_Peneira", "Floresta", "Peneira", etc.
    pick_col <- function(group) {
      hits <- grep(group, cn, ignore.case = TRUE, value = TRUE)
      if (length(hits) == 0) return(NA_character_)
      # prefer columns containing "zero" or "ind"
      pref <- hits[grep("zero|ind", hits, ignore.case = TRUE)]
      if (length(pref) > 0) return(pref[1])
      hits[1]
    }
    
    col_flo <- pick_col("Floresta")
    col_pen <- pick_col("Peneira")
    
    if (is.na(col_flo) || is.na(col_pen)) {
      message(">>> Step 10: zero_ind found but couldn't identify Floresta/Peneira indicator columns. Skipping zero_ind sampling.")
    } else {
      # Coerce to integer-ish 0/1 (TRUE/FALSE also ok)
      z <- copy(zero_dt)
      z[, flo := as.integer(as.logical(get(col_flo)) | get(col_flo) == 1)]
      z[, pen := as.integer(as.logical(get(col_pen)) | get(col_pen) == 1)]
      
      # present in Peneira not in Floresta  => pen==1 & flo==0
      pen_only <- z[pen == 1 & flo == 0]
      # present in Floresta not in Peneira => flo==1 & pen==0
      flo_only <- z[flo == 1 & pen == 0]
      
      set.seed(seed)
      if (nrow(pen_only) > 0) pen_only <- pen_only[sample(.N, min(zero_n_each_dir, .N))]
      if (nrow(flo_only) > 0) flo_only <- flo_only[sample(.N, min(zero_n_each_dir, .N))]
      
      zero_pick <- data.table::rbindlist(list(
        pen_only[, .(genus, zero_class = "Peneira_only")],
        flo_only[, .(genus, zero_class = "Floresta_only")]
      ), use.names = TRUE, fill = TRUE)
      
      zero_pick <- unique(zero_pick, by = "genus")
      zero_pick[, source := "zero_ind"]
      zero_pick[, icon_type := "binary"]
    }
  }
  
  # ---------- 3) union of selected genera ----------
  keep <- unique(c(top_res$genus, zero_pick$genus))
  keep <- keep[is.finite(match(keep, keep))]  # drop NA if any
  keep <- keep[keep != ""]
  keep <- keep[!is.na(keep)]
  
  if (length(keep) < min_unique_tips) {
    message(">>> Step 10: too few genera selected for phylogeny (n=", length(keep), "). Skipping.")
    return(invisible(NULL))
  }
  
  # ---------- 4) resolve to OTL IDs and induced subtree ----------
  message(">>> Step 10: resolving genus names in OpenTree (tnrs_match_names)...")
  m <- rotl::tnrs_match_names(keep, do_approximate_matching = TRUE)
  
  # rotl returns a data.frame-like object; keep ott_id and unique name mapping
  mdt <- data.table::as.data.table(m)
  if (!all(c("search_string", "ott_id") %in% names(mdt))) {
    stop(">>> Step 10: unexpected tnrs_match_names output; missing search_string/ott_id columns.")
  }
  
  #filter out NA + flagged OTT IDs (pruned/deprecated/etc.), then best hit per name
  mdt <- filter_otl_matches(mdt)
  
  if (nrow(mdt) == 0) {
    stop(">>> Step 10: none of the selected genera could be resolved in OpenTree (after filtering flags).")
  }
  
  resolved_keep <- mdt$search_string
  ott_ids <- mdt$ott_id
  
  if (length(ott_ids) < min_unique_tips) {
    message(">>> Step 10: too few resolved genera in OpenTree (n=", length(ott_ids), "). Skipping.")
    return(invisible(NULL))
  }
  
  #use safe wrapper (drops missing ott ids and retries)
  message(">>> Step 10: fetching induced subtree from OpenTree (tol_induced_subtree; safe)...")
  tr <- safe_tol_induced_subtree(ott_ids = ott_ids)
  
  # The tree tips are ott IDs by default; replace with your genus labels
  tip_map <- data.table::data.table(
    tip_label = tr$tip.label
  )
  tip_map[, ott_id := suppressWarnings(as.integer(gsub("^ott", "", tip_label)))]
  tip_map <- merge(tip_map, mdt[, .(search_string, ott_id)], by = "ott_id", all.x = TRUE)
  
  new_labels <- ifelse(is.na(tip_map$search_string), tip_map$tip_label, tip_map$search_string)
  tr$tip.label <- new_labels
  
  # ---------- 5) build tip icon colors ----------
  # Build a unified table with per-genus attributes
  meta_tip <- data.table::data.table(genus = tr$tip.label)
  
  # LFC from res summary (for continuous coloring)
  res_lfc <- res_dt[, .(genus, lfc, q_val)]
  meta_tip <- merge(meta_tip, res_lfc, by = "genus", all.x = TRUE)
  
  # zero_ind class (binary)
  if (nrow(zero_pick) > 0) {
    meta_tip <- merge(meta_tip, zero_pick[, .(genus, zero_class)], by = "genus", all.x = TRUE)
  } else {
    meta_tip[, zero_class := NA_character_]
  }
  
  # Define tip icon types:
  # - if genus has lfc: continuous (filled circle colored by lfc)
  # - else if genus has zero_class: binary (filled circle with two colors)
  # - else: grey
  meta_tip[, icon_group := fifelse(!is.na(lfc), "lfc",
                                   fifelse(!is.na(zero_class), "binary", "other"))]
  
  # continuous palette for lfc
  lfc_vals <- meta_tip[icon_group == "lfc", lfc]
  if (length(lfc_vals) == 0) {
    lfc_min <- -1; lfc_max <- 1
  } else {
    lfc_min <- min(lfc_vals, na.rm = TRUE)
    lfc_max <- max(lfc_vals, na.rm = TRUE)
    if (!is.finite(lfc_min) || !is.finite(lfc_max) || lfc_min == lfc_max) {
      lfc_min <- -1; lfc_max <- 1
    }
  }
  
  # Map lfc to colors (no dependencies)
  pal <- grDevices::colorRampPalette(c("blue3", "white", "red3"))
  ncol <- 101
  cols_lfc <- pal(ncol)
  lfc_to_col <- function(x) {
    if (!is.finite(x)) return("grey70")
    t <- (x - lfc_min) / (lfc_max - lfc_min)
    t <- max(0, min(1, t))
    idx <- 1 + floor(t * (ncol - 1))
    cols_lfc[idx]
  }
  
  tip_cols <- rep("grey80", nrow(meta_tip))
  tip_cols[meta_tip$icon_group == "lfc"] <- vapply(meta_tip[meta_tip$icon_group == "lfc", lfc], lfc_to_col, character(1))
  
  # binary colors for zero_ind classes
  # (choose fixed colors; you can tweak later)
  tip_cols[meta_tip$icon_group == "binary" & meta_tip$zero_class == "Peneira_only"]  <- "darkorange2"
  tip_cols[meta_tip$icon_group == "binary" & meta_tip$zero_class == "Floresta_only"] <- "deepskyblue3"
  
  # ---------- 6) plot tree + tip icons ----------
  out_png <- file.path(PLOTS_DIR, paste0(prefix, "_otol_tree_lfc_zeroind.png"))
  out_pdf <- file.path(PLOTS_DIR, paste0(prefix, "_otol_tree_lfc_zeroind.pdf"))
  
  # Save a metadata table used for the plot (reproducibility)
  data.table::fwrite(meta_tip, file.path(TABLES_DIR, paste0(prefix, "_otol_tree_tip_metadata.tsv")), sep = "\t")
  
  png(out_png, width = 1800, height = 1400, res = 150)
  par(mar = c(2, 2, 2, 10))
  ape::plot.phylo(tr, cex = 0.7, label.offset = 0.01, no.margin = TRUE)
  ape::tiplabels(pch = 16, col = tip_cols, cex = 1.1)
  title(main = paste0("OpenTree induced genus phylogeny | tip icons: LFC (res) + presence/absence (zero_ind)\n",
                      "LFC range: [", signif(lfc_min, 3), ", ", signif(lfc_max, 3), "]"))
  # small legend
  legend("topleft",
         legend = c("res (LFC<0 -> Floresta)", "res (LFC>0 -> Peneira)", "zero_ind: Peneira-only", "zero_ind: Floresta-only", "other"),
         pch = 16,
         col = c("blue3", "red3", "darkorange2", "deepskyblue3", "grey80"),
         bty = "n", cex = 0.8)
  dev.off()
  
  pdf(out_pdf, width = 14, height = 10)
  par(mar = c(2, 2, 2, 10))
  ape::plot.phylo(tr, cex = 0.7, label.offset = 0.01, no.margin = TRUE)
  ape::tiplabels(pch = 16, col = tip_cols, cex = 1.1)
  title(main = paste0("OpenTree induced genus phylogeny | tip icons: LFC (res) + presence/absence (zero_ind)\n",
                      "LFC range: [", signif(lfc_min, 3), ", ", signif(lfc_max, 3), "]"))
  legend("topleft",
         legend = c("res (LFC<0 -> Floresta)", "res (LFC>0 -> Peneira)", "zero_ind: Peneira-only", "zero_ind: Floresta-only", "other"),
         pch = 16,
         col = c("blue3", "red3", "darkorange2", "deepskyblue3", "grey80"),
         bty = "n", cex = 0.8)
  dev.off()
  
  message(">>> Step 10: phylogeny written to: ", out_png, " and ", out_pdf)
  invisible(list(tree = tr, tip_meta = meta_tip, tip_colors = tip_cols, resolved = mdt))
}

# ==========================================================
# Calling functions
# ==========================================================
merge_batches_if_needed(infile)

obj0 <- step0_read_stage(infile, outdir, prefix, USE_COUNTS_0_4, USE_COUNTS_5, TABLES_DIR)
dt_raw <- obj0$dt_raw

#step1_stacked_bars(dt_raw, VAL_0_4, outdir, prefix, PLOTS_DIR)
#obj2 <- step2_alpha(dt_raw, VAL_0_4, outdir, prefix, TABLES_DIR, PLOTS_DIR)
#rel  <- obj2$rel
#meta <- obj2$meta
#obj3 <- step3_beta_pcoa(rel, meta, outdir, prefix, TABLES_DIR, PLOTS_DIR)
#bray <- obj3$bray
#step4_beta_stats(bray, meta, outdir, prefix, TABLES_DIR)
#step5_concordance(clr_obj, outdir, prefix, MODE, TABLES_DIR, PLOTS_DIR)

clr_obj <- obj0$clr_obj

# pooled scatter across all codes (each dot = code × replicate × genus)
step6c_allcodes_clr_concordance(
  clr_obj    = clr_obj,
  outdir     = outdir,
  prefix     = prefix,
  MODE       = MODE,
  TABLES_DIR = TABLES_DIR,
  PLOTS_DIR  = PLOTS_DIR,
  color_by   = "pairing_code"   # or "replicate" or "none"
)

#Capture Step 7 output; use USE_COUNTS_ANCOM
#res7 <- step7_ancombc2(dt_raw, outdir, prefix, USE_COUNTS_ANCOM, TABLES_DIR, PLOTS_DIR)

#keep raw ANCOM-BC2 object path for zero_ind retrieval
#ancom_raw_rds <- file.path(TABLES_DIR, paste0(prefix, "_ancombc2_raw.rds"))

#Heatmap based on Step 7 significant genera
#step9_heatmap_ancom_sig(dt_raw, res7, outdir, prefix, USE_COUNTS_ANCOM, TABLES_DIR, PLOTS_DIR)

# Step 8 uses USE_COUNTS_ANCOM
#step8_ancombc2_by_code(dt_raw, outdir, prefix, USE_COUNTS_ANCOM, TABLES_DIR, PLOTS_DIR)

#build true phylogeny from OTL + overlay LFC and zero_ind
# step10_otl_phylogeny_lfc_zeroind(
#   res7_summary    = res7,
#   ancombc2_raw_rds = ancom_raw_rds,
#   outdir          = outdir,
#   prefix          = prefix,
#   PLOTS_DIR       = PLOTS_DIR,
#   TABLES_DIR      = TABLES_DIR
# )

message(">>> Done. Outputs in:")
message(">>>   tables: ", TABLES_DIR)
message(">>>   plots : ", PLOTS_DIR)
