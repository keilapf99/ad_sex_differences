#!/usr/bin/env Rscript

# Title : null_model_dea.R
# Author: Keila Adonai Pérez-Fuentes
# Date: August 18th 2026
# Objective: Build a permutation-based null model for each cohort/region/sex
#            comparison that matters downstream (AD vs control, split by sex):
#            shuffle the diagnosis/ceradsc label 1000 times among the same
#            patients (same two group sizes preserved), run a full DESeq2 fit
#            with the exact design formula used in the base script for each
#            cohort, and keep the Wald "stat" column of each run. Output: one
#            genes x 1000 matrix per combination, saved under
#            results/07_master_regulators/null_model_dea/<COHORT>_<REGION>_<sex>_null_model_DEA.tsv

time0 <- Sys.time()

cat("# Load libraries\n")
suppressPackageStartupMessages({
  library(dplyr)
  library(parallel)
  library(DESeq2)
})

# 0. Config ---------------------------------------------------------------------

SEED     <- 10
N_PERM   <- as.integer(Sys.getenv("NULL_MODEL_NPERM", "1000"))
N_CORES  <- as.integer(Sys.getenv("NULL_MODEL_CORES", "60"))

source(file.path(Sys.getenv("ADSEX_ROOT", "."), "config", "paths.R"))
mod_dir      <- module_results("07_master_regulators")
mayo_dir   <- file.path(mod_dir, "0_preprocess")
rosmap_dir <- file.path(mod_dir, "0_preprocess_rosmap")
out_dir    <- file.path(mod_dir, "null_model_dea")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("N_PERM =", N_PERM, " N_CORES =", N_CORES, " SEED =", SEED, "\n")

# 1. Build the list of 10 combinations (cohort x region x sex) ------------------

combos <- list()

for (region in c("CRB", "TC")) {
  for (sex in c("male", "female")) {
    combos[[length(combos) + 1]] <- list(
      cohort         = "MAYO",
      region         = region,
      sex            = sex,
      sex_col        = "sex",
      sex_val        = sex,
      cond_col       = "diagnosis",
      cond_vals      = c("Alzheimer Disease", "control"),
      design_formula = "~ specimenIdSource + pmi + diagnosis",
      factor_cols    = c("specimenIdSource"),
      counts_file    = file.path(mayo_dir, paste0("MAYO_", region, "_counts_preprocessed.rds")),
      meta_file      = file.path(mayo_dir, paste0("MAYO_", region, "_metadata_preprocessed.txt"))
    )
  }
}

for (region in c("DLPFC", "HCN", "PCC")) {
  for (sexval in c(1, 0)) {
    sexlabel <- ifelse(sexval == 1, "male", "female")
    combos[[length(combos) + 1]] <- list(
      cohort         = "ROSMAP",
      region         = region,
      sex            = sexlabel,
      sex_col        = "msex",
      sex_val        = sexval,
      cond_col       = "ceradsc",
      cond_vals      = c("1", "4"),
      design_formula = "~ sequencingBatch + pmi + RIN + ceradsc",
      factor_cols    = c("sequencingBatch"),
      counts_file    = file.path(rosmap_dir, paste0("ROSMAP_", region, "_counts_preprocessed.rds")),
      meta_file      = file.path(rosmap_dir, paste0("ROSMAP_", region, "_metadata_preprocessed.txt"))
    )
  }
}

cat("Total combinaciones:", length(combos), "\n")

# 2. Functions --------------------------------------------------------------------

process_combo <- function(combo, n_perm, n_cores, out_dir) {
  t0 <- Sys.time()
  cat(sprintf("\n==== %s %s %s ==== (%s)\n", combo$cohort, combo$region, combo$sex, format(t0)))

  counts   <- readRDS(combo$counts_file)
  metadata <- read.delim(combo$meta_file)

  sub_meta <- metadata[metadata[[combo$sex_col]] == combo$sex_val &
                        metadata[[combo$cond_col]] %in% combo$cond_vals, ]
  sub_counts <- counts[, sub_meta$specimenID]

  for (fc in combo$factor_cols) {
    sub_meta[[fc]] <- droplevels(as.factor(sub_meta[[fc]]))
  }
  sub_meta[[combo$cond_col]] <- droplevels(factor(as.character(sub_meta[[combo$cond_col]]),
                                                    levels = combo$cond_vals))

  cnt_tab <- table(sub_meta[[combo$cond_col]])
  cat("n =", nrow(sub_meta), " | ", paste(names(cnt_tab), cnt_tab, sep = "=", collapse = ", "), "\n")

  real_labels <- sub_meta[[combo$cond_col]]

  # Pre-generate all N_PERM permutations sequentially (single RNG stream) so the
  # exact set of permutations used is fully reproducible regardless of n_cores.
  perm_labels <- vector("list", n_perm)
  for (i in seq_len(n_perm)) {
    perm_labels[[i]] <- sample(real_labels)
  }

  design_formula <- as.formula(combo$design_formula)
  contrast <- c(combo$cond_col, combo$cond_vals[1], combo$cond_vals[2])
  counts_mat <- as.matrix(sub_counts)

  run_one <- function(i) {
    cd <- sub_meta
    cd[[combo$cond_col]] <- perm_labels[[i]]

    tryCatch({
      suppressMessages(suppressWarnings({
        dds <- DESeqDataSetFromMatrix(countData = counts_mat, colData = cd, design = design_formula)
        dds <- DESeq(dds, quiet = TRUE)
        res <- results(dds, contrast = contrast)
      }))
      as.numeric(res$stat)
    }, error = function(e) {
      message(sprintf("[%s %s %s] Permutacion %d fallo: %s",
                       combo$cohort, combo$region, combo$sex, i, conditionMessage(e)))
      rep(NA_real_, nrow(counts_mat))
    })
  }

  stat_list <- mclapply(seq_len(n_perm), run_one, mc.cores = n_cores, mc.preschedule = TRUE)

  stat_mat <- do.call(cbind, stat_list)
  rownames(stat_mat) <- rownames(counts_mat)
  colnames(stat_mat) <- paste0("perm_", seq_len(n_perm))

  out_file <- file.path(out_dir, sprintf("%s_%s_%s_null_model_DEA.tsv", combo$cohort, combo$region, combo$sex))
  out_df <- data.frame(gene = rownames(stat_mat), stat_mat, check.names = FALSE)
  write.table(out_df, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

  t1 <- Sys.time()
  n_failed <- sum(vapply(stat_list, function(x) all(is.na(x)), logical(1)))
  cat(sprintf("Guardado: %s | tiempo: %.1f min | permutaciones fallidas: %d\n",
              out_file, as.numeric(t1 - t0, units = "mins"), n_failed))
}

###### MAIN #####

set.seed(SEED)

for (combo in combos) {
  process_combo(combo, n_perm = N_PERM, n_cores = N_CORES, out_dir = out_dir)
}

cat("\nEND ---- \n")
cat("Elapsed time total: ", as.numeric(Sys.time() - time0, units = "hours"), "h\n")
