#!/usr/bin/env Rscript

# Title : MRA_by_sex.R
# Author: Keila Adonai Pérez-Fuentes
# Date: August 19th 2026
# Objective: Master Regulator Analysis (VIPER/msviper) for each cohort/region/sex
#            AD-only combination. Adapted from mra_signature.R:
#            - Expression matrix for the regulon: TMM AD-only counts (data/counts_tmm)
#            - Molecular signature: full differential expression table
#              (data/full_diff_expr), same "stat" column as before
#            - Null model: the permutation-based null model built in
#              null_model_dea.R (results/null_model_dea), replacing the ad-hoc
#              name-shuffled nullmodel of the original script
#            - Regulon: TF-filtered TRN (data/trn_networks), not the raw
#              ARACNe coexpression network

time0 <- Sys.time()

cat("# Load libraries\n")
suppressPackageStartupMessages(library(viper))

# 0. Config ---------------------------------------------------------------------

CORES <- as.integer(Sys.getenv("MRA_CORES", "10"))

base_dir     <- path.expand("~/ad_sex_differences")
counts_dir   <- file.path(base_dir, "data", "counts_tmm", "expression_alzheimer")
de_dir       <- file.path(base_dir, "data", "full_diff_expr")
trn_dir      <- file.path(base_dir, "data", "trn_networks")
nullmodel_dir <- file.path(base_dir, "results", "null_model_dea")
out_dir      <- file.path(base_dir, "results", "MRA_by_sex")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("CORES =", CORES, "\n")

# 1. Build the list of 10 combinations (cohort x region x sex) ------------------

combos <- list()

for (region in c("CRB", "TC")) {
  for (sex in c("male", "female")) {
    combos[[length(combos) + 1]] <- list(
      cohort      = "MAYO",
      region      = region,
      sex         = sex,
      counts_file = file.path(counts_dir, sprintf("ROSMAP_RNAseq_%s_AD_counts_%s.tsv", sex, region)),
      de_file     = file.path(de_dir, sprintf("Mayo_%s_DE_full_gene_list_control_%s_vs_AD_%s.txt", region, sex, sex)),
      trn_file    = file.path(trn_dir, sprintf("%s_AD_%s_trn.txt", sex, region)),
      null_file   = file.path(nullmodel_dir, sprintf("MAYO_%s_%s_null_model_DEA.tsv", region, sex))
    )
  }
}

for (region in c("DLPFC", "HCN", "PCC")) {
  for (sex in c("male", "female")) {
    msex <- ifelse(sex == "male", 1, 0)
    combos[[length(combos) + 1]] <- list(
      cohort      = "ROSMAP",
      region      = region,
      sex         = sex,
      counts_file = file.path(counts_dir, sprintf("ROSMAP_RNAseq_msex%d_cerad_1_counts_%s.tsv", msex, region)),
      de_file     = file.path(de_dir, sprintf("%s_ROSMAP_DE_full_gene_list_control_%s_vs_AD_%s.txt", region, sex, sex)),
      trn_file    = file.path(trn_dir, sprintf("%s_AD_%s_trn.txt", sex, region)),
      null_file   = file.path(nullmodel_dir, sprintf("ROSMAP_%s_%s_null_model_DEA.tsv", region, sex))
    )
  }
}

test_n <- Sys.getenv("MRA_TEST_N", "")
if (nzchar(test_n)) combos <- combos[seq_len(as.integer(test_n))]

cat("Total combinaciones:", length(combos), "\n")

# 2. Functions --------------------------------------------------------------------

load_counts <- function(path) {
  counts <- read.table(path, header = TRUE, check.names = FALSE, sep = "\t")
  rownames(counts) <- counts$feature
  counts$feature <- NULL
  as.matrix(counts)
}

load_signature <- function(path) {
  expr_diff <- read.delim(path, header = TRUE, stringsAsFactors = FALSE)
  signature <- expr_diff$stat
  names(signature) <- expr_diff$ensembl_gene_id

  signature <- signature[!is.na(signature)]
  signature <- signature[names(signature) != ""]
  signature <- signature[!duplicated(names(signature))]
  sort(signature, decreasing = TRUE)
}

load_nullmodel <- function(path) {
  nm <- read.delim(path, header = TRUE, check.names = FALSE)
  mat <- as.matrix(nm[, -1])
  rownames(mat) <- nm$gene
  mat
}

process_combo <- function(combo, cores, out_dir) {
  label <- sprintf("%s_%s_%s", combo$cohort, combo$region, combo$sex)
  t0 <- Sys.time()
  cat(sprintf("\n==== %s ==== (%s)\n", label, format(t0)))

  missing <- Filter(function(f) !file.exists(f),
                     c(combo$counts_file, combo$de_file, combo$trn_file, combo$null_file))
  if (length(missing) > 0) {
    warning(sprintf("[%s] Faltan archivos, se omite: %s", label, paste(missing, collapse = "; ")))
    return(invisible(NULL))
  }

  alzheimer <- load_counts(combo$counts_file)
  storage.mode(alzheimer) <- "double"

  regulon <- aracne2regulon(combo$trn_file, alzheimer, verbose = FALSE)

  signature <- load_signature(combo$de_file)

  nullmodel <- load_nullmodel(combo$null_file)

  mrs <- msviper(signature, regulon, nullmodel, verbose = FALSE, cores = cores)

  n_reg <- length(mrs$es$nes)
  results_table <- summary(mrs, n_reg)

  saveRDS(mrs, file.path(out_dir, paste0(label, "_mrs.rds")))
  write.table(results_table, file.path(out_dir, paste0(label, "_MRA_results.txt")),
              sep = "\t", quote = FALSE, row.names = FALSE)

  pdf(file.path(out_dir, paste0(label, "_mrs_plot.pdf")))
  plot(mrs, mrs = min(25, n_reg))
  dev.off()

  t1 <- Sys.time()
  cat(sprintf("Guardado: %s | reguladores probados: %d | tiempo: %.1f min\n",
              file.path(out_dir, paste0(label, "_MRA_results.txt")), n_reg,
              as.numeric(t1 - t0, units = "mins")))
}

###### MAIN #####

for (combo in combos) {
  process_combo(combo, cores = CORES, out_dir = out_dir)
}

cat("\nEND ---- \n")
cat("Elapsed time total: ", as.numeric(Sys.time() - time0, units = "hours"), "h\n")
