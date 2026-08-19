#!/usr/bin/env Rscript

# Title : 0_preprocess.R
# Author: Keila Adonai Pérez-Fuentes
# Date: August 17th 2026
# Objective: Adapt the MAYO raw count matrices (data/raw_counts/MAYO_rawcounts)
#            into the same format used in expr_diff_base_mayo.R right before
#            differential expression, and save the resulting counts + metadata
#            in results/0_preprocess for downstream analysis.

time <- Sys.time()

cat("# Load libraries\n")
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
})

# 0. Paths --------------------------------------------------------------------

raw_dir <- path.expand("~/ad_sex_differences/data/raw_counts/MAYO_rawcounts")
out_dir <- path.expand("~/ad_sex_differences/results/0_preprocess")

condition_column <- "diagnosis"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Define functions -----------------------------------------------------------

load_data <- function(counts_path, metadata_path) {
  counts <- readRDS(counts_path) %>% as.data.frame()
  metadata <- read.delim(metadata_path)
  return(list(counts = counts, metadata = metadata))
}

# Function to delete duplicated features (same criterion as expr_diff_base_mayo.R)
del_dupl <- function(counts) {
  repeated_values <- counts %>%
    group_by(feature) %>%
    filter(n() > 1) %>%
    distinct(feature) %>%
    pull(feature)

  repeated_rows <- counts[counts$feature %in% repeated_values, ]
  repeated_rows <- repeated_rows[order(repeated_rows$feature), ]
  repeated_rows <- repeated_rows %>%
    group_by(feature) %>%
    summarize(across(everything(), \(x) median(x, na.rm = TRUE)))

  counts <- counts %>% filter(!feature %in% repeated_rows$feature)
  counts <- bind_rows(counts, repeated_rows)
  counts <- counts %>% mutate(across(-feature, as.integer))

  return(counts)
}

# Pre-process data (same logic as preprocess_data() in expr_diff_base_mayo.R)
preprocess_data <- function(counts, metadata, condition_column) {
  cat("Raw dims -> metadata:", dim(metadata), " counts:", dim(counts), "\n")

  # Drop STAR summary rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
  counts <- counts[-c(1:4), ]

  counts <- counts %>%
    mutate(feature = str_remove(feature, "\\..*$")) %>%
    del_dupl() %>%
    column_to_rownames(var = "feature")

  cat("# Filter metadata based on conditions\n")
  if ("exclude" %in% colnames(metadata)) {
    metadata <- metadata[metadata$exclude != TRUE | is.na(metadata$exclude), ]
  }
  metadata <- metadata %>% filter(!is.na(.data[[condition_column]]))

  cat("# Filter NAs in sequencingBatch\n")
  if ("sequencingBatch" %in% colnames(metadata)) {
    metadata <- metadata[!is.na(metadata$sequencingBatch), ]
  }

  specimenIDs_in_counts <- colnames(counts)
  metadata <- metadata %>%
    filter(specimenID %in% specimenIDs_in_counts,
           !is.na(.data[[condition_column]]),
           !is.na(pmi))

  counts <- counts %>% dplyr::select(any_of(metadata$specimenID))

  cat("Final dims -> metadata:", dim(metadata), " counts:", dim(counts), "\n")

  return(list(counts = counts, metadata = metadata))
}

###### MAIN #####

region_dirs <- list.dirs(raw_dir, recursive = FALSE)

if (length(region_dirs) == 0) {
  stop("No region subdirectories found in ", raw_dir, call. = FALSE)
}

for (region_dir in region_dirs) {

  rds_file <- list.files(region_dir, pattern = "_rawcounts_.*\\.rds$", full.names = TRUE)
  meta_file <- list.files(region_dir, pattern = "_metadata_.*\\.txt$", full.names = TRUE)

  if (length(rds_file) != 1 || length(meta_file) != 1) {
    warning("Skipping ", region_dir,
            ": expected exactly one counts .rds and one metadata .txt file")
    next
  }

  region <- sub(".*_rawcounts_(.*)\\.rds$", "\\1", basename(rds_file))

  cat("\n==== Processing region:", region, "====\n")
  cat("Counts file: ", rds_file, "\n")
  cat("Metadata file: ", meta_file, "\n")

  data <- load_data(rds_file, meta_file)

  preprocessed_data <- preprocess_data(data$counts,
                                       metadata = data$metadata,
                                       condition_column = condition_column)

  counts_out <- file.path(out_dir, paste0("MAYO_", region, "_counts_preprocessed.rds"))
  metadata_out <- file.path(out_dir, paste0("MAYO_", region, "_metadata_preprocessed.txt"))

  saveRDS(preprocessed_data$counts, counts_out)
  write.table(preprocessed_data$metadata, metadata_out,
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("Saved: ", counts_out, "\n")
  cat("Saved: ", metadata_out, "\n")
}

cat("\nEND ---- \n")
cat("Elapsed time: ", Sys.time() - time, "\n")
