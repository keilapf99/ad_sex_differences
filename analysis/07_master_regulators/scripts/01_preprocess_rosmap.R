#!/usr/bin/env Rscript

# Title : 0_preprocess_rosmap.R
# Author: Keila Adonai Pérez-Fuentes
# Date: August 17th 2026
# Objective: Adapt the ROSMAP raw count matrices (data/raw_counts/ROSMAP_counts)
#            into the same format used in expr_diff_base_rosmap.R right before
#            differential expression, and save the resulting counts + metadata
#            in results/07_master_regulators/0_preprocess_rosmap for downstream analysis.

time <- Sys.time()

cat("# Load libraries\n")
suppressPackageStartupMessages({
  library(dplyr)
  library(stringr)
  library(tibble)
  library(tidyr)
})

# 0. Paths --------------------------------------------------------------------

source(file.path(Sys.getenv("ADSEX_ROOT", "."), "config", "paths.R"))
mod_dir      <- module_results("07_master_regulators")
raw_dir <- file.path(DATA_DIR, "raw_counts", "ROSMAP_counts")
out_dir <- file.path(mod_dir, "0_preprocess_rosmap")

condition_column <- "ceradsc"

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Define functions -----------------------------------------------------------

load_data <- function(counts_path, metadata_path) {
  counts <- readRDS(counts_path) %>% as.data.frame()
  metadata <- read.delim(metadata_path)
  return(list(counts = counts, metadata = metadata))
}

# Function to delete duplicated features (same criterion as expr_diff_base_rosmap.R)
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
    summarise(across(where(is.numeric), ~ median(.x, na.rm = TRUE)), .groups = "drop")

  counts <- counts %>% filter(!feature %in% repeated_rows$feature)
  counts <- bind_rows(counts, repeated_rows)
  counts <- counts %>% mutate(across(-feature, as.integer))

  return(counts)
}

# Pre-process data (same logic as preprocess_data() in expr_diff_base_rosmap.R)
preprocess_data <- function(counts, metadata, condition_column) {
  cat("Raw dims -> metadata:", dim(metadata), " counts:", dim(counts), "\n")

  # Drop STAR summary rows (N_unmapped, N_multimapping, N_noFeature, N_ambiguous)
  counts <- counts[-c(1:4), ]

  counts <- counts %>% mutate(feature = str_remove(feature, "\\..*$"))
  counts <- del_dupl(counts)
  rownames(counts) <- NULL
  counts <- counts %>% column_to_rownames(var = "feature")

  cat("# Filter metadata based on conditions\n")
  if ("exclude" %in% colnames(metadata)) {
    metadata <- metadata[metadata$exclude != TRUE | is.na(metadata$exclude), ]
  }
  metadata <- metadata %>% filter(!is.na(.data[[condition_column]]))

  cat("#Filter NAs in sequencingBatch\n")
  if ("sequencingBatch" %in% colnames(metadata)) {
    metadata <- metadata[!is.na(metadata$sequencingBatch), ]
  }

  #Remove value '0, 6, 7' if it exists
  if ("0, 6, 7" %in% metadata$sequencingBatch) {
    metadata <- metadata %>% filter(sequencingBatch != "0, 6, 7")
  }

  specimenIDs_in_counts <- colnames(counts)  #Assuming that the first column is not a specimenID
  metadata <- metadata %>%
    dplyr::filter(specimenID %in% specimenIDs_in_counts)

  specimenIDs <- metadata$specimenID
  counts <- counts %>% dplyr::select(any_of(specimenIDs))

  specimenIDs_in_counts <- colnames(counts)
  metadata <- metadata %>% dplyr::filter(specimenID %in% specimenIDs_in_counts)

  #Drop NA in metadata
  metadata <- metadata %>% drop_na(pmi)
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
  meta_file <- list.files(region_dir, pattern = "_metadata_.*\\.csv$", full.names = TRUE)

  if (length(rds_file) != 1 || length(meta_file) != 1) {
    warning("Skipping ", region_dir,
            ": expected exactly one counts .rds and one metadata .csv file")
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

  counts_out <- file.path(out_dir, paste0("ROSMAP_", region, "_counts_preprocessed.rds"))
  metadata_out <- file.path(out_dir, paste0("ROSMAP_", region, "_metadata_preprocessed.txt"))

  saveRDS(preprocessed_data$counts, counts_out)
  write.table(preprocessed_data$metadata, metadata_out,
              sep = "\t", quote = FALSE, row.names = FALSE)

  cat("Saved: ", counts_out, "\n")
  cat("Saved: ", metadata_out, "\n")
}

cat("\nEND ---- \n")
cat("Elapsed time: ", Sys.time() - time, "\n")
