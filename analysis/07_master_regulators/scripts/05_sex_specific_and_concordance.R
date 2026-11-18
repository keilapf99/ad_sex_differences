#!/usr/bin/env Rscript

# Title : sex_specific_and_concordance.R
# Author: Keila Adonai Pérez-Fuentes
# Date: October 20th 2026
# Objective: Build the final MR (master regulator) comparison between sexes,
#            per cohort/region, following this flow (decided with the user):
#
#   1. Compare the full universe of regulators tested by msviper (male vs
#      female MRA_results.txt) per region -> shared / male_only / female_only.
#      Regulators not considered by msviper in one sex were run under the
#      exact same ARACNe/msviper parameters as the other sex, so their
#      absence from one sex's result table is itself taken as evidence of
#      sex-specific regulatory architecture (not treated as "untestable").
#   2. Regulators found in only one sex -> "sex_specific_master_regulators.txt"
#      (one consolidated table, all regions), annotated with gene symbol and
#      a column stating whether it's male_only or female_only, per region.
#   3. Regulators shared between both sexes (what's left after step 2) ->
#      plot the |NES| distribution per region/sex (restricted to the shared
#      set) and keep only the top quartile (Q3-Q4, |NES| >= Q3).
#   4. Build the final shared/top-quartile table: a regulator qualifies if it
#      is in the top quartile of AT LEAST ONE of the two sexes (union); report
#      NES/p.value/FDR for BOTH sexes side by side, gene symbol, and a
#      "concordant" (same direction of NES in both sexes) vs "discordant"
#      (opposite direction) column.
#
# Output:
#   - results/07_master_regulators/MRA_by_sex/sex_specific_MRs/sex_specific_master_regulators.txt
#   - results/07_master_regulators/MRA_by_sex/shared_NES_distribution_plots/<region>_<sex>_shared_NES_distribution.{pdf,png}
#   - results/07_master_regulators/MRA_by_sex/top_quartile_NES_shared/<region>_<sex>_top_quartile_NES_shared.txt
#   - results/07_master_regulators/MRA_by_sex/top_quartile_NES_shared/summary_quartile_filter_shared.txt
#   - results/07_master_regulators/MRA_by_sex/top_quartile_NES_shared/shared_top_quartile_concordance.txt

suppressPackageStartupMessages(library(ggplot2))

source(file.path(Sys.getenv("ADSEX_ROOT", "."), "config", "paths.R"))
mod_dir      <- module_results("07_master_regulators")
mra_dir    <- file.path(mod_dir, "MRA_by_sex")
annot_file <- file.path(DATA_DIR, "gene_id_to_symbol.tsv")

specific_dir <- file.path(mra_dir, "sex_specific_MRs")
plot_dir     <- file.path(mra_dir, "shared_NES_distribution_plots")
shared_dir   <- file.path(mra_dir, "top_quartile_NES_shared")

dir.create(specific_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plot_dir,     recursive = TRUE, showWarnings = FALSE)
dir.create(shared_dir,   recursive = TRUE, showWarnings = FALSE)

# 0. Annotation -------------------------------------------------------------

annot <- read.delim(annot_file, header = TRUE, stringsAsFactors = FALSE)
annot <- annot[!duplicated(annot$ensembl_gene_id), ]

add_symbol <- function(df, id_col = "Regulon") {
  df$gene_symbol <- annot$gene_symbol[match(df[[id_col]], annot$ensembl_gene_id)]
  df
}

regions <- c("MAYO_CRB", "MAYO_TC", "ROSMAP_DLPFC", "ROSMAP_HCN", "ROSMAP_PCC")

read_region <- function(region_label, sex) {
  f <- file.path(mra_dir, sprintf("%s_%s_MRA_results.txt", region_label, sex))
  read.delim(f, header = TRUE, stringsAsFactors = FALSE)
}

sex_specific_list   <- list()
shared_quartile_list <- list()
concordance_list     <- list()
summary_rows         <- list()

for (region_label in regions) {

  cat(sprintf("\n==== %s ====\n", region_label))

  male_df   <- read_region(region_label, "male")
  female_df <- read_region(region_label, "female")

  male_ids   <- male_df$Regulon
  female_ids <- female_df$Regulon

  shared_ids     <- intersect(male_ids, female_ids)
  male_only_ids  <- setdiff(male_ids, female_ids)
  female_only_ids <- setdiff(female_ids, male_ids)

  cat(sprintf("Total male=%d female=%d | shared=%d | male_only=%d | female_only=%d\n",
              length(male_ids), length(female_ids), length(shared_ids),
              length(male_only_ids), length(female_only_ids)))

  # --- 1-2. Sex-specific regulators (not shared) ---------------------------

  male_only_df <- male_df[male_df$Regulon %in% male_only_ids, ]
  male_only_df$region <- region_label
  male_only_df$sex_specific <- "male_only"

  female_only_df <- female_df[female_df$Regulon %in% female_only_ids, ]
  female_only_df$region <- region_label
  female_only_df$sex_specific <- "female_only"

  sex_specific_list[[region_label]] <- rbind(male_only_df, female_only_df)

  # --- 3. Shared subset: |NES| distribution + top quartile -----------------

  male_shared_df   <- male_df[male_df$Regulon %in% shared_ids, ]
  female_shared_df <- female_df[female_df$Regulon %in% shared_ids, ]

  for (sex in c("male", "female")) {
    df <- if (sex == "male") male_shared_df else female_shared_df
    df$abs_NES <- abs(df$NES)
    q <- quantile(df$abs_NES, probs = c(0.25, 0.50, 0.75))
    q3 <- q[["75%"]]

    label <- sprintf("%s_%s", region_label, sex)
    cat(sprintf("  [shared %s] n=%d Q1=%.3f Q2=%.3f Q3=%.3f\n",
                sex, nrow(df), q[["25%"]], q[["50%"]], q3))

    p <- ggplot(df, aes(x = abs_NES)) +
      geom_histogram(aes(y = after_stat(density)), bins = 40,
                      fill = "#4C72B0", color = "white", linewidth = 0.2, alpha = 0.85) +
      geom_density(color = "#2B2B2B", linewidth = 0.6) +
      geom_vline(xintercept = q3, color = "#C44E52", linetype = "dashed", linewidth = 0.7) +
      annotate("text", x = q3, y = Inf, label = sprintf("Q3 = %.2f", q3),
               color = "#C44E52", hjust = -0.1, vjust = 1.5, size = 3.5) +
      labs(title = sprintf("Distribucion de |NES| (solo compartidos) - %s", label),
           subtitle = sprintf("n = %d reguladores compartidos entre sexos | top cuartil pasa al siguiente paso", nrow(df)),
           x = "|NES|", y = "Densidad") +
      theme_minimal(base_size = 12) +
      theme(plot.title = element_text(face = "bold"))

    ggsave(file.path(plot_dir, sprintf("%s_shared_NES_distribution.pdf", label)), p, width = 7, height = 5)
    ggsave(file.path(plot_dir, sprintf("%s_shared_NES_distribution.png", label)), p, width = 7, height = 5, dpi = 150)

    df_top <- df[df$abs_NES >= q3, ]
    df_top <- df_top[order(-df_top$abs_NES), ]
    df_top$abs_NES <- NULL
    write.table(df_top, file.path(shared_dir, sprintf("%s_top_quartile_NES_shared.txt", label)),
                sep = "\t", quote = FALSE, row.names = FALSE)

    summary_rows[[label]] <- data.frame(
      region = region_label, sex = sex, n_shared_total = nrow(df),
      Q1 = q[["25%"]], Q2 = q[["50%"]], Q3 = q3, n_retained = nrow(df_top)
    )

    if (sex == "male") male_top_ids <- df_top$Regulon else female_top_ids <- df_top$Regulon
  }

  # --- 4. Final table: union of top-quartile regulators, both-sex NES, concordance ---

  union_ids <- union(male_top_ids, female_top_ids)

  m <- male_shared_df[match(union_ids, male_shared_df$Regulon),
                       c("Regulon", "Size", "NES", "p.value", "FDR")]
  names(m) <- c("Regulon", "Size_male", "NES_male", "p.value_male", "FDR_male")

  f <- female_shared_df[match(union_ids, female_shared_df$Regulon),
                         c("Size", "NES", "p.value", "FDR")]
  names(f) <- c("Size_female", "NES_female", "p.value_female", "FDR_female")

  combined <- cbind(region = region_label, m, f)
  combined$in_top_quartile_male   <- combined$Regulon %in% male_top_ids
  combined$in_top_quartile_female <- combined$Regulon %in% female_top_ids
  combined$concordance <- ifelse(sign(combined$NES_male) == sign(combined$NES_female),
                                  "concordant", "discordant")

  concordance_list[[region_label]] <- combined
}

# --- Save consolidated outputs ---------------------------------------------

sex_specific_all <- do.call(rbind, sex_specific_list)
sex_specific_all <- add_symbol(sex_specific_all, "Regulon")
sex_specific_all <- sex_specific_all[, c("region", "sex_specific", "Regulon", "gene_symbol",
                                          "Size", "NES", "p.value", "FDR")]
write.table(sex_specific_all,
            file.path(specific_dir, "sex_specific_master_regulators.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, file.path(shared_dir, "summary_quartile_filter_shared.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

concordance_all <- do.call(rbind, concordance_list)
concordance_all <- add_symbol(concordance_all, "Regulon")
concordance_all <- concordance_all[, c("region", "Regulon", "gene_symbol",
                                        "Size_male", "NES_male", "p.value_male", "FDR_male",
                                        "Size_female", "NES_female", "p.value_female", "FDR_female",
                                        "in_top_quartile_male", "in_top_quartile_female", "concordance")]
concordance_all <- concordance_all[order(concordance_all$region, -abs(concordance_all$NES_male)), ]
write.table(concordance_all,
            file.path(shared_dir, "shared_top_quartile_concordance.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n==== Resumen general ====\n")
cat(sprintf("Reguladores sexo-especificos (total, todas las regiones): %d\n", nrow(sex_specific_all)))
cat(sprintf("  male_only: %d | female_only: %d\n",
            sum(sex_specific_all$sex_specific == "male_only"),
            sum(sex_specific_all$sex_specific == "female_only")))
cat(sprintf("Reguladores compartidos en top cuartil (union male u female), todas las regiones: %d\n",
            nrow(concordance_all)))
cat(sprintf("  concordantes: %d | discordantes: %d\n",
            sum(concordance_all$concordance == "concordant"),
            sum(concordance_all$concordance == "discordant")))

cat("\nEND ----\n")
