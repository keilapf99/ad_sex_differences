#!/usr/bin/env Rscript

# Title : nes_quartile_filter.R
# Author: Keila Adonai Pérez-Fuentes
# Date: October 20th 2026
# Objective: For each cohort/region/sex MRA combination, plot the distribution
#            of |NES| across all regulators tested by msviper, and keep only
#            the regulators in the top quartile (Q3-Q4, i.e. the 25% with the
#            highest |NES|) for downstream analysis. This replaces a fixed
#            FDR<0.05 cutoff (too strict given how few combinations reach it)
#            and a fixed top-N cutoff (unfair across combinations that tested
#            very different numbers of regulators, from ~400 to ~713).
#            Output:
#            - results/07_master_regulators/MRA_by_sex/nes_distribution_plots/<label>_NES_distribution.pdf
#            - results/07_master_regulators/MRA_by_sex/top_quartile_NES/<label>_top_quartile_NES.txt

suppressPackageStartupMessages(library(ggplot2))

source(file.path(Sys.getenv("ADSEX_ROOT", "."), "config", "paths.R"))
mod_dir      <- module_results("07_master_regulators")
mra_dir   <- file.path(mod_dir, "MRA_by_sex")
plot_dir  <- file.path(mra_dir, "nes_distribution_plots")
out_dir   <- file.path(mra_dir, "top_quartile_NES")

dir.create(plot_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)

result_files <- list.files(mra_dir, pattern = "_MRA_results\\.txt$", full.names = TRUE)
cat("Combinaciones encontradas:", length(result_files), "\n")

summary_rows <- list()

for (f in result_files) {
  label <- sub("_MRA_results\\.txt$", "", basename(f))
  cat("\n====", label, "====\n")

  df <- read.delim(f, header = TRUE, stringsAsFactors = FALSE, check.names = FALSE)
  df$abs_NES <- abs(df$NES)

  q <- quantile(df$abs_NES, probs = c(0.25, 0.50, 0.75))
  q3 <- q[["75%"]]

  cat(sprintf("n reguladores = %d | Q1=%.3f Q2=%.3f Q3=%.3f\n",
              nrow(df), q[["25%"]], q[["50%"]], q3))

  # --- Plot: distribución de |NES| con línea de corte en Q3 ---------------
  p <- ggplot(df, aes(x = abs_NES)) +
    geom_histogram(aes(y = after_stat(density)), bins = 40,
                    fill = "#4C72B0", color = "white", linewidth = 0.2, alpha = 0.85) +
    geom_density(color = "#2B2B2B", linewidth = 0.6) +
    geom_vline(xintercept = q3, color = "#C44E52", linetype = "dashed", linewidth = 0.7) +
    annotate("text", x = q3, y = Inf, label = sprintf("Q3 = %.2f", q3),
             color = "#C44E52", hjust = -0.1, vjust = 1.5, size = 3.5) +
    labs(title = sprintf("Distribucion de |NES| - %s", label),
         subtitle = sprintf("n = %d reguladores | top cuartil (|NES| >= Q3) pasa al siguiente paso", nrow(df)),
         x = "|NES|", y = "Densidad") +
    theme_minimal(base_size = 12) +
    theme(plot.title = element_text(face = "bold"))

  ggsave(file.path(plot_dir, sprintf("%s_NES_distribution.pdf", label)), p, width = 7, height = 5)
  ggsave(file.path(plot_dir, sprintf("%s_NES_distribution.png", label)), p, width = 7, height = 5, dpi = 150)

  # --- Filtro: top cuartil (75-100%) por |NES| -----------------------------
  df_top <- df[df$abs_NES >= q3, ]
  df_top <- df_top[order(-df_top$abs_NES), ]
  df_top$abs_NES <- NULL

  out_file <- file.path(out_dir, sprintf("%s_top_quartile_NES.txt", label))
  write.table(df_top, out_file, sep = "\t", quote = FALSE, row.names = FALSE)

  cat(sprintf("Guardado: %s | reguladores retenidos: %d / %d (%.1f%%)\n",
              out_file, nrow(df_top), nrow(df), 100 * nrow(df_top) / nrow(df)))

  summary_rows[[label]] <- data.frame(
    combination = label,
    n_total = nrow(df),
    Q1 = q[["25%"]], Q2 = q[["50%"]], Q3 = q3,
    n_retained = nrow(df_top)
  )
}

summary_df <- do.call(rbind, summary_rows)
write.table(summary_df, file.path(out_dir, "summary_quartile_filter.txt"),
            sep = "\t", quote = FALSE, row.names = FALSE)

cat("\n==== Resumen ====\n")
print(summary_df, row.names = FALSE)
cat("\nEND ----\n")
