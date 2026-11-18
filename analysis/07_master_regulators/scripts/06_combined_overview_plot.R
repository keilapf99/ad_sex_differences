#!/usr/bin/env Rscript

# Title : combined_overview_plot.R
# Author: Keila Adonai Pérez-Fuentes
# Date: October 21th 2026
# Objective: Single combined figure presenting both sex-comparison stories at
#            once, per region:
#              - Shared master regulators -> scatter of NES_male (x) vs
#                NES_female (y). Same-sign quadrants (I/III) = concordant,
#                opposite-sign quadrants (II/IV) = discordant.
#              - Sex-specific master regulators (no value in the other sex,
#                so they can't be placed inside the scatter) -> marginal rug
#                ticks: male-specific along the top (x axis), female-specific
#                along the right (y axis).
#            One unified categorical legend (4 categories), using the
#            dataviz-validated all-pairs-safe slots (blue/orange/aqua) plus
#            neutral gray for the de-emphasized "concordant" majority.
# Output: results/07_master_regulators/MRA_by_sex/combined_overview.{pdf,png}

suppressPackageStartupMessages(library(ggplot2))

source(file.path(Sys.getenv("ADSEX_ROOT", "."), "config", "paths.R"))
mod_dir      <- module_results("07_master_regulators")
mra_dir      <- file.path(mod_dir, "MRA_by_sex")
shared_file  <- file.path(mra_dir, "top_quartile_NES_shared", "shared_top_quartile_concordance.txt")
specific_file <- file.path(mra_dir, "sex_specific_MRs", "sex_specific_master_regulators.txt")

shared   <- read.delim(shared_file, header = TRUE, stringsAsFactors = FALSE)
specific <- read.delim(specific_file, header = TRUE, stringsAsFactors = FALSE)

region_order <- c("MAYO_CRB", "MAYO_TC", "ROSMAP_DLPFC", "ROSMAP_HCN", "ROSMAP_PCC")
shared$region   <- factor(shared$region,   levels = region_order)
specific$region <- factor(specific$region, levels = region_order)

shared$category <- ifelse(shared$concordance == "discordant", "Discordant (shared)", "Concordant (shared)")
specific$category <- ifelse(specific$sex_specific == "male_only", "Male-specific", "Female-specific")

cat_levels <- c("Concordant (shared)", "Discordant (shared)", "Male-specific", "Female-specific")
color_pal  <- c("Concordant (shared)" = "#c3c2b7",   # muted gray - de-emphasized majority
                 "Discordant (shared)" = "#1baf7a",   # aqua - slot 3 (all-pairs safe with blue/orange)
                 "Male-specific"       = "#2a78d6",   # blue - slot 1
                 "Female-specific"     = "#eb6834")   # orange - slot 2

shared$category   <- factor(shared$category,   levels = cat_levels)
specific$category <- factor(specific$category, levels = cat_levels)

male_only_df   <- specific[specific$sex_specific == "male_only", ]
female_only_df <- specific[specific$sex_specific == "female_only", ]

p <- ggplot() +
  geom_vline(xintercept = 0, color = "#c3c2b7", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "#c3c2b7", linewidth = 0.4) +
  geom_point(data = shared, aes(x = NES_male, y = NES_female, color = category),
             alpha = 0.75, size = 1.8) +
  geom_rug(data = male_only_df, aes(x = NES, color = category),
           sides = "t", length = unit(0.05, "npc"), alpha = 0.8) +
  geom_rug(data = female_only_df, aes(y = NES, color = category),
           sides = "r", length = unit(0.05, "npc"), alpha = 0.8) +
  facet_wrap(~region, ncol = 3, scales = "free") +
  scale_color_manual(values = color_pal, name = NULL, drop = FALSE) +
  labs(title = "Panorama combinado: especificidad sexual y concordancia de master regulators",
       subtitle = paste(
         "Centro: reguladores compartidos (NES hombres vs mujeres) - mismo cuadrante de signo = concordante, cuadrante opuesto = discordante.",
         "Margenes: reguladores sexo-especificos (arriba = solo en hombres, derecha = solo en mujeres).", sep = "\n"),
       x = "NES (hombres)", y = "NES (mujeres)") +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        plot.subtitle = element_text(size = 9, color = "#52514e"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#e1e0d9"),
        legend.position = "top",
        strip.text = element_text(face = "bold"))

ggsave(file.path(mra_dir, "combined_overview.pdf"), p, width = 11, height = 8)
ggsave(file.path(mra_dir, "combined_overview.png"), p, width = 11, height = 8, dpi = 150)

cat("Guardado:", file.path(mra_dir, "combined_overview.pdf"), "\n")
cat("Guardado:", file.path(mra_dir, "combined_overview.png"), "\n")
