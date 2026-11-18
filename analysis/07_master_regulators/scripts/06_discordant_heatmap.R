#!/usr/bin/env Rscript

# Title : discordant_heatmap.R
# Author: Keila Adonai Pérez-Fuentes
# Date: October 21th 2026
# Objective: Heatmap of the "discordant" shared master regulators (same
#            Regulon tested and in the top-NES-quartile in at least one sex,
#            but with opposite direction of NES between male and female),
#            from shared_top_quartile_concordance.txt. Two columns (Male NES,
#            Female NES) per regulator, rows grouped by region, diverging
#            color scale (blue = activated / positive NES, red = repressed /
#            negative NES, neutral gray midpoint at 0) per the dataviz
#            palette's diverging pair.
# Output: results/07_master_regulators/MRA_by_sex/top_quartile_NES_shared/discordant_heatmap.{pdf,png}

suppressPackageStartupMessages(library(ggplot2))

source(file.path(Sys.getenv("ADSEX_ROOT", "."), "config", "paths.R"))
mod_dir      <- module_results("07_master_regulators")
shared_dir <- file.path(mod_dir, "MRA_by_sex", "top_quartile_NES_shared")
in_file    <- file.path(shared_dir, "shared_top_quartile_concordance.txt")

df <- read.delim(in_file, header = TRUE, stringsAsFactors = FALSE)
disc <- df[df$concordance == "discordant", ]
cat("Reguladores discordantes:", nrow(disc), "\n")

region_order <- c("MAYO_CRB", "MAYO_TC", "ROSMAP_DLPFC", "ROSMAP_HCN", "ROSMAP_PCC")
disc$region <- factor(disc$region, levels = intersect(region_order, unique(disc$region)))

# Order rows within region by magnitude of discordance (|NES_male - NES_female|)
disc$delta <- abs(disc$NES_male - disc$NES_female)
disc <- disc[order(disc$region, -disc$delta), ]
disc$gene_symbol[is.na(disc$gene_symbol) | disc$gene_symbol == ""] <- disc$Regulon[is.na(disc$gene_symbol) | disc$gene_symbol == ""]

# Row order preserved within each region facet (top-to-bottom = largest divergence first)
gene_levels <- unique(disc$gene_symbol[order(disc$region, disc$delta)])
disc$gene_symbol <- factor(disc$gene_symbol, levels = gene_levels)

long <- rbind(
  data.frame(region = disc$region, gene_symbol = disc$gene_symbol, sex = "Male",   NES = disc$NES_male),
  data.frame(region = disc$region, gene_symbol = disc$gene_symbol, sex = "Female", NES = disc$NES_female)
)
long$sex <- factor(long$sex, levels = c("Male", "Female"))

p <- ggplot(long, aes(x = sex, y = gene_symbol, fill = NES)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.2f", NES)), size = 2.6, color = "#0b0b0b") +
  facet_grid(region ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = "#e34948", mid = "#f0efec", high = "#2a78d6", midpoint = 0,
                        name = "NES") +
  labs(title = "Master regulators discordantes entre sexos",
       subtitle = sprintf("n = %d reguladores compartidos con direccion opuesta de NES (top cuartil en al menos un sexo)", nrow(disc)),
       x = NULL, y = NULL) +
  theme_minimal(base_size = 11) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid = element_blank(),
        strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 0),
        strip.placement = "outside",
        axis.text.y = element_text(size = 8))

n_row <- nrow(disc)
fig_height <- max(5, 0.17 * n_row + 1.2)

ggsave(file.path(shared_dir, "discordant_heatmap.pdf"), p, width = 6, height = fig_height, limitsize = FALSE)
ggsave(file.path(shared_dir, "discordant_heatmap.png"), p, width = 6, height = fig_height, dpi = 150, limitsize = FALSE)

cat("Guardado:", file.path(shared_dir, "discordant_heatmap.pdf"), "\n")
cat("Guardado:", file.path(shared_dir, "discordant_heatmap.png"), "\n")
