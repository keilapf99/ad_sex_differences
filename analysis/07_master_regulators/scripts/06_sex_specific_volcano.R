#!/usr/bin/env Rscript

# Title : sex_specific_volcano.R
# Author: Keila Adonai Pérez-Fuentes
# Date: October 21th 2026
# Objective: Volcano plot (NES vs -log10 p.value) of the sex-specific master
#            regulators (male_only / female_only per region, from
#            sex_specific_and_concordance.R), faceted by region. Colored by
#            specificity category (2-series categorical: blue = male-specific,
#            orange = female-specific — slots 1-2 of the validated dataviz
#            palette, pre-validated all-pairs safe for up to 3 series).
#            The most significant regulators per region/category are labeled
#            with their gene symbol.
# Output: results/07_master_regulators/MRA_by_sex/sex_specific_MRs/sex_specific_volcano.{pdf,png}

suppressPackageStartupMessages(library(ggplot2))

source(file.path(Sys.getenv("ADSEX_ROOT", "."), "config", "paths.R"))
mod_dir      <- module_results("07_master_regulators")
specific_dir <- file.path(mod_dir, "MRA_by_sex", "sex_specific_MRs")
in_file      <- file.path(specific_dir, "sex_specific_master_regulators.txt")

df <- read.delim(in_file, header = TRUE, stringsAsFactors = FALSE)
df$neglog10p <- -log10(df$p.value)

region_order <- c("MAYO_CRB", "MAYO_TC", "ROSMAP_DLPFC", "ROSMAP_HCN", "ROSMAP_PCC")
df$region <- factor(df$region, levels = region_order)

label_pal  <- c(male_only = "Male-specific", female_only = "Female-specific")
color_pal  <- c(male_only = "#2a78d6", female_only = "#eb6834")   # dataviz palette slots 1-2

# Top 3 most significant per region x category -> gene_symbol label
df <- df[order(df$region, df$sex_specific, -df$neglog10p), ]
df$rank <- ave(-df$neglog10p, df$region, df$sex_specific, FUN = rank)
df$label <- ifelse(df$rank <= 3 & !is.na(df$gene_symbol) & df$gene_symbol != "",
                    df$gene_symbol, NA)

p <- ggplot(df, aes(x = NES, y = neglog10p, color = sex_specific)) +
  geom_vline(xintercept = 0, color = "#c3c2b7", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#898781", linewidth = 0.4) +
  geom_point(alpha = 0.75, size = 1.8) +
  geom_text(aes(label = label), color = "#0b0b0b", size = 2.8,
            vjust = -0.7, na.rm = TRUE, check_overlap = TRUE) +
  facet_wrap(~region, ncol = 3, scales = "free") +
  scale_color_manual(values = color_pal, labels = label_pal, name = "Especificidad") +
  labs(title = "Master regulators sexo-especificos por region",
       subtitle = "Regulones probados por msviper en un solo sexo (misma red TRN, mismos parametros ARACNe/VIPER) | linea punteada: p = 0.05",
       x = "NES", y = expression(-log[10](p-value))) +
  theme_minimal(base_size = 12) +
  theme(plot.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "#e1e0d9"),
        legend.position = "top",
        strip.text = element_text(face = "bold"))

ggsave(file.path(specific_dir, "sex_specific_volcano.pdf"), p, width = 11, height = 7)
ggsave(file.path(specific_dir, "sex_specific_volcano.png"), p, width = 11, height = 7, dpi = 150)

cat("Guardado:", file.path(specific_dir, "sex_specific_volcano.pdf"), "\n")
cat("Guardado:", file.path(specific_dir, "sex_specific_volcano.png"), "\n")
