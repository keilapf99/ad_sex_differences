#!/usr/bin/env Rscript

# Title : three_panel_summary.R
# Author: Keila Adonai Pérez-Fuentes
# Date: October 26th 2026 (restructured October 27th 2026)
# Objective: Summary figure(s) of the sex-comparison MRA. Currently composing
#            only Panel A + Panel B (Panel C kept defined below, not drawn,
#            for when it's added back):
#   Panel A - Composition of the FULL (unfiltered) master regulator universe
#             per region: male-only / shared / female-only, as a DIVERGING
#             stacked bar centered on zero (population-pyramid style): shared
#             sits centered (neutral gray), male-only extends left (blue),
#             female-only extends right (orange). Chosen instead of a plain
#             stacked bar so the male/female asymmetry per region reads
#             immediately from bar length on each side, not just segment width.
#   Panel B - Direct Spearman correlation between the female-AD NES profile
#             and the male-AD NES profile, one value per region, computed on
#             the regulators shared/tested in both sexes for that region (the
#             full, unfiltered set). Spearman because NES profiles aren't
#             guaranteed linear/normal and what matters is rank consistency of
#             regulator activity between sexes.
#   Panel C - (not drawn for now) Combined scatter, NES_male vs NES_female,
#             faceted by region; sex-exclusive regulators on the axis of the
#             sex that has them.
#
# No patchwork/gridExtra/cowplot available in this environment -> composed
# manually with base grid viewports.
#
# Output: results/07_master_regulators/MRA_by_sex/three_panel_summary.{pdf,png}

suppressPackageStartupMessages({
  library(ggplot2)
  library(grid)
})

source(file.path(Sys.getenv("ADSEX_ROOT", "."), "config", "paths.R"))
mod_dir      <- module_results("07_master_regulators")
mra_dir       <- file.path(mod_dir, "MRA_by_sex")
shared_file   <- file.path(mra_dir, "top_quartile_NES_shared", "shared_top_quartile_concordance.txt")
specific_file <- file.path(mra_dir, "sex_specific_MRs", "sex_specific_master_regulators.txt")

region_order <- c("MAYO_CRB", "MAYO_TC", "ROSMAP_DLPFC", "ROSMAP_HCN", "ROSMAP_PCC")

read_region <- function(region_label, sex) {
  f <- file.path(mra_dir, sprintf("%s_%s_MRA_results.txt", region_label, sex))
  df <- read.delim(f, header = TRUE, stringsAsFactors = FALSE)
  df$region <- region_label
  df
}

all_male   <- do.call(rbind, lapply(region_order, read_region, sex = "male"))
all_female <- do.call(rbind, lapply(region_order, read_region, sex = "female"))

# ============================================================================
# Panel A - Composition per region, diverging stacked bar (population-pyramid)
# ============================================================================

comp_rows <- lapply(region_order, function(region_label) {
  m_ids <- all_male$Regulon[all_male$region == region_label]
  f_ids <- all_female$Regulon[all_female$region == region_label]
  data.frame(region = region_label,
             male_only = length(setdiff(m_ids, f_ids)),
             shared = length(intersect(m_ids, f_ids)),
             female_only = length(setdiff(f_ids, m_ids)))
})
comp_df <- do.call(rbind, comp_rows)
comp_df$y <- length(region_order):1   # numeric y so geom_rect can center-stack per row

# diverging stacked layout: shared centered on 0, male_only flanks left, female_only flanks right
comp_df$shared_half  <- comp_df$shared / 2
comp_df$male_xmax    <- -comp_df$shared_half
comp_df$male_xmin    <- comp_df$male_xmax - comp_df$male_only
comp_df$female_xmin  <- comp_df$shared_half
comp_df$female_xmax  <- comp_df$female_xmin + comp_df$female_only

rect_df <- rbind(
  data.frame(y = comp_df$y, xmin = comp_df$male_xmin,   xmax = comp_df$male_xmax,
             category = "Male-only",   label = comp_df$male_only,
             label_x = comp_df$male_xmin - 15),
  data.frame(y = comp_df$y, xmin = -comp_df$shared_half, xmax = comp_df$shared_half,
             category = "Shared",      label = comp_df$shared,
             label_x = 0),
  data.frame(y = comp_df$y, xmin = comp_df$female_xmin,  xmax = comp_df$female_xmax,
             category = "Female-only", label = comp_df$female_only,
             label_x = comp_df$female_xmax + 15)
)
rect_df$category <- factor(rect_df$category, levels = c("Male-only", "Shared", "Female-only"))

comp_colors <- c("Male-only" = "#2a78d6", "Shared" = "#c3c2b7", "Female-only" = "#eb6834")

panel_A <- ggplot() +
  geom_vline(xintercept = 0, color = "#898781", linewidth = 0.4) +
  geom_rect(data = rect_df, aes(xmin = xmin, xmax = xmax, ymin = y - 0.35, ymax = y + 0.35, fill = category),
            color = "white", linewidth = 0.6) +
  geom_text(data = rect_df, aes(x = label_x, y = y, label = label), size = 2.8, color = "#0b0b0b") +
  scale_y_continuous(breaks = comp_df$y, labels = comp_df$region, expand = expansion(add = 0.6)) +
  scale_fill_manual(values = comp_colors, name = NULL) +
  labs(tag = "A",
       title = "Composicion del universo total de MR por region",
       subtitle = "Centrado en 'shared' - male-only a la izquierda, female-only a la derecha",
       x = "Numero de reguladores (relativo al centro)", y = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "#52514e"),
        plot.tag = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank(),
        axis.text.x = element_blank(),
        legend.position = "top")

# ============================================================================
# Panel B - Direct Spearman correlation: female NES vs male NES, per region
# ============================================================================

corr_rows <- lapply(region_order, function(region_label) {
  m <- all_male[all_male$region == region_label, ]
  f <- all_female[all_female$region == region_label, ]
  shared_ids <- intersect(m$Regulon, f$Regulon)
  m_nes <- m$NES[match(shared_ids, m$Regulon)]
  f_nes <- f$NES[match(shared_ids, f$Regulon)]
  data.frame(region = region_label, n_shared = length(shared_ids),
             rho = cor(m_nes, f_nes, method = "spearman"))
})
corr_df <- do.call(rbind, corr_rows)
corr_df$region <- factor(corr_df$region, levels = rev(region_order))

panel_B <- ggplot(corr_df, aes(x = rho, y = region, fill = rho)) +
  geom_col(width = 0.6) +
  geom_text(aes(label = sprintf("%.2f (n=%d)", rho, n_shared)),
            hjust = -0.08, size = 3, color = "#0b0b0b") +
  scale_fill_gradient(low = "#cde2fb", high = "#104281", limits = c(0, 1), guide = "none") +
  scale_x_continuous(limits = c(0, 1.15), expand = c(0, 0)) +
  labs(tag = "B",
       title = "Correlacion de Spearman: NES mujeres vs NES hombres, por region",
       x = expression(rho), y = NULL) +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.tag = element_text(face = "bold", size = 14),
        axis.text.y = element_text(size = 8),
        panel.grid.minor = element_blank(),
        panel.grid.major.y = element_blank())

# ============================================================================
# Panel C - NES_male vs NES_female, by region (kept for later, not drawn now)
# ============================================================================

shared   <- read.delim(shared_file, header = TRUE, stringsAsFactors = FALSE)
specific <- read.delim(specific_file, header = TRUE, stringsAsFactors = FALSE)

cat_levels <- c("Concordant (shared)", "Discordant (shared)", "Male-specific", "Female-specific")
color_pal  <- c("Concordant (shared)" = "#c3c2b7",
                 "Discordant (shared)" = "#1baf7a",
                 "Male-specific"       = "#2a78d6",
                 "Female-specific"     = "#eb6834")

pts_shared <- data.frame(
  region = shared$region,
  x = shared$NES_male, y = shared$NES_female,
  category = ifelse(shared$concordance == "discordant", "Discordant (shared)", "Concordant (shared)")
)
male_only_df   <- specific[specific$sex_specific == "male_only", ]
female_only_df <- specific[specific$sex_specific == "female_only", ]
pts_male_only   <- data.frame(region = male_only_df$region,   x = male_only_df$NES,   y = 0, category = "Male-specific")
pts_female_only <- data.frame(region = female_only_df$region, x = 0, y = female_only_df$NES, category = "Female-specific")

pts_all <- rbind(pts_shared, pts_male_only, pts_female_only)
pts_all$category <- factor(pts_all$category, levels = cat_levels)
pts_all$region <- factor(pts_all$region, levels = region_order)
pts_all <- pts_all[order(pts_all$category != "Concordant (shared)"), ]

panel_C <- ggplot(pts_all, aes(x = x, y = y, color = category)) +
  geom_vline(xintercept = 0, color = "#c3c2b7", linewidth = 0.4) +
  geom_hline(yintercept = 0, color = "#c3c2b7", linewidth = 0.4) +
  geom_point(alpha = 0.6, size = 1.4) +
  facet_wrap(~region, ncol = 3, scales = "free") +
  scale_color_manual(values = color_pal, name = NULL, drop = FALSE) +
  labs(tag = "C",
       title = "NES hombres vs NES mujeres, por region",
       subtitle = "Sexo-exclusivos sobre el eje del sexo donde si se probaron (no tienen NES en el otro sexo)",
       x = "NES (hombres)", y = "NES (mujeres)") +
  theme_minimal(base_size = 10) +
  theme(plot.title = element_text(size = 10, face = "bold"),
        plot.subtitle = element_text(size = 8, color = "#52514e"),
        plot.tag = element_text(face = "bold", size = 14),
        panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold"),
        legend.position = "top")

# ============================================================================
# Compose A + B only, for now
# ============================================================================

draw_composite <- function() {
  grid.newpage()
  grid.text("Comparacion de master regulators entre sexos y regiones cerebrales",
            y = unit(0.985, "npc"), gp = gpar(fontface = "bold", fontsize = 13))
  layout <- grid.layout(nrow = 1, ncol = 2)
  pushViewport(viewport(y = unit(0.47, "npc"), height = unit(0.92, "npc"), layout = layout))
  print(panel_A, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
  print(panel_B, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
  popViewport()
}

out_pdf <- file.path(mra_dir, "three_panel_summary.pdf")
out_png <- file.path(mra_dir, "three_panel_summary.png")

pdf(out_pdf, width = 12, height = 6)
draw_composite()
dev.off()

png(out_png, width = 12, height = 6, units = "in", res = 150)
draw_composite()
dev.off()

cat("Guardado:", out_pdf, "\n")
cat("Guardado:", out_png, "\n")
