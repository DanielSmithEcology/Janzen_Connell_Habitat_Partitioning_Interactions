#!/usr/bin/env Rscript
# ==========================================================
# Script:Figure_4_heatmaps_scatter.R
#
# Description:
#   Generates Fig. 4 (1×3 panel) for the JC–HP manuscript:
#     • Left: species richness heatmap
#     • Middle: scaled covariance heatmap
#     • Right: scatterplot of richness vs. covariance
#
# Inputs:
#   - CSV: HP_JCE_Sims/Outputs/Fig_4/simulation_outputs_Fig_4/Fig_4AB/Summary_Fig_4AB.csv
#
# Outputs:
#   - Figures/Fig4_1x3_Richness_Covariance_Scatter.svg
#   - Figures/Fig4_1x3_Richness_Covariance_Scatter.png
#
# Usage:
#   Rscript Code_For_Figures/Figure_4_heatmaps_scatter.R
#
# Notes:
#   Paths are resolved relative to the script’s location,
#   so the script is portable across machines once the repo
#   structure is preserved.
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(grid)
  library(scales)
})

# ---- Portable path helpers (minimal + robust) ----
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  i <- grep("^--file=", args)
  if (length(i)) return(dirname(normalizePath(sub("^--file=", "", args[i]))))
  if (interactive() &&
      requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()                           # e.g., .../Hab_Par_New/Code_For_Figures
repo_dir   <- normalizePath(file.path(script_dir, "..")) # .../Hab_Par_New
fig4_dir   <- file.path(repo_dir, "HP_JCE_Sims", "Outputs", "Fig_4", "simulation_outputs_Fig_4", "Fig_4AB")
figures_dir <- file.path(repo_dir, "Figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ---- Load data from file  ----
df <- read.csv(
  file.path(fig4_dir, "Summary_Fig_4AB.csv"),
  stringsAsFactors = FALSE
)

# Expect columns: species_richness, Sigma_str, range_str, scaled_cov
# ---- Clean + order factors ----
df <- df %>%
  mutate(
    Sigma = as.numeric(gsub("_", ".", Sigma_str)),
    Range = as.numeric(gsub("_", ".", range_str))
  )

# y (Sigma): top→bottom = high→low; x (Range): left→right = low→high
sigma_levels <- sort(unique(df$Sigma), decreasing = TRUE)
range_levels <- sort(unique(df$Range), decreasing = FALSE)

df <- df %>%
  mutate(
    Sigma_f = factor(Sigma, levels = sigma_levels),
    Range_f = factor(Range, levels = range_levels)
  )

# ---- Common axis labels ----
x_lab <- "Range of habitat variation"
y_lab <- "Nugget (noise) of variation"

# ---- Shared legend style helper (vertical, title on left) ----
legend_style <- list(
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 12, angle = -90), # vertical, reads downward
    legend.text       = element_text(size = 11),
    legend.key.height = unit(2, "cm"),
    legend.key.width  = unit(.35, "cm")
  ),
  guides(
    fill = guide_colorbar(
      title.position = "left",
      title.vjust    = 1,
      label.position = "right",
      direction      = "vertical",
      ticks.colour   = NA
    )
  )
)

# ===================================================================
# Left heatmap: Species richness
# ===================================================================
r_min    <- min(df$species_richness, na.rm = TRUE)
r_max    <- max(df$species_richness, na.rm = TRUE)
r_breaks <- pretty(c(r_min, r_max), n = 4)

p_rich <- ggplot(df, aes(x = Range_f, y = Sigma_f, fill = species_richness)) +
  geom_tile(color = "black") +
  geom_text(aes(label = species_richness),
            color = "grey70", fontface = "bold", size = 3.6) +
  scale_fill_viridis(option = "D", limits = c(r_min, r_max), breaks = r_breaks,
                     name = "Species richness") +
  labs(x = x_lab, y = y_lab) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border     = element_rect(colour = "black", fill = NA, size = 2),
    axis.text.x      = element_text(color = "black", size = 12),
    axis.text.y      = element_text(color = "black", size = 12),
    axis.title.x     = element_text(size = 16, margin = margin(t = 8)),
    axis.title.y     = element_text(size = 16, margin = margin(r = 8))
  ) +
  legend_style

# ===================================================================
# Middle heatmap: Scaled covariance
# ===================================================================
c_min    <- min(df$scaled_cov, na.rm = TRUE)
c_max    <- max(df$scaled_cov, na.rm = TRUE)
c_breaks <- pretty(c(c_min, c_max), n = 4)

p_cov <- ggplot(df, aes(x = Range_f, y = Sigma_f, fill = scaled_cov)) +
  geom_tile(color = "black") +
  geom_text(aes(label = sprintf("%.2f", scaled_cov)),
            color = "grey70", fontface = "bold", size = 3.2) +
  scale_fill_viridis(option = "D", limits = c(c_min, c_max), breaks = c_breaks,
                     name = expression(bar(Cov)[x]~bgroup("(", J(x)~","~~H(x)/mu[H(x)], ")"))) +
  labs(x = x_lab, y = y_lab) +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  theme(
    panel.background = element_rect(fill = "white"),
    panel.border     = element_rect(colour = "black", fill = NA, size = 2),
    axis.text.x      = element_text(color = "black", size = 12),
    axis.text.y      = element_text(color = "black", size = 12),
    axis.title.x     = element_text(size = 16, margin = margin(t = 8)),
    axis.title.y     = element_text(size = 16, margin = margin(r = 8))
  ) +
  legend_style

# ===================================================================
# Right panel: scatter richness vs Cov
# ===================================================================
p_scatter <- ggplot(df, aes(x = scaled_cov, y = species_richness)) +
  geom_point(color = "black", fill = "grey40", shape = 21, size = 4.5, alpha = 1.0) +
  labs(
    x = expression(bar(Cov)[x]~bgroup("(", J(x)~","~~H(x)/mu[H(x)], ")")),
    y = "Species richness"
  ) +
  theme_classic(base_size = 14) +
  theme(
    axis.title.x = element_text(size = 18, margin = margin(t = 8)),
    axis.title.y = element_text(size = 18, margin = margin(r = 8)),
    axis.text.x  = element_text(color = "black", size = 16),
    axis.text.y  = element_text(color = "black", size = 16),
    axis.line    = element_line(color = "black", linewidth = 1.5),
    axis.ticks   = element_line(color = "black", linewidth = 1.5)
  ) +
  scale_x_continuous(breaks = pretty_breaks(n = 3)) +
  scale_y_continuous(breaks = pretty_breaks(n = 3))

# ---- Arrange 1×3 (keep both legends from heatmaps) ----
final_row <- p_rich + p_cov + p_scatter + plot_layout(ncol = 3, guides = "keep")

# ---- Save combined figure to top-level Figures ----
svg_out <- file.path(figures_dir, "Fig4_1x3_Richness_Covariance_Scatter.svg")
png_out <- file.path(figures_dir, "Fig4_1x3_Richness_Covariance_Scatter.png")

ggsave(svg_out, plot = final_row, width = 15, height = 5, units = "in", dpi = 300)
ggsave(png_out, plot = final_row, width = 15, height = 5, units = "in", dpi = 300)

message("Saved: ", svg_out, " and ", png_out)
