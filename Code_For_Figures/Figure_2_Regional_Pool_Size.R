#!/usr/bin/env Rscript
# ==========================================================
# Script: Figure_2_Regional_Pool_Size.R
#
# Description:
#   Generates Fig. 2: Local species richness vs regional pool size,
#   across cases (JCEs only, HP only, both with varying autocorr).
#
# Input:
#   HP_JCE_Sims/Outputs/Fig_SpeciesPoolSize/simulation_outputs_Fig_Nspecies/Richness_vs_nspecies.csv
#
# Outputs:
#   Figures/Fig2_Richness_vs_nspecies.svg
#   Figures/Fig2_Richness_vs_nspecies.pdf
#
# Usage:
#   Rscript Code_For_Figures/Figure_2_Regional_Pool_Size.R
# ==========================================================

library(tidyverse)
library(grid)  # unit(), margin()

# ---- Portable path helpers ----
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

script_dir  <- get_script_dir()                           # .../Hab_Par_New/Code_For_Figures
repo_dir    <- normalizePath(file.path(script_dir, "..")) # .../Hab_Par_New
data_dir    <- file.path(repo_dir, "HP_JCE_Sims", "Outputs", "Fig_SpeciesPoolSize", "simulation_outputs_Fig_Nspecies")
figures_dir <- file.path(repo_dir, "Figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# === File path (relative) ===
f <- file.path(data_dir, "Richness_vs_nspecies.csv")

# === Load & prep ===
df <- read_csv(f, show_col_types = FALSE) %>%
  mutate(
    nsp_str = as.numeric(nsp_str),
    case_pretty = recode(
      case_label,
      "No_HP_sigma_M3"            = "JCEs only",
      "sigma_0p02_none"           = "HP only",
      "sigma_0p02_M3_lam_10p0"    = "Both, low autocorr.",
      "sigma_0p02_M3_lam_0p1"     = "Both, med autocorr.",
      "sigma_0p02_M3_lam_0p0"     = "Both, high autocorr.",
      .default = case_label
    )
  )

# Desired legend order
legend_order <- c("Both, high autocorr.",
                  "Both, med autocorr.",
                  "Both, low autocorr.",
                  "JCEs only",
                  "HP only")

df <- df %>%
  mutate(case_pretty = factor(case_pretty, levels = legend_order))

# Define custom colors (left as provided)
color_vals <- c(
  "Both, high autocorr." = "#7E4794",
  "Both, med autocorr."  = "#d95f02",
  "Both, low autocorr."  = "#66a61e",
  "JCEs only"            = "#1a80bb",
  "HP only"              = "#f2c45f"
)

# Plot
p <- ggplot(df, aes(x = nsp_str, y = species_richness,
                    color = case_pretty, group = case_pretty)) +
  geom_line(linewidth = 1, alpha = .7) +
  geom_point(size = 3.5) +
  scale_color_manual(values = color_vals) +
  labs(
    x = "Regional species pool size",
    y = "Local species richness (N)",
    color = NULL
  ) +
  theme_minimal() +
  theme(
    panel.grid        = element_blank(),
    axis.line.x       = element_line(color = "black", linewidth = 2),
    axis.line.y       = element_line(color = "black", linewidth = 2),
    axis.ticks        = element_line(color = "black", linewidth = 1.6),
    axis.ticks.length = unit(0.35, "cm"),
    axis.text         = element_text(size = 20, color = "black"),
    axis.title        = element_text(size = 24, color = "black"),
    plot.margin       = margin(10, 14, 10, 10),
    legend.position   = c(0.02, 0.98),
    legend.justification = c(0, 1),
    legend.background = element_rect(fill = scales::alpha("white", 0.85),
                                     color = "black"),
    legend.key.width  = unit(1.7, "lines"),
    legend.key.height = unit(1.2, "lines"),
    legend.text       = element_text(size = 13, color = "black")
  )

print(p)

# Save to Figures/
svg_out <- file.path(figures_dir, "Fig2_Richness_vs_nspecies.svg")
pdf_out <- file.path(figures_dir, "Fig2_Richness_vs_nspecies.pdf")

ggsave(svg_out, p, width = 8 * .8, height = 10 * .7, dpi = 300)
ggsave(pdf_out, p, width = 7, height = 7)

message("Saved: ", svg_out, " and ", pdf_out)
