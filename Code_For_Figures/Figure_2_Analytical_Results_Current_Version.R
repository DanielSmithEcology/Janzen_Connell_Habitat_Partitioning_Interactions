#!/usr/bin/env Rscript
# ==========================================================
# Script: Figure_2_Analytical_Results.R
#
# Description:
#   Generates analytical results for Figure 2.
#   Plots stacked bar charts of stabilizing effects:
#     - CNDD (Janzen–Connell effects),
#     - HP (habitat partitioning),
#     - COV (their covariance upper bound).
#
# Inputs: none (fully analytical, parameterized within script)
#
# Outputs:
#   .../Hab_Par_New/Figures/Fig2_AnalyticalResults.svg
#
# Usage:
#   Rscript Figure_2_Analytical_Results.R
# ==========================================================

# === Libraries ===
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)
library(purrr)

# === Functions ===
Mean_CNDD <- function(S, a, M) {
  (1 - exp(-a)) * (M^2) / S
}

Mean_HP <- function(S, sig) {
  rho <- 2 * sig * sqrt(pi) 
  (1/S) * ((1 / rho) - 1)
}

COV_Upper_Bound <- function(S, a, M, sig) {
  rho <- sig * sqrt(2 * pi)
  CNDD_Var <- (1 - exp(-a * M^2)) * sqrt((1/S) - (1/S^2))
  HP_Var   <- sqrt((1 / (sqrt(2) * rho)) - 1)
  CNDD_Var * HP_Var
}

# === Parameters ===
S_vals  <- c(50, 100, 250, 500, 1000)
M_vals  <- c(3, 5)
sig_vals <- c(0.075, 0.05)
a <- 0.5

# === Create parameter dataframe ===
param_grid <- expand.grid(S = S_vals, M = M_vals, sig = sig_vals) %>%
  arrange(M, sig, S) %>%
  mutate(
    CNDD  = mapply(Mean_CNDD, S, a = a, M = M),
    HP    = mapply(Mean_HP,   S, sig),
    COV   = mapply(COV_Upper_Bound, S, a = a, M = M, sig = sig),
    Total = CNDD + HP + COV
  ) %>%
  mutate(sig = factor(sig, levels = c(0.075, 0.05)))

# === Tidy data for plotting ===
plot_data <- param_grid %>%
  pivot_longer(cols = c("CNDD", "HP", "COV"),
               names_to = "mechanism", values_to = "effect") %>%
  mutate(mechanism = factor(mechanism, levels = c("HP", "CNDD", "COV")))

# --- 1) Compute global max stack height ---
global_max <- plot_data %>%
  group_by(M, sig, S) %>%
  summarise(total = sum(effect), .groups = "drop_last") %>%
  summarise(panel_max = max(total), .groups = "drop") %>%
  summarise(global_max = max(panel_max), .groups = "drop") %>%
  pull(global_max)

# --- 2) Panel plotting function ---
make_panel_plot <- function(df, global_max, top_pad_frac = 0.15) {
  y_top     <- global_max * (1 + top_pad_frac)
  rect_ymin <- global_max * (1 + 0.06)
  rect_ymax <- global_max * (1 + 0.12)
  text_y    <- global_max * (1 + 0.09)
  
  x_vals <- sort(unique(df$S))
  x_mid  <- mean(seq_along(x_vals))
  x_min  <- x_mid - 0.7
  x_max  <- x_mid + 0.7
  
  label_text <- deparse(
    bquote(M == .(unique(df$M)) * ", " * sigma[h] * " = " * .(as.numeric(as.character(unique(df$sig)))))
  )
  
  ggplot(df, aes(x = factor(S), y = effect, fill = mechanism)) +
    geom_col(position = position_stack(reverse = TRUE), color = "grey65") +
    scale_fill_manual(values = c("HP" = "#f2c45f", "CNDD" = "#1a80bb", "COV" = "#7E4794")) +
    scale_y_continuous(
      limits = c(-0.008, y_top * 1.05),
      breaks = pretty(c(0, global_max), n = 4),
      expand = c(0, 0)
    ) +
    scale_x_discrete(expand = c(0.18, 0.18)) +
    annotate("rect",
             xmin = x_min - 1.5, xmax = x_max + 1.75,
             ymin = rect_ymin - .015, ymax = rect_ymax + .015,
             fill = "grey85", color = "black", linewidth = 1.5) +
    annotate("text", x = x_mid + 0.1, y = text_y, label = label_text,
             size = 5.8, fontface = "bold", color = "black", parse = TRUE) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.grid = element_blank(),
      axis.title = element_blank(),
      plot.background = element_rect(fill = "white"),
      axis.text.y = element_text(size = 21, color = "black"),
      axis.text.x = element_text(size = 21, color = "black", angle = -40, vjust = -0.1),
      axis.line = element_line(colour = "black", linewidth = 2),
      axis.ticks = element_line(color = "black", linewidth = 2),
      axis.ticks.length = unit(0.3, "cm"),
      legend.position = "none"
    )
}

# --- 3) Build panels with shared scale ---
plots <- split(plot_data, interaction(plot_data$M, plot_data$sig)) %>%
  purrr::map(~ make_panel_plot(.x, global_max))

combined_plot <- wrap_plots(plots, nrow = 1)
print(combined_plot)

# --- Save into unified Figures/ folder ---
script_dir <- normalizePath(dirname(rstudioapi::getSourceEditorContext()$path), winslash = "/")
repo_dir   <- normalizePath(file.path(script_dir, ".."), winslash = "/")
figures_dir <- file.path(repo_dir, "Figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

svg_out <- file.path(figures_dir, "Fig2_AnalyticalResults.svg")
png_out <- file.path(figures_dir, "Fig2_AnalyticalResults.png")

ggsave(svg_out, plot = combined_plot, width = 12.5, height = 5.6, dpi = 400)
ggsave(png_out, plot = combined_plot, width = 12.5, height = 5.6, dpi = 400)

message("Saved: ", svg_out, " and ", png_out)
