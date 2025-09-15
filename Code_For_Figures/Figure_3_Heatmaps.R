#!/usr/bin/env Rscript
# ==========================================================
# Script: Figure_3_Heatmaps.R
#
# Description:
#   Generates Fig. 3 as a 2×2 heatmap panel of species richness.
#   Columns: Low vs High autocorrelation
#   Rows:    Global dispersal (σ_D = ∞) vs Dispersal-limited (σ_D = 1.5)
#
# Inputs (relative to repo root):
#   HP_JCE_Sims/Outputs/Fig_3/simulation_outputs_Fig_3A/Richness_M_sigmah_No_Disp_Lim_Low_Autocorr.csv
#   HP_JCE_Sims/Outputs/Fig_3/simulation_outputs_Fig_3B/Richness_M_sigmah_No_Disp_Lim_High_Autocorr.csv
#   HP_JCE_Sims/Outputs/Fig_3/simulation_outputs_Fig_3C/Richness_M_sigmah_Disp_Lim_sigmaD_1_5_Low_Autocorr.csv
#   HP_JCE_Sims/Outputs/Fig_3/simulation_outputs_Fig_3D/Richness_M_sigmah_Disp_Lim_sigmaD_1_5_High_Autocorr.csv
#
# Outputs:
#   Figures/Fig3_2x2_RichnessHeatmaps.svg
#   Figures/Fig3_2x2_RichnessHeatmaps.png
#
# Usage:
#   Rscript Code_For_Figures/Figure_3_Heatmaps.R
#   # or source() from RStudio; paths resolve relative to this script
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(patchwork)  # layout + shared legend
  library(grid)       # unit(), margins
})

# ---- File paths (portable & simple) ----
get_script_dir <- function() {
  # 1) Rscript: capture --file=...
  args <- commandArgs(trailingOnly = FALSE)
  i <- grep("^--file=", args)
  if (length(i)) return(dirname(normalizePath(sub("^--file=", "", args[i]))))
  # 2) RStudio: active editor path
  if (interactive() &&
      requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    p <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(p)) return(dirname(normalizePath(p)))
  }
  # 3) Fallback: current working directory
  normalizePath(getwd())
}

script_dir <- get_script_dir()                           # .../Hab_Par_New/Code_For_Figures
repo_dir   <- normalizePath(file.path(script_dir, "..")) # .../Hab_Par_New
fig3_dir   <- file.path(repo_dir, "HP_JCE_Sims", "Outputs", "Fig_3")

files <- c(
  A = file.path(fig3_dir, "simulation_outputs_Fig_3A", "Richness_M_sigmah_No_Disp_Lim_Low_Autocorr.csv"),
  C = file.path(fig3_dir, "simulation_outputs_Fig_3B", "Richness_M_sigmah_No_Disp_Lim_High_Autocorr.csv"),
  G = file.path(fig3_dir, "simulation_outputs_Fig_3C", "Richness_M_sigmah_Disp_Lim_sigmaD_1_5_Low_Autocorr.csv"),
  I = file.path(fig3_dir, "simulation_outputs_Fig_3D", "Richness_M_sigmah_Disp_Lim_sigmaD_1_5_High_Autocorr.csv")
)

# ---- Reader ----
read_sim_csv <- function(path) {
  df_try <- try(read.csv(path, header = TRUE, stringsAsFactors = FALSE), silent = TRUE)
  df <- if (inherits(df_try, "try-error")) read.csv(path, header = FALSE, stringsAsFactors = FALSE) else df_try
  expected <- c("species_richness","sigma_str","lambda_str","neigh_str")
  if (!all(expected %in% names(df))) {
    if (ncol(df) >= 4) names(df)[1:4] <- expected
  }
  df[, expected, drop = FALSE]
}

# ---- Cleaning & ordering ----
clean_labels <- function(df) {
  df <- df %>%
    mutate(
      neigh_str = as.character(neigh_str),
      sigma_str = as.character(sigma_str),
      
      # neigh_str mapping
      neigh_str = dplyr::case_when(
        neigh_str == "none" ~ "No JCEs",
        neigh_str == "0"    ~ "1",
        neigh_str == "1"    ~ "3",
        neigh_str == "2"    ~ "5",
        neigh_str == "3"    ~ "7",
        TRUE ~ neigh_str
      ),
      
      # sigma_str: map to labels
      sigma_label = ifelse(sigma_str == "9_999999e6", "No HP", gsub("_", ".", sigma_str))
    )
  
  # x-axis order
  neigh_levels <- c("No JCEs", "1", "3", "5", "7")
  df$neigh_str <- factor(df$neigh_str, levels = neigh_levels)
  
  # y-axis order: top→bottom = "No HP", 0.2, 0.1, 0.05, 0.02
  # ggplot interprets factor levels bottom→top; so set: 0.02, 0.05, 0.1, 0.2, "No HP"
  target_bottom_to_top <- rev(c("0.02", "0.05", "0.1", "0.2", "No HP"))
  present_levels <- target_bottom_to_top[target_bottom_to_top %in% unique(df$sigma_label)]
  df$sigma_label <- factor(df$sigma_label, levels = present_levels)
  
  df
}

# ---- Build a single panel ----
plot_heatmap <- function(df) {
  ggplot(df, aes(x = neigh_str, y = sigma_label, fill = species_richness)) +
    geom_tile(color = "black") +
    geom_text(aes(label = round(species_richness, 0)),
              color = "grey70", fontface = "bold", size = 3.8) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "JC-effect Moore nighborhood (M)",
         y = expression("Habitat specialization" ~ (sigma[h]))) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(colour = "black", fill = NA, size = 2),
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(size = 14, color = "black"),
      plot.background = element_rect(fill = "white"),
      legend.key.height = unit(2.2, "cm"),
      legend.key.width  = unit(.35, "cm"),
      axis.title.x = element_text(size = 18, margin = margin(t = 8)),
      axis.title.y = element_text(size = 18),
      legend.text = element_text(size = 18)
    )
}

# ---- Load, clean, global scale ----
dfs <- lapply(files, function(f) clean_labels(read_sim_csv(f)))
rich_min <- min(sapply(dfs, function(d) min(d$species_richness, na.rm = TRUE)))
rich_max <- max(sapply(dfs, function(d) max(d$species_richness, na.rm = TRUE)))
breaks   <- pretty(c(rich_min, rich_max), n = 4)

fill_scale <- scale_fill_viridis(
  option = "D",
  limits = c(rich_min, rich_max),
  breaks = breaks,
  name   = "Species richness"  # legend title (will be vertical on LEFT, reading downward)
)

pA <- plot_heatmap(dfs$A) + fill_scale
pC <- plot_heatmap(dfs$C) + fill_scale
pG <- plot_heatmap(dfs$G) + fill_scale
pI <- plot_heatmap(dfs$I) + fill_scale
pA <- pA + theme(axis.title.x = element_blank())
pC <- pC + theme(axis.title.x = element_blank())
pC <- pC + theme(axis.title.y = element_blank())
pI <- pI + theme(axis.title.y = element_blank())

# ---- Column titles (top), Row labels (RIGHT side) ----
col_titles <- c("Low autocorrelation", "High autocorrelation")
row_titles <- c(bquote("Global dispersal: "~sigma[D] == infinity), bquote("Disperal limited: "~sigma[D] == 1.5))

lab_plot <- function(txt, angle = 0, parse = FALSE,
                     size = 7, margin = margin(),
                     fill = "grey80", width = 1.9, height = 1.9,
                     x = 0.5, y = 0.5, hjust = 0.5) {
  ggplot() +
    theme_void() +
    theme(plot.margin = margin) +
    annotate("label", x = x, y = y, label = txt,
             size = size, angle = angle, parse = parse,
             fill = fill, label.size = NA,
             label.padding = unit(0.6, "lines"),
             hjust = hjust)
}

top_low   <- lab_plot(col_titles[1], margin = margin(b = -6), x = 0.5, hjust = .575)
top_high  <- lab_plot(col_titles[2], margin = margin(b = -6), x = 0.5, hjust = .4)
right_row1 <- lab_plot(deparse(row_titles[1]), angle = -90, parse = TRUE,
                       margin = margin(l = -18), size = 5.5)
right_row2 <- lab_plot(deparse(row_titles[2]), angle = -90, parse = TRUE,
                       margin = margin(l = -18), size = 5.5)

# ---- Assemble: columns = [pA, pC, row_label], rows stacked; header over first two cols ----
header_row <- top_low + top_high + plot_spacer() + plot_layout(widths = c(1, 1, 0.12))
row1 <- pA + pC + right_row1 + plot_layout(widths = c(1, 1, 0.12))
row2 <- pG + pI + right_row2 + plot_layout(widths = c(1, 1, 0.12))

final_plot <-
  header_row / row1 / row2 +
  plot_layout(heights = c(0.08, 1, 1), guides = "collect") &
  theme(
    legend.position = "right",
    # Legend title on LEFT of bar, rotated to read downward
    legend.title = element_text(angle = -90, vjust = 0.5, hjust = 0.5,size=18)
  )

final_plot <- final_plot & guides(
  fill = guide_colorbar(
    title.position = "left",
    title.vjust    = 1,
    label.position = "right",
    direction      = "vertical", 
    ticks.colour = NA
  )
)

# View
print(final_plot)

# Save
figures_dir <- file.path(repo_dir, "Figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

svg_out <- file.path(figures_dir, "Fig3_2x2_RichnessHeatmaps.svg")
png_out <- file.path(figures_dir, "Fig3_2x2_RichnessHeatmaps.png")
ggsave(svg_out, plot = final_plot, width = 10, height = 8.5, units = "in", dpi = 300)
ggsave(png_out, plot = final_plot, width = 10, height = 8.5, units = "in", dpi = 300)

message("Saved: ", svg_out, " and ", png_out)
