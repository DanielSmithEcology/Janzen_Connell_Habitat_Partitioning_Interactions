#!/usr/bin/env Rscript
# ==========================================================
# Script: Figure_3_All.R
#
# Produces:
#   1) Figures/Fig3_2x2_RichnessHeatmaps.{svg,png}
#   2) Figures/Fig3_2x4_Heatmaps_and_Ratios.{svg,png}
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(scales)
  library(stringr)
  library(ggplot2)
  library(viridis)
  library(patchwork)  # layout + shared legend
  library(grid)       # unit(), margins
})

# ---- File paths (portable & simple) ----
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

script_dir <- get_script_dir()                           # .../Code_For_Figures
repo_dir   <- normalizePath(file.path(script_dir, "..")) # repo root
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
      
      # neigh_str mapping to display levels & later numeric M
      neigh_str = dplyr::case_when(
        neigh_str == "none" ~ "No JCEs",
        neigh_str == "0"    ~ "1",
        neigh_str == "1"    ~ "3",
        neigh_str == "2"    ~ "5",
        neigh_str == "3"    ~ "7",
        TRUE ~ neigh_str
      ),
      
      # sigma label for display
      sigma_label = ifelse(sigma_str == "9_999999e6", "No HP", gsub("_", ".", sigma_str))
    )
  
  # x-axis order
  neigh_levels <- c("No JCEs", "1", "3", "5", "7")
  df$neigh_str <- factor(df$neigh_str, levels = neigh_levels)
  
  # y-axis order: top→bottom = "No HP", 0.2, 0.1, 0.05, 0.02
  target_bottom_to_top <- rev(c("0.02", "0.05", "0.1", "0.2", "No HP"))
  present_levels <- target_bottom_to_top[target_bottom_to_top %in% unique(df$sigma_label)]
  df$sigma_label <- factor(df$sigma_label, levels = present_levels)
  
  df
}

# ---- Heatmap panel ----
plot_heatmap <- function(df) {
  ggplot(df, aes(x = neigh_str, y = sigma_label, fill = species_richness)) +
    geom_tile(color = "black") +
    geom_text(aes(label = round(species_richness, 0)),
              color = "grey70", fontface = "bold", size = 3.8) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(x = "JC-effect spatial scale (M)",
         y = expression("Habitat response breadth" ~ (sigma[h]))) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.border = element_rect(colour = "black", fill = NA, size = 2),
      axis.text.x = element_text(color = "black", size = 14),
      axis.text.y = element_text(size = 14, color = "black"),
      plot.background = element_rect(fill = "white"),
      legend.key.height = unit(1.8, "cm"),
      legend.key.width  = unit(.3, "cm"),
      axis.title.x = element_text(size = 14, margin = margin(t = 8)),  # ↓ matched to ratio panels
      axis.title.y = element_text(size = 17),
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
  name   = "Species richness"
)

# Build per-dataset heatmaps (consistent axis titles)
pA <- plot_heatmap(dfs$A) + fill_scale
pC <- plot_heatmap(dfs$C) + fill_scale + theme(axis.title.y = element_blank()) # 2) remove y title on 2nd heatmap
pG <- plot_heatmap(dfs$G) + fill_scale + theme(axis.title.y = element_blank())
pI <- plot_heatmap(dfs$I) + fill_scale + theme(axis.title.y = element_blank())

# ---- Assemble simple 2×2 richness heatmaps (legend title on right side of bar) ----
final_heatmaps <-
  ((pA | pC) / (pG | pI)) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(angle = -90, vjust = 0.5, hjust = 0.5, size = 16) # vertical title
  )

final_heatmaps <- final_heatmaps & guides(
  fill = guide_colorbar(
    title.position = "right",  # 4) title on the *right* side of the vertical bar
    title.vjust    = 0.5,
    label.position = "right",
    direction      = "vertical",
    ticks.colour   = NA
  )
)

# ---- Save 2×2 heatmaps ----
figures_dir <- file.path(repo_dir, "Figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

svg_out <- file.path(figures_dir, "Fig3_2x2_RichnessHeatmaps.svg")
png_out <- file.path(figures_dir, "Fig3_2x2_RichnessHeatmaps.png")
ggsave(svg_out, plot = final_heatmaps, width = 10, height = 8.5, units = "in", dpi = 300)
ggsave(png_out, plot = final_heatmaps, width = 10, height = 8.5, units = "in", dpi = 300)
message("Saved: ", svg_out, " and ", png_out)

# ==========================================================
# SUM OF PARTS vs BOTH  +  RATIO PANELS (2×4)
# ==========================================================

# ---- Prep for Sum-of-Parts vs Both ---------------------------------------
prep_df <- function(d, dataset_name){
  d %>%
    rename(richness = species_richness) %>%
    mutate(
      dataset      = dataset_name,
      no_JCEs_flag = neigh_str == "No JCEs",
      no_HP_flag   = sigma_label == "No HP" | sigma_str == "9_999999e6",
      M            = suppressWarnings(as.numeric(as.character(neigh_str)))  # "1","3","5","7" -> 1,3,5,7 ; NA for "No JCEs"
    )
}

sum_vs_both_one <- function(dd){
  # No HP richness by M
  nohp_M <- dd %>%
    filter(no_HP_flag, !is.na(M)) %>%
    distinct(M, .keep_all = TRUE) %>%
    select(M, richness_nohp = richness)
  
  # No JCEs richness by sigma (one per sigma)
  nojce_sigma <- dd %>%
    filter(no_JCEs_flag, !no_HP_flag) %>%
    distinct(sigma_label, .keep_all = TRUE) %>%
    select(sigma_label, richness_nojce = richness)
  
  # Both operating by (sigma, M)
  both <- dd %>%
    filter(!no_HP_flag, !no_JCEs_flag, !is.na(M)) %>%
    select(sigma_label, M, richness_both = richness)
  
  # Build wide table at the (sigma, M) combos that exist in 'both'
  out <- both %>%
    left_join(nohp_M, by = "M") %>%
    left_join(nojce_sigma, by = "sigma_label") %>%
    mutate(
      richness_sum_parts = richness_nohp + richness_nojce,
      dataset = unique(dd$dataset)
    ) %>%
    select(dataset, sigma_label, M,
           richness_nohp, richness_nojce,
           richness_sum_parts, richness_both)
  out
}

# Construct cmp_wide for all datasets
dfs_prepped <- purrr::imap(dfs, prep_df)
cmp_wide    <- purrr::map_dfr(dfs_prepped, sum_vs_both_one)

# ---- Ratio table + global y-limits ---------------------------------------
cmp_ratio <- cmp_wide %>%
  mutate(ratio = richness_both / richness_sum_parts) %>%
  filter(is.finite(ratio), !is.na(ratio))

y_limits <- range(cmp_ratio$ratio, na.rm = TRUE)
# Optional gentle clipping:
# q <- quantile(cmp_ratio$ratio, c(0.01, 0.99), na.rm = TRUE); y_limits <- c(q[1], q[2])

sigma_keep <- c("0.02","0.05","0.1","0.2")  # which σ lines to draw in ratio panels

ratio_theme <- theme(
  panel.background = element_rect(fill = "white"),
  panel.border = element_rect(colour = "black", fill = NA, size = 2),
  plot.background  = element_rect(fill = "white"),
  axis.text.x = element_text(color = "black", size = 14),
  axis.text.y = element_text(color = "black", size = 14),
  axis.title.x = element_text(size = 14, margin = margin(t = 6)),
  axis.title.y = element_text(size = 18),
  legend.position = "right"
)

build_ratio_plot <- function(ds, show_y_title = FALSE){
  ggplot(
    filter(cmp_ratio, dataset == ds, sigma_label %in% sigma_keep),
    aes(M, ratio, color = sigma_label, group = sigma_label)
  ) +
    geom_hline(yintercept = 1, linetype = "dashed",size=1) +
    geom_point(size=2.5) +
    geom_line(size=.8) +
    scale_x_continuous(breaks = c(1,3,5,7), limits = c(1,7)) +
    coord_cartesian(ylim = y_limits) +
    scale_y_continuous(breaks = scales::pretty_breaks(3)) +
  #  scale_y_log10(
  #    breaks = c(0.5, 1, 2, 4),
  #    labels = c(".5", "1", "2", "4")
  #  ) +
 #   coord_cartesian(ylim = c(0.5, 4)) +
    labs(
      x = "JC-effect spatial scale (M)",
      y = if (show_y_title) "Synergy" else NULL,
      color = "σ"
    ) +
    ratio_theme
}

# Only the left-most ratio plot shows the y-title
prA <- build_ratio_plot("A", show_y_title = TRUE)
prC <- build_ratio_plot("C", show_y_title = FALSE)
prG <- build_ratio_plot("G", show_y_title = FALSE)
prI <- build_ratio_plot("I", show_y_title = FALSE)

# ---- Assemble 2×4 (top heatmaps, bottom ratios), no header/row boxes ----
heat_row <-
  (pA | pC | pG | pI) +
  plot_layout(guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(angle = -90, vjust = 0.5, hjust = 0.5, size = 16)
  )

# Heatmap legend title right side of bar
heat_row <- heat_row & guides(
  fill = guide_colorbar(
    title.position = "right",
    title.vjust    = 0.5,
    label.position = "right",
    direction      = "vertical",
    ticks.colour   = NA
  )
)

# Replace current ratio_row block with this:

ratio_row <-
  (prA | prC | prG | prI) +
  plot_layout(guides = "collect")

# Style ONLY the color legend (second legend, for σ)
ratio_row <- ratio_row &
  guides(
    color = guide_legend(
      title       = expression(sigma[h]),
      title.position = "top",   # title above swatches
      title.hjust    = 0.5,     # centered title
      ncol        = 1,          # 1=vertical, 2=two columns, etc.
      byrow       = TRUE,
      override.aes = list(      # make legend lines clean & thicker
        linetype = 1,
        shape    = NA,          # hide points in the legend
        size     = 1.6,
        alpha    = 1
      ),
      keywidth   = unit(1.3, "lines"),
      keyheight  = unit(0.9, "lines")
    )
  ) &
  theme(
    legend.position   = "right",
    legend.title      = element_text(size = 16),
    legend.text       = element_text(size = 14),
    legend.key.width  = unit(1.3, "lines"),
    legend.key.height = unit(0.9, "lines")
  )

final_combo <- heat_row / ratio_row + plot_layout(heights = c(1, 1))

# Keep σ legend order
final_combo <- final_combo & scale_color_discrete(limits = sigma_keep)

# ---- Save 2×4 combined figure --------------------------------------------
out_svg <- file.path(figures_dir, "Fig3_2x4_Heatmaps_and_Ratios.svg")
out_png <- file.path(figures_dir, "Fig3_2x4_Heatmaps_and_Ratios.png")
ggsave(out_svg, plot = final_combo, width = 16.5, height = 8.5, units = "in", dpi = 600)
ggsave(out_png, plot = final_combo, width = 16.5, height = 8.5, units = "in", dpi = 600)
message("Saved: ", out_svg, " and ", out_png)

# ------------------------------- END ---------------------------------------
