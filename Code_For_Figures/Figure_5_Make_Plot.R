#!/usr/bin/env Rscript
# ==========================================================
# Script: Figure_5_Make_Plot.R
#
# Description:
#   Generates Fig. 5 (2×4 panel) for the JC–HP manuscript:
#     • Columns: different parameter regimes
#     • Rows: richness vs. parameter (top), richness vs. κ (bottom)
#
# Inputs:
#   - CSVs in: HP_JCE_Sims/Outputs/Fig_5/
#       * SimSet1_JCstrength/Summary_SimSet1_JCstrength.csv
#       * SimSet2_Dispersal/Summary_SimSet2_Dispersal.csv
#       * SimSet3_HP_only/Summary_SimSet3_HP_only.csv
#       * SimSet4_HP_plus_JCE/Summary_SimSet4_HP_plus_JCE.csv
#
# Outputs:
#   - Figures/Fig5_2x4.svg
#
# Usage:
#   Rscript Code_For_Figures/Figure_5_Make_Plot.R
#
# Notes:
#   Paths are resolved relative to the script’s location,
#   so the script is portable across machines as long as
#   the repo structure is preserved.
# ==========================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(grid)
  library(scales)
  library(viridis)
})

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
fig5_dir    <- file.path(repo_dir, "HP_JCE_Sims", "Outputs", "Fig_5")
figures_dir <- file.path(repo_dir, "Figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- Load & clean helpers ----------
load_and_clean <- function(path) {
  readr::read_csv(path, show_col_types = FALSE) %>%
    mutate(
      param_str = trimws(param_str),
      param_raw = sub("^[^_]*_", "", param_str),
      param     = suppressWarnings(as.numeric(gsub("_", ".", param_raw)))
    )
}

add_param_labels <- function(df) {
  df %>%
    mutate(
      is_global   = tolower(param_raw) == "none" | is.na(param),
      param_lab   = ifelse(is_global, "Global",
                           ifelse(!is.na(param),
                                  format(param, trim = TRUE, scientific = FALSE),
                                  param_raw)),
      .num_order  = ifelse(is_global, Inf, suppressWarnings(as.numeric(param))),
      param_label = factor(param_lab, levels = unique(param_lab[order(.num_order, na.last = TRUE)]))
    ) %>%
    select(-.num_order)
}

make_palette <- function(levels_vec) {
  cols <- viridisLite::magma(length(levels_vec))
  names(cols) <- levels_vec
  cols
}

# ---------- Load data ----------
df1 <- load_and_clean(file.path(fig5_dir, "SimSet1_JCstrength",  "Summary_SimSet1_JCstrength.csv"))
df2 <- load_and_clean(file.path(fig5_dir, "SimSet2_Dispersal",   "Summary_SimSet2_Dispersal.csv"))
df3 <- load_and_clean(file.path(fig5_dir, "SimSet3_HP_only",     "Summary_SimSet3_HP_only.csv"))
df4 <- load_and_clean(file.path(fig5_dir, "SimSet4_HP_plus_JCE", "Summary_SimSet4_HP_plus_JCE.csv"))

# Color-coded versions + palettes
df1_c <- add_param_labels(df1); pal1 <- make_palette(levels(df1_c$param_label))
df2_c <- add_param_labels(df2); pal2 <- make_palette(levels(df2_c$param_label))
df3_c <- add_param_labels(df3); pal3 <- make_palette(levels(df3_c$param_label))
df4_c <- add_param_labels(df4); pal4 <- make_palette(levels(df4_c$param_label))

# ---------- Common y-limits ----------
y_max_all <- max(c(df1$species_richness, df2$species_richness,
                   df3$species_richness, df4$species_richness), na.rm = TRUE)
y_limits <- c(0, 1.1 * y_max_all)

# ---------- Theme ----------
theme_axes_heavy <- function(show_legend = FALSE) {
  theme_minimal() +
    theme(
      legend.position = if (show_legend) "right" else "none",
      panel.grid = element_blank(),
      panel.border = element_blank(),
      axis.line.x = element_line(color = "black", linewidth = 2),
      axis.line.y = element_line(color = "black", linewidth = 2),
      axis.ticks = element_line(color = "black", linewidth = 2),
      axis.ticks.length = unit(0.3, "cm"),
      axis.text  = element_text(size = 20, color = "black"),
      axis.title = element_blank()
    )
}

# ---------- Plot helpers ----------
p_rich_vs_param <- function(df, palette, log_x = FALSE, show_legend = FALSE) {
  d <- df %>% arrange(param)
  p <- ggplot(d, aes(param, species_richness)) +
    geom_line(linetype = "dotted", linewidth = 1, color = "grey40") +
    geom_point(aes(fill = param_label), size = 4.5, shape = 21, stroke = 1.2, color = "black") +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_y_continuous(limits = y_limits, breaks = c(0, 100, 200),
                       expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = NULL, y = NULL) +
    theme_axes_heavy(show_legend)
  if (log_x) p <- p + scale_x_log10()
  p
}

p_rich_vs_kappa <- function(df, breaksx, palette, show_legend = FALSE) {
  ggplot(df %>% arrange(mean_kappa),
         aes(mean_kappa, species_richness)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "red3", size = 1) +
    geom_line(linetype = "dotted", linewidth = 1, color = "grey40") +
    geom_point(aes(fill = param_label), size = 4.5, shape = 21, stroke = 1.2, color = "black") +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_y_continuous(limits = y_limits, breaks = c(0, 100, 200),
                       expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = NULL, y = NULL) +
    theme_axes_heavy(show_legend) +
    scale_x_log10(breaks = breaksx,
                  limits = c(min(min(df$mean_kappa, na.rm = TRUE), .99),
                             max(df$mean_kappa, na.rm = TRUE)))
}

# SimSet2: "Global" as special labeled point with axis break
make_p2A_with_simple_break <- function(df, palette, show_legend = FALSE) {
  finite <- df %>% filter(!is.na(param) & is.finite(param))
  finite_max <- if (nrow(finite) > 0) max(finite$param) else 1
  finite_min <- if (nrow(finite) > 0) min(finite$param) else 0
  span <- max(1e-8, finite_max - finite_min)
  global_val <- finite_max + span / 2
  
  d <- df %>%
    mutate(
      is_global  = tolower(param_raw) == "none" | is.na(param),
      param_plot = ifelse(is_global, global_val, param)
    )
  
  p <- ggplot() +
    geom_line(
      data = d %>% filter(!is_global) %>% arrange(param_plot),
      aes(param_plot, species_richness),
      linetype = "dotted", linewidth = 1, color = "grey40"
    ) +
    geom_point(
      data = d %>% arrange(param_plot),
      aes(param_plot, species_richness, fill = param_label),
      size = 4.5, shape = 21, stroke = 1.2, color = "black"
    ) +
    scale_fill_manual(values = palette, drop = FALSE) +
    scale_x_continuous(
      breaks = c(scales::pretty_breaks(4)(c(finite_min, finite_max)), global_val),
      labels = function(x) ifelse(x == global_val, "Global", x),
      expand = expansion(mult = c(0.06, 0.08))
    ) +
    scale_y_continuous(limits = y_limits, breaks = c(0, 100, 200),
                       expand = expansion(mult = c(0.02, 0.05))) +
    labs(x = NULL, y = NULL) +
    theme_axes_heavy(show_legend) +
    coord_cartesian(clip = "off")
  
  y_pos <- y_limits[1] + 0.02 * diff(y_limits)
  x_break <- (finite_max + global_val) / 2
  
  p + annotate("text", x = x_break, y = y_pos, label = "//", size = 10, vjust = 1, hjust = .7)
}

# ---------- Build panels ----------
p1A <- p_rich_vs_param(df1_c, palette = pal1, log_x = TRUE)
p1B <- p_rich_vs_kappa(df1_c, breaksx = c(1, 1.4, 2), palette = pal1)

p2A <- make_p2A_with_simple_break(df2_c, palette = pal2)
p2B <- p_rich_vs_kappa(df2_c, breaksx = c(1.0, 1.4, 2), palette = pal2)

p3A <- p_rich_vs_param(df3_c, palette = pal3, log_x = TRUE)
p3B <- p_rich_vs_kappa(df3_c, breaksx = c(1, 3, 10), palette = pal3)

p4A <- p_rich_vs_param(df4_c, palette = pal4, log_x = TRUE)
p4B <- p_rich_vs_kappa(df4_c, breaksx = c(1, 1.2, 1.4), palette = pal4)

# ---------- Axis labels ----------
axis_title_theme <- theme(axis.title = element_text(size = 18, color = "black"))

p1A_l <- p1A + labs(x = expression("JC-effect strength " * italic(a)), y = "Species richness") + axis_title_theme
p2A_l <- p2A + labs(x = expression("Dispersal distance " * sigma[D]), y = NULL) + axis_title_theme
p4A_l <- p4A + labs(x = expression("Noise (nugget) of variation"), y = NULL) + axis_title_theme
p3A_l <- p3A + labs(x = expression("Noise (nugget) of variation"), y = NULL) + axis_title_theme

agg_x_lab <- expression("Aggregation " * bar(kappa)[M])
p1B_l <- p1B + labs(x = agg_x_lab, y = "Species richness") + axis_title_theme
p2B_l <- p2B + labs(x = agg_x_lab, y = NULL) + axis_title_theme
p4B_l <- p4B + labs(x = agg_x_lab, y = NULL) + axis_title_theme
p3B_l <- p3B + labs(x = agg_x_lab, y = NULL) + axis_title_theme

# ---------- Column titles ----------
col_title <- function(lbl_expr,
                      fill = "grey80",
                      border = "black",
                      txt_size = 18,
                      pad_w = grid::unit(6, "pt"),
                      pad_h = grid::unit(5, "pt"),
                      radius = grid::unit(2, "pt")) {
  tg <- grid::textGrob(
    lbl_expr, x = 0.5, y = 0.5,
    gp = grid::gpar(col = "black", fontsize = txt_size)
  )
  w <- grid::grobWidth(tg)  + 2 * pad_w
  h <- grid::grobHeight(tg) + 2 * pad_h
  rg <- grid::roundrectGrob(
    x = 0.5, y = 0.5, width = w, height = h, r = radius,
    gp = grid::gpar(fill = fill, col = border, lwd = 1.9)
  )
  patchwork::wrap_elements(grid::grobTree(rg, tg))
}

t1 <- col_title(expression("Only JC-effects; " * italic(a) * " varies"))
t2 <- col_title(expression("Only JC-effects; " * sigma[D] * " varies"))
t3 <- col_title("JC-effect and HP")
t4 <- col_title("Only HP")

# ---------- Assemble 2×4 ----------
fig_5_2x4 <- (t1 | t2 | t3 | t4) /
  (p1A_l | p2A_l | p4A_l | p3A_l) /
  (p1B_l | p2B_l | p4B_l | p3B_l) +
  patchwork::plot_layout(heights = c(0.10, 1, 1)) &
  theme(plot.margin = margin(6, 6, 6, 6))

# ---------- Save to top-level Figures ----------
out_svg <- file.path(figures_dir, "Fig5_2x4.svg")
ggsave(out_svg, fig_5_2x4, width = 16, height = 10, dpi = 300)

message("Saved: ", out_svg)
