#!/usr/bin/env Rscript
# ==========================================================
# Script: Figure_6_Make_Plot.R
#
# Description:
#   Generates Fig. 6 (2×2), JC-only vs JC+HP, with:
#     - JC-only richness panel overlaid with Dispersal-Limited curve (red triangles)
#     - JC-only functional panel comparing TWO κ-adjusted functions:
#         σ_D = ∞   (no dispersal limit; uses Summary_JCE_Only_Kappa_Values_no_disp_lim.csv)
#         σ_D = 1.0 (dispersal-limited; uses Summary_JCE_Only_Disp_Lim_Kappa_Values.csv)
#       Global (M=99) is forced to 25 and any negative values are removed (NA)
#     - HP panels grouped by range_str
#
# Inputs (relative to repo root):
#   HP_JCE_Sims/Outputs/Fig_6/simulation_outputs_Fig_6/Fig_6_JCE_Only/Summary_JCE_Only.csv
#   HP_JCE_Sims/Outputs/Fig_6/simulation_outputs_Fig_6/Fig_6_JCE_HP/Summary_JCE_HP.csv
#   HP_JCE_Sims/Outputs/Fig_6/simulation_outputs_Fig_6/Fig_6_JCE_Only_Disp_Lim/Summary_JCE_Only_Disp_Lim.csv
#   HP_JCE_Sims/Outputs/Fig_6/simulation_outputs_Fig_6/Fig_6_JCE_Only_Disp_Lim/Summary_JCE_Only_Disp_Lim_Kappa_Values.csv
#   HP_JCE_Sims/Outputs/Fig_6/simulation_outputs_Fig_6/Fig_6_JCE_Only/Summary_JCE_Only_Kappa_Values_no_disp_lim.csv
#
# Outputs:
#   Figures/Fig6_2x2_grouped.svg
#   Figures/Fig6_2x2_grouped.png
# ==========================================================

# Packages
suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(grid)
  library(scales)
})

# ---------- Portable paths ----------
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
fig6_base   <- file.path(repo_dir, "HP_JCE_Sims", "Outputs", "Fig_6", "simulation_outputs_Fig_6")
figures_dir <- file.path(repo_dir, "Figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ---------- File paths ----------
f_only            <- file.path(fig6_base, "Fig_6_JCE_Only",          "Summary_JCE_Only.csv")
f_hp              <- file.path(fig6_base, "Fig_6_JCE_HP",            "Summary_JCE_HP.csv")
f_only_DispLim    <- file.path(fig6_base, "Fig_6_JCE_Only_Disp_Lim", "Summary_JCE_Only_Disp_Lim.csv")
f_kappa_disp      <- file.path(fig6_base, "Fig_6_JCE_Only_Disp_Lim", "Summary_JCE_Only_Disp_Lim_Kappa_Values.csv")
f_kappa_nodisp    <- file.path(fig6_base, "Fig_6_JCE_Only",          "Summary_JCE_Only_Kappa_Values_no_disp_lim.csv")

# ---------- Load data ----------
d_only <- readr::read_csv(f_only, show_col_types = FALSE)
d_hp   <- readr::read_csv(f_hp,   show_col_types = FALSE) %>%
  mutate(range_str = factor(range_str, levels = c("0_5", "0_1", "0_05")))

# Dispersal-limited JC-only series (for richness overlay)
d_only_disp <- readr::read_csv(f_only_DispLim, show_col_types = FALSE)

# κ values (M, Kappa): dispersal-limited and no-dispersal-limit
d_kappa_disp   <- readr::read_csv(f_kappa_disp,   show_col_types = FALSE) %>% mutate(M = as.integer(M))
d_kappa_nodisp <- readr::read_csv(f_kappa_nodisp, show_col_types = FALSE) %>% mutate(M = as.integer(M))

# ---------- Axis break helper ----------
prep_global <- function(df, global_offset_factor = 0.35) {
  finite <- df %>% filter(M != 99)
  finite_max <- max(finite$M, na.rm = TRUE)
  finite_min <- min(finite$M, na.rm = TRUE)
  span <- max(1e-6, finite_max - finite_min)
  global_val <- finite_max + global_offset_factor * span
  
  df_aug <- df %>%
    mutate(
      is_global  = M == 99,
      M_plot     = ifelse(is_global, global_val, M),
      M_label    = ifelse(is_global, "Global", as.character(M))
    ) %>%
    arrange(M_plot)
  
  list(
    global_val = global_val,
    finite_min = finite_min,
    finite_max = finite_max,
    df = df_aug
  )
}

# ---------- Theme ----------
theme_axes_heavy <- function() {
  theme_minimal() +
    theme(
      panel.grid        = element_blank(),
      axis.line.x       = element_line(color = "black", linewidth = 2),
      axis.line.y       = element_line(color = "black", linewidth = 2),
      axis.ticks        = element_line(color = "black", linewidth = 1.6),
      axis.ticks.length = unit(0.35, "cm"),
      axis.text         = element_text(size = 20, color = "black"),
      axis.title        = element_text(size = 24, color = "black"),
      plot.margin       = margin(10, 14, 10, 10)
    )
}

# ---------- Plot with global break (supports grouping) ----------
plot_with_global <- function(df_info, yvar, ylab,
                             slash_pos = 0.55,
                             gap_start_frac = 0.475,
                             gap_end_frac   = 0.61,
                             x_breaks = c(1, 7, 13, 19),
                             group_var = NULL,
                             color_values = NULL,
                             shape_values = NULL,
                             legend_title = NULL,   # can be character: e.g. "sigma[D]"
                             bridge_start = NULL,
                             bridge_end   = NULL,
                             labels_map   = NULL) { # named CHAR vector: c(key="sigma[D]==infinity", ...)
  
  df <- df_info$df
  global_val <- df_info$global_val
  finite_max <- df_info$finite_max
  
  ysym <- rlang::sym(yvar)
  brks <- unique(c(x_breaks, global_val))
  p <- ggplot()
  
  # Lines (finite only)
  if (!is.null(group_var)) {
    p <- p +
      geom_line(
        data = df %>% filter(!is_global),
        aes(x = M_plot, y = !!ysym, color = .data[[group_var]], group = .data[[group_var]]),
        linetype = "dotted", linewidth = 0.8
      )
  } else {
    p <- p +
      geom_line(
        data = df %>% filter(!is_global),
        aes(x = M_plot, y = !!ysym),
        linetype = "dotted", linewidth = 0.8, color = "grey35"
      )
  }
  
  # Scales/theme
  p <- p +
    scale_x_continuous(
      breaks = brks,
      labels = function(x) ifelse(x == global_val, "Global", x),
      expand = expansion(mult = c(0.06, 0.08))
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.02))) +
    labs(x = "Neighborhood size  M", y = ylab) +
    theme_axes_heavy() +
    theme(axis.line.x = element_blank()) +
    coord_cartesian(clip = "off")
  
  # Legend
  if (!is.null(group_var)) {
    if (is.null(color_values)) {
      color_values <- c("0_5"  = "#F59E0B",
                        "0_1"  = "#10B981",
                        "0_05" = "#3B82F6")
    }
    if (is.null(shape_values)) {
      shape_values <- c("0_5" = 15, "0_1" = 17, "0_05" = 16)
    }
    
    lvls <- levels(df[[group_var]])
    if (is.null(lvls)) lvls <- unique(df[[group_var]])
    
    if (!is.null(labels_map)) {
      if (is.character(labels_map)) {
        lbl_chr <- unname(labels_map[lvls]); lbls <- parse(text = lbl_chr)
      } else if (is.list(labels_map)) {
        lbls <- as.expression(labels_map[lvls])
      } else if (is.expression(labels_map)) {
        lbls <- labels_map[lvls]
      } else stop("labels_map must be character, list, or expression.")
    } else {
      lbls <- gsub("_", ".", lvls, fixed = TRUE)
    }
    
    p <- p +
      scale_color_manual(values = color_values, limits = lvls, labels = lbls, name = NULL) +
      scale_shape_manual(values = shape_values, limits = lvls, labels = lbls, name = NULL) +
      guides(color = guide_legend(override.aes = list(linetype = 0, size = 4.2)),
             shape = guide_legend(override.aes = list(linetype = 0, size = 4.2))) +
      theme(
        legend.position      = c(0.98, 0.98),
        legend.justification = c("right", "top"),
        legend.background    = element_rect(fill = alpha("white", 0.75), color = NA),
        legend.title         = element_text(size = 16),
        legend.text          = element_text(size = 14)
      )
    
    if (!is.null(legend_title)) {
      lt <- if (is.character(legend_title)) parse(text = legend_title)[[1]] else legend_title
      p <- p + labs(color = lt, shape = lt)
    }
  }
  
  # Gap geometry & connectors
  yrng <- range(df[[yvar]], na.rm = TRUE)
  yrng_span <- diff(yrng); if (!is.finite(yrng_span) || yrng_span == 0) yrng_span <- max(1e-6, abs(yrng[1]))
  y_line <- yrng[1] - 0.06 * yrng_span
  
  gap_span <- (global_val - finite_max)
  x_gap_start <- finite_max + gap_start_frac * gap_span
  x_gap_end   <- finite_max + gap_end_frac   * gap_span
  
  if (is.null(bridge_start)) bridge_start <- gap_start_frac
  if (is.null(bridge_end))   bridge_end   <- gap_end_frac
  bridge_start <- max(0, min(1, bridge_start))
  bridge_end   <- max(0, min(1, bridge_end))
  x_bridge_start <- finite_max + bridge_start * gap_span
  x_bridge_end   <- finite_max + bridge_end   * gap_span
  
  if (!is.null(group_var)) {
    last_pts <- df %>% filter(!is_global) %>%
      group_by(.data[[group_var]]) %>%
      filter(M_plot == max(M_plot, na.rm = TRUE)) %>%
      slice_tail(n = 1) %>%
      ungroup() %>%
      transmute(group = .data[[group_var]], x1 = M_plot, y1 = !!ysym)
    
    glob_pts <- df %>% filter(is_global) %>%
      group_by(.data[[group_var]]) %>%
      summarize(group = first(.data[[group_var]]), x2 = global_val, y2 = first(!!ysym), .groups = "drop")
    
    seg <- inner_join(last_pts, glob_pts, by = "group") %>%
      mutate(
        x_start = x_bridge_start, x_end = x_bridge_end,
        y_start = y1 + (y2 - y1) * ((x_start - x1) / (x2 - x1)),
        y_end   = y1 + (y2 - y1) * ((x_end   - x1) / (x2 - x1))
      ) %>%
      filter(is.finite(y_start), is.finite(y_end))
    seg[[group_var]] <- seg$group
    
    p <- p + geom_segment(
      data = seg,
      aes(x = x_start, y = y_start, xend = x_end, yend = y_end, color = .data[[group_var]]),
      linetype = "dotted", linewidth = 0.9
    )
  } else {
    last_pt <- df %>% filter(!is_global) %>% arrange(desc(M_plot)) %>% slice(1) %>% transmute(x1 = M_plot, y1 = !!ysym)
    glob_pt <- df %>% filter(is_global) %>% transmute(x2 = global_val, y2 = !!ysym) %>% slice(1)
    if (nrow(last_pt) == 1 && nrow(glob_pt) == 1 && (glob_pt$x2 - last_pt$x1) != 0) {
      y_start <- last_pt$y1 + (glob_pt$y2 - last_pt$y1) * ((x_bridge_start - last_pt$x1) / (glob_pt$x2 - last_pt$x1))
      y_end   <- last_pt$y1 + (glob_pt$y2 - last_pt$y1) * ((x_bridge_end   - last_pt$x1) / (glob_pt$x2 - last_pt$x1))
      p <- p + geom_segment(aes(x = x_bridge_start, y = y_start, xend = x_bridge_end, yend = y_end),
                            linetype = "dotted", linewidth = 0.9, color = "grey35")
    }
  }
  
  # Points on top
  if (!is.null(group_var)) {
    size_values <- setNames(rep(4.0, length(lvls)), lvls)
    if ("0_5" %in% names(size_values)) size_values["0_5"] <- 5.0
    for (lvl in lvls) {
      df_sub <- df %>% filter(.data[[group_var]] == lvl)
      p <- p +
        geom_point(
          data = df_sub,
          aes(x = M_plot, y = !!ysym, color = .data[[group_var]], shape = .data[[group_var]]),
          size = size_values[[lvl]], stroke = 1.0
        )
    }
  } else {
    p <- p +
      geom_point(
        data = df,
        aes(x = M_plot, y = !!ysym),
        size = 4.6, shape = 21, stroke = 1.3, color = "black", fill = "grey85"
      )
  }
  
  # Draw broken axis line & slashes
  p +
    annotate("segment", x = -Inf, xend = x_gap_start, y = y_line, yend = y_line,
             linewidth = 2.2, color = "black", lineend = "butt") +
    annotate("segment", x = x_gap_end, xend = Inf, y = y_line, yend = y_line,
             linewidth = 2.2, color = "black", lineend = "butt") +
    annotate("text", x = finite_max + slash_pos * (global_val - finite_max),
             y = y_line, label = "//", size = 11, vjust = 0.5, hjust = 0.5)
}

# ---------- Prep datasets ----------
d_only_info <- prep_global(d_only, global_offset_factor = 0.35)
d_hp_info   <- prep_global(d_hp,   global_offset_factor = 0.35)

# JC-only richness OVERLAY: σ_D ∞ (grey) vs 1.0 (red)
d_only$disp_str      <- factor("infty", levels = c("infty", "1_0"))
d_only_disp$disp_str <- factor("1_0",  levels = c("infty", "1_0"))
d_only_both <- bind_rows(d_only, d_only_disp) %>%
  mutate(fun_val = (1 - exp(-a)) * (M^2))
d_only_both_info <- prep_global(d_only_both, global_offset_factor = 0.35)

# JC-only: baseline fun_val and Global override for reference (used for y-label only)
d_only_info$df <- d_only_info$df %>%
  mutate(fun_val = (1 - exp(-a)) * (M^2))
d_only_info$df$fun_val[d_only_info$df$M == 99] <- 25

# ---------- Build panels ----------
# (A) JC-only richness with σ_D overlay
p_only_rich <- plot_with_global(
  d_only_both_info, "species_richness", "Species richness",
  gap_start_frac = 0.47, gap_end_frac = 0.61,
  bridge_start   = 0.02, bridge_end   = 0.98,
  group_var      = "disp_str",
  color_values   = c("infty" = "grey35", "1_0" = "#E11D48"),
  shape_values   = c("infty" = 16,       "1_0" = 17),
  legend_title   = "sigma[D]",
  labels_map     = c("infty" = "sigma[D]==infinity", "1_0" = "sigma[D]==1.0")
) +
  theme(
    legend.position      = c(0.98, 0.86),
    legend.justification = c("right", "top")
  )

# ---------- Functional panel: TWO κ-adjusted curves (no baseline line)
# Compute fun_plot = a * M^2 * (1 - (a/2) * Kappa) for both σ_D variants; force Global=25; drop negatives
fun_infty <- d_only_info$df %>%                     # contains M, a, M_plot, etc.
  left_join(d_kappa_nodisp, by = "M") %>%           # Kappa from no-disp-limit run
  transmute(M, M_plot, is_global, a,
            disp_str = factor("infty", levels = c("infty", "1_0")),
            fun_plot = a * (M^2) * (1 - (a/2) * Kappa))

fun_1_0 <- d_only_info$df %>%
  left_join(d_kappa_disp, by = "M") %>%             # Kappa from disp-limit run
  transmute(M, M_plot, is_global, a,
            disp_str = factor("1_0",  levels = c("infty", "1_0")),
            fun_plot = a * (M^2) * (1 - (a/2) * Kappa))

fun_two_info <- d_only_info
fun_two_info$df <- bind_rows(fun_infty, fun_1_0) %>%
  mutate(
    fun_plot = ifelse(M == 99, 25, fun_plot),       # Global override
    fun_plot = ifelse(fun_plot < 0, NA_real_, fun_plot)
  )

p_only_fun <- plot_with_global(
  fun_two_info,
  "fun_plot",
  # y-axis label kept as canonical form for continuity
  expression( bar(E[x] * "[" * J(x) * "]")/N),
  gap_start_frac = 0.47, gap_end_frac = 0.61,
  bridge_start   = 0.02,  bridge_end   = 0.98,
  group_var      = "disp_str",
  color_values   = c("infty" = "grey35", "1_0" = "#E11D48"),
  shape_values   = c("infty" = 16,       "1_0" = 17),
  legend_title   = "sigma[D]",
  labels_map     = c("infty" = "sigma[D]==infinity", "1_0" = "sigma[D]==1.0")
) +
  theme(
    legend.position      = c(0.98, 0.86),
    legend.justification = c("right", "top")
  )

# (C) HP richness (grouped by range)
p_hp_rich <- plot_with_global(
  d_hp_info, "species_richness", "Species richness",
  group_var = "range_str",
  gap_start_frac = 0.47, gap_end_frac = 0.61,
  bridge_start   = 0.02, bridge_end   = 0.98,
  legend_title = "Range"
)

# (D) HP covariance (grouped by range)
p_hp_cov <- plot_with_global(
  d_hp_info, "scaled_cov",
  expression(bar(Cov)[x]~bgroup("(", J(x)~","~~H(x)/mu[H(x)], ")")),
  group_var = "range_str",
  gap_start_frac = 0.47, gap_end_frac = 0.61,
  bridge_start   = 0.02, bridge_end   = 0.98,
  legend_title = "Range"
)

# ---------- Assemble 2×2 ----------
fig_6 <- (p_only_rich | p_only_fun) / (p_hp_rich | p_hp_cov)

# ---------- Save to top-level Figures ----------
out_svg <- file.path(figures_dir, "Fig6_2x2_grouped.svg")
out_png <- file.path(figures_dir, "Fig6_2x2_grouped.png")

ggsave(out_svg, fig_6, width = 11, height = 10, dpi = 300)
ggsave(out_png, fig_6, width = 11, height = 10, dpi = 300)

message("Saved: ", out_svg, " and ", out_png)
