#!/usr/bin/env Rscript
# ==========================================================
# Script: pSI_Figure_Heatmaps_Tol_Fec.R
# Purpose: 2 panels (Low vs High autocorr), No dispersal limitation
# Inputs (relative to repo root):
#   - HP_JCE_Sims/Outputs/Appendix_SST/Low_Auto_No_Disp_Lim/Richness_M_Z_No_Disp_Lim_Low_Autocorr.csv
#   - HP_JCE_Sims/Outputs/Appendix_SST/High_Auto_No_Disp_Lim/Richness_M_Z_No_Disp_Lim_High_Autocorr.csv
# Outputs:
#   - Figures/Fig_SST_2x1_RichnessHeatmaps.svg
#   - Figures/Fig_SST_2x1_RichnessHeatmaps.png
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(viridis)
  library(patchwork)
  library(grid)
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
appx_base   <- file.path(repo_dir, "HP_JCE_Sims", "Outputs", "Appendix_SST")
figures_dir <- file.path(repo_dir, "Figures")
if (!dir.exists(figures_dir)) dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

# ---- File paths (relative) ----
files <- c(
  Low  = file.path(appx_base, "Low_Auto_No_Disp_Lim",  "Richness_M_Z_No_Disp_Lim_Low_Autocorr.csv"),
  High = file.path(appx_base, "High_Auto_No_Disp_Lim", "Richness_M_Z_No_Disp_Lim_High_Autocorr.csv")
)

# ---- Reader (robust to missing/odd headers) ----
read_sim_csv <- function(path) {
  df_try <- try(read.csv(path, header = TRUE, stringsAsFactors = FALSE), silent = TRUE)
  df <- if (inherits(df_try, "try-error")) read.csv(path, header = FALSE, stringsAsFactors = FALSE) else df_try
  
  nm <- tolower(names(df))
  richness_col <- which(nm %in% c("species_richness","richness","sr"))
  neigh_col    <- which(nm %in% c("neigh_str","m","neigh","moore","moore_neigh","moore_n"))
  z_col        <- which(nm %in% c("z_str","z","sigma_str","sigmah","sigma_h"))
  
  if (length(richness_col) == 0 && ncol(df) >= 1) richness_col <- 1
  if (length(z_col)        == 0 && ncol(df) >= 2) z_col        <- 2
  if (length(neigh_col)    == 0 && ncol(df) >= 3) neigh_col    <- 3
  
  out <- data.frame(
    species_richness = as.numeric(df[[richness_col]]),
    Z_str            = as.character(df[[z_col]]),
    neigh_str        = as.character(df[[neigh_col]]),
    stringsAsFactors = FALSE
  )
  names(out) <- c("species_richness","Z_str","neigh_str")
  out
}

# ---- Cleaning & ordering ----
clean_labels <- function(df) {
  df <- df %>%
    mutate(
      neigh_str = as.character(neigh_str),
      Z_str     = as.character(Z_str),
      
      neigh_str = dplyr::case_when(
        tolower(neigh_str) %in% c("none","no_jces","no jces") ~ "No JCEs",
        neigh_str == "0" ~ "1",
        neigh_str == "1" ~ "3",
        neigh_str == "2" ~ "5",
        neigh_str == "3" ~ "7",
        TRUE ~ neigh_str
      ),
      
      Z_label = dplyr::case_when(
        tolower(Z_str) %in% c("inf","infinity","infinite") ~ "No HP",
        Z_str %in% c("9_999999e6","9.999999e6","9999999","1e9") ~ "No HP",
        TRUE ~ gsub("_","\\.", Z_str)
      )
    )
  
  neigh_levels <- c("No JCEs","1","3","5","7")
  df$neigh_str <- factor(df$neigh_str, levels = neigh_levels[neigh_levels %in% df$neigh_str])
  
  uniqZ <- unique(df$Z_label)
  has_nohp <- "No HP" %in% uniqZ
  z_numeric <- suppressWarnings(as.numeric(uniqZ))
  if (all(is.na(z_numeric))) {
    lev <- setdiff(sort(uniqZ), "No HP")
    if (has_nohp) lev <- c(lev, "No HP")
  } else {
    z_vals <- suppressWarnings(as.numeric(gsub("[^0-9.]+","", uniqZ)))
    tmp <- data.frame(lbl = uniqZ, val = z_vals)
    tmp <- tmp[order(tmp$val), , drop = FALSE]
    lev <- tmp$lbl[!is.na(tmp$val)]
    if (has_nohp && !("No HP" %in% lev)) lev <- c(lev, "No HP")
  }
  df$Z_label <- factor(df$Z_label, levels = lev)
  
  df
}

# ---- Single heatmap panel ----
plot_heatmap <- function(df) {
  ggplot(df, aes(x = neigh_str, y = Z_label, fill = species_richness)) +
    geom_tile(color = "black", linewidth = 0.6) +
    geom_text(aes(label = round(species_richness, 0)),
              color = "grey70", fontface = "bold", size = 3.8) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0)) +
    labs(
      x = "JC-effect Moore neighborhood (M)",
      y = expression(italic(Z))
    ) +
    theme(
      panel.background = element_rect(fill = "white"),
      panel.border     = element_rect(colour = "black", fill = NA, size = 2),
      plot.background  = element_rect(fill = "white"),
      axis.text.x      = element_text(color = "black", size = 14),
      axis.text.y      = element_text(color = "black", size = 14),
      axis.title.x     = element_text(size = 18, margin = margin(t = 8)),
      axis.title.y     = element_text(size = 18),
      legend.text      = element_text(size = 18),
      legend.key.height = unit(2.2, "cm"),
      legend.key.width  = unit(.35, "cm")
    )
}

# ---- Load & prepare ----
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

p_low  <- plot_heatmap(dfs$Low)  + fill_scale + labs(title = "Low autocorrelation")
p_high <- plot_heatmap(dfs$High) + fill_scale + labs(title = "High autocorrelation")

title_theme <- theme(
  plot.title = element_text(hjust = 0.5, size = 16, face = "bold", margin = margin(b = 4))
)

p_low  <- p_low  + title_theme
p_high <- p_high + title_theme

# ---- Assemble 1×2, shared legend on right ----
final_plot <-
  (p_low + p_high) +
  plot_layout(guides = "collect", widths = c(1, 1)) &
  theme(
    legend.position = "right",
    legend.title = element_text(angle = -90, vjust = 0.5, hjust = 0.5, size = 18)
  )

final_plot <- final_plot & guides(
  fill = guide_colorbar(
    title.position = "left",
    title.vjust    = 1,
    label.position = "right",
    direction      = "vertical",
    ticks.colour   = NA
  )
)

# View
print(final_plot)

# ---- Save to top-level Figures ----
svg_out <- file.path(figures_dir, "Fig_SST_2x1_RichnessHeatmaps.svg")
png_out <- file.path(figures_dir, "Fig_SST_2x1_RichnessHeatmaps.png")
ggsave(svg_out, plot = final_plot, width = 10, height = 5.0, units = "in", dpi = 300)
ggsave(png_out, plot = final_plot, width = 10, height = 5.0, units = "in", dpi = 300)

message("Saved: ", svg_out, " and ", png_out)
