# ==========================================================
# Script: presence_L_2x4_centers.R
# Purpose: 2x4 panels = presence maps + Besag's L(r)-r
# Assumption: each individual is at the CENTER of its occupied grid cell.
# Geometry: Euclidean distances on a unit square grid (cell_size = 1).
# ==========================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
  library(spatstat)   # meta-package; provides Kest/Lest (or Kest + transform)
})

`%||%` <- function(a, b) if (is.null(a)) b else a

# ---------- INPUT DATA ----------
paths <- c(
  "C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims/Outputs/Fig_5/SimSet1_JCstrength/comm_JC_100.0.csv",
  "C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims/Outputs/Fig_5/SimSet2_Dispersal/comm_Disp_2.75.csv",
  "C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims/Outputs/Fig_5/SimSet4_HP_plus_JCE/comm_Rough_1000_0.csv",
  "C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/HP_JCE_Sims/Outputs/Fig_5/SimSet4_HP_plus_JCE/comm_Rough_0_015848931924611134.csv"
)

labels <- c("Only JC (a = 100)",
            "Dispersal (σ[D] = 2.75)",
            "HP+JC (rough = 1000)",
            "HP+JC (rough ≈ 0.016)")

# ---------- THEME ----------
theme_fig5_child <- function() {
  theme_classic(base_size = 13) +
    theme(
      panel.border   = element_rect(color = "black", fill = NA, linewidth = 2.0),
      axis.line      = element_line(color = "black", linewidth = 1.2),
      axis.ticks     = element_line(color = "black", linewidth = 1.2),
      axis.text      = element_text(size = 12, color = "black"),
      plot.title     = element_text(size = 13, face = "bold", hjust = 0),
      panel.grid     = element_blank(),
      legend.position= "none",
      plot.margin    = margin(6,6,6,6)
    )
}

# ---------- HELPERS ----------
# Read grid; drop first row (labels) and optional first column (row labels)
load_species_grid <- function(path) {
  raw <- read.csv(path, header = FALSE, stringsAsFactors = FALSE,
                  check.names = FALSE, na.strings = c("", "NA"))
  stopifnot(nrow(raw) >= 2, ncol(raw) >= 2)
  raw <- raw[-1, , drop = FALSE]  # drop CSV header row
  
  # If first column is non-numeric (row labels), drop it
  first_col_num <- suppressWarnings(as.numeric(gsub(",", "", raw[[1]])))
  if (any(is.na(first_col_num))) raw <- raw[, -1, drop = FALSE]
  
  X <- as.matrix(apply(raw, 2, function(x) as.numeric(gsub(",", "", x))))
  storage.mode(X) <- "integer"
  X
}

# Convert grid to ppp of cell CENTERS; marks = species IDs
# Center of cell (row i, col j) -> (x = j-0.5, y = i-0.5) * cell_size
grid_centers_ppp <- function(X, cell_size = 1) {
  n <- nrow(X); stopifnot(ncol(X) == n)
  ij <- which(!is.na(X), arr.ind = TRUE)
  if (nrow(ij) == 0) stop("No occupied cells found.")
  x <- (ij[,2] - 0.5) * cell_size
  y <- (ij[,1] - 0.5) * cell_size
  W <- owin(c(0, n * cell_size), c(0, n * cell_size))
  marks <- factor(X[ij], levels = sort(unique(na.omit(as.integer(X)))))
  ppp(x, y, window = W, marks = marks)
}

# Species-specific Besag's L (reported as L(r) - r) at integer cell widths (>= 1)
# We compute Kest with translation correction, transform to L = sqrt(K/pi),
# then plot L - r. r must start at 0 for spatstat; we drop r < 1 cell afterwards.
L_species_centers <- function(X, species_id, rmax_cells = 40, cell_size = 1) {
  n <- nrow(X)
  rmax <- min(rmax_cells * cell_size, (n/2 - 1) * cell_size)
  r_eval_full <- c(0, seq(from = 1 * cell_size, to = rmax, by = 1 * cell_size))
  
  pp <- grid_centers_ppp(X, cell_size = cell_size)
  if (!(as.character(species_id) %in% levels(pp$marks)))
    stop("Species not present in this grid.")
  sub <- subset(pp, marks == as.character(species_id))
  if (sub$n < 2) stop("Too few individuals for this species to compute L(r).")
  
  kest <- Kest(sub, r = r_eval_full, correction = "translation")
  # transform: L(r) = sqrt(K(r)/pi); use the translation-corrected column
  Ltrans <- sqrt(kest$trans / pi)
  
  # drop r < 1 cell
  keep <- kest$r >= 1 * cell_size
  data.frame(r = kest$r[keep],
             L = Ltrans[keep],
             Lminusr = (Ltrans - kest$r)[keep])
}

# Build presence + L(r)-r panels for one dataset
make_panels <- function(csv_path, mean_abund = 300, tol = 0.10,
                        rmax_cells = 40, cell_size = 1,
                        map_lab = NULL, line_lab = NULL) {
  X <- load_species_grid(csv_path)
  n <- nrow(X)
  
  # --- choose focal species: closest to target mean (with optional ±tol) ---
  species_ids    <- sort(unique(as.integer(X)))
  species_counts <- as.integer(table(factor(as.integer(X), levels = species_ids)))
  diffs      <- abs(species_counts - mean_abund)
  candidates <- which(diffs == min(diffs))
  #within_tol <- which(abs(species_counts - mean_abund) <= tol * mean_abund)
  #if (length(within_tol) > 0) candidates <- within_tol
  set.seed(123)
  chosen_idx     <- candidates[sample(1:length(candidates), 1)]
  chosen_species <- species_ids[chosen_idx]
  
  # --- presence map (origin lower-left) ---
  pres_mat <- (X == chosen_species) * 1L
  df_pres  <- expand.grid(x = 1:n, y = 1:n)
  df_pres$y    <- n + 1 - df_pres$y
  df_pres$pres <- factor(as.vector(pres_mat), levels = c(0,1))
  p_map <- ggplot(df_pres, aes(x, y, fill = pres)) +
    geom_raster() +
    coord_equal(expand = FALSE) +
    scale_fill_manual(values = c("0" = "#f0f0f0", "1" = "#2c7fb8")) +
    labs(x = NULL, y = NULL,
         title = map_lab %||% sprintf("Presence (sp %s)", chosen_species)) +
    theme_fig5_child()
  
  # --- Besag's L(r)-r from centers, at integer cell widths ---
  df_L <- L_species_centers(X, chosen_species, rmax_cells = rmax_cells, cell_size = cell_size)
  
  p_L <- ggplot(df_L, aes(r / cell_size, Lminusr / cell_size)) +  # in cell-width units
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_point(size = 2.2) +
    geom_line(linewidth = 0.8) +
    labs(x = "Distance (cell widths)",
         y = expression(L(r) - r),
         title = line_lab %||% sprintf("Besag's L(r) - r (sp %s)", chosen_species)) +
    theme_fig5_child()
  
  list(p_map = p_map, p_gr = p_L)  # keep list names consistent with caller
}

# ---------- BUILD ALL PANELS ----------
panels <- Map(function(pth, lab) {
  make_panels(
    csv_path   = pth,
    mean_abund = 375, tol = 0.10,         # your current target
    rmax_cells = 40, cell_size = 1,
    map_lab  = paste0(lab, "\nPresence"),
    line_lab = paste0(lab, "\nBesag's L(r) - r")
  )
}, paths, labels)

# ---------- ASSEMBLE 2×4 ----------
row1 <- panels[[1]]$p_map | panels[[2]]$p_map | panels[[3]]$p_map | panels[[4]]$p_map
row2 <- panels[[1]]$p_gr  | panels[[2]]$p_gr  | panels[[3]]$p_gr  | panels[[4]]$p_gr
fig_2x4 <- row1 / row2 + plot_layout(heights = c(1, 0.6))

# ---------- VIEW / SAVE ----------
fig_2x4
# ggsave("presence_L_2x4_centers.png", fig_2x4, width = 16, height = 8, dpi = 300)
