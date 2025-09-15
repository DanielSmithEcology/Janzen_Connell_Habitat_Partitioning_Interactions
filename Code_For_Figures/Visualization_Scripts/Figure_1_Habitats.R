# === Libraries ===
library(ggplot2)
library(viridis)
library(dplyr)
library(tidyr)
library(readr)
library(purrr)
library(stringr)
library(grid)
library(ggpattern)

# === Output directory ===
dir.create("PNG_Exports", showWarnings = FALSE)

# === Parameters ===
L <- 9  # size of the square subsample
all_targets <- paste0("Habb", 1:9)
sub_targets <- paste0("Habb", 7:9)

# === Helper function: extract right-center L × L region ===
extract_right_center <- function(gdat, L) {
  max_x <- max(gdat$x) - 30
  mid_y <- median(unique(gdat$y)) +5
  
  x_vals <- (max_x - L + 1):max_x
  y_vals <- round(mid_y - L/2 + 1):(round(mid_y + L/2))
  
  gdat %>%
    dplyr::filter(x %in% x_vals, y %in% y_vals)
}

# === Texture background for subplots ===
texture <- rasterGrob(matrix(runif(10000), nrow = 200, ncol = 200), 
                      width = unit(1, "npc"), height = unit(1, "npc"), 
                      interpolate = TRUE)

# === Plot full grid with black box around subsample ===
make_full_plot_with_box <- function(full_gdat, sub_gdat) {
  #full_gdat$habitat <- 1 - full_gdat$habitat
  sub_box <- data.frame(
    xmin = min(sub_gdat$x) - 0.5,
    xmax = max(sub_gdat$x) + 0.5,
    ymin = min(sub_gdat$y) - 0.5,
    ymax = max(sub_gdat$y) + 0.5
  )
  
  ggplot(full_gdat, aes(x, y)) +
#    geom_tile(aes(fill = habitat), color = NA) +
    geom_tile(aes(fill = habitat), color = "grey12", size = 1/174, alpha = 1.0) +
    geom_rect(data = sub_box,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE,
              color = "red3", linewidth = 3, fill = NA) +
    
    geom_rect(aes(xmin = min(full_gdat$x) - 0.65,
                  xmax = max(full_gdat$x) + 0.65,
                  ymin = min(full_gdat$y) - 0.65,
                  ymax = max(full_gdat$y) + 0.65),
              inherit.aes = FALSE,
              color = "black", fill = NA, linewidth = 4.5) +
    
    scale_fill_viridis_c(option = "D") +
    coord_fixed(expand = FALSE) +
    theme_void() +
    theme(legend.position = "none")
}
make_textured_plot <- function(gdat, m = 5) {
  require(ggpattern)
  
  # Center and Moore neighborhood
  center_x <- median(unique(gdat$x))
  center_y <- median(unique(gdat$y))
  
  moore_neigh <- gdat %>%
    filter(abs(x - center_x) <= floor(m / 2),
           abs(y - center_y) <= floor(m / 2))
  
  ggplot() +
    # Texture background (goes under everything)
    annotation_custom(texture,
                      xmin = min(gdat$x) - 0.5,
                      xmax = max(gdat$x) + 0.5,
                      ymin = min(gdat$y) - 0.5,
                      ymax = max(gdat$y) + 0.5) +
    
    # Base fill (for entire grid)
    geom_tile(data = gdat,
              aes(x = x, y = y, fill = habitat),
              color = "grey12", size = 1.1, alpha = 0.8) +
    
    # Pattern stripes
    geom_tile_pattern(data = moore_neigh,
                      aes(x = x, y = y),
                      pattern = "pch",
                      pattern_shape = 3,
                      pattern_angle = 45,
                      pattern_colour = "red3",
                      pattern_fill = "red3",
                      pattern_density = 0.35,
                      pattern_spacing = 0.075,
                      pattern_size = 2.0,
                      fill = NA,
                      color = NA,
                      size = 1.1,
                      alpha = 0.0) +
    
    # Outline Moore neighborhood
    annotate("rect",
             xmin = min(moore_neigh$x) - 0.5,
             xmax = max(moore_neigh$x) + 0.5,
             ymin = min(moore_neigh$y) - 0.5,
             ymax = max(moore_neigh$y) + 0.5,
             color = "darkred", fill = NA, linewidth = 4) +
    
    annotate("rect",
             xmin = min(moore_neigh$x) - 0.5,
             xmax = max(moore_neigh$x) + 0.5,
             ymin = min(moore_neigh$y) - 0.5,
             ymax = max(moore_neigh$y) + 0.5,
             color = "#FF4756", fill = NA, linewidth = 2) +
    
    # Plot border
    geom_rect(aes(xmin = min(gdat$x) - 0.5,
                  xmax = max(gdat$x) + 0.5,
                  ymin = min(gdat$y) - 0.5,
                  ymax = max(gdat$y) + 0.5),
              inherit.aes = FALSE,
              color = "black", fill = NA, linewidth = 4.5) +
    
    scale_fill_viridis_c(option = "D", limits = c(0, 1)) +
    coord_fixed(expand = FALSE) +
    theme_void() +
    theme(
      plot.margin = margin(5, 5, 5, 5, "mm"),
      legend.position = "none"
    )
}





# === Load and process ALL Habb1–Habb9 ===
full_grids_all <- map(all_targets, ~ {
  df <- read_csv(file.path("C:/Users/smith/OneDrive/Desktop/CNDD New Ideas/Hab_Par_New/Code_For_Figures", paste0(.x, ".csv")),
                 col_names = FALSE, show_col_types = FALSE)
  colnames(df) <- paste0("V", seq_len(ncol(df)))
  
  df <- df %>%
    dplyr::mutate(y = nrow(.) - dplyr::row_number() + 1) %>%
    pivot_longer(cols = starts_with("V"), names_to = "x", values_to = "habitat") %>%
    dplyr::mutate(
      x = as.integer(str_extract(x, "\\d+")),
      habitat = as.numeric(habitat)
    )
  
  df
}) %>% set_names(all_targets)

# === Extract subset from Habb7–Habb9 ===
full_grids <- full_grids_all[sub_targets]
sub_grids <- map(full_grids, extract_right_center, L = L)

# === Save full plots with box ===
walk(names(full_grids), function(name) {
  p <- make_full_plot_with_box(full_grids[[name]], sub_grids[[name]])
  ggsave(file.path("PNG_Exports", paste0(name, "_full_with_box.png")),
         plot = p, width = 8, height = 8, dpi = 600)
})

# === Save zoomed-in textured plots ===
walk2(sub_grids, names(sub_grids), function(sub_data, name) {
  p <- make_textured_plot(sub_data)
  ggsave(file.path("PNG_Exports", paste0(name, "_subsample.png")),
         plot = p, width = 8, height = 8, dpi = 600)
})


# === Save full plots for Habb1–Habb6 (no subsample box needed) ===
# === Save full plots for Habb1–Habb6 (no subsample box needed) ===
#walk(names(full_grids_all)[1:6], function(name) {
#  gdat <- full_grids_all[[name]]
  # gdat$habitat <- 1 - gdat$habitat
  
#  p <- ggplot(gdat, aes(x, y)) +
#    geom_tile(aes(fill = habitat), color = "grey12", size = 1/174, alpha = 1.0) +

#    geom_rect(aes(xmin = min(gdat$x) - 0.65,
#                  xmax = max(gdat$x) + 0.65,
#                  ymin = min(gdat$y) - 0.65,
#                  ymax = max(gdat$y) + 0.65),
#              inherit.aes = FALSE,
#              color = "black", fill = NA, linewidth = 4.5)+ 
#    scale_fill_viridis_c(option = "D") +
#    coord_fixed() +
#    theme_void() +
#    theme(legend.position = "none") 
  
  # Save the plot (outside the ggplot chain!)
#  ggsave(file.path("PNG_Exports", paste0(name, "_full.png")),
#         plot = p, width = 8, height = 8, dpi = 600)
#})




walk(names(full_grids_all)[1:6], function(name) {
  gdat <- full_grids_all[[name]]
  sub_gdat <- gdat[1:L^2, ]  # dummy to reuse function
  p <- make_full_plot_with_box(gdat, sub_gdat[0, ])  # pass empty subset to skip red box
  ggsave(file.path("PNG_Exports", paste0(name, "_full.png")),
         plot = p, width = 8, height = 8, dpi = 600)
})


