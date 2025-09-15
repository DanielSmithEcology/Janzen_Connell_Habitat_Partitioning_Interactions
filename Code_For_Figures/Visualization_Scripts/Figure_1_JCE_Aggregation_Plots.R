library(ggplot2)
library(dplyr)
library(grid)
library(ggpattern)

# === Parameters ===
grid_size <- 19
M <- 5
radius <- floor(M / 2)

# === Texture background ===
texture <- rasterGrob(
  matrix(runif(10000), nrow = 100, ncol = 100),
  width = unit(1, "npc"), height = unit(1, "npc"),
  interpolate = TRUE
)

# === Helper to compute Moore neighborhoods ===
compute_moore_neigh <- function(grid, trees, radius) {
  grid %>%
    rowwise() %>%
    mutate(
      in_neigh = any(abs(trees$x - X) <= radius & abs(trees$y - Y) <= radius)
    ) %>%
    ungroup()
}

# === Function to generate annotations for all Moore boxes ===
moore_box_annotations <- function(trees, radius) {
  boxes <- lapply(1:nrow(trees), function(i) {
    center_x <- trees$x[i]
    center_y <- trees$y[i]
    list(
      annotate("rect",
               xmin = center_x - radius - 0.5,
               xmax = center_x + radius + 0.5,
               ymin = center_y - radius - 0.5,
               ymax = center_y + radius + 0.5,
               color = "darkred", fill = NA, linewidth = 4),
      annotate("rect",
               xmin = center_x - radius - 0.5,
               xmax = center_x + radius + 0.5,
               ymin = center_y - radius - 0.5,
               ymax = center_y + radius + 0.5,
               color = "#FF4756", fill = NA, linewidth = 2)
    )
  })
  do.call(list, unlist(boxes, recursive = FALSE))
}

# === Plot Function ===
make_plot <- function(grid, trees, filename) {
  grid_marked <- compute_moore_neigh(grid, trees, radius)
  moore_neigh <- grid_marked %>% filter(in_neigh)
  
  p <- ggplot(grid_marked, aes(X, Y)) +
    annotation_custom(
      grob = texture,
      xmin = 0.5, xmax = grid_size + 0.5,
      ymin = 0.5, ymax = grid_size + 0.5
    ) +
    geom_tile(fill = "#6B3500", color = "grey12", alpha = 0.9, size = .6) +
    geom_tile_pattern(data = moore_neigh,
                      aes(x = X, y = Y),
                      pattern = "pch",
                      pattern_shape = 3,
                      pattern_angle = 45,
                      pattern_colour = "red3",
                      pattern_fill = "red3",
                      pattern_density = 0.35,
                      pattern_spacing = 0.075,
                      pattern_size = 2,
                      fill = NA,
                      color = NA, size = 1.1, alpha = 0.0) +
    moore_box_annotations(trees, radius) +
    geom_point(data = trees, aes(x = x, y = y),
               color = "black", size = 1.2, inherit.aes = FALSE) +
    annotate("rect",
             xmin = min(grid$X) - 0.5, xmax = max(grid$X) + 0.5,
             ymin = min(grid$Y) - 0.5, ymax = max(grid$Y) + 0.5,
             color = "black", fill = NA, linewidth = 2) +
    coord_fixed() +
    theme_void() +
    theme(legend.position = "none")
  
  ggsave(filename, plot = p, width = 8, height = 8, dpi = 300)
}

# === Grid ===
X_vals <- 1:grid_size
Y_vals <- 1:grid_size
grid <- expand.grid(X = X_vals, Y = Y_vals)

# === Plot 1: Evenly spaced trees ===
even_trees <- data.frame(
  x = c(4, 5, 10, 16, 15),
  y = c(6, 14, 7, 4, 16)
)
make_plot(grid, even_trees, "Moore_Even_Spacing.svg")

# === Plot 2: Clustered trees near center ===
clustered_trees <- data.frame(
  x = c(8, 10, 11, 11, 13),
  y = c(10, 9, 11, 8, 10)
)
make_plot(grid, clustered_trees, "Moore_Clustered_Spacing.svg")
