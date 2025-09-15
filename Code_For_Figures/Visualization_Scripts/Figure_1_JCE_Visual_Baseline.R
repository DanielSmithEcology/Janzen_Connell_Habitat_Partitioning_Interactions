library(ggplot2)
library(dplyr)
library(grid)
library(ggpattern)  # <- required for pattern fill

# Grid and tree positions
grid_size <- 9
X_vals <- 1:grid_size
Y_vals <- 1:grid_size
grid <- expand.grid(X = X_vals, Y = Y_vals)

# Tree locations
trees <- data.frame(
  x = c(5),
  y = c(5) 
)

# Moore neighborhood definition (M = 5)
M <- 5
radius <- floor(M / 2)

# Identify Moore neighborhood tiles
grid$in_neigh <- apply(grid, 1, function(row) {
  x <- row["X"]
  y <- row["Y"]
  any(abs(trees$x - x) <= radius & abs(trees$y - y) <= radius)
})

# Subset grid into affected tiles (for stripes)
moore_neigh <- grid %>% filter(in_neigh)

# === Texture background ===
texture <- rasterGrob(
  matrix(runif(10000), nrow = 100, ncol = 100),
  width = unit(1, "npc"), height = unit(1, "npc"),
  interpolate = TRUE
)

# === Plot ===
# === Plot ===
p <- ggplot(grid, aes(X, Y)) +
  # Texture
  annotation_custom(
    grob = texture,
    xmin = 0.5, xmax = grid_size + 0.5,
    ymin = 0.5, ymax = grid_size + 0.5
  ) +
  
  # Base fill (brown)
  geom_tile(fill = "#6B3500", color = "grey12", alpha = 0.9, size = 1.1) +
  
  # Light red background in Moore neighborhood
  #geom_tile(data = moore_neigh,
  #          aes(x = X, y = Y),
  #          fill = "darkred", alpha = 0.06, size=4
  #          color = NA, inherit.aes = FALSE) +

  #geom_tile(data = moore_neigh,
  #          aes(x = X, y = Y),
  #          fill = "#FF4756", alpha = 0.06
  #          color = NA, inherit.aes = FALSE) +
  
  
    
  # Pattern stripes
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
  
  # Outline Moore neighborhood
  annotate("rect",
           xmin = min(moore_neigh$X) - 0.5,
           xmax = max(moore_neigh$X) + 0.5,
           ymin = min(moore_neigh$Y) - 0.5,
           ymax = max(moore_neigh$Y) + 0.5,
           color = "darkred", fill = NA, linewidth = 4) +
  
  annotate("rect",
           xmin = min(moore_neigh$X) - 0.5,
           xmax = max(moore_neigh$X) + 0.5,
           ymin = min(moore_neigh$Y) - 0.5,
           ymax = max(moore_neigh$Y) + 0.5,
           color = "#FF4756", fill = NA, linewidth = 2) +
  
  
  # Tree
  geom_point(data = trees, aes(x = x, y = y),
             color = "black", size = 1.2, inherit.aes = FALSE) +
  
  # Final layout
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none") +
  annotate("rect",
           xmin = min(grid$X) - 0.5, xmax = max(grid$X) + 0.5,
           ymin = min(grid$Y) - 0.5, ymax = max(grid$Y) + 0.5,
           color = "black", fill = NA, linewidth = 2)

p
# Save output
ggsave("JC_effect_pattern.svg", plot = p, width = 8, height = 8, dpi = 300)
