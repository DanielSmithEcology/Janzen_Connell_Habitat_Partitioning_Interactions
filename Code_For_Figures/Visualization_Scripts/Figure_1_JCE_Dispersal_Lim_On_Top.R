# === Libraries ===
library(ggplot2)
library(dplyr)
library(grid)
library(ggpattern)

# === Grid and tree position ===
grid_size <- 11
X_vals <- 1:grid_size
Y_vals <- 1:grid_size
grid <- expand.grid(X = X_vals, Y = Y_vals)

# === Tree location ===
trees <- data.frame(
  x = 6,
  y = 6
)

# === Dispersal kernel (Gaussian) ===
sigma <- 3.5
tree_x <- trees$x[1]
tree_y <- trees$y[1]

grid <- grid %>%
  mutate(
    dist2 = (X - tree_x)^2 + (Y - tree_y)^2,
    dispersal = exp(-dist2 / (2 * sigma^2)),
    dispersal_scaled = dispersal / max(dispersal)  # normalize to [0,1]
  )

# === Moore neighborhood (M = 5) ===
M <- 5
radius <- floor(M / 2)

grid$in_neigh <- apply(grid, 1, function(row) {
  x <- row["X"]
  y <- row["Y"]
  abs(tree_x - x) <= radius & abs(tree_y - y) <= radius
})

moore_neigh <- grid %>% filter(in_neigh)

# === Texture background ===
texture <- rasterGrob(
  matrix(runif(10000), nrow = 100, ncol = 100),
  width = unit(1, "npc"), height = unit(1, "npc"),
  interpolate = TRUE
)

# === Combined plot ===
p <- ggplot(grid, aes(X, Y)) +
  # Texture
  annotation_custom(texture,
                    xmin = 0.5, xmax = grid_size + 0.5,
                    ymin = 0.5, ymax = grid_size + 0.5) +
  
  # Base brown fill
  geom_tile(fill = "#6B3500", color = "grey12", alpha = 0.9, size = 1.1) +
  
  # Dispersal kernel (orange gradient)
  geom_tile(aes(alpha = dispersal_scaled),
            fill = "orange", color = NA) +
  
  # Pattern stripes in Moore neighborhood
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
  
  geom_tile(fill = NA, color = "grey12", alpha = 0.9, size = .75) +
  
  
  # Moore neighborhood outline
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
             color = "black", size = .75, inherit.aes = FALSE) +
  
  # Outer border
  annotate("rect",
           xmin = min(grid$X) - 0.5, xmax = max(grid$X) + 0.5,
           ymin = min(grid$Y) - 0.5, ymax = max(grid$Y) + 0.5,
           color = "black", fill = NA, linewidth = 2) +
  
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none") +
  scale_alpha(range = c(0, 0.8))

# === Save or show ===
p
ggsave("Dispersal_JC_Overlay.svg", plot = p, width = 8, height = 8, dpi = 300)
