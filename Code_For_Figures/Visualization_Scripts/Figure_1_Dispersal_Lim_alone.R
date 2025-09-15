# === Libraries ===
library(ggplot2)
library(dplyr)
library(grid)
library(ggpattern)  # for patterns if needed

# === Grid and tree positions ===
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
sigma <- 3.5  # controls the spread of dispersal
tree_x <- trees$x[1]
tree_y <- trees$y[1]

grid <- grid %>%
  mutate(
    dist2 = (X - tree_x)^2 + (Y - tree_y)^2,
    dispersal = exp(-dist2 / (2 * sigma^2)),
    dispersal_scaled = dispersal / max(dispersal)  # normalize to [0,1]
  )

# === Texture background ===
texture <- rasterGrob(
  matrix(runif(10000), nrow = 100, ncol = 100),
  width = unit(1, "npc"), height = unit(1, "npc"),
  interpolate = TRUE
)

# === Plot with dispersal overlay ===
p <- ggplot(grid, aes(X, Y)) +
  # Texture
  annotation_custom(
    grob = texture,
    xmin = 0.5, xmax = grid_size + 0.5,
    ymin = 0.5, ymax = grid_size + 0.5
  ) +
  
  # Base brown fill
  geom_tile(fill = "#6B3500", color = "grey12", alpha = 0.9, size = 0.75) +
  
  # Dispersal kernel overlay (blue)
  geom_tile(aes(alpha = dispersal_scaled), 
            fill = "orange", color = NA) +
  
  # Tree location
  geom_point(data = trees, aes(x = x, y = y),
             color = "black", size = 1.2, inherit.aes = FALSE) +
  geom_tile(fill = NA, color = "grey12", alpha = 0.9, size = .75) +
  
  # Border
  annotate("rect",
           xmin = min(grid$X) - 0.5, xmax = max(grid$X) + 0.5,
           ymin = min(grid$Y) - 0.5, ymax = max(grid$Y) + 0.5,
           color = "black", fill = NA, linewidth = 2) +
  
  # Layout
  coord_fixed() +
  theme_void() +
  theme(legend.position = "none") +
  scale_alpha(range = c(0, 0.8))

# === Show plot ===
p

# === Save output ===
ggsave("Dispersal_Kernel_Plot.svg", plot = p, width = 8, height = 8, dpi = 300)
