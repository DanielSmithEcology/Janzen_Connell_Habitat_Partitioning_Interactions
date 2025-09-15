# === Libraries ===
library(ggplot2)
library(grid)
library(png)
library(RandomFields)

# === Habitat Grid Generator with Optional Gaussian Clump ===
generate_habitat_grid_with_center_clump <- function(grid_size = 13,
                                                    scale = 2,
                                                    sigma_mean = 2,
                                                    seed = NULL,
                                                    transform_uniform = TRUE) {
  # Create grid
  grid_data <- expand.grid(x = 1:grid_size, y = 1:grid_size)
  x_coords <- unique(grid_data$x)
  y_coords <- unique(grid_data$y)
  
  # Set or reset seed
  if (!is.null(seed)) {
    RFoptions(seed = seed, spConform = TRUE)
  } else if (exists(".Random.seed", envir = globalenv())) {
    rm(.Random.seed, envir = globalenv())
  }
  
  # Gaussian bump centered in the grid
  cx <- mean(x_coords)
  cy <- mean(y_coords)
  bump <- outer(
    1:grid_size, 1:grid_size,
    function(i, j) exp(-((i - cx)^2 + (j - cy)^2) / (2 * sigma_mean^2))
  )
  mean_field <- as.vector(t(bump))  # transpose to match expand.grid order
  
  # Simulate spatial autocorrelation with Gaussian RF
  model <- RMgauss(var = 1, scale = scale)
  sim <- RFsimulate(model, x = x_coords, y = y_coords, grid = TRUE)
  
  # Add spatial autocorrelation and central bump
  grid_data$habitat <- as.vector(sim@data[[1]]) + mean_field
  
  # Normalize to [0, 1]
  if (transform_uniform) {
    grid_data$habitat <- rank(grid_data$habitat, ties.method = "average") / nrow(grid_data)
  } else {
    grid_data$habitat <- (grid_data$habitat - min(grid_data$habitat)) / 
      (max(grid_data$habitat) - min(grid_data$habitat))
  }
  
  return(grid_data)
}

# === Center Highlight Rectangle Function ===
add_center_rectangle <- function(center = c(7, 7), size = 5, 
                                 fill = NA, color = "#FF4756", thickness = 2) {
  half <- floor(size / 2)
  xmin <- center[1] - half + .5 - 1 - .2
  xmax <- center[1] + half + .5 + .2
  ymin <- center[2] - half + .5 - 1 - .2
  ymax <- center[2] + half + .5 + .2
  
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = fill,
    color = color,
    size = thickness
  )
}

# === Set Seed and Generate Grid ===
SEEDSET <- 8641
grid_data <- generate_habitat_grid_with_center_clump(
  grid_size = 13,
  scale = 10,
  sigma_mean = 2,
  seed = SEEDSET, 
  transform_uniform = TRUE
)

# === Reverse habitat for plotting aesthetic ===
grid_data$habitat <- 1 - grid_data$habitat

# === Background Texture Layer ===
texture <- rasterGrob(
  matrix(runif(10000), nrow = 200, ncol = 200), 
  width = unit(1, "npc"), height = unit(1, "npc"), 
  interpolate = TRUE
)

# === Create Plot ===
p <- ggplot(grid_data, aes(x, y)) +
  annotation_custom(texture, xmin = 0.5, xmax = 13.5, ymin = 0.5, ymax = 13.5) +
  geom_tile(aes(fill = habitat), color = "grey12", size = 1.1, alpha = 0.8) +
  geom_rect(
    aes(xmin = 0.35, xmax = 13.65, ymin = 0.35, ymax = 13.65),
    inherit.aes = FALSE,
    color = "black",
    fill = NA,
    linewidth = 7
  ) +
  scale_fill_viridis_c(option = "D", name = "Habitat") +
  coord_fixed(expand = FALSE) +
  add_center_rectangle(center = c(7, 7), size = 5, color = "darkred", thickness = 4) +
  add_center_rectangle(center = c(7, 7), size = 5, color = "#FF4756", thickness = 2) +
  theme_void() +
  theme(
    plot.margin = margin(5, 5, 5, 5, "mm"),
    legend.position = "none"
  )

# === Save and Draw Plot ===
ggsave("Test_Agg.svg", p, width = 8, height = 8, dpi = 400)
grid.newpage()
grid.draw(p)
