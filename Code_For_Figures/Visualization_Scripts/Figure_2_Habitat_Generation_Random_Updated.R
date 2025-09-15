# === Libraries ===
library(ggplot2)
library(grid)
library(png)

# === Function: Basic Habitat Grid Generator ===
generate_habitat_grid <- function(grid_size = 13, scale = 2, seed = 123, transform_uniform = TRUE) {
  if (!requireNamespace("RandomFields", quietly = TRUE)) {
    stop("Please install the 'RandomFields' package.")
  }
  library(RandomFields)
  
  grid_data <- expand.grid(x = 1:grid_size, y = 1:grid_size)
  RFoptions(seed = seed)
  model <- RMgauss(var = 1, scale = scale)
  sim <- RFsimulate(model, x = 1:grid_size, y = 1:grid_size)
  grid_data$habitat <- as.vector(sim$variable1)
  
  if (transform_uniform) {
    grid_data$habitat <- rank(grid_data$habitat, ties.method = "average") / nrow(grid_data)
  } else {
    grid_data$habitat <- (grid_data$habitat - min(grid_data$habitat)) / 
      (max(grid_data$habitat) - min(grid_data$habitat))
  }
  
  return(grid_data)
}

# === Function: Habitat Grid with Central Gaussian Clump ===
generate_habitat_grid_with_center_clump <- function(grid_size = 13,
                                                    scale = 2,
                                                    sigma_mean = 2,
                                                    seed = NULL,
                                                    transform_uniform = TRUE) {
  if (!requireNamespace("RandomFields", quietly = TRUE)) {
    stop("Please install the 'RandomFields' package.")
  }
  library(RandomFields)
  
  grid_data <- expand.grid(x = 1:grid_size, y = 1:grid_size)
  x_coords <- unique(grid_data$x)
  y_coords <- unique(grid_data$y)
  
  if (!is.null(seed)) {
    RFoptions(seed = seed)
  } else if (exists(".Random.seed", envir = globalenv())) {
    rm(.Random.seed, envir = globalenv())
  }
  
  cx <- mean(x_coords)
  cy <- mean(y_coords)
  bump <- matrix(0, nrow = grid_size, ncol = grid_size)
  for (i in 1:grid_size) {
    for (j in 1:grid_size) {
      bump[i, j] <- exp(-((i - cx)^2 + (j - cy)^2) / (2 * sigma_mean^2))
    }
  }
  
  mean_field <- as.vector(t(bump))
  model <- RMgauss(var = 1, scale = scale)
  sim <- RFsimulate(model, x = x_coords, y = y_coords)
  grid_data$habitat <- as.vector(sim$variable1) + mean_field
  
  if (transform_uniform) {
    grid_data$habitat <- rank(grid_data$habitat, ties.method = "average") / nrow(grid_data)
  } else {
    grid_data$habitat <- (grid_data$habitat - min(grid_data$habitat)) / 
      (max(grid_data$habitat) - min(grid_data$habitat))
  }
  
  return(grid_data)
}

# === Function: Add Highlight Rectangle ===
add_center_rectangle <- function(center = c(7, 7), size = 5, fill = NA, color = "#FF4756") {
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
    size = 2
  )
}

add_center_rectangle_black <- function(center = c(7, 7), size = 5, fill = NA, color = "darkred") {
  half <- floor(size / 2)
  xmin <- center[1] - half + .5 - 1 - .2
  xmax <- center[1] + half + .5 + .2
  ymin <- center[2] - half + .5 - 1 - .2
  ymax <- center[2] + half + .5 + .2
  
  geom_rect(
    aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = fill,
    color = "darkred",
    size = 4
  )
}

# === Seed and Grid Generation ===
SEEDSET <- 8641

grid_data <- generate_habitat_grid_with_center_clump(
  grid_size = 13,
  scale = 10,
  sigma_mean = 2,
  seed = SEEDSET, 
  transform_uniform = TRUE
)

# === Background Texture ===
texture <- rasterGrob(
  matrix(runif(10000), nrow = 200, ncol = 200), 
  width = unit(1, "npc"), height = unit(1, "npc"), 
  interpolate = TRUE
)

# === Reverse Habitat Values ===
grid_data$habitat <- -1 * (grid_data$habitat) + 1

# === Save Value at Center, Then Randomize ===
Save_Center <- grid_data$habitat[which(grid_data$x == 7 & grid_data$y == 7)]
set.seed(1575)
grid_data$habitat <- sample(grid_data$habitat, replace = FALSE)

# === Set Center to 0 ===
grid_data$habitat[which(grid_data$x == 7 & grid_data$y == 7)] <- 0

# === Plot ===
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
  add_center_rectangle_black(center = c(7, 7), size = 5) +
  add_center_rectangle(center = c(7, 7), size = 5) +
  theme_void() +
  theme(
    plot.margin = margin(5, 5, 5, 5, "mm"),
    legend.position = "none"
  )

# === Save and Draw Plot ===
ggsave("Test_Agg_RAND.svg", p, width = 8, height = 8, dpi = 400)
grid.newpage()
grid.draw(p)
