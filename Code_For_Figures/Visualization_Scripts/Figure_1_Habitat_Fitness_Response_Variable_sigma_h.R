library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# === Parameters ===
hopt_vals <- c(0.1, 0.4, 0.6, 0.8)
sigma_vals <- c(0.02, 0.2, 0.1, 0.05)

# === Generate data ===
df <- purrr::map2_dfr(hopt_vals, sigma_vals, function(hopt, sigma_h) {
  x_local <- seq(hopt - 5 * sigma_h, hopt + 5 * sigma_h, length.out = 300)
  x_local <- x_local[x_local >= 0 & x_local <= 1]
  data.frame(
    x = x_local,
    hopt = hopt,
    sigma_h = sigma_h,
    fitness = exp(-((x_local - hopt)^2) / (2 * sigma_h^2))
  )
})

# === Plot ===
ggplot(df, aes(x = x, y = fitness, group = hopt, color = hopt)) +
  geom_line(linewidth = 3.0) +
  scale_color_viridis_c(option = "D") +#limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1.05), expand = c(0, 0)) +
  labs(x = " ", y = " ") +
  theme_minimal(base_size = 53) +
  theme(
    legend.position = "none",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    axis.line.x = element_line(color = "black", linewidth = 3),
    axis.line.y = element_line(color = "black", linewidth = 3),
    axis.ticks = element_line(color = "black", linewidth = 3),
    axis.ticks.length = unit(0.3, "cm"),
    axis.text = element_text(size = 31, color = "black"),
    axis.title = element_text(size = 32, color = "black")
  )

# === Save as SVG ===
ggsave("habitat_partitioning_varsigmah.svg", width = 5*3.5/1.65, height = 5*1.5/1.65, dpi = 300)
