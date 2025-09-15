# === Libraries ===
library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)

# === Parameters ===
sigma_h <- 0.05
hopt_vals <- seq(0.1, 0.9, length.out = 6)

# === Generate data (truncated around each hopt) ===
df <- lapply(hopt_vals, function(hopt) {
  x_local <- seq(hopt - 5 * sigma_h, hopt + 5 * sigma_h, length.out = 300)
  x_local <- x_local[x_local >= 0 & x_local <= 1]
  data.frame(
    x = x_local,
    hopt = hopt,
    fitness = exp(-((x_local - hopt)^2) / (2 * sigma_h^2))
  )
}) %>% bind_rows()

# === Plot ===
ggplot(df, aes(x = x, y = fitness, group = hopt, color = hopt)) +
  geom_line(linewidth = 3.0) +
  scale_color_viridis_c(option = "D", limits = c(0, 1)) +
  scale_x_continuous(breaks = c(0, 0.5, 1), limits = c(0, 1),expand=c(0,0)) +
  scale_y_continuous(breaks = c(0, 1), limits = c(0, 1.05),expand=c(0,0)) +
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


# Save as SVG
ggsave("habitat_partitioning.svg", width = 5*3.5/1.65, height = 5*1.5/1.65, dpi = 300)

