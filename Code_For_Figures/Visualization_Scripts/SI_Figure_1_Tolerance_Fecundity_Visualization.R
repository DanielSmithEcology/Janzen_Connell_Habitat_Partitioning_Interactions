library(ggplot2)
library(dplyr)
library(tidyr)
library(viridis)
library(purrr)

# === Parameters ===
params <- tibble(
  Ri_over_S = c(0.1, 0.3, 0.5, 0.7),
  Z         = rev(c(15,   25,   50,  150))
)
b0 <- 0  # baseline threshold location; set to 0 so threshold ~ Ri_over_S

# === Generate data for H(x) ===
df <- pmap_dfr(params, function(Ri_over_S, Z) {
  x <- seq(0, 1, length.out = 400)
  tibble(
    x = x,
    Ri_over_S = Ri_over_S,
    Z = Z,
    H = 1 / (1 + exp(-Z * (b0 + Ri_over_S - x)))
  )
})
df <- df %>% arrange(desc(Ri_over_S))
# === Plot ===
ggplot(df, aes(x = x, y = H, group = rev(Ri_over_S), color = Ri_over_S)) +
  geom_line(linewidth = 3.0) +
  scale_color_viridis_c(option = "D") +
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
ggsave("tolerance_fecundity_Hx.svg", width = 5*3.5/1.65, height = 5*1.5/1.65, dpi = 300)
