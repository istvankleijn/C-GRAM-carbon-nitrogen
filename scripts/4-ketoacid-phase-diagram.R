source(here::here("scripts", "0-preprocessing.R"))

fig_4A <- df %>%
  filter(
    experiment == "k_Kex=20.0_k_N=20.0"
  ) %>%
  ggplot(aes(x = C_per_N, y = k_Kre)) +
  geom_raster(aes(fill = zero_alloc)) +
  geom_hline(yintercept = c(5.0, 10.0), linetype = "dashed") +
  scale_y_log10(
    breaks = c(0.5, 1, 2, 3, 5, 10, 20, 30, 50),
    labels = c("0.5", "1", "2", "3", "5", "10", "20", "30", "50")
  ) +
  scale_fill_manual(values = palette_zeros) +
  coord_cartesian(expand = FALSE) +
  annotation_logticks(
    sides = "lr",
    mid = unit(0.1, "cm"),
    long = unit(0.2, "cm"),
    size = 0.25,
    colour = "grey20"
  ) +
  labs(
    x = "Nutrient C-to-N ratio",
    y = expression("Recycling efficiency" ~ k[Kre]),
    fill = "Zero allocation"
  ) +
  guides(
    fill = guide_legend(title.position = "left")
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical"
  )

fig_4B <- df %>%
  filter(
    experiment == "k_Kre=10.0_k_N=20.0"
  ) %>%
  ggplot(aes(x = C_per_N, y = k_Kex)) +
  geom_raster(aes(fill = zero_alloc), show.legend = FALSE) +
  geom_hline(yintercept = 20.0, linetype = "dashed") +
  scale_y_log10(
    breaks = c(2, 3, 5, 10, 20, 30, 50, 100, 200),
    labels = c("2", "3", "5", "10", "20", "30", "50", "100", "200")
  ) +
  scale_fill_manual(values = palette_zeros) +
  coord_cartesian(expand = FALSE) +
  annotation_logticks(
    sides = "lr",
    mid = unit(0.1, "cm"),
    long = unit(0.2, "cm"),
    size = 0.25,
    colour = "grey20"
  ) +
  labs(
    x = "Nutrient C-to-N ratio",
    y = expression("Excretion efficiency" ~ k[Kex])
  )

fig_4D <- df %>%
  filter(
    experiment == "k_Kre=5.0_k_N=20.0"
  ) %>%
  ggplot(aes(x = C_per_N, y = k_Kex)) +
  geom_raster(aes(fill = zero_alloc), show.legend = FALSE) +
  geom_hline(yintercept = 20.0, linetype = "dashed") +
  scale_y_log10(
    breaks = c(2, 3, 5, 10, 20, 30, 50, 100, 200),
    labels = c("2", "3", "5", "10", "20", "30", "50", "100", "200")
  ) +
  scale_fill_manual(values = palette_zeros) +
  coord_cartesian(expand = FALSE) +
  annotation_logticks(
    sides = "lr",
    mid = unit(0.1, "cm"),
    long = unit(0.2, "cm"),
    size = 0.25,
    colour = "grey20"
  ) +
  labs(
    x = "Nutrient C-to-N ratio",
    y = expression("Excretion efficiency" ~ k[Kex])
  )

fig_4C <- df %>%
  filter(
    experiment == "k_Kex=20.0_k_Kre=10.0"
  ) %>%
  ggplot(aes(x = C_per_N, y = k_N)) +
  geom_raster(aes(fill = zero_alloc), show.legend = FALSE) +
  geom_hline(yintercept = 20.0, linetype = "dashed") +
  scale_y_log10(
    breaks = c(0.5, 1, 2, 3, 5, 10, 20, 30, 50),
    labels = c("0.5", "1", "2", "3", "5", "10", "20", "30", "50")
  ) +
  scale_fill_manual(values = palette_zeros) +
  coord_cartesian(expand = FALSE) +
  annotation_logticks(
    sides = "lr",
    mid = unit(0.1, "cm"),
    long = unit(0.2, "cm"),
    size = 0.25,
    colour = "grey20"
  ) +
  labs(
    x = "Nutrient C-to-N ratio",
    y = expression("Uptake efficiency" ~ k[N])
  )

fig_4E <- df %>%
  filter(
    experiment == "k_Kex=20.0_k_Kre=5.0"
  ) %>%
  ggplot(aes(x = C_per_N, y = k_N)) +
  geom_raster(aes(fill = zero_alloc), show.legend = FALSE) +
  geom_hline(yintercept = 20.0, linetype = "dashed") +
  scale_y_log10(
    breaks = c(0.5, 1, 2, 3, 5, 10, 20, 30, 50),
    labels = c("0.5", "1", "2", "3", "5", "10", "20", "30", "50")
  ) +
  scale_fill_manual(values = palette_zeros) +
  coord_cartesian(expand = FALSE) +
  annotation_logticks(
    sides = "lr",
    mid = unit(0.1, "cm"),
    long = unit(0.2, "cm"),
    size = 0.25,
    colour = "grey20"
  ) +
  labs(
    x = "Nutrient C-to-N ratio",
    y = expression("Uptake efficiency" ~ k[N])
  )

fig_4BCDE <- plot_grid(
  fig_4B, fig_4C,
  fig_4D, fig_4E,
  nrow = 2,
  labels = c("B", "C", "D", "E"),
  label_fontfamily = sansfamily,
  label_size = 10
)
fig_4 <- plot_grid(
  fig_4A, fig_4BCDE,
  ncol = 2,
  rel_widths = c(1, 2),
  labels = c("A", ""),
  label_fontfamily = sansfamily,
  label_size = 8
)
ggsave(
  "4-ketoacid-phase-diagram.svg",
  fig_4,
  path = plotdir,
  width = 17.5,
  height = 10,
  units = "cm"
)
ggsave(
  "4-ketoacid-phase-diagram.pdf",
  fig_4,
  device = cairo_pdf,
  path = plotdir,
  width = 17.5,
  height = 12.5,
  units = "cm"
)
