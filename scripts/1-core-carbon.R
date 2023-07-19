source(here::here("scripts", "0-preprocessing.R"))

fig_1B <- df %>%
  filter(
    experiment == "sweep_kC_nores"
  ) %>%
  ggplot(aes(x = k_C, y = mu)) +
  geom_point(colour = "#555555", size = 1) +
  scale_x_log10(
    n.breaks = 6,
    labels = scales::label_log(base = 10, digits = 1)
  ) +
  coord_cartesian(
    xlim = c(10^-2.2, 10^3.2),
    ylim = c(0, 3.0),
    expand = FALSE
  ) +
  labs(
    x = expression("Carbon uptake rate" ~ italic(k)[C] ~ (h^{
      -1
    })),
    y = expression("Growth rate" ~ italic(mu) ~ (h^{
      -1
    }))
  )
fig_1B

fig_1C <- df_alloc %>%
  filter(
    experiment == "sweep_kC_nores"
  ) %>%
  ggplot(aes(x = mu, y = allocation)) +
  geom_point(aes(colour = fraction), size = 1) +
  scale_colour_manual(
    values = palette_fractions,
    labels = c(
      "f_C" = expression("carbon uptake" ~ E[C]),
      "f_N" = expression("nitrogen uptake" ~ E[N]),
      "f_Z" = expression("housekeeping" ~ Z),
      "f_R" = expression("ribosome" ~ R),
      "f_Ef" = expression("amino acid synthesis" ~ E[A])
    ),
    breaks = c("f_R", "f_Ef", "f_Z", "f_N", "f_C")
  ) +
  guides(
    colour = guide_legend(title.position = "left")
  ) +
  coord_cartesian(xlim = c(0, 3.0), ylim = c(0, 0.8), expand = FALSE) +
  labs(
    x = expression("Growth rate" ~ italic(mu) ~ (h^{
      -1
    })),
    y = "Allocation fraction",
    colour = "Protein"
  )
fig_1C

fig_1D <- df_mass %>%
  filter(
    experiment == "sweep_kC_nores"
  ) %>%
  filter(state_var %in% c("c", "a", "n", "total_metabolite")) %>%
  ggplot(aes(x = mu, y = mass_fraction)) +
  geom_point(aes(colour = state_var), size = 1) +
  scale_colour_manual(
    values = palette_mass,
    labels = c(
      "a" = expression("amino acid" ~ italic(a)),
      "c" = expression("carbon" ~ italic(c)),
      "n" = expression("nitrogen" ~ italic(n)),
      "total_metabolite" = "total"
    ),
    breaks = c("total_metabolite", "c", "a", "n")
  ) +
  guides(
    colour = guide_legend(title.position = "left")
  ) +
  coord_cartesian(xlim = c(0, 3.0), ylim = c(0, 0.45), expand = FALSE) +
  labs(
    x = expression("Growth rate" ~ italic(mu) ~ (h^{
      -1
    })),
    y = "Mass fraction",
    colour = "Metabolite"
  )
fig_1D

fig_1AB <- plot_grid(
  NULL, fig_1B,
  nrow = 1,
  rel_widths = c(7, 5),
  labels = c("A", "B"),
  label_size = 10,
  label_fontfamily = sansfamily
)
fig_1CD <- plot_grid(
  fig_1C, fig_1D,
  nrow = 1,
  rel_widths = c(6, 6),
  align = "h",
  labels = c("C", "D"),
  label_size = 10,
  label_fontfamily = sansfamily
)
fig_1CD

fig_1 <- plot_grid(
  fig_1AB, fig_1CD,
  nrow = 2,
  rel_heights = c(4, 6)
)

ggsave(
  "1-core-carbon.pdf",
  fig_1,
  path = plotdir,
  device = cairo_pdf,
  width = 12,
  height = 10,
  units = "cm"
)
ggsave(
  "1-core-carbon.svg",
  fig_1,
  path = plotdir,
  width = 12,
  height = 10,
  units = "cm"
)
