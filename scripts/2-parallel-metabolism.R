source(here::here("scripts", "0-preprocessing.R"))

fig_2B <- df %>%
  filter(
    experiment %in% c("sweep_kC", "sweep_kN")
  ) %>%
  mutate(
    varied_rate = case_when(
      experiment == "sweep_kC" ~ k_C,
      experiment == "sweep_kN" ~ k_N,
      TRUE ~ NA_real_
    )
  ) %>%
  ggplot(aes(x = varied_rate, y = mu)) +
  geom_point(aes(colour = experiment), size = 1) +
  scale_colour_manual(
    values = palette_experiment,
    labels = c(
      "sweep_kC" = expression("carbon uptake" ~ italic(k)[C]),
      "sweep_kN" = expression("nitrogen uptake" ~ italic(k)[N])
    ),
    breaks = c("sweep_kC", "sweep_kN")
  ) +
  scale_x_log10(
    n.breaks = 6,
    labels = scales::label_log(base = 10, digits = 1)
  ) +
  coord_cartesian(
    xlim = c(10^-2.2, 10^3.2),
    ylim = c(0, 3.0),
    expand = FALSE
  ) +
  guides(
    colour = guide_legend(title.position = "left")
  ) +
  labs(
    x = expression("Uptake rate" ~ (h^{
      -1
    })),
    y = expression("Growth rate" ~ italic(mu) ~ (h^{
      -1
    })),
    colour = "Parameter\nvaried"
  )
fig_2B

fig_2C <- df_alloc %>%
  filter(
    experiment %in% c("sweep_kC", "sweep_kN"),
    fraction %in% c("f_C", "f_N", "f_Ef", "f_Er")
  ) %>%
  mutate(
    facet = case_when(
      experiment == "sweep_kC" ~ "italic(k)[C]~varied",
      experiment == "sweep_kN" ~ "italic(k)[N]~varied"
    )
  ) %>%
  {
    ggplot(., aes(x = mu, y = allocation)) +
      geom_vline(
        data = filter(., experiment == "sweep_kC"),
        aes(xintercept = 0.54),
        linetype = "22",
      ) +
      geom_point(aes(colour = fraction)) +
      facet_wrap(~facet, ncol = 1, labeller = "label_parsed") +
      scale_colour_manual(
        values = palette_fractions,
        labels = c(
          "f_C" = expression("carbon uptake" ~ E[C]),
          "f_N" = expression("nitrogen uptake" ~ E[N]),
          "f_Ef" = expression("respirofermentative enzyme" ~ E[Af]),
          "f_Er" = expression("purely respiratory enzyme" ~ E[Ar])
        ),
        breaks = c("f_Ef", "f_Er", "f_N", "f_C")
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
  }
fig_2C

fig_2D <- df_mass %>%
  filter(
    experiment %in% c("sweep_kC", "sweep_kN")
  ) %>%
  mutate(
    facet = case_when(
      experiment == "sweep_kC" ~ "italic(k)[C]~varied",
      experiment == "sweep_kN" ~ "italic(k)[N]~varied"
    )
  ) %>%
  filter(state_var %in% c("c", "a", "n", "total_metabolite")) %>%
  {
    ggplot(., aes(x = mu, y = mass_fraction)) +
      geom_vline(
        data = filter(., experiment == "sweep_kC"),
        aes(xintercept = 0.54),
        linetype = "22",
      ) +
      geom_point(aes(colour = state_var), size = 1) +
      facet_wrap(~facet, ncol = 1, labeller = "label_parsed") +
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
  }
fig_2D

fig_2AB <- plot_grid(
  NULL, fig_2B,
  nrow = 1,
  rel_widths = c(7, 5),
  labels = c("A", "B"),
  label_size = 10,
  label_fontfamily = sansfamily
)
fig_2CD <- plot_grid(
  fig_2C, fig_2D,
  nrow = 1,
  rel_widths = c(6, 6),
  align = "h",
  labels = c("C", "D"),
  label_size = 10,
  label_fontfamily = sansfamily
)

fig_2 <- plot_grid(
  fig_2AB, fig_2CD,
  nrow = 2,
  rel_heights = c(6, 8)
)
ggsave(
  "2-parallel-metabolism.pdf",
  fig_2,
  path = plotdir,
  device = cairo_pdf,
  width = 12,
  height = 14,
  units = "cm"
)
ggsave(
  "2-parallel-metabolism.svg",
  fig_2,
  path = plotdir,
  width = 12,
  height = 14,
  units = "cm"
)
