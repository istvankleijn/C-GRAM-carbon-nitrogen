source(here::here("scripts", "0-preprocessing.R"))

fig_3B <- df %>%
  filter(
    experiment %in% c("k_Kre=10.0_k_N=20.0", "k_Kre=5.0_k_N=20.0"),
    19.0 < k_Kex, k_Kex < 21.0
  ) %>%
  ggplot(aes(x = C_per_N, y = mu)) +
  geom_point(aes(colour = experiment), size = 1) +
  scale_colour_manual(
    values = palette_experiment,
    labels = c(
      "k_Kre=10.0_k_N=20.0" = expression(italic(k)[Kre] == 10),
      "k_Kre=5.0_k_N=20.0" = expression(italic(k)[Kre] == 5)
    ),
    breaks = c("k_Kre=10.0_k_N=20.0", "k_Kre=5.0_k_N=20.0")
  ) +
  coord_cartesian(
    ylim = c(0, 2.0),
    expand = FALSE
  ) +
  guides(
    colour = guide_legend(title.position = "left")
  ) +
  labs(
    x = "Nutrient C-to-N ratio",
    y = expression("Growth rate" ~ italic(mu) ~ (h^{
      -1
    })),
    colour = "Recycling\nrate"
  ) +
  theme(
    legend.text.align = 0
  )
fig_3B


fig_3C <- df_alloc %>%
  filter(
    experiment %in% c("k_Kre=10.0_k_N=20.0", "k_Kre=5.0_k_N=20.0"),
    19.0 < k_Kex, k_Kex < 21.0,
    fraction %in% c("f_C", "f_N", "f_Ef", "f_R", "f_Kex", "f_Kre")
  ) %>%
  mutate(
    facet = case_when(
      experiment == "k_Kre=10.0_k_N=20.0" ~ "italic(k)[Kre] == 10",
      experiment == "k_Kre=5.0_k_N=20.0" ~ "italic(k)[Kre] == 5"
    )
  ) %>%
  {
    ggplot(., aes(x = mu, y = allocation)) +
      geom_vline(
        data = filter(., experiment == "k_Kre=10.0_k_N=20.0"),
        aes(xintercept = 1.44),
        linetype = "22"
      ) +
      geom_vline(
        data = filter(., experiment == "k_Kre=10.0_k_N=20.0"),
        aes(xintercept = 1.56),
        linetype = "22"
      ) +
      geom_vline(
        data = filter(., experiment == "k_Kre=5.0_k_N=20.0"),
        aes(xintercept = 1.49),
        linetype = "22"
      ) +
      geom_point(aes(colour = fraction), size = 1) +
      facet_wrap(~facet, ncol = 1, labeller = "label_parsed") +
      scale_colour_manual(
        values = palette_fractions,
        labels = c(
          "f_C" = expression("carbon uptake" ~ E[C]),
          "f_N" = expression("nitrogen uptake" ~ E[N]),
          "f_R" = expression("ribosome" ~ R),
          "f_Ef" = expression("amino acid synthesis" ~ E[Af]),
          "f_Kre" = expression("ketoacid recycling" ~ E[Kre]),
          "f_Kex" = expression("ketoacid excretion" ~ E[Kex])
        ),
        breaks = c("f_R", "f_Ef", "f_C", "f_Kex", "f_Kre", "f_N")
      ) +
      guides(
        colour = guide_legend(title.position = "left")
      ) +
      coord_cartesian(xlim = c(1.0, 2.0), ylim = c(0, 0.4), expand = FALSE) +
      labs(
        x = expression("Growth rate" ~ italic(mu) ~ (h^{
          -1
        })),
        y = "Allocation fraction",
        colour = "Protein"
      )
  }
fig_3C


fig_3D <- df_mass %>%
  filter(
    experiment %in% c("k_Kre=10.0_k_N=20.0", "k_Kre=5.0_k_N=20.0"),
    19.0 < k_Kex, k_Kex < 21.0,
    state_var %in% c("k", "c", "a", "n")
  ) %>%
  mutate(
    facet = case_when(
      experiment == "k_Kre=10.0_k_N=20.0" ~ "italic(k)[Kre] == 10",
      experiment == "k_Kre=5.0_k_N=20.0" ~ "italic(k)[Kre] == 5"
    )
  ) %>%
  {
    ggplot(., aes(x = mu, y = mass_fraction)) +
      geom_vline(
        data = filter(., experiment == "k_Kre=10.0_k_N=20.0"),
        aes(xintercept = 1.44),
        linetype = "22"
      ) +
      geom_vline(
        data = filter(., experiment == "k_Kre=10.0_k_N=20.0"),
        aes(xintercept = 1.56),
        linetype = "22"
      ) +
      geom_vline(
        data = filter(., experiment == "k_Kre=5.0_k_N=20.0"),
        aes(xintercept = 1.49),
        linetype = "22"
      ) +
      geom_point(aes(colour = state_var), size = 1) +
      facet_wrap(~facet, ncol = 1, labeller = "label_parsed") +
      scale_colour_manual(
        values = palette_mass,
        labels = c(
          "k" = expression("ketoacid" ~ italic(k)),
          "a" = expression("amino acid" ~ italic(a)),
          "c" = expression("carbon" ~ italic(c)),
          "n" = expression("nitrogen" ~ italic(n))
        ),
        breaks = c("k", "c", "a", "n")
      ) +
      guides(
        colour = guide_legend(title.position = "left")
      ) +
      coord_cartesian(xlim = c(1.0, 2.0), ylim = c(0, 0.5), expand = FALSE) +
      labs(
        x = expression("Growth rate" ~ italic(mu) ~ (h^{
          -1
        })),
        y = "Mass fraction",
        colour = "Metabolite"
      )
  }
fig_3D

fig_3AB <- plot_grid(
  NULL, fig_3B,
  nrow = 1,
  rel_widths = c(7, 5),
  labels = c("A", "B"),
  label_size = base_size,
  label_fontfamily = sansfamily
)
fig_3CD <- plot_grid(
  fig_3C, fig_3D,
  nrow = 1,
  rel_widths = c(6, 6),
  align = "h",
  labels = c("C", "D"),
  label_size = base_size,
  label_fontfamily = sansfamily
)

fig_3 <- plot_grid(
  fig_3AB, fig_3CD,
  nrow = 2,
  rel_heights = c(5, 9)
)
ggsave(
  "3-ketoacid-sweeps.pdf",
  fig_3,
  path = plotdir,
  device = cairo_pdf,
  width = 12,
  height = 14,
  units = "cm"
)
ggsave(
  "3-ketoacid-sweeps.svg",
  fig_3,
  path = plotdir,
  width = 12,
  height = 14,
  units = "cm"
)
