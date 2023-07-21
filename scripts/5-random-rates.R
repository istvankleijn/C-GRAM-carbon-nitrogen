source(here::here("scripts", "0-preprocessing.R"))

fig_5A <- df_alloc %>%
  filter(
    str_detect(experiment, "random_kCkEfkN"),
    k_Kre == 10.0, γ_K == 0
  ) %>%
  ggplot(aes(x = mu, y = allocation)) +
  geom_point(aes(colour = fraction), alpha = 0.3, shape = 1) +
  scale_colour_manual(
    values = palette_fractions,
    labels = c(
      "f_C" = expression("carbon uptake" ~ E[C]),
      "f_N" = expression("nitrogen uptake" ~ E[N]),
      "f_Ef" = expression("respirofermentative enzyme" ~ E[Af]),
      "f_Er" = expression("purely respiratory enzyme" ~ E[Ar]),
      "f_R" = expression("ribosome R")
    ),
    breaks = c("f_R", "f_Ef", "f_Er", "f_N", "f_C")
  ) +
  guides(
    colour = guide_legend(title.position = "left")
  ) +
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 0.800001), expand = FALSE) +
  labs(
    x = expression("Growth rate" ~ italic(mu) ~ (h^{
      -1
    })),
    y = "Allocation fraction",
    colour = "Protein"
  )
fig_5A

fig_5B <- df_mass %>%
  filter(
    str_detect(experiment, "random_kCkEfkN"),
    state_var %in% c("c", "a", "n", "total_metabolite"),
    k_Kre == 10.0, γ_K == 0
  ) %>%
  ggplot(aes(x = mu, y = mass_fraction)) +
  geom_point(aes(colour = state_var), alpha = 0.3, shape = 1) +
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
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, NA), expand = FALSE) +
  labs(
    x = expression("Growth rate" ~ italic(mu) ~ (h^{
      -1
    })),
    y = "Mass fraction",
    colour = "Metabolite"
  )
fig_5B


df_5C <- df %>%
  filter(
    str_detect(experiment, "random_kCkEfkN"),
    γ_K == 0, k_Kre == 10.0
  )

model_5C <- lm(
  f_R ~ a,
  data = df_5C
)
delta <- coef(model_5C)["a"]
asat <- -coef(model_5C)[["(Intercept)"]] / delta
label <- substitute(
  italic(f)[R] == delta %*% (italic(a) - asat),
  list(
    delta = format(delta, digits = 3),
    asat = format(asat, digits = 3)
  )
)
fig_5C <- df_5C %>%
  ggplot(aes(x = a, y = f_R)) +
  geom_point(aes(colour = energy_type, shape = energy_type), alpha = 0.3) +
  geom_abline(
    slope = delta,
    intercept = -delta * asat,
    colour = "#555555",
    linewidth = 0.25,
    linetype = "21"
  ) +
  scale_colour_manual(values = palette_fermresp) +
  scale_shape_manual(values = shapes_fermresp) +
  coord_cartesian(xlim = c(0, 0.12), ylim = c(0, 0.45), expand = FALSE) +
  labs(
    x = "Amino acid mass fraction",
    y = "Ribosomal allocation fraction",
    colour = "Energy\nmetabolism",
    shape = "Energy\nmetabolism"
  ) +
  annotate(
    "text",
    label = label,
    x = 0.115, y = 0.04,
    size = 8 / .pt,
    family = sansfamily,
    colour = "#555555",
    hjust = 1
  ) +
  theme(
    legend.position = c(0.004, 0.995),
    legend.justification = c(0, 1),
    legend.box.background = element_rect(colour = "grey20")
  )
fig_5C


df_5D <- df %>%
  filter(
    str_detect(experiment, "random_kCkEfkN"),
    C_per_N != 0
  ) %>%
  arrange(γ_K) %>%
  mutate(
    C_per_N_parsed = str_c("C/N == ", C_per_N) %>%
      as_factor(),
    k_Kre_parsed = str_c("k[Kre] == ", k_Kre)
  )
fig_5D <- df_5D %>%
  ggplot(aes(x = a)) +
  geom_point(
    aes(
      y = f_R,
      shape = energy_type,
      colour = ketoacid_type
    ),
    alpha = 0.3,
    shape = 1
  ) +
  geom_abline(
    slope = delta,
    intercept = -delta * asat,
    colour = "#555555",
    linewidth = 0.25,
    linetype = "21"
  ) +
  facet_grid(k_Kre_parsed ~ C_per_N_parsed, labeller = "label_parsed") +
  scale_x_continuous(breaks = seq(0, 0.09, 0.03)) +
  scale_colour_manual(values = palette_ketoacid) +
  scale_shape_manual(values = shapes_fermresp) +
  coord_cartesian(xlim = c(0, 0.09), ylim = c(0, 0.4), expand = FALSE) +
  labs(
    x = "Amino acid mass fraction",
    y = "Ribosomal allocation fraction",
    colour = "Ketoacid\nmetabolism"
  ) +
  guides(
    colour = guide_legend(override.aes = list(shape = 1))
  ) +
  theme(
    legend.position = "right"
  )
fig_5D


fig_5ABC <- plot_grid(
  fig_5A, fig_5B, fig_5C,
  nrow = 1,
  labels = c("A", "B", "C"),
  label_fontfamily = sansfamily,
  label_size = 10
)
fig_5 <- plot_grid(
  fig_5ABC, fig_5D,
  ncol = 1,
  rel_heights = c(6, 9),
  labels = c("", "D"),
  label_fontfamily = sansfamily,
  label_size = 10
)
ggsave(
  "5-random-rates.svg",
  fig_5,
  path = plotdir,
  width = 17.5,
  height = 15,
  units = "cm"
)
ggsave(
  "5-random-rates.pdf",
  fig_5,
  device = cairo_pdf,
  path = plotdir,
  width = 17.5,
  height = 12.5,
  units = "cm"
)
