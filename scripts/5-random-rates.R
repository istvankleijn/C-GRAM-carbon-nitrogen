source(here::here("scripts", "0-preprocessing.R"))

datadir <- here("data", "pomballoc")

cultures <- read_csv(file.path(datadir, "TableS1.csv"))
proteome_fractions <- read_csv(file.path(datadir, "TableS4.csv"))
protein_annotations <- read_csv(file.path(datadir, "TableS16.csv"))

df_5A <- df_alloc %>%
  filter(
    str_detect(experiment, "random_kCkEfkN"),
    k_Kre == 10.0, γ_K == 0,
    fraction %in% c("f_R", "f_Ef", "f_Er", "f_N", "f_C")
  )
fig_5A <- df_5A %>%
  ggplot(aes(x = mu, y = allocation, colour = fraction)) +
  geom_smooth(
    data = filter(df_5A, fraction == "f_R"),
    alpha = 0.3,
    linetype = "21",
    linewidth = 0.5,
    show.legend = FALSE,
    method = "lm",
    se = FALSE,
    fullrange = TRUE
  ) +
  geom_point(alpha = 0.3, shape = 1) +
  scale_x_continuous(limits = c(0, 2.5)) +
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
    x = expression("Growth rate" ~ italic(µ) ~ (h^{
      -1
    })),
    y = "Allocation fraction",
    colour = "Protein"
  ) +
  theme(
    axis.title.y = element_text(
      margin = margin(l = base_size / 2, r = base_size / 4, unit = "pt")
    )
  )
fig_5A

df_5B <- proteome_fractions %>%
  left_join(cultures, by = c("medium", "replicate")) %>%
  right_join(protein_annotations, by = "PomBaseIDs") %>%
  mutate(
    coarse_grained_equivalent = case_when(
      class == "0: Translation/RiBi" ~ "f_R",
      class == "1: Glycolysis" ~ "f_C",
      class == "2: Precursors/Energy" ~ "f_Ef",
      class == "3: Amino acids" ~ "f_Ef",
      TRUE ~ NA_character_
    )
  ) %>%
  group_by(coarse_grained_equivalent, medium, growth_rate) %>%
  summarise(
    proteome_fraction = sum(proteome_fraction, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  drop_na()
fig_5B <- df_5B %>%
  ggplot(
    aes(
      x = growth_rate,
      y = proteome_fraction,
      colour = coarse_grained_equivalent
    )
  ) +
  geom_smooth(
    data = filter(df_5B, coarse_grained_equivalent == "f_R"),
    alpha = 0.3,
    linetype = "21",
    linewidth = 0.5,
    show.legend = FALSE,
    method = "lm",
    se = FALSE,
    fullrange = TRUE
  ) +
  geom_point(alpha = 0.7, shape = 0) +
  scale_x_continuous(limits = c(0, 0.3)) +
  scale_colour_manual(
    values = palette_fractions,
    labels = c(
      "f_C" = expression("carbon uptake" ~ E[C]),
      "f_Ef" = expression("amino acid synthesis" ~ E[A]),
      "f_R" = expression("ribosome R")
    ),
    breaks = c("f_R", "f_Ef", "f_C")
  ) +
  guides(
    colour = guide_legend(title.position = "left")
  ) +
  coord_cartesian(xlim = c(0, 0.30001), ylim = c(0, 0.4), expand = FALSE) +
  guides(
    colour = guide_legend(title.position = "left")
  ) +
  labs(
    x = expression("Growth rate" ~ italic(µ) ~ (h^{
      -1
    })),
    y = "Proteome fraction",
    colour = "Coarse-grained\nequivalent"
  ) +
  theme(
    axis.title.y = element_text(
      margin = margin(l = base_size / 2, r = base_size / 4, unit = "pt")
    )
  )
fig_5B

fig_5C <- df_mass %>%
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
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 0.7), expand = FALSE) +
  labs(
    x = expression("Growth rate" ~ italic(µ) ~ (h^{
      -1
    })),
    y = "Mass fraction",
    colour = "Metabolite"
  )
fig_5C


df_5D <- df %>%
  filter(
    str_detect(experiment, "random_kCkEfkN"),
    γ_K == 0, k_Kre == 10.0
  )

model_5D <- lm(
  f_R ~ a,
  data = df_5D
)
delta <- coef(model_5D)["a"]
asat <- -coef(model_5D)[["(Intercept)"]] / delta
label <- substitute(
  italic(f)[R] == delta %*% (italic(a) - asat),
  list(
    delta = format(delta, digits = 3),
    asat = format(asat, digits = 3)
  )
)
fig_5D <- df_5D %>%
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
fig_5D


df_5E <- df %>%
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
fig_5E <- df_5E %>%
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
fig_5E


fig_5ABC <- plot_grid(
  fig_5A, fig_5B, fig_5C,
  nrow = 1,
  align = "h",
  labels = c("A", "B", "C"),
  label_fontfamily = sansfamily,
  label_size = 10
)
fig_5DE <- plot_grid(
  fig_5D, fig_5E,
  nrow = 1,
  rel_widths = c(1, 2),
  labels = c("D", "E"),
  label_fontfamily = sansfamily,
  label_size = 10
)
fig_5 <- plot_grid(
  fig_5ABC, fig_5DE,
  ncol = 1,
  rel_heights = c(6, 8),
  labels = c("", "D"),
  label_fontfamily = sansfamily,
  label_size = 10
)
ggsave(
  "5-random-rates.svg",
  fig_5,
  path = plotdir,
  width = 17.5,
  height = 14,
  units = "cm"
)
ggsave(
  "5-random-rates.pdf",
  fig_5,
  device = cairo_pdf,
  path = plotdir,
  width = 17.5,
  height = 14,
  units = "cm"
)
