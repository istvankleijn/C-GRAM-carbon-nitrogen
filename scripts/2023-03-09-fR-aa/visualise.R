library(here)
source(here("src", "plot_setup.R"))

modelname <- "ketoacid-excretion"
date <- "2023-03-09"
datadir <- here("data", modelname, date)
plotdir <- here("plots", "2023-03-09-fR-aa")

cutoff <- 1e-8

df <- read_csv(file.path(datadir, "results.csv")) %>%
  rename(
    Delta_CN = `<U+0394>_CN`,
    gamma_K = `<U+03B3>_K`,
    alpha_Cf = `a_Cf`, # String marks do not work for \alpha
    alpha_Cr = `a_Cr`, # but backticks do not work for other Greek characters
  ) %>%
  mutate(
    # zero_alloc = case_when(
    #   f_C <= cutoff & f_Kex <= cutoff & f_Kre <= cutoff ~ "C,Kex,Kre",
    #   f_C <= cutoff & f_Kex <= cutoff                   ~ "C,Kex",
    #   f_C <= cutoff &                   f_Kre <= cutoff ~ "C,Kre",
    #                   f_Kex <= cutoff & f_Kre <= cutoff ~ "Kex,Kre",
    #                   f_Kex <= cutoff                   ~ "Kex",
    #                                     f_Kre <= cutoff ~ "Kre",
    #   f_C <= cutoff                                     ~ "C",
    #   TRUE                                              ~ "none"
    # ),
    experiment = str_replace_all(experiment, "\u03B3", "gamma")
  )
df_alloc <- df %>%
  pivot_longer(
    cols = starts_with("f_"),
    names_to = "fraction",
    values_to = "allocation"
  )
state_vars <- c(
  "k", "c", "a", "n",
  "e_C", "e_Af", "e_Ar", "e_N", "e_Kre", "e_Kex", "r", "z"
)
df_mass <- df %>%
  pivot_longer(
    cols = all_of(state_vars),
    names_to = "state_var",
    values_to = "mass_fraction"
  ) %>%
  mutate(
    var_type = if_else(
      state_var %in% c("k", "c", "a", "n"),
      "metabolite",
      "protein"
    )
  )



labels_proteome <- tribble(
  ~fraction, ~label, ~colour,
  "f_Ef", "E[Af]", "#009E73",
  "f_Er", "E[Ar]", "#F0E442",
  "f_C", "E[C]", "#56B4E9",
  "f_Kex", "E[Kex]", "#D55E00",
  "f_Kre", "E[Kre]", "#630a75",
  "f_N", "E[N]", "#CC79A7",
  # "f_Z",      "Z",      "#555555",
  "f_R", "R", "#E69F00"
)

df_rf <- df_alloc %>%
  mutate(
    metabolism_type = case_when(
      e_Af > cutoff & e_Ar <= cutoff ~ "fermentation",
      e_Af <= cutoff & e_Ar > cutoff ~ "respiration",
      TRUE ~ "mixed"
    ),
    ketoacid_type = case_when(
      e_Kre > cutoff & e_Kex <= cutoff ~ "recycling",
      e_Kre <= cutoff & e_Kex > cutoff ~ "excretion",
      e_Kre <= cutoff & e_Kex <= cutoff ~ "none",
      TRUE ~ "mixed",
    )
  ) %>%
  # filter(
  #     experiment == "random_kCkEfkN"
  # ) %>%
  filter(fraction == "f_R")


df_alloc %>%
  mutate(
    metabolism_type = case_when(
      e_Af > cutoff & e_Ar <= cutoff ~ "fermentation",
      e_Af <= cutoff & e_Ar > cutoff ~ "respiration",
      TRUE ~ "mixed"
    ),
    ketoacid_type = case_when(
      e_Kre > cutoff & e_Kex <= cutoff ~ "recycling",
      e_Kre <= cutoff & e_Kex > cutoff ~ "excretion",
      e_Kre <= cutoff & e_Kex <= cutoff ~ "none",
      TRUE ~ "mixed",
    )
  ) %>%
  filter(metabolism_type == "mixed")

ribosome_fraction_fig <- df_rf %>%
  ggplot(aes(x = mu, y = allocation)) +
  geom_point(aes(colour = ketoacid_type, shape = gamma_K), alpha = 0.1) +
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 0.5), expand = FALSE) +
  labs(
    x = expression("Growth rate" ~ italic(mu) ~ (h^{
      -1
    })),
    y = "Ribosome allocation fraction",
    shape = expression(italic(gamma)[K]),
    colour = "Ketoacid\nmetabolism"
  ) +
  theme(
    legend.text.align = 0
  )
ribosome_fraction_fig
ggsave(
  "ribosome_fraction_ketoacid.svg",
  ribosome_fraction_fig,
  path = plotdir,
  width = 8.75,
  height = 8.75,
  units = "cm"
)

ribosome_fraction_metab_fig <- df_rf %>%
  ggplot(aes(x = mu, y = allocation)) +
  geom_point(aes(colour = metabolism_type, shape = gamma_K, alpha = metabolism_type)) +
  scale_alpha_manual(values = c(0.1, 1.0, 0.1)) +
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 0.5), expand = FALSE) +
  labs(
    x = expression("Growth rate" ~ italic(mu) ~ (h^{
      -1
    })),
    y = "Ribosome allocation fraction",
    shape = expression(italic(gamma)[K]),
    colour = "Amino acid\nmetabolism"
  ) +
  theme(
    legend.text.align = 0
  )
ribosome_fraction_metab_fig
ggsave(
  "ribosome_fraction_aa_metab.svg",
  ribosome_fraction_metab_fig,
  path = plotdir,
  width = 8.75,
  height = 8.75,
  units = "cm"
)


r_TER <- df %>%
  mutate(
    energy_type = case_when(
      f_Ef > cutoff & f_Er <= cutoff ~ "fermentation",
      f_Ef <= cutoff & f_Er > cutoff ~ "respiration",
      TRUE ~ "mixed"
    ),
    ketoacid_type = case_when(
      f_Kre > cutoff & f_Kex <= cutoff ~ "recycling",
      f_Kre <= cutoff & f_Kex > cutoff ~ "excretion",
      f_Kre <= cutoff & f_Kex <= cutoff ~ "none",
      TRUE ~ "mixed"
    ),
    metabolism_type = str_c(ketoacid_type, energy_type, sep = ":")
  ) %>%
  ggplot(aes(x = r)) +
  geom_point(aes(y = a / (a + a_sat), colour = metabolism_type), alpha = 0.5) +
  coord_cartesian(xlim = c(0, 0.5), ylim = c(0, 1), expand = FALSE) +
  # geom_line(aes(y = pred), colour = "black") +
  # facet_grid(gamma_K ~ k_Kre) +
  labs(
    x = "Ribosome mass fraction",
    y = "relative translation elongation rate",
    colour = "Metabolism"
  )
r_TER
ggsave(
  "r-TER.svg",
  r_TER,
  path = plotdir,
  width = 17.5,
  height = 12,
  units = "cm"
)


labels_5A <- tribble(
  ~fraction, ~label, ~colour,
  "f_Ef", "E[Af]", "#009E73",
  "f_Er", "E[Ar]", "#F0E442",
  "f_C", "E[C]", "#56B4E9",
  "f_N", "E[N]", "#CC79A7",
  # "f_Z",      "Z",      "#555555",
  "f_R", "R", "#E69F00",
)
fig_5A <- df_alloc %>%
  filter(
    k_Kre == 10.0, gamma_K == 0
  ) %>%
  left_join(labels_5A, by = "fraction") %>%
  drop_na(label) %>%
  ggplot(aes(x = mu, y = allocation)) +
  geom_point(aes(colour = label), alpha = 0.3, shape = 1) +
  scale_colour_manual(
    values = deframe(select(labels_5A, label, colour)),
    labels = scales::parse_format()
  ) +
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 0.800001), expand = FALSE) +
  labs(
    x = expression("Growth rate" ~ italic(mu) ~ (h^{
      -1
    })),
    y = "Allocation fraction",
    colour = "Protein"
  ) +
  theme(
    legend.text.align = 0
  )
# fig_5A
ggsave(
  "fig_5A.svg",
  fig_5A,
  path = plotdir,
  width = 6,
  height = 4.5,
  units = "cm"
)


palette_5B <- c(
  "a"    = "#D55E00", # vermillion
  "c"    = "#0072B2", # blue
  "n"    = "#630a75" # purple
)
fig_5B <- df_mass %>%
  filter(
    state_var %in% c("c", "a", "n"),
    k_Kre == 10.0, gamma_K == 0
  ) %>%
  ggplot(aes(x = mu, y = mass_fraction)) +
  geom_point(aes(colour = state_var), alpha = 0.3, shape = 1) +
  scale_colour_manual(values = palette_5B) +
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 0.300001), expand = FALSE) +
  labs(
    x = expression("Growth rate" ~ italic(mu) ~ (h^{
      -1
    })),
    y = "Mass fraction",
    colour = "Metabolite"
  )
# fig_5B
ggsave(
  "fig_5B.svg",
  fig_5B,
  path = plotdir,
  width = 6,
  height = 4.5,
  units = "cm"
)


palette_5C <- c(
  "fermentation" = "#555555", # grey
  "respiration"  = "#E69F00" # orange
)
shapes_5C <- c(
  "fermentation" = 1, # circle
  "respiration" = 0 # square
)
df_5C <- df %>%
  filter(
    gamma_K == 0, k_Kre == 10.0
  ) %>%
  mutate(
    energy_type = case_when(
      f_Ef > cutoff & f_Er <= cutoff ~ "fermentation",
      f_Ef <= cutoff & f_Er > cutoff ~ "respiration",
      TRUE ~ "mixed"
    ),
    ketoacid_type = case_when(
      f_Kre <= cutoff & f_Kex <= cutoff ~ "none",
      f_Kre > cutoff & f_Kex <= cutoff ~ "recycling",
      f_Kre <= cutoff & f_Kex > cutoff ~ "excretion",
      TRUE ~ "mixed"
    ),
    metabolism_type = str_c(ketoacid_type, energy_type, sep = ":")
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
    size = 0.25,
    linetype = "21"
  ) +
  scale_colour_manual(values = palette_5C) +
  scale_shape_manual(values = shapes_5C) +
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
    size = 6 / .pt,
    family = sansfamily,
    colour = "#555555",
    hjust = 1
  ) +
  theme(
    legend.position = c(0.004, 0.995),
    legend.justification = c(0, 1),
    legend.box.background = element_rect(colour = "grey20")
  )
ggsave(
  "fig_5C.svg",
  fig_5C,
  path = plotdir,
  width = 6,
  height = 4.5,
  units = "cm"
)

palette_5D <- c(
  "none"      = "#555555", # grey
  "excretion" = "#0072B2", # blue
  "recycling" = "#D55E00", # vermillion
  "mixed"     = "#009E73" # green
)
df_5D <- df %>%
  filter(
    gamma_K != 0
  ) %>%
  mutate(
    energy_type = case_when(
      f_Ef > cutoff & f_Er <= cutoff ~ "fermentation",
      f_Ef <= cutoff & f_Er > cutoff ~ "respiration",
      TRUE ~ "mixed"
    ),
    ketoacid_type = case_when(
      f_Kre <= cutoff & f_Kex <= cutoff ~ "none",
      f_Kre > cutoff & f_Kex <= cutoff ~ "recycling",
      f_Kre <= cutoff & f_Kex > cutoff ~ "excretion",
      TRUE ~ "mixed"
    ),
    metabolism_type = str_c(ketoacid_type, energy_type, sep = ":"),
    k_Kre_parsed = str_c("k[Kre]==", k_Kre)
  ) %>%
  separate(
    C_per_N,
    into = c("C_per_N_numerator", "C_per_N_denominator"),
    sep = "//",
    remove = FALSE,
    convert = TRUE
  ) %>%
  arrange(gamma_K) %>%
  mutate(
    C_per_N_number = C_per_N_numerator / C_per_N_denominator,
    C_per_N_parsed = str_c("C/N==", C_per_N_number) %>%
      as_factor()
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
    size = 0.25,
    linetype = "21"
  ) +
  facet_grid(k_Kre_parsed ~ C_per_N_parsed, labeller = "label_parsed") +
  scale_x_continuous(breaks = seq(0, 0.09, 0.03)) +
  scale_colour_manual(values = palette_5D) +
  scale_shape_manual(values = shapes_5C) +
  coord_cartesian(xlim = c(0, 0.09), ylim = c(0, 0.4), expand = FALSE) +
  labs(
    x = "Amino acid mass fraction",
    y = "Ribosomal allocation fraction",
    colour = "Ketoacid\nmetabolism"
  ) +
  guides(
    colour = guide_legend(override.aes = list(shape = 1))
  )
ggsave(
  "fig_5D.svg",
  fig_5D,
  path = plotdir,
  width = 18,
  height = 9,
  units = "cm"
)


fig_5ABC <- plot_grid(
  fig_5A, fig_5B, fig_5C,
  nrow = 1,
  labels = c("A", "B", "C"),
  label_fontfamily = sansfamily,
  label_size = 8
)
fig_5ABCD <- plot_grid(
  fig_5ABC, fig_5D,
  ncol = 1,
  rel_heights = c(4.5, 8),
  labels = c("", "D"),
  label_fontfamily = sansfamily,
  label_size = 8
)
ggsave(
  "5-amino_acid-ribosome.svg",
  fig_5ABCD,
  path = plotdir,
  width = 17.5,
  height = 12.5,
  units = "cm"
)
ggsave(
  "5-amino_acid-ribosome.pdf",
  fig_5ABCD,
  device = cairo_pdf,
  path = plotdir,
  width = 17.5,
  height = 12.5,
  units = "cm"
)