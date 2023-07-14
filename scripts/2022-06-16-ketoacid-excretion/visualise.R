library(here)
source(here("src", "plot_setup.R"))

modelname <- "ketoacid-excretion"
date <- "2022-06-16"
datadir <- here("data", modelname, date)
plotdir <- here("plots", modelname, date)

old_df <- read_csv(here("data", modelname, "2021-12-15", "results.csv"))

# cutoff <- 1e-6

df <- read_csv(file.path(datadir, "results.csv")) %>%
    rename(
      Delta_CN = "Δ_CN",
      gamma_K = "γ_K",
      alpha_Cf = `α_Cf`,  # String marks do not work for \alpha
      alpha_Cr = `α_Cr`,  # but backticks do not work for other Greek characters
    )
  #   ) %>%
  # mutate(
  #   f_zeros = case_when(
  #     f_C <= cutoff & f_Kex <= cutoff & f_Kre <= cutoff ~ "C,Kex,Kre",
  #     f_C <= cutoff & f_Kex <= cutoff                   ~ "C,Kex",
  #     f_C <= cutoff &                   f_Kre <= cutoff ~ "C,Kre",
  #                     f_Kex <= cutoff & f_Kre <= cutoff ~ "Kex,Kre",
  #                     f_Kex <= cutoff                   ~ "Kex",
  #                                       f_Kre <= cutoff ~ "Kre",
  #     f_C <= cutoff                                     ~ "C",
  #     TRUE                                              ~ "none"
  #   )
  # )
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

fig_1B <- df %>%
  filter(
      experiment == "sweep_kC_nores"
  ) %>%
  ggplot(aes(x = k_C, y = mu)) +
  geom_point(colour = "#555555") +
  scale_x_log10(n.breaks = 6) +
  coord_cartesian(
    xlim = c(10^-2.2, 10^3.2),
    ylim = c(0, 3.0),
    expand = FALSE
  ) +
  labs(
    x = expression("Carbon uptake rate"~italic(k)[C]~(h^{-1})),
    y = expression("Growth rate"~italic(mu)~(h^{-1}))
  )
fig_1B

labels_1C <- tribble(
  ~fraction,  ~label,  ~colour,
  "f_C",      "E[C]",  "#56B4E9",
  "f_N",      "E[N]",  "#CC79A7",
  "f_Z",      "Z",     "#555555",
  "f_R",      "R",     "#E69F00",
  "f_Ef",     "E[A]",  "#009E73",
)
fig_1C <- df_alloc %>%
  filter(
      experiment == "sweep_kC_nores"
  ) %>%
  left_join(labels_1C, by = "fraction") %>%
  drop_na(label) %>%
  ggplot(aes(x = mu, y = allocation)) +
  geom_point(aes(colour = label)) +
  scale_colour_manual(
    values = deframe(select(labels_1C, label, colour)),
    labels = scales::parse_format(),
  ) +
  coord_cartesian(xlim = c(0, 3.0), ylim = c(0, 0.8), expand = FALSE) +
  labs(
    x = expression("Growth rate"~italic(mu)~(h^{-1})),
    y = "Allocation fraction",
    colour = "Protein"
  ) +
  theme(
    legend.text.align = 0
  )
fig_1C
# model <-  lm(
#   f_R ~ a,
#   data = df
# )
# label <- substitute(
#   italic(f)[R] == delta %*%~~(italic(a) - asat),
#   list(
#     delta = format(unname(coef(model)["a"]), digits = 3),
#     asat = format(
#       -unname(coef(model)["(Intercept)"]/coef(model)["a"]),
#       digits = 3
#     )
#   )
# )
# fig_C <- df %>%
#   filter(
#     experiment == "sweep_kC"
#   ) %>%
#   modelr::add_predictions(model) %>%
#   ggplot(aes(x = a)) +
#   geom_point(aes(y = f_R, colour = k_C), shape = 16) +
#   geom_line(aes(y = pred), colour = "black") +
#   scale_colour_continuous(
#     trans = "log10",
#     limits = c(1e-2, 1e3),
#     n.breaks = 6
#   ) +
#   coord_cartesian(xlim = c(0, 0.2), ylim = c(0, 0.6), expand = FALSE) +
#   annotate(
#     "text",
#     label = label,
#     x = 0.02, y = 0.52,
#     size = 8 / .pt,
#     family = seriffamily,
#     hjust = 0
#   ) +
#   labs(
#     x = "Amino acid mass fraction",
#     y = "Ribosomal allocation fraction",
#     colour = "Transporter rate"
#   )

fig_1D <- df_mass %>%
  filter(
      experiment == "sweep_kC_nores"
  ) %>%
  filter(state_var %in% c("c", "a", "n")) %>%
  ggplot(aes(x = mu, y = mass_fraction)) +
  geom_point(aes(colour = state_var)) +
  scale_colour_manual(values = palette_mass) +
  coord_cartesian(xlim = c(0, 3.0), ylim = c(0, 0.2), expand = FALSE) +
  labs(
    x = expression("Growth rate"~italic(mu)~(h^{-1})),
    y = "Mass fraction",
    colour = "Metabolite"
  )
fig_1D


fig_1 <- plot_grid(
  NULL, fig_1B,
  fig_1C, fig_1D,
  nrow = 2,
  rel_widths = c(7, 7),
  align = "v",
  axis = "lr",
  labels = "AUTO",
  label_size = 10,
  label_fontfamily = sansfamily
)
ggsave(
  "1-sweep_kC_nores.pdf",
  fig_1,
  path = plotdir,
  device = cairo_pdf,
  width = 14,
  height = 10,
  units = "cm"
)
ggsave(
  "1-sweep_kC_nores.svg",
  fig_1,
  path = plotdir,
  width = 14,
  height = 10,
  units = "cm"
)

labels_2B <- tribble(
  ~experiment,  ~label,  ~colour,
  "sweep_kC",   "k[C]",  "#56B4E9",
  "sweep_kN",   "k[N]",  "#CC79A7"
)
fig_2B <- df %>%
  filter(
      experiment %in% c("sweep_kC", "sweep_kN")
  ) %>%
  mutate(
    varied_rate = case_when(
      experiment == "sweep_kC" ~ k_C,
      experiment == "sweep_kN" ~ k_N,
      TRUE                     ~ NA_real_
    )
  ) %>%
  left_join(labels_2B, by = "experiment") %>%
  ggplot(aes(x = varied_rate, y = mu)) +
  geom_point(aes(colour = label)) +
  scale_colour_manual(
    values = deframe(select(labels_2B, label, colour)),
    labels = scales::parse_format(),
  ) +
  scale_x_log10(n.breaks = 6) +
  coord_cartesian(
    xlim = c(10^-2.2, 10^3.2),
    ylim = c(0, 3.0),
    expand = FALSE
  ) +
  labs(
    x = expression("Uptake rate"~(h^{-1})),
    y = expression("Growth rate"~italic(mu)~(h^{-1})),
    colour = "Rate\nvaried"
  )
fig_2B

labels_2C <- tribble(
  ~fraction,  ~label,  ~colour,
  "f_C",      "E[C]",   "#56B4E9",
  "f_N",      "E[N]",   "#CC79A7",
  # "f_Z",      "Z",      "#555555",
  # "f_R",      "R",      "#E69F00",
  "f_Ef",     "E[Af]",  "#009E73",
  "f_Er",     "E[Ar]",  "#F0E442",
)
fig_2C <- df_alloc %>%
  filter(
      experiment %in% c("sweep_kC", "sweep_kN")
  ) %>%
  mutate(
    facet = case_when(
      experiment == "sweep_kC" ~ "italic(k)[C]~varied",
      experiment == "sweep_kN" ~ "italic(k)[N]~varied"
    )
  ) %>%
  left_join(labels_2C, by = "fraction") %>%
  drop_na(label) %>%
  ggplot(aes(x = mu, y = allocation)) +
  geom_point(aes(colour = label)) +
  facet_wrap(~facet, ncol = 1, labeller = "label_parsed") +
  scale_colour_manual(
    values = deframe(select(labels_2C, label, colour)),
    labels = scales::parse_format(),
  ) +
  coord_cartesian(xlim = c(0, 3.0), ylim = c(0, 0.8), expand = FALSE) +
  labs(
    x = expression("Growth rate"~italic(mu)~(h^{-1})),
    y = "Allocation fraction",
    colour = "Protein"
  ) +
  theme(
    legend.text.align = 0
  )
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
  filter(state_var %in% c("c", "a", "n")) %>%
  ggplot(aes(x = mu, y = mass_fraction)) +
  geom_point(aes(colour = state_var)) +
  facet_wrap(~facet, ncol = 1, labeller = "label_parsed") +
  scale_colour_manual(values = palette_mass) +
  coord_cartesian(xlim = c(0, 3.0), ylim = c(0, 0.2), expand = FALSE) +
  labs(
    x = expression("Growth rate"~italic(mu)~(h^{-1})),
    y = "Mass fraction",
    colour = "Metabolite"
  )
fig_2D


fig_2 <- plot_grid(
  NULL, fig_2B,
  fig_2C, fig_2D,
  nrow = 2,
  rel_widths = c(7, 7),
  rel_heights = c(5, 9),
  align = "v",
  axis = "lr",
  labels = "AUTO",
  label_size = 10,
  label_fontfamily = sansfamily
)
ggsave(
  "2-sweep_kCkN_fermres.pdf",
  fig_2,
  path = plotdir,
  device = cairo_pdf,
  width = 14,
  height = 14,
  units = "cm"
)
ggsave(
  "2-sweep_kCkN_fermres.svg",
  fig_2,
  path = plotdir,
  width = 14,
  height = 14,
  units = "cm"
)




labels_5A <- tribble(
  ~fraction,  ~label,  ~colour,
  "f_C",      "E[C]",   "#56B4E9",
  "f_N",      "E[N]",   "#CC79A7",
  # "f_Z",      "Z",      "#555555",
  "f_R",      "R",      "#E69F00",
  "f_Ef",     "E[Af]",  "#009E73",
  "f_Er",     "E[Ar]",  "#F0E442",
)
fig_5A <- df_alloc %>%
  filter(
      experiment == "random_kCkEfkN"
  ) %>%
  left_join(labels_5A, by = "fraction") %>%
  drop_na(label) %>%
  ggplot(aes(x = mu, y = allocation)) +
  geom_point(aes(colour = label), alpha = 0.3) +
  scale_colour_manual(
    values = deframe(select(labels_5A, label, colour)),
    labels = scales::parse_format(),
  ) +
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 0.8), expand = FALSE) +
  labs(
    x = expression("Growth rate"~italic(mu)~(h^{-1})),
    y = "Allocation fraction",
    colour = "Protein"
  ) +
  theme(
    legend.text.align = 0
  )
fig_5A

fig_5B <- df_mass %>%
  filter(
      experiment == "random_kCkEfkN"
  ) %>%
  filter(state_var %in% c("c", "a", "n")) %>%
  ggplot(aes(x = mu, y = mass_fraction)) +
  geom_point(aes(colour = state_var), alpha = 0.3) +
  scale_colour_manual(values = palette_mass) +
  coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 0.15), expand = FALSE) +
  labs(
    x = expression("Growth rate"~italic(mu)~(h^{-1})),
    y = "Mass fraction",
    colour = "Metabolite"
  )
fig_5B


# model <-  lm(
#   f_R ~ a,
#   data = df
# )
# label <- substitute(
#   italic(f)[R] == delta %*%~~(italic(a) - asat),
#   list(
#     delta = format(unname(coef(model)["a"]), digits = 3),
#     asat = format(
#       -unname(coef(model)["(Intercept)"]/coef(model)["a"]),
#       digits = 3
#     )
#   )
# )
delta <- 5.0
asat <- 0.0167
label <- substitute(
  italic(f)[R] == delta %*% (italic(a) - asat),
  list(
    delta = format(delta, digits = 3),
    asat = format(asat, digits = 3)
  )
)
fig_5C <- df %>%
  filter(
      experiment == "random_kCkEfkN"
  ) %>%
  mutate(
    mean_k = (k_C*k_Ef*k_N)^(1/3)
  ) %>%
  ggplot(aes(x = a)) +
  geom_point(aes(y = f_R, colour = mean_k), alpha = 0.5) +
  # geom_line(aes(y = pred), colour = "black") +
  geom_abline(slope = delta, intercept = -delta*asat, colour = "#555555") +
  scale_colour_continuous(limits = c(0, 20)) +
  coord_cartesian(xlim = c(0, 0.12), ylim = c(0, 0.5), expand = FALSE) +
  annotate(
    "text",
    label = label,
    colour = "#555555",
    x = 0.01, y = 0.44,
    size = 8 / .pt,
    family = seriffamily,
    hjust = 0
  ) +
  labs(
    x = "Amino acid mass fraction",
    y = "Ribosomal allocation fraction",
    colour = "Geometric\nmean of rates"
  )

fig_5 <- plot_grid(
  fig_5A, fig_5B, 
  fig_5C, NULL,
  nrow = 2,
  rel_widths = c(7, 7),
  align = "v",
  axis = "lr",
  labels = "AUTO",
  label_size = 10,
  label_fontfamily = sansfamily
)
ggsave(
  "5-random_kCkEfkN.pdf",
  fig_5,
  path = plotdir,
  device = cairo_pdf,
  width = 14,
  height = 10,
  units = "cm"
)
ggsave(
  "5-random_kCkEfkN.svg",
  fig_5,
  path = plotdir,
  width = 14,
  height = 10,
  units = "cm"
)













# fig_A <- df_alloc %>%
#   filter(
#       experiment == "random_kCkEfkN"
#   ) %>%
#   filter(!(protein %in% c("Er", "Q", "K"))) %>%
#   ggplot(aes(x = mu, y = allocation)) +
#   geom_point(aes(colour = protein), alpha = 0.5) +
#   scale_colour_manual(values = palette_proteins) +
#   coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 0.8), expand = FALSE) +
#   labs(
#     x = "Growth rate",
#     y = "Allocation fraction",
#     colour = "Protein"
#   )


# fig_B <- df_mass %>%
#   filter(
#       experiment == "random_kCkEfkN"
#   ) %>%
#   filter(state_var %in% c("c", "a", "n")) %>%
#   ggplot(aes(x = mu, y = mass_fraction)) +
#   geom_point(aes(colour = state_var), alpha = 0.5) +
#   scale_colour_manual(values = palette_mass) +
#   coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 0.32), expand = FALSE) +
#   labs(
#     x = "Growth rate",
#     y = "Mass fraction",
#     colour = "Metabolite"
#   )

# # model <-  lm(
# #   f_R ~ a,
# #   data = df
# # )
# # label <- substitute(
# #   italic(f)[R] == delta %*%~~(italic(a) - asat),
# #   list(
# #     delta = format(unname(coef(model)["a"]), digits = 3),
# #     asat = format(
# #       -unname(coef(model)["(Intercept)"]/coef(model)["a"]),
# #       digits = 3
# #     )
# #   )
# # )
# delta <- 5.0
# asat <- 0.0167
# label <- substitute(
#   italic(f)[R] == delta %*%~~(italic(a) - asat),
#   list(
#     delta = format(delta, digits = 3),
#     asat = format(asat, digits = 3)
#   )
# )
# fig_C <- df %>%
#   filter(
#       experiment == "random_kCkEfkN"
#   ) %>%
#   mutate(
#     mean_k = (k_C*k_Ef*k_N)^(1/3)
#   ) %>%
#   ggplot(aes(x = a)) +
#   geom_point(aes(y = f_R, colour = mean_k), alpha = 0.5) +
#   # geom_line(aes(y = pred), colour = "black") +
#   geom_abline(slope = delta, intercept = -delta*asat, colour = "#555555") +
#   scale_colour_continuous(limits = c(0, 20)) +
#   coord_cartesian(xlim = c(0, 0.12), ylim = c(0, 0.5), expand = FALSE) +
#   annotate(
#     "text",
#     label = label,
#     x = 0.01, y = 0.44,
#     size = 8 / .pt,
#     family = seriffamily,
#     hjust = 0
#   ) +
#   labs(
#     x = "Amino acid mass fraction",
#     y = "Ribosomal allocation fraction",
#     colour = "Geometric\nmean of rates"
#   )

# fig <- plot_grid(
#   fig_A, fig_B, fig_C,
#   nrow = 1,
#   labels = "AUTO",
#   label_size = 10,
#   label_fontfamily = sansfamily
# )
# ggsave(
#   "random_kCkEfkN.pdf",
#   fig,
#   path = plotdir,
#   device = cairo_pdf,
#   width = 18,
#   height = 5,
#   units = "cm"
# )
