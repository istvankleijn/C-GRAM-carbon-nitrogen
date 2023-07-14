library(here)
source(here("src", "plot_setup.R"))

modelname <- "ketoacid-excretion"
date <- "2022-06-21"
datadir <- here("data", modelname, date)
plotdir <- here("plots", modelname, date)

cutoff <- 1e-8

df <- read_csv(file.path(datadir, "results.csv")) %>%
  rename(
    Delta_CN = "Δ_CN",
    gamma_K = "γ_K",
    alpha_Cf = `α_Cf`,  # String marks do not work for \alpha
    alpha_Cr = `α_Cr`,  # but backticks do not work for other Greek characters
  ) %>%
  mutate(
    zero_alloc = case_when(
      f_C <= cutoff & f_Kex <= cutoff & f_Kre <= cutoff ~ "C,Kex,Kre",
      f_C <= cutoff & f_Kex <= cutoff                   ~ "C,Kex",
      f_C <= cutoff &                   f_Kre <= cutoff ~ "C,Kre",
                      f_Kex <= cutoff & f_Kre <= cutoff ~ "Kex,Kre",
                      f_Kex <= cutoff                   ~ "Kex",
                                        f_Kre <= cutoff ~ "Kre",
      f_C <= cutoff                                     ~ "C",
      TRUE                                              ~ "none"
    )
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




labels_3B <- tribble(
  ~experiment,  ~label,  ~colour,
  "k_Kre=10.0_k_N=20.0",  "k[Kre]==10.0",  "#56B4E9",
  "k_Kre=5.0_k_N=20.0",   "k[Kre]==5.0",   "#CC79A7"
)
fig_3B <- df %>%
  filter(
    experiment %in% c("k_Kre=10.0_k_N=20.0", "k_Kre=5.0_k_N=20.0"),
    19.0 < k_Kex, k_Kex < 21.0
  ) %>%
  left_join(labels_3B, by = "experiment") %>%
  ggplot(aes(x = C_per_N, y = mu)) +
  geom_point(aes(colour = label)) +
  scale_colour_manual(
    values = deframe(select(labels_3B, label, colour)),
    labels = scales::parse_format(),
  ) +
  coord_cartesian(
    # xlim = c(0, 1.0),
    ylim = c(0, 2.0),
    expand = FALSE
  ) +
  labs(
    x = "Nutrient C-to-N ratio",
    y = expression("Growth rate"~italic(mu)~(h^{-1})),
    colour = "Recycling\nrate"
  ) +
  theme(
    legend.text.align = 0
  )
fig_3B


labels_3C <- tribble(
  ~fraction,  ~label,  ~colour,
  "f_C",      "E[C]",   "#56B4E9",
  "f_N",      "E[N]",   "#CC79A7",
  # "f_Z",      "Z",      "#555555",
  "f_R",      "R",      "#E69F00",
  "f_Ef",     "E[Af]",  "#009E73",
  # "f_Er",     "E[Ar]",  "#F0E442",
  "f_Kex",    "E[Kex]", "#D55E00",
  "f_Kre",    "E[Kre]", "#630a75",
)
fig_3C <- df_alloc %>%
  filter(
    experiment %in% c("k_Kre=10.0_k_N=20.0", "k_Kre=5.0_k_N=20.0"),
    19.0 < k_Kex, k_Kex < 21.0
  ) %>%
  left_join(labels_3C, by = "fraction") %>%
  left_join(select(labels_3B, experiment, facet = label), by = "experiment") %>%
  drop_na(label) %>%
  ggplot(aes(x = mu, y = allocation)) +
  geom_point(aes(colour = label)) +
  facet_wrap(~facet, ncol = 1, labeller = "label_parsed") +
  scale_colour_manual(
    values = deframe(select(labels_3C, label, colour)),
    labels = scales::parse_format(),
  ) +
  coord_cartesian(xlim = c(1.0, 2.0), ylim = c(0, 0.4), expand = FALSE) +
  labs(
    x = expression("Growth rate"~italic(mu)~(h^{-1})),
    y = "Allocation fraction",
    colour = "Protein"
  ) +
  theme(
    legend.text.align = 0
  )
fig_3C


fig_3D <- df_mass %>%
  filter(
    experiment %in% c("k_Kre=10.0_k_N=20.0", "k_Kre=5.0_k_N=20.0"),
    19.0 < k_Kex, k_Kex < 21.0,
    state_var %in% c("k", "c", "a", "n")
  ) %>%
  left_join(select(labels_3B, experiment, facet = label), by = "experiment") %>%
  ggplot(aes(x = mu, y = mass_fraction)) +
  geom_point(aes(colour = state_var)) +
  facet_wrap(~facet, ncol = 1, labeller = "label_parsed") +
  scale_colour_manual(values = palette_mass) +
  coord_cartesian(xlim = c(1.0, 2.0), ylim = c(0, 0.5), expand = FALSE) +
  labs(
    x = expression("Growth rate"~italic(mu)~(h^{-1})),
    y = "Mass fraction",
    colour = "Metabolite"
  )
fig_3D


fig_3 <- plot_grid(
  NULL, fig_3B,
  fig_3C, fig_3D,
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
  "3-CtoN-allocs-nores.pdf",
  fig_3,
  path = plotdir,
  device = cairo_pdf,
  width = 14,
  height = 14,
  units = "cm"
)
ggsave(
  "3-CtoN-allocs-nores.svg",
  fig_3,
  path = plotdir,
  width = 14,
  height = 14,
  units = "cm"
)


palette_zeros <- c(
  "C" = "#555555",
  "Kex" = "#D55E00",
  "Kre" = "#0072B2",
  "Kex,Kre" = "#CC79A7",
  "C,Kex" = "#E69F00"
)



zeros_kKre <- df %>%
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
    y = expression("Recycling efficiency"~k[Kre]),
    fill = "Zero allocation"
  ) +
  guides(
    fill = guide_legend(title.position = "left")
  ) +
  theme(
    legend.position = "bottom",
    legend.direction = "vertical"
  )
# ggsave(
#   "zeros_kKre.svg",
#   zeros_kKre,
#   path = plotdir,
#   width = 5,
#   height = 10,
#   units = "cm"
# )

zeros_kKex_10 <- df %>%
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
    y = expression("Excretion efficiency"~k[Kex])
  )
# ggsave(
#   "zeros_kKex_10.svg",
#   zeros_kKex_10,
#   path = plotdir,
#   width = 5,
#   height = 5,
#   units = "cm"
# )

zeros_kKex_5 <- df %>%
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
    y = expression("Excretion efficiency"~k[Kex])
  )
# ggsave(
#   "zeros_kKex_5.svg",
#   zeros_kKex_5,
#   path = plotdir,
#   width = 5,
#   height = 5,
#   units = "cm"
# )

zeros_kN_10 <- df %>%
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
    y = expression("Uptake efficiency"~k[N])
  )
# ggsave(
#   "zeros_kN_10.svg",
#   zeros_kN_10,
#   path = plotdir,
#   width = 5,
#   height = 5,
#   units = "cm"
# )

zeros_kN_5 <- df %>%
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
    y = expression("Uptake efficiency"~k[N])
  )
# ggsave(
#   "zeros_kN_5.svg",
#   zeros_kN_5,
#   path = plotdir,
#   width = 5,
#   height = 5,
#   units = "cm"
# )

zeros_BCDE <- plot_grid(
  zeros_kKex_10, zeros_kN_10,
  zeros_kKex_5, zeros_kN_5,
  nrow = 2,
  labels = c("B", "C", "D", "E"),
  label_fontfamily = sansfamily,
  label_size = 8
)
zeros_ABCDE <- plot_grid(
  zeros_kKre, zeros_BCDE,
  ncol = 2,
  rel_widths = c(1, 2),
  labels = c("A", ""),
  label_fontfamily = sansfamily,
  label_size = 8
)
ggsave(
  "4-zeros_rates_gamma.svg",
  zeros_ABCDE,
  path = plotdir,
  width = 15,
  height = 10,
  units = "cm"
)
ggsave(
  "4-zeros_rates_gamma.pdf",
  zeros_ABCDE,
  device = cairo_pdf,
  path = plotdir,
  width = 15,
  height = 10,
  units = "cm"
)
