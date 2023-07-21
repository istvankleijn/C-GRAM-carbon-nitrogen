source(here::here("scripts", "0-preprocessing.R"))

datadir <- here("data", "pomballoc")

cultures <- read_csv(file.path(datadir, "TableS1.csv"))
proteome_fractions <- read_csv(file.path(datadir, "TableS4.csv"))
protein_annotations <- read_csv(file.path(datadir, "TableS16.csv"))

fig_6 <- proteome_fractions %>%
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
  drop_na() %>%
  ggplot(aes(x = growth_rate, y = proteome_fraction)) +
  geom_point(aes(colour = coarse_grained_equivalent), alpha = 0.7, shape = 1) +
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
    x = expression("Growth rate" ~ italic(mu) ~ (h^{
      -1
    })),
    y = "Allocation fraction",
    colour = "Coarse-grained\nequivalent"
  ) +
  theme(
    legend.justification = c(1, 1)
  )
fig_6
ggsave(
  "6-compare-experiments.svg",
  fig_6,
  path = plotdir,
  width = 6,
  height = 6,
  units = "cm"
)
ggsave(
  "6-compare-experiments.pdf",
  fig_6,
  path = plotdir,
  width = 6,
  height = 6,
  units = "cm",
  device = cairo_pdf
)

