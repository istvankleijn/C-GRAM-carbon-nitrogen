library(tidyverse)
library(here)
library(cowplot)
library(extrafont)

# Needs to run once, to import fonts into extrafont database
# and register with PDF output device
if (is.null(fonts())) {
  font_import()
  loadfonts()
}

seriffamily <- "Sitka Text"
sansfamily <- "Gill Sans MT"

okabe_grey <- c(
  "#555555", # grey
  "#E69F00", # orange
  "#56B4E9", # cyan
  "#009E73", # green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7" # pink
)
options(
  ggplot2.discrete.colour = okabe_grey,
  ggplot2.discrete.fill = okabe_grey
)


theme_set(
  theme_bw(
    base_size = 10,
    base_family = sansfamily,
    base_line_size = 0.25,
    base_rect_size = 0.5
  ) +
    theme(
      legend.key.size = unit(8, "pt"),
      legend.spacing = unit(0, "pt"),
      legend.position = "bottom",
      legend.direction = "vertical",
      legend.box.margin = margin(),
      legend.text.align = 0
    )
)

palette_proteins <- c(
  "Q"  = "#555555",
  "R"  = "#E69F00",
  "C"  = "#56B4E9",
  "Ef" = "#009E73",
  "Er" = "#F0E442",
  "K"  = "#555555",
  "N"  = "#CC79A7",
  "S"  = "#D55E00",
  "Z"  = "#555555"
)
palette_mass <- c(
  "q" = "#555555", # grey
  "r" = "#E69F00", # orange
  "t_C" = "#56B4E9", # sky blue
  "e_Af" = "#009E73", # green
  "e_Ar" = "#F0E442", # yellow
  "t_N" = "#CC79A7", # pink
  "c" = "#0072B2", # blue
  "a" = "#D55E00", # vermillion
  "n" = "#630a75", # purple
  "k" = "#555555",
  "l" = "#555555",
  "total_metabolite" = "#555555",
  "total_protein" = "#555555"
)
palette_fractions <- c(
  "f_Q" = "#555555",
  "f_R" = "#E69F00",
  "f_C" = "#56B4E9",
  "f_Ef" = "#009E73",
  "f_Er" = "#F0E442",
  "f_Kex" = "#D55E00",
  "f_Kre" = "#630a75",
  "f_N" = "#CC79A7",
  "f_Z" = "#555555"
)
palette_experiment <- c(
  "sweep_kC" = "#56B4E9",
  "sweep_kN" = "#CC79A7",
  "k_Kre=10.0_k_N=20.0" = "#D55E00",
  "k_Kre=5.0_k_N=20.0" = "#0072B2"
)
palette_zeros <- c(
  "C" = "#555555",
  "Kex" = "#D55E00",
  "Kre" = "#0072B2",
  "Kex,Kre" = "#CC79A7",
  "C,Kex" = "#E69F00"
)
palette_fermresp <- c(
  "fermentation" = "#555555", # grey
  "respiration"  = "#E69F00" # orange
)
palette_ketoacid <- c(
  "none"      = "#555555", # grey
  "excretion" = "#0072B2", # blue
  "recycling" = "#D55E00", # vermillion
  "mixed"     = "#009E73" # green
)

shapes_fermresp <- c(
  "fermentation" = 1, # circle
  "respiration" = 0 # square
)
