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

plotdir <- here("plots")

cutoff <- 1e-8

df <- read_csv(
  here("data", "results.csv"),
  col_types = cols(
    k = col_double(),
    c = col_double(),
    a = col_double(),
    n = col_double(),
    e_Kre = col_double(),
    e_Kex = col_double(),
    e_C = col_double(),
    e_Af = col_double(),
    e_Ar = col_double(),
    e_N = col_double(),
    r = col_double(),
    z = col_double(),
    k_Kre = col_double(),
    k_Kex = col_double(),
    k_C = col_double(),
    k_Ef = col_double(),
    k_Er = col_double(),
    k_N = col_double(),
    k_R = col_double(),
    k_sat = col_double(),
    c_sat = col_double(),
    a_sat = col_double(),
    n_sat = col_double(),
    Δ_CN = col_double(),
    f_Kre = col_double(),
    f_Kex = col_double(),
    f_C = col_double(),
    f_Ef = col_double(),
    f_Er = col_double(),
    f_N = col_double(),
    f_R = col_double(),
    f_Z = col_double(),
    γ_K = col_character(),
    α_Cr = col_character(),
    α_Cf = col_character(),
    mu = col_double(),
    C_per_N = col_character(),
    experiment = col_character()
  )
) %>%
  mutate(
    C_per_N = as.double(str_remove_all(C_per_N, "//1")),
    zero_alloc = case_when(
      f_C <= cutoff & f_Kex <= cutoff & f_Kre <= cutoff ~ "C,Kex,Kre",
      f_C <= cutoff & f_Kex <= cutoff ~ "C,Kex",
      f_C <= cutoff & f_Kre <= cutoff ~ "C,Kre",
      f_Kex <= cutoff & f_Kre <= cutoff ~ "Kex,Kre",
      f_Kex <= cutoff ~ "Kex",
      f_Kre <= cutoff ~ "Kre",
      f_C <= cutoff ~ "C",
      TRUE ~ "none"
    ),
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
  mutate(
    total_metabolite = k + c + a + n,
    total_protein = e_Kre + e_Kex + e_C + e_Af + e_Ar + e_N + r + z
  ) %>%
  pivot_longer(
    cols = all_of(c(state_vars, "total_metabolite", "total_protein")),
    names_to = "state_var",
    values_to = "mass_fraction"
  ) %>%
  mutate(
    var_type = if_else(
      state_var %in% c("k", "c", "a", "n", "total_metabolite"),
      "metabolite",
      "protein"
    )
  )

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
      legend.justification = c(1, 1),
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
