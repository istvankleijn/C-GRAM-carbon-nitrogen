library(here)
source(here("src", "plot_setup.R"))

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
