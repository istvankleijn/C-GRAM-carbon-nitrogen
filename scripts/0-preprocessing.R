library(here)
source(here("src", "plot_setup.R"))

plotdir <- here("plots")

df <- read_csv(here("data", "results.csv"))
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
