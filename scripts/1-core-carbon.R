library(here)
source(here("src", "plot_setup.R"))

modelname <- "ketoacid-excretion"
date <- "2022-06-16"
datadir <- here("data", modelname, date)
plotdir <- here("plots")

df <- read_csv(file.path(datadir, "results.csv"))
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

fig_1B <- df %>%
	filter(
		experiment == "sweep_kC_nores"
	) %>%
	ggplot(aes(x = k_C, y = mu)) +
	geom_point(colour = "#555555", size = 1) +
	scale_x_log10(
	  n.breaks = 6,
	  labels = scales::label_log(base = 10, digits = 1)
	) +
	coord_cartesian(
		xlim = c(10^-2.2, 10^3.2),
		ylim = c(0, 3.0),
		expand = FALSE
	) +
	labs(
		x = expression("Carbon uptake rate" ~ italic(k)[C] ~ (h^{
			-1
		})),
		y = expression("Growth rate" ~ italic(mu) ~ (h^{
			-1
		}))
	)
fig_1B

labels_1C <- tribble(
	~fraction,  ~label,  ~colour,
	"f_C",      "paste('carbon uptake ', E[C])",  "#56B4E9",
	"f_N",      "paste('nitrogen uptake ', E[N])",  "#CC79A7",
	"f_Z",      "paste('housekeeping ', Z)", "#555555",
	"f_R",      "paste('ribosome ', R)",     "#E69F00",
	"f_Ef",     "paste('amino acid synthesis ', E[A])",  "#009E73",
)
fig_1C <- df_alloc %>%
	filter(
		experiment == "sweep_kC_nores"
	) %>%
	left_join(labels_1C, by = "fraction") %>%
	drop_na(label) %>%
	ggplot(aes(x = mu, y = allocation)) +
	geom_point(aes(colour = fraction), size = 1) +
	scale_colour_manual(
	  values = palette_fractions,
	  labels = c(
	    "f_C" = expression("carbon uptake" ~ E[C]),
	    "f_N" = expression("nitrogen uptake" ~ E[N]),
	    "f_Z" = expression("housekeeping" ~ Z),
	    "f_R" = expression("ribosome" ~ R),
	    "f_Ef" = expression("amino acid synthesis" ~ E[A])
	  ),
	  breaks = c("f_R", "f_Ef", "f_Z", "f_N", "f_C")
	) +
  guides(
    colour = guide_legend(title.position = "left")
  ) +
	coord_cartesian(xlim = c(0, 3.0), ylim = c(0, 0.8), expand = FALSE) +
	labs(
		x = expression("Growth rate" ~ italic(mu) ~ (h^{-1})),
		y = "Allocation fraction",
		colour = "Protein"
	)
fig_1C

fig_1D <- df_mass %>%
	filter(
		experiment == "sweep_kC_nores"
	) %>%
	filter(state_var %in% c("c", "a", "n", "total_metabolite")) %>%
	ggplot(aes(x = mu, y = mass_fraction)) +
	geom_point(aes(colour = state_var), size = 1) +
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
	coord_cartesian(xlim = c(0, 3.0), ylim = c(0, 0.45), expand = FALSE) +
	labs(
		x = expression("Growth rate" ~ italic(mu) ~ (h^{
			-1
		})),
		y = "Mass fraction",
		colour = "Metabolite",
		size = "Metabolite"
	)
fig_1D

fig_1 <- plot_grid(
  NULL, fig_1B,
  fig_1C, fig_1D,
  nrow = 2,
  rel_widths = c(7, 5),
  rel_heights = c(4, 6),
  labels = "AUTO",
  label_size = 10,
  label_fontfamily = sansfamily
)
fig_1

fig_1AB <- plot_grid(
  NULL, fig_1B,
  nrow = 1,
  rel_widths = c(7, 5),
  labels = c("A", "B"),
  label_size = 10,
  label_fontfamily = sansfamily
)
fig_1CD <- plot_grid(
  fig_1C, fig_1D,
  nrow = 1,
  rel_widths = c(6, 6),
  align = "h",
  labels = c("C", "D"),
  label_size = 10,
  label_fontfamily = sansfamily
)
fig_1CD

fig_1 <- plot_grid(
  fig_1AB, fig_1CD,
  nrow = 2,
  rel_heights = c(4, 6)
)

ggsave(
	"1-sweep_kC_nores.pdf",
	fig_1,
	path = plotdir,
	device = cairo_pdf,
	width = 12,
	height = 10,
	units = "cm"
)
ggsave(
	"1-sweep_kC_nores.svg",
	fig_1,
	path = plotdir,
	width = 12,
	height = 10,
	units = "cm"
)
