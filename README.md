# A coarse-grained resource allocation model of carbon and nitrogen metabolism in unicellular microbes

This repository contains code to simulate microbial growth on carbon and nitrogen using a coarse-grained resource allocation model.
It also includes code to generate the figures for our associated manuscript (Kleijn *et al.*, bioRxiv 2023).

### Folder structure

General Julia code for simulating the models is provided in `src/`.
The Julia code under `scripts/` was used to generate simulation data in `data/ketoacid-excretion/`, which were aggregated manually as `data/results.csv`.
The R code in `scripts/` was used to generate the figures in `plots/`.
For Figure 5B, supplementary data was taken from Kleijn, Martínez-Segura *et al.* (Life Science Alliance 2022), it is included here in `data/pomballoc/`.

### Citation

Kleijn, I.T., Marguerat, S. and Shahrezaei, V. (2023) ‘A coarse-grained resource allocation model of carbon and nitrogen metabolism in unicellular microbes’. bioRxiv, p. 2023.04.04.535571. Available at: https://doi.org/10.1101/2023.04.04.535571.

### Reference

Kleijn, I.T., Martínez-Segura A., et al. (2022) ‘Growth-rate-dependent and nutrient-specific gene expression resource allocation in fission yeast’, Life Science Alliance, 5(5), p. e202101223. Available at: https://doi.org/10.26508/lsa.202101223.




