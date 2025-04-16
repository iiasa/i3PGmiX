# i3PGmiX
i3PGmiX forest growth model is an enhanced version of 3PGmix, built upon the [r3PG package](https://github.com/trotsiuk/r3PG) (Trotsiuk et al. 2020). i3PGmiX expands the representation of GPP, NPP, soil processes and natural disturbances, improving the climate sensitivity and the representation of carbon and water balances in the model. The model integrates the P-model (Stocker et al. 2020) and the P-hydro models (Joshi et al. 2022) of GPP. The latter models are optimality-based Farquhar–von Caemmerer–Berry (FvBC) models of photosynthesis that compute optimal leaf internal to ambient CO2 ratios to balance the costs of maintaining carboxylation and transpiration. i3PGmiX also includes different soil models, allowing to account for the interactions and feedbacks between soil and vegetation. Specifically, i3PGmiX includes the Yasso20 and ICBM/2N soil carbon models  (Xenakis et al. 2008) are available. 

i3PGmiX can simulate water, carbon and nutrient cycles in vegetation and soil. The model provides estimates of the impacts of climate change, natural disturbances and forest management drivers on forest ecosystems, allowing to derive parameters relevant for management planning and ecosystem services indicators. 

# Usage
The following code snippet illustrates a basic usage of the i3PGmiX package for running a model simulation, setting up the climate data, and defining various model settings with specific simulation flags. 
```r
library(i3PGmiX)
    out <- run_3PG(
      site = d_site,
      species = d_species,
      climate = d_climate,
      thinning = d_thinning,
      parameters = d_parameters,
      size_dist = d_sizeDist,
      settings = list(light_model = 2,
                      transp_model = 1, 
                      phys_model = 1,
                      soil_model=2,
                      method_jmaxlim=3,
                      nut_lim=0,
                      wind_dist=1, 
                      beetle_dist=1,
                      dist_start=1,
                      correct_bias = 1,
                      calculate_d13c = 0,
                      canopy_cond=1,
                      gpp_model=2, 
                      maint_resp=2),
      check_input = TRUE,
      df_out = FALSE)
```
# Installation
For the installation, users are required to install directly from the GitHub repository using the R package devtools. The following commands can be used to install the model in R:
```r
if(!require(devtools)){install.packages(devtools)}
devtools::install_github( "andreyaugustynczik/i3PGmiX", build_vignettes = TRUE )
library(i3PGmiX)
```
