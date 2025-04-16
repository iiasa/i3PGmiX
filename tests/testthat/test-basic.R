library(i3PGmiX)
library(testthat)

context("Basic Model runs work")

test_that("basic model run", {
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
                      correct_bias = 1,
                      calculate_d13c = 0,
                      canopy_cond=1,
                      gpp_model=2,
                      maint_resp=2),
      check_input = TRUE, df_out = FALSE)

    testthat::expect_true(class(out) == "array")
})


