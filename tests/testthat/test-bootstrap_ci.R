test_that("multiplication works", {
  W <- simulate_paper_data(n =3,
                           J = 5,
                           B_multiplier = 1,
                           distrib = "NB",
                           seed = 0)


  fitted_model <- fit_simulation_model(W,"reweighted_Poisson")

  cis <- bootstrap_ci(W = W,
                      fitted_model = fitted_model,
                      n_boot = 10,
                      m = NULL,
                      alpha = 0.05,
                      parallelize = TRUE,
                      ncores = 5,
                      seed = 3,
                      return_models = FALSE

  )

  expect_type(cis,"list")


})
