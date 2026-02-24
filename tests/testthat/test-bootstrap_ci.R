test_that("bootstrap_ci works, and reweighting works", {
  W <- tinyvamp:::simulate_paper_data(n =3,
                                      J = 5,
                                      B_multiplier = 1,
                                      distrib = "NB",
                                      seed = 0)
  
  fitted_model_rw <- tinyvamp:::fit_simulation_model(W,"reweighted_Poisson")
  fitted_model_poi <- tinyvamp:::fit_simulation_model(W,"Poisson")
  
  expect_true(
    mean(abs(fitted_model_poi$varying$value - fitted_model_rw$varying$value))
    >0.01)
  
  cis <- bootstrap_ci(W = W,
                      fitted_model = fitted_model_rw,
                      n_boot = 3,
                      m = NULL,
                      alpha = 0.05,
                      parallelize = FALSE,
                      seed = 3,
                      return_models = FALSE
                      
  )
  
  expect_type(cis,"list")
  
  
})
