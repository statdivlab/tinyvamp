test_that("log penalty is correct", {
  require(logsum)

  P <- matrix(1/7,ncol = 7, nrow = 1)
  P_tilde <- matrix(1/7,ncol = 7, nrow = 1)

  param_dfs <- parameters_to_dataframes(P = P,
                           P_fixed_indices = matrix(FALSE,ncol = 7, nrow = 1),
                           P_tilde = P_tilde,
                           P_tilde_fixed_indices = matrix(FALSE,ncol = 7, nrow = 1),
                           B = matrix(0, ncol = 7, nrow = 1),
                           B_fixed_indices = matrix(FALSE,ncol = 7, nrow = 1),
                           gammas = 0,
                           gammas_fixed_indices = FALSE,
                           gamma_tilde = matrix(0,ncol = 1, nrow = 1),
                           gamma_tilde_fixed_indices = FALSE
                           )

  varying_lr_df <- ra_to_lr(param_dfs$varying)

  penalty_from_function <- calculate_log_penalty(varying_lr_df = varying_lr_df,
                        fixed_df = param_dfs$fixed,
                        barrier_t = 1)

  direct_penalty <- -log(1/7)*14

  expect_equal(penalty_from_function, direct_penalty)
})
