test_that("P derivatives work in simple case", {

  X <- matrix(1)
  X_tilde <- matrix(1)
  Z <- matrix(1)
  Z_tilde <- matrix(1)
  Z_tilde_gamma_cols <- 1

  P <- matrix(rep(1/8, 8), nrow = 1)
  P_fixed_indices <- matrix(rep(FALSE,8), nrow = 1)
  P_tilde <- matrix((1:8)/sum(1:8), nrow = 1)
  P_tilde_fixed_indices <- matrix(rep(TRUE,8), nrow = 1)
  B <- matrix(c(-3:3,0),nrow = 1)
  B_fixed_indices <- matrix(rep(TRUE,8), nrow = 1)
  gammas <- 8
  gammas_fixed_indices <- TRUE
  gamma_tilde <- matrix(log(100))
  gamma_tilde_fixed_indices <- matrix(TRUE)

  dfs <- parameters_to_dataframes(P,
                                  P_fixed_indices,
                                  P_tilde,
                                  P_tilde_fixed_indices,
                                  B,
                                  B_fixed_indices,
                                  gammas,
                                  gammas_fixed_indices,
                                  gamma_tilde,
                                  gamma_tilde_fixed_indices)

  params <- dataframes_to_parameters(dfs$fixed,
                                     dfs$varying)

  fixed_status <- dfs$fixed
  fixed_status$value <- 0
  varying_status <- dfs$varying
  varying_status$value <- 1

  param_status <- dataframes_to_parameters(fixed_status,
                                           varying_status)

  function_output <- par_to_jacobian_row(params,
                      param_status,
                      i = 1,
                      j = 1,
                      X,
                      Z,
                      X_tilde,
                      Z_tilde,
                      Z_tilde_gamma_cols)
  theoretical_output <- c(exp(gammas + B[1]),rep(0,7))

  expect_equal(function_output, theoretical_output)
})


test_that("More mean derivatives work in simple case", {

  X <- matrix(1)
  X_tilde <- matrix(1)
  Z <- matrix(1)
  Z_tilde <- matrix(1)
  Z_tilde_gamma_cols <- 1

  P <- matrix(rep(1/8, 8), nrow = 1)
  P_fixed_indices <- matrix(rep(FALSE,8), nrow = 1)
  P_tilde <- matrix((1:8)/sum(1:8), nrow = 1)
  P_tilde_fixed_indices <- matrix(rep(TRUE,8), nrow = 1)
  B <- matrix(c(-3:3,0),nrow = 1)
  B_fixed_indices <- matrix(rep(FALSE,8), nrow = 1)
  gammas <- 8
  gammas_fixed_indices <- FALSE
  gamma_tilde <- matrix(log(100))
  gamma_tilde_fixed_indices <- matrix(FALSE)

  dfs <- parameters_to_dataframes(P,
                                  P_fixed_indices,
                                  P_tilde,
                                  P_tilde_fixed_indices,
                                  B,
                                  B_fixed_indices,
                                  gammas,
                                  gammas_fixed_indices,
                                  gamma_tilde,
                                  gamma_tilde_fixed_indices)

  params <- dataframes_to_parameters(dfs$fixed,
                                     dfs$varying)

  fixed_status <- dfs$fixed
  if(nrow(fixed_status)>0){
  fixed_status$value <- 0}
  varying_status <- dfs$varying
  varying_status$value <- 1

  param_status <- dataframes_to_parameters(fixed_status,
                                           varying_status)

  function_output <- par_to_jacobian_row(params,
                                         param_status,
                                         i = 1,
                                         j = 1,
                                         X,
                                         Z,
                                         X_tilde,
                                         Z_tilde,
                                         Z_tilde_gamma_cols)
  theoretical_output <- c(exp(gammas + B[1]),rep(0,7),
                          c(P[,1]*exp(gammas + B[,1]) +
                          P_tilde[,1]*exp(gammas + gamma_tilde + B[,1]),
                          rep(0,6)),
                          P[,1]*exp(gammas + B[,1]) +
                            P_tilde[,1]*exp(gammas + gamma_tilde + B[,1]),
                            P_tilde[,1]*exp(gammas + gamma_tilde + B[,1]))




  expect_equal(function_output, theoretical_output)
})
