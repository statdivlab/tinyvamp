test_that("get_A_tilde_k_list works when  P_tilde entries fixed", {

  X <- matrix(1)
  X_tilde <- matrix(1)
  Z <- matrix(1)
  Z_tilde <- matrix(1)
  Z_tilde_gamma_cols <- 1

  P <- matrix(rep(1/8, 8), nrow = 1)
  P_fixed_indices <- matrix(rep(FALSE,8), nrow = 1)
  P_tilde <- matrix((1:8)/sum(1:8), nrow = 1)
  P_tilde_fixed_indices <- matrix(rep(FALSE,8), nrow = 1)
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

  varying_lr_df <- ra_to_lr(dfs$varying)

  expect_equal(get_A_tilde_k_list(fixed_df = dfs$fixed,
                     varying_df = dfs$varying,
                     varying_lr_df = varying_lr_df),
               list(diag(8)))

})

test_that("get_A_tilde_k_list works when some P_tilde entries fixed", {

  X <- matrix(1)
  X_tilde <- matrix(1)
  Z <- matrix(1)
  Z_tilde <- matrix(1)
  Z_tilde_gamma_cols <- 1

  P <- matrix(rep(1/8, 8), nrow = 1)
  P_fixed_indices <- matrix(rep(FALSE,8), nrow = 1)
  P_tilde <- matrix((1:8)/sum(1:8), nrow = 1)
  P_tilde_fixed_indices <- matrix(c(FALSE, TRUE, TRUE, rep(FALSE,5)), nrow = 1)
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

  varying_lr_df <- ra_to_lr(dfs$varying)

  expect_equal(get_A_tilde_k_list(fixed_df = dfs$fixed,
                                  varying_df = dfs$varying,
                                  varying_lr_df = varying_lr_df),
               list(rbind(diag(6)[1,],
                          0*diag(6)[1:2,],
                          diag(6)[2:6,])))

})
