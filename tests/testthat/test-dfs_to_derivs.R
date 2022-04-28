#
#
# test_that("We get output at all",
#           {
#
#             set.seed(0)
#             W <- matrix(sapply(1:10,function(x) rpois(1,1000)),
#                         ncol = 5)
#             X <- matrix(1,ncol = 1, nrow = 2)
#             Z <- matrix(1,nrow = 2, ncol = 1)
#             Z_tilde <- matrix(0,nrow = 2, ncol = 1)
#             Z_tilde_gamma_cols <- 1
#             gammas <- apply(W,1,function(x) log(sum(x)))
#             gammas_fixed_indices <- rep(F,2)
#             P <- matrix(1/5, nrow = 1, ncol = 5)
#             P_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 5)
#             B <- matrix(0,ncol = 5, nrow = 1)
#             B_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
#             X_tilde <- matrix(0,ncol = 1, nrow = 1)
#             P_tilde <- matrix(1/5,ncol = 5, nrow = 1)
#             P_tilde_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
#             gamma_tilde <- matrix(0,nrow = 1, ncol = 1)
#             gamma_tilde_fixed_indices <- matrix(TRUE, nrow = 1, ncol = 1)
#
#
#             parameter_dfs <- parameters_to_dataframes(P,
#                                                       P_fixed_indices,
#                                                       P_tilde,
#                                                       P_tilde_fixed_indices,
#                                                       B,
#                                                       B_fixed_indices,
#                                                       gammas,
#                                                       gammas_fixed_indices,
#                                                       gamma_tilde,
#                                                       gamma_tilde_fixed_indices)
#
#
#             dfs_to_derivs(varying_df = parameter_dfs$varying,
#                                       varying_lr_df = NULL,
#                                       fixed_df = parameter_dfs$fixed_df,
#                                       X = X,
#                                       Z = Z,
#                                       X_tilde = X_tilde,
#                                       Z_tilde = Z_tilde,
#                                       Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                                       criterion = "Poisson",
#                                       gmm_inv_wts = NULL)
#
#           })
