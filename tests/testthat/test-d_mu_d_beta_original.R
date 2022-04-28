# test_that("beta ", {
#   load("data/karstens_beta_problems.rdata")
#
#   # mean_function <- function(x,i,j){
#     temp <- varying_lr_df
#     temp <- lr_to_ra(fixed_df,varying_lr_df,
#                      varying_df)
#     temp$value[temp$param == "B"] <- x
#     temp_params <- dataframes_to_parameters(fixed_df,
#                                             temp)
#
#     means <- meaninate(gammas = temp_params$gammas,
#                        B = temp_params$B,
#                        X = X,
#                        Z = Z,
#                        P = temp_params$P,
#                        X_tilde = X_tilde,
#                        Z_tilde = Z_tilde,
#                        Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                        P_tilde = temp_params$P_tilde,
#                        gamma_tilde = temp_params$gamma_tilde)
#     return(means[i,j])
#   }
#
#   numeric_deriv <- numDeriv::grad(func = function(x) mean_function(x, 1,247),
#                               varying_lr_df$value[varying_lr_df$param == "B"])[7]
#
#   varying_df <- lr_to_ra(fixed_df,varying_lr_df,
#                    varying_df)
#   params <- dataframes_to_parameters(fixed_df,varying_df)
#   analytic_deriv <- mu_d_beta(i = 1,
#                                j = 247,
#                                q = 1,
#                                gammas = params$gammas,
#                                B = params$B,
#                                X = X,
#                                Z = Z,
#                                P = params$P,
#                                X_tilde = X_tilde,
#                                Z_tilde = Z_tilde,
#                                Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                                P_tilde = params$P_tilde,
#                                gamma_tilde = params$gamma_tilde)
#
#   expect_equal(numeric_deriv,as.numeric(analytic_deriv))
# })
