
do_one_boot <- function(W,
                        fitted_model,
                        m = NULL,
                        seed = NULL,
                        boot_weights = NULL){
  
  
  n <- nrow(W)
  if(is.null(m)){
    m <- sqrt(n)
  }
  
  
  if(is.null(boot_weights)){
    boot_weights <- rgamma(n,m/n)
    boot_weights <- boot_weights/sum(boot_weights)
    boot_weights <- rep(boot_weights, each = J)
  }
  
  if(!is.null(fitted_model$Z_tilde_list)){
    fitted_model$Z_tilde <- NULL
  }
  
  boot_model <- estimate_parameters(W = W,
                                    X = fitted_model$X,
                                    Z = fitted_model$Z,
                                    Z_tilde = fitted_model$Z_tilde,
                                    Z_tilde_gamma_cols =
                                      fitted_model$Z_tilde_gamma_cols,
                                    gammas =
                                      fitted_model$gammas,
                                    gammas_fixed_indices =
                                      fitted_model$gammas_fixed_indices,
                                    P = fitted_model$P,
                                    P_fixed_indices =
                                      fitted_model$P_fixed_indices,
                                    B = fitted_model$B,
                                    B_fixed_indices =
                                      fitted_model$B_fixed_indices,
                                    X_tilde = fitted_model$X_tilde,
                                    P_tilde = fitted_model$P_tilde,
                                    P_tilde_fixed_indices =
                                      fitted_model$P_tilde_fixed_indices,
                                    gamma_tilde = fitted_model$gamma_tilde,
                                    gamma_tilde_fixed_indices =
                                      fitted_model$gamma_tilde_fixed_indices,
                                    alpha_tilde =
                                      fitted_model$alpha_tilde,
                                    Z_tilde_list =
                                      fitted_model$Z_tilde_list,
                                    barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                                    barrier_scale = 10, #increments for value of barrier penalty
                                    max_barrier = 1e12, #maximum value of barrier_t
                                    initial_conv_tol = 1000,
                                    final_conv_tol = 0.1,
                                    constraint_tolerance = 1e-10,
                                    hessian_regularization = 0.01,
                                    criterion = "Poisson",
                                    profile_P = TRUE,
                                    wts = boot_weights,
                                    profiling_maxit = 25)
  
  return(boot_model)
  
}
