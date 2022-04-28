
fit_simple_model <- function(W,
                             B_fixed_at_zero = FALSE,
                             reweight = FALSE){

  n <- nrow(W)
  J <- ncol(W)

  if(B_fixed_at_zero){
    B_fixed_indices <-  matrix(TRUE,nrow = 1, ncol = J)
  }else{
    B_fixed_indices <-  matrix(
      c(rep(FALSE, J - 1),TRUE),nrow = 1, ncol = J)
  }

  if(reweight){
    criterion <- "reweighted_Poisson"
  } else{
    criterion <- "Poisson"
  }



fitted_model <- estimate_parameters(W = W,
                      X = matrix(1, nrow = n,ncol = 1),
                      Z = matrix(1, nrow = n, ncol = 1),
                      P = matrix(1/J,nrow = 1, ncol = J),
                      P_fixed_indices = matrix(TRUE, nrow = 1, ncol = J),
                      X_tilde = matrix(0,nrow = 1, ncol = 1),
                      Z_tilde = matrix(0,nrow = n, ncol = 1),
                      Z_tilde_gamma_cols = 1,
                      P_tilde = matrix(1/J,nrow = 1, ncol = J),
                      P_tilde_fixed_indices = matrix(TRUE, nrow = 1, ncol = J),
                      gammas = apply(W,1,function(x) log(sum(x))),
                      gammas_fixed_indices = rep(FALSE, n),
                      B = matrix(0, ncol = J, nrow = 1),
                      B_fixed_indices = B_fixed_indices,
                      gamma_tilde = matrix(0,ncol = 1, nrow = 1),
                      gamma_tilde_fixed_indices = matrix(TRUE,
                                                         ncol = 1, nrow = 1),
                      alpha_tilde = NULL,
                      Z_tilde_list = NULL,
                      barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                      barrier_scale = 10, #increments for value of barrier penalty
                      max_barrier = 1e12, #maximum value of barrier_t
                      initial_conv_tol = 1000,
                      final_conv_tol = 0.1,
                      final_f = 1e-6,
                      constraint_tolerance = 1e-10,
                      hessian_regularization = 0.01,
                      criterion = criterion,
                      subproblem_method = "Newton",
                      profile_P = TRUE,
                      profiling_maxit = 25,
                      wts = NULL,
                      verbose = FALSE,
                      bootstrap_failure_cutoff = NULL)

return(fitted_model)

}
