
do_one_lrt <- function(W0,
                       full_model,
                       null_model,#null model specification
                       boot_method = "bayesian_subsample",
                       m = NULL,
                       seed = NULL,
                       boot_weights = NULL,
                       return_models = FALSE){


  n <- nrow(W0)
  if(is.null(m)){
    m <- sqrt(n)
  }


  if(is.null(boot_weights)){
    stop("Bootstrapping weights boot_weights must be provided.")
  }

  boot_full <- estimate_parameters(W = W0,
                                   X = full_model$X,
                                   Z = full_model$Z,
                                   Z_tilde = full_model$Z_tilde,
                                   Z_tilde_gamma_cols =
                                     full_model$Z_tilde_gamma_cols,
                                   gammas =
                                     full_model$gammas,
                                   gammas_fixed_indices =
                                     full_model$gammas_fixed_indices,
                                   P = full_model$P,
                                   P_fixed_indices =
                                     full_model$P_fixed_indices,
                                   B = full_model$B,
                                   B_fixed_indices =
                                     full_model$B_fixed_indices,
                                   X_tilde = full_model$X_tilde,
                                   P_tilde = full_model$P_tilde,
                                   P_tilde_fixed_indices =
                                     full_model$P_tilde_fixed_indices,
                                   gamma_tilde = full_model$gamma_tilde,
                                   gamma_tilde_fixed_indices =
                                     full_model$gamma_tilde_fixed_indices,
                                   alpha_tilde =
                                     full_model$alpha_tilde,
                                   Z_tilde_list =
                                     full_model$Z_tilde_list,
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
                                   verbose = TRUE,
                                   profiling_maxit = 25)

  print(paste("full model weights sum to "), sum(boot_weights),
        sep = "", collapse = "")

  # if(boot_full$criterion == "reweighted_Poisson"){
  #   null_model$criterion <- "Poisson"
  #   boot_weights <- boot_full$weights
  # }

  boot_null <- estimate_parameters(W = W0,
                                   X = null_model$X,
                                   Z = null_model$Z,
                                   Z_tilde = null_model$Z_tilde,
                                   Z_tilde_gamma_cols =
                                     null_model$Z_tilde_gamma_cols,
                                   gammas =
                                     null_model$gammas,
                                   gammas_fixed_indices =
                                     null_model$gammas_fixed_indices,
                                   P = null_model$P,
                                   P_fixed_indices =
                                     null_model$P_fixed_indices,
                                   B = null_model$B,
                                   B_fixed_indices =
                                     null_model$B_fixed_indices,
                                   X_tilde = null_model$X_tilde,
                                   P_tilde = null_model$P_tilde,
                                   P_tilde_fixed_indices =
                                     null_model$P_tilde_fixed_indices,
                                   gamma_tilde = null_model$gamma_tilde,
                                   gamma_tilde_fixed_indices =
                                     null_model$gamma_tilde_fixed_indices,
                                   alpha_tilde =
                                     null_model$alpha_tilde,
                                   Z_tilde_list =
                                     null_model$Z_tilde_list,
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
                                   verbose = TRUE,
                                   profiling_maxit = 25)

  print(paste("null model weights sum to "), sum(boot_weights),
        sep = "", collapse = "")


  lr_stat <- 2*(boot_null$objective -
                                            boot_full$objective)


  if(!return_models){
  return(lr_stat)
  } else{
      return(list("lr_stat" = lr_stat,
                  "full_model" = boot_full_model,
                  "null_model" = boot_null_model))
    }

}
