

refit_model <- function(fitted_model,
                        m,
                        seed = NULL,
                        bootstrap_method = "bayesian_subsample",
                        bootstrap_failure_cutoff = -1e4){

  n <- nrow(fitted_model$W)
  J <- ncol(fitted_model$W)

  if(!is.null(seed)){
    set.seed(seed)
  }

  if(!(boot_method %in% c("bayesian_subsample",
                        "subsample"))){
    stop("Argument boot_method must be equal to `bayesian_subsample`
or `subsample`.")
  }

if(boot_method == "bayesian_subsample"){
    bootstrap_weights <- rgamma(n_effective,shape = m/n_effective)
    bootstrap_weights <- bootstrap_weights/sum(bootstrap_weights)
    }
if(boot_method == "subsample"){
  bootstrap_weights <- rmultinom(1,m,rep(1/n,n))
  bootstrap_weights <- bootstrap_weights/sum(bootstrap_weights)
}

  bootstrap_weights <- rep(bootstrap_weights, each = J)

  refit <-
    estimate_parameters(W = fitted_model$W,
                        X = fitted_model$X,
                        Z = fitted_model$Z,
                        Z_tilde = fitted_model$Z_tilde,
                        Z_tilde_gamma_cols = fitted_model$Z_tilde_gamma_cols,
                        gammas = fitted_model$gammas,
                        gammas_fixed_indices =
                          fitted_model$gammas_fixed_indices,
                        P = fitted_model$P + .01,
                        P_fixed_indices = fitted_model$P_fixed_indices,
                        B = fitted_model$B,
                        B_fixed_indices = fitted_model$B_fixed_indices,
                        X_tilde = fitted_model$X_tilde,
                        P_tilde = fitted_model$P_tilde + .01,
                        P_tilde_fixed_indices =
                          fitted_model$P_tilde_fixed_indices,
                        gamma_tilde = fitted_model$gamma_tilde,
                        gamma_tilde_fixed_indices =
                          fitted_model$gamma_tilde_fixed_indices,
                        alpha_tilde = fitted_model$alpha_tilde,
                        Z_tilde_list = fitted_model$Z_tilde_list,
                        barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                        barrier_scale = 10, #increments for value of barrier penalty
                        max_barrier = 1e10, #maximum value of barrier_t
                        initial_conv_tol = 1000,
                        final_conv_tol = 0.1,
                        constraint_tolerance = 1e-10,
                        hessian_regularization = .01,
                        criterion = fitted_model$criterion,
                        profile_P = TRUE,
                        profiling_maxit = 25,
                        wts = bootstrap_weights,
                        verbose = FALSE,
                        bootstrap_failure_cutoff = bootstrap_failure_cutoff)

  return(list("varying" = refit$varying,
              "objective" = refit$objective))
}
