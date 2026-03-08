## March 8 2026
## AW: I was trying to clean up redundancy in bootstrap_ci and bootstrap_lrt
## and ChatGPT suggested the following internal help of functions as a replacement to
## previous do_one_lrt and do_one_boot. 
## Revert at the first sign of trouble. 

.tvamp_default_boot_m <- function(n, m) {
  if (is.null(m)) sqrt(n) else m
}

.tvamp_default_seed <- function(seed) {
  if (is.null(seed)) 0L else as.integer(seed)
}

.tvamp_base_weights <- function(fitted_model, n, J) {
  wts <- fitted_model$wts
  if (is.null(wts)) {
    rep(1, n * J)
  } else {
    wts
  }
}

.tvamp_generate_boot_weights <- function(n, J, m, base_wts, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }

  bwts <- rgamma(n, shape = m / n)
  bwts <- rep(bwts, each = J)
  bwts <- bwts * base_wts
  J * bwts / sum(bwts)
}

.tvamp_estimate_from_model <- function(
  W,
  model_spec,
  wts = NULL,
  verbose = FALSE,
  criterion = "Poisson"
) {
  args <- list(
    W = W,
    X = model_spec$X,
    Z = model_spec$Z,
    Z_tilde = model_spec$Z_tilde,
    Z_tilde_gamma_cols = model_spec$Z_tilde_gamma_cols,
    gammas = model_spec$gammas,
    gammas_fixed_indices = model_spec$gammas_fixed_indices,
    P = model_spec$P,
    P_fixed_indices = model_spec$P_fixed_indices,
    B = model_spec$B,
    B_fixed_indices = model_spec$B_fixed_indices,
    X_tilde = model_spec$X_tilde,
    P_tilde = model_spec$P_tilde,
    P_tilde_fixed_indices = model_spec$P_tilde_fixed_indices,
    gamma_tilde = model_spec$gamma_tilde,
    gamma_tilde_fixed_indices = model_spec$gamma_tilde_fixed_indices,
    alpha_tilde = model_spec$alpha_tilde,
    Z_tilde_list = model_spec$Z_tilde_list,
    barrier_t = 1,
    barrier_scale = 10,
    max_barrier = 1e12,
    initial_conv_tol = 1000,
    final_conv_tol = 0.1,
    constraint_tolerance = 1e-10,
    hessian_regularization = 0.01,
    criterion = criterion,
    profile_P = TRUE,
    verbose = verbose,
    wts = wts,
    profiling_maxit = 25
  )

  do.call(estimate_parameters, args)
}

.tvamp_boot_apply <- function(
  n_boot,
  worker,
  parallelize = FALSE,
  ncores = 1
) {
  idx <- seq_len(n_boot)

  if (parallelize) {
    parallel::mclapply(
      X = idx,
      FUN = worker,
      mc.cores = ncores,
      mc.set.seed = FALSE
    )
  } else {
    lapply(idx, worker)
  }
}

.tvamp_boot_matrix <- function(boot_results, theta_hat, m) {
  vapply(
    boot_results,
    function(res) sqrt(m) * (res$varying$value - theta_hat),
    numeric(length(theta_hat))
  )
}

.tvamp_clip_probability_ci <- function(summary_df) {
  prob_idx <- summary_df$param %in% c("P", "P_tilde")

  summary_df$lower_ci[prob_idx] <- pmax(summary_df$lower_ci[prob_idx], 0)
  summary_df$upper_ci[prob_idx] <- pmin(summary_df$upper_ci[prob_idx], 1)

  summary_df
}

.tvamp_boot_quantiles <- function(boot_matrix, alpha) {
  list(
    lower = apply(boot_matrix, 1, quantile, probs = alpha / 2, names = FALSE),
    upper = apply(boot_matrix, 1, quantile, probs = 1 - alpha / 2, names = FALSE)
  )
}

.tvamp_lrt_boot_criterion_and_weights <- function(fitted_model) {
  if (identical(fitted_model$criterion, "reweighted_Poisson")) {
    list(
      criterion = "Poisson",
      weights = fitted_model$weights
    )
  } else {
    list(
      criterion = fitted_model$criterion,
      weights = NULL
    )
  }
}