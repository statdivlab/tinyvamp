#' Apply a likelihod ratio test
#'
#' @param W An \eqn{n \times J} matrix of numeric HTS output (e.g., read counts, coverages, etc.)
#' @param fitted_model The model fitted under the alternative
#' @param null_param The model fitted under the null. This could be a modified copy of the model fit under the alternative.
#' @param n_boot The number of bootstrap resamples to perform
#' @param m description
#' @param recalculate_W0 Boolean. Should rescaling W. TODO. Defaults to FALSE.
#' @param parallelize Boolean. do you want to parallelize it?
#' @param ncores the number of cores to parallelize over
#' @param save_models TODO
#' @param verbose Do you want to know what I'm doing? Defaults to FALSE.
#'
#' @export
bootstrap_lrt <- function(
  W,
  fitted_model,
  null_param,
  n_boot,
  m = NULL,
  recalculate_W0 = FALSE,
  parallelize = FALSE,
  ncores = 1,
  save_models = FALSE,
  verbose = FALSE
) {
  n <- nrow(W)
  J <- ncol(W)

  m <- .tvamp_default_boot_m(n, m)

  if (recalculate_W0) {
    stop("recalculate_W0 = TRUE is not implemented")
  }
  if (save_models) {
    stop("save_models = TRUE is not implemented")
  }

  fit_info <- .tvamp_lrt_boot_criterion_and_weights(fitted_model)
  criterion <- fit_info$criterion
  null_weights <- fit_info$weights
  base_wts <- if (is.null(null_weights)) rep(1, n * J) else null_weights

  fitted_null <- .tvamp_estimate_from_model(
    W = W,
    model_spec = null_param,
    wts = null_weights,
    verbose = verbose,
    criterion = criterion
  )

  W0 <- nullify(W, full_model = fitted_model, null_model = fitted_null)

  worker <- function(i) {
    if (!parallelize && verbose) {
      message("Bootstrap iteration ", i)
    }

    boot_wts <- .tvamp_generate_boot_weights(
      n = n,
      J = J,
      m = m,
      base_wts = base_wts
    )

    boot_full <- .tvamp_estimate_from_model(
      W = W0,
      model_spec = fitted_model,
      wts = boot_wts,
      verbose = verbose,
      criterion = criterion
    )

    boot_null <- .tvamp_estimate_from_model(
      W = W0,
      model_spec = null_param,
      wts = boot_wts,
      verbose = verbose,
      criterion = criterion
    )

    list(
      lr_stat = 2 * (boot_null$objective - boot_full$objective),
      full_model = if (save_models) boot_full else NULL,
      null_model = if (save_models) boot_null else NULL
    )
  }

  boot_results <- .tvamp_boot_apply(
    n_boot = n_boot,
    worker = worker,
    parallelize = parallelize,
    ncores = ncores
  )

  observed_lr_stat <- 2 * (fitted_null$objective - fitted_model$objective)

  boot_lr_stats <- 2 * m * vapply(
    boot_results,
    function(x) x$lr_stat,
    numeric(1)
  )

  list(
    observed_lr_stat = observed_lr_stat,
    boot_lr_stats = boot_lr_stats,
    boot_pval = mean(boot_lr_stats >= observed_lr_stat),
    n_boot = n_boot,
    criterion = fitted_model$criterion,
    recalculate_W0 = recalculate_W0,
    boot_models = if (save_models) boot_results else NULL
  )
}

# bootstrap_lrt <- function(
#   W,
#   fitted_model,
#   null_param,
#   n_boot,
#   m = NULL,
#   recalculate_W0 = FALSE,
#   parallelize = FALSE,
#   ncores = 1,
#   save_models = FALSE,
#   verbose = FALSE
# ) {
#   n <- nrow(W)
#   J <- ncol(W)

#   m <- if (is.null(m)) sqrt(n) else m

#   if (recalculate_W0) {
#     stop("recalculate_W0 = TRUE is not implemented")
#   }
#   if (save_models) {
#     stop("save_models = TRUE is not implemented")
#   }

#   if (fitted_model$criterion == "reweighted_Poisson") {
#     criterion <- "Poisson"
#     null_weights <- fitted_model$weights
#   } else {
#     criterion <- fitted_model$criterion
#     null_weights <- NULL
#   }

#   fit_model <- function(W, model_spec, wts = NULL, verbose = FALSE) {
#     args <- list(
#       W = W,
#       X = model_spec$X,
#       Z = model_spec$Z,
#       Z_tilde = model_spec$Z_tilde,
#       Z_tilde_gamma_cols = model_spec$Z_tilde_gamma_cols,
#       gammas = model_spec$gammas,
#       gammas_fixed_indices = model_spec$gammas_fixed_indices,
#       P = model_spec$P,
#       P_fixed_indices = model_spec$P_fixed_indices,
#       B = model_spec$B,
#       B_fixed_indices = model_spec$B_fixed_indices,
#       X_tilde = model_spec$X_tilde,
#       P_tilde = model_spec$P_tilde,
#       P_tilde_fixed_indices = model_spec$P_tilde_fixed_indices,
#       gamma_tilde = model_spec$gamma_tilde,
#       gamma_tilde_fixed_indices = model_spec$gamma_tilde_fixed_indices,
#       alpha_tilde = model_spec$alpha_tilde,
#       Z_tilde_list = model_spec$Z_tilde_list,
#       barrier_t = 1,
#       barrier_scale = 10,
#       max_barrier = 1e12,
#       initial_conv_tol = 1000,
#       final_conv_tol = 0.1,
#       constraint_tolerance = 1e-10,
#       hessian_regularization = 0.01,
#       criterion = criterion,
#       profile_P = TRUE,
#       wts = wts,
#       verbose = verbose,
#       profiling_maxit = 25
#     )

#     do.call(estimate_parameters, args)
#   }

#   fitted_null <- fit_model(
#     W = W,
#     model_spec = null_param,
#     wts = null_weights,
#     verbose = verbose
#   )

#   base_wts <- if (is.null(null_weights)) rep(1, n * J) else null_weights

#   make_boot_weights <- function() {
#     bwts <- rgamma(n, shape = m / n)
#     bwts <- rep(bwts, each = J)
#     bwts <- bwts * base_wts
#     J * bwts / sum(bwts)
#   }

#   W0 <- nullify(W, full_model = fitted_model, null_model = fitted_null)

#   run_boot <- function(i) {
#     if (!parallelize && verbose) {
#       message("Bootstrap iteration ", i)
#     }

#     boot_wts <- make_boot_weights()

#     boot_full <- fit_model(
#       W = W0,
#       model_spec = fitted_model,
#       wts = boot_wts,
#       verbose = verbose
#     )

#     boot_null <- fit_model(
#       W = W0,
#       model_spec = null_param,
#       wts = boot_wts,
#       verbose = verbose
#     )

#     list(
#       lr_stat = 2 * (boot_null$objective - boot_full$objective)
#     )
#   }

#   boot_results <- if (parallelize) {
#     parallel::mclapply(
#       X = seq_len(n_boot),
#       FUN = run_boot,
#       mc.cores = ncores,
#       mc.set.seed = TRUE
#     )
#   } else {
#     lapply(seq_len(n_boot), run_boot)
#   }

#   observed_lr_stat <- 2 * (fitted_null$objective - fitted_model$objective)

#   boot_lr_stats <- 2 * m * vapply(
#     boot_results,
#     function(x) x$lr_stat,
#     numeric(1)
#   )

#   list(
#     observed_lr_stat = observed_lr_stat,
#     boot_lr_stats = boot_lr_stats,
#     boot_pval = mean(boot_lr_stats >= observed_lr_stat),
#     n_boot = n_boot,
#     criterion = fitted_model$criterion,
#     recalculate_W0 = recalculate_W0
#   )
# }
# bootstrap_lrt <- function(
#   W,
#   fitted_model,
#   null_param,
#   n_boot,
#   m = NULL,
#   recalculate_W0 = FALSE,
#   parallelize = FALSE,
#   ncores = 1,
#   save_models = FALSE,
#   verbose = FALSE
# ) {
#   n <- nrow(W)
#   J <- ncol(W)

#   if (is.null(m)) {
#     m <- sqrt(n)
#   }

#   if (fitted_model$criterion == "reweighted_Poisson") {
#     criterion <- "Poisson"
#     null_weights <- fitted_model$weights
#   } else {
#     criterion <- fitted_model$criterion
#     null_weights <- NULL
#   }

#   fitted_null <- estimate_parameters(
#     W = W,
#     X = null_param$X,
#     Z = null_param$Z,
#     Z_tilde = null_param$Z_tilde,
#     Z_tilde_gamma_cols = null_param$Z_tilde_gamma_cols,
#     gammas = null_param$gammas,
#     gammas_fixed_indices = null_param$gammas_fixed_indices,
#     P = null_param$P,
#     P_fixed_indices = null_param$P_fixed_indices,
#     B = null_param$B,
#     B_fixed_indices = null_param$B_fixed_indices,
#     X_tilde = null_param$X_tilde,
#     P_tilde = null_param$P_tilde,
#     P_tilde_fixed_indices = null_param$P_tilde_fixed_indices,
#     gamma_tilde = null_param$gamma_tilde,
#     gamma_tilde_fixed_indices = null_param$gamma_tilde_fixed_indices,
#     alpha_tilde = null_param$alpha_tilde,
#     Z_tilde_list = null_param$Z_tilde_list,
#     barrier_t = 1, #starting value of reciprocal barrier penalty coef.
#     barrier_scale = 10, #increments for value of barrier penalty
#     max_barrier = 1e12, #maximum value of barrier_t
#     initial_conv_tol = 1000,
#     final_conv_tol = 0.1,

#     constraint_tolerance = 1e-10,
#     hessian_regularization = 0.01,
#     criterion = "Poisson",

#     profile_P = TRUE,
#     wts = null_weights,
#     verbose = verbose,
#     profiling_maxit = 25
#   )

#   if (is.null(null_weights)) {
#     wts <- rep(1, n * J)
#   } else {
#     wts <- null_weights
#   }

#   if (!recalculate_W0) {
#     W0 <- nullify(W, full_model = fitted_model, null_model = fitted_null)
#   }

#   boot_results <- vector(n_boot, mode = "list")

#   boot_weights <- lapply(1:n_boot, function(x) {
#     # set.seed(x)
#     bwts <- rgamma(n, m / n)
#     bwts <- rep(bwts, each = J)
#     bwts <- bwts * wts
#     bwts <- J * bwts / sum(bwts)
#     return(bwts)
#   })

#   if (!parallelize) {
#     for (boot_iter in 1:n_boot) {
#       if (recalculate_W0) {
#         stop("David you haven't implemented recalculation of W0")
#       }
#       message(paste("Bootstrap iteration ", boot_iter, sep = "", collapse = ""))
#       # boot_weights <- rgamma(n,m/n)
#       # boot_weights <- boot_weights/sum(boot_weights)
#       # boot_weights <- rep(boot_weights, each = J)

#       boot_full <- estimate_parameters(
#         W = W0,
#         X = fitted_model$X,
#         Z = fitted_model$Z,
#         Z_tilde = fitted_model$Z_tilde,
#         Z_tilde_gamma_cols = fitted_model$Z_tilde_gamma_cols,
#         gammas = fitted_model$gammas,
#         gammas_fixed_indices = fitted_model$gammas_fixed_indices,
#         P = fitted_model$P,
#         P_fixed_indices = fitted_model$P_fixed_indices,
#         B = fitted_model$B,
#         B_fixed_indices = fitted_model$B_fixed_indices,
#         X_tilde = fitted_model$X_tilde,
#         P_tilde = fitted_model$P_tilde,
#         P_tilde_fixed_indices = fitted_model$P_tilde_fixed_indices,
#         gamma_tilde = fitted_model$gamma_tilde,
#         gamma_tilde_fixed_indices = fitted_model$gamma_tilde_fixed_indices,
#         alpha_tilde = fitted_model$alpha_tilde,
#         Z_tilde_list = fitted_model$Z_tilde_list,
#         barrier_t = 1, #starting value of reciprocal barrier penalty coef.
#         barrier_scale = 10, #increments for value of barrier penalty
#         max_barrier = 1e12, #maximum value of barrier_t
#         initial_conv_tol = 1000,
#         final_conv_tol = 0.1,

#         constraint_tolerance = 1e-10,
#         hessian_regularization = 0.01,
#         criterion = "Poisson",

#         profile_P = TRUE,
#         verbose = verbose,
#         wts = boot_weights[[boot_iter]],
#         profiling_maxit = 25
#       )

#       boot_null <- estimate_parameters(
#         W = W0,
#         X = null_param$X,
#         Z = null_param$Z,
#         Z_tilde = null_param$Z_tilde,
#         Z_tilde_gamma_cols = null_param$Z_tilde_gamma_cols,
#         gammas = null_param$gammas,
#         gammas_fixed_indices = null_param$gammas_fixed_indices,
#         P = null_param$P,
#         P_fixed_indices = null_param$P_fixed_indices,
#         B = null_param$B,
#         B_fixed_indices = null_param$B_fixed_indices,
#         X_tilde = null_param$X_tilde,
#         P_tilde = null_param$P_tilde,
#         P_tilde_fixed_indices = null_param$P_tilde_fixed_indices,
#         gamma_tilde = null_param$gamma_tilde,
#         gamma_tilde_fixed_indices = null_param$gamma_tilde_fixed_indices,
#         alpha_tilde = null_param$alpha_tilde,
#         Z_tilde_list = null_param$Z_tilde_list,
#         barrier_t = 1, #starting value of reciprocal barrier penalty coef.
#         barrier_scale = 10, #increments for value of barrier penalty
#         max_barrier = 1e12, #maximum value of barrier_t
#         initial_conv_tol = 1000,
#         final_conv_tol = 0.1,

#         constraint_tolerance = 1e-10,
#         hessian_regularization = 0.01,
#         criterion = "Poisson",

#         profile_P = TRUE,
#         verbose = verbose,
#         wts = boot_weights[[boot_iter]],
#         profiling_maxit = 25
#       )

#       boot_results[[boot_iter]]$lr_stat <- 2 *
#         (boot_null$objective -
#           boot_full$objective)
#       boot_results[[boot_iter]]$weights <- boot_weights

#       if (save_models) {
#         stop("David you haven't implemented save_models")
#       }
#     }
#   }
#   if (parallelize) {
#     boot_results <- parallel::mclapply(
#       1:n_boot,
#       function(k) {
#         do_one_lrt(
#           W0 = W0,
#           full_model = fitted_model,
#           null_model = fitted_null, #null model specification
#           boot_method = "bayesian_subsample",
#           boot_weights = boot_weights[[k]],
#           return_models = FALSE
#         )
#       },
#       mc.cores = ncores,
#       mc.set.seed = TRUE
#     )
#   }
#   observed_lr_stat <- 2 * (fitted_null$objective - fitted_model$objective)
#   if (!parallelize) {
#     boot_lr_stats <- 2 *
#       m *
#       (sapply(1:n_boot, function(k) {
#         boot_results[[k]]$lr_stat
#       }))
#   } else {
#     boot_lr_stats <- 2 *
#       m *
#       (sapply(1:n_boot, function(k) {
#         boot_results[[k]]
#       }))
#   }

#   return(list(
#     "observed_lr_stat" = observed_lr_stat,
#     "boot_lr_stats" = boot_lr_stats,
#     "boot_pval" = mean(boot_lr_stats >= observed_lr_stat),
#     "n_boot" = n_boot,
#     "criterion" = fitted_model$criterion,
#     "recalculate_W0" = recalculate_W0
#   ))
# }
