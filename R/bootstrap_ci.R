#' Apply the Bayesian subsampled bootstrap to a fitted tinyvamp model
#'
#' Compute bootstrap confidence intervals for the varying parameters of a
#' fitted tinyvamp model using a Bayesian subsampled bootstrap procedure.
#' Each bootstrap replicate refits the model using randomly generated
#' observation weights, and confidence intervals are formed from the
#' empirical quantiles of the bootstrap distribution.
#'
#' @param W A numeric matrix of observed counts with \eqn{n} rows
#'   (samples) and \eqn{J} columns (taxa).
#' @param fitted_model A fitted tinyvamp model object, typically returned by
#'   \code{estimate_parameters()}, containing the original fit and all inputs
#'   needed for refitting during bootstrap iterations.
#' @param n_boot Integer; the number of bootstrap replicates.
#' @param m Optional numeric subsample size parameter for the Bayesian
#'   subsampled bootstrap. If \code{NULL}, defaults to \code{sqrt(n)}, where
#'   \eqn{n} is the number of rows of \code{W}.
#' @param alpha Significance level used to construct two-sided
#'   confidence intervals. The default \code{0.05} yields 95% intervals.
#' @param parallelize Logical; if \code{TRUE}, bootstrap replicates are fit in
#'   parallel using \code{parallel::mclapply()}. If \code{FALSE}, replicates are
#'   fit sequentially.
#' @param ncores Integer; number of cores to use when
#'   \code{parallelize = TRUE}.
#' @param seed Optional integer random seed. If \code{NULL}, a seed value of
#'   \code{0} is used.
#' @param return_models Logical; if \code{TRUE}, return the fitted bootstrap
#'   model objects in addition to the confidence interval summary.
#' @param verbose Logical; if \code{TRUE}, print progress messages during
#'   sequential bootstrap fitting and pass verbose output to the internal model
#'   fitting routine.
#' @param adjust Logical; if \code{TRUE}, apply a finite-sample adjustment to
#'   the bootstrap deviations based on the number of non-\code{"gamma"}
#'   varying parameters.
#'
#' @details
#' For each bootstrap replicate, the function generates a vector of Bayesian
#' bootstrap weights from a gamma distribution, rescales them, and uses them as
#' observation weights in a call to \code{estimate_parameters()}. The resulting
#' bootstrap distribution is centered at the original fitted values and scaled
#' by \eqn{\sqrt{m}}. Confidence intervals are then obtained by inverting the
#' bootstrap quantiles and scaling by \eqn{1 / \sqrt{n}}.
#'
#' When \code{adjust = TRUE}, the bootstrap deviations are multiplied by a
#' correction factor
#' \deqn{\sqrt{\frac{nJ}{nJ - p}}}
#' where \eqn{p} is the number of varying parameters whose \code{param} field is
#' not equal to \code{"gamma"}.
#'
#'
#' @return
#' A list with components \code{ci} and \code{bootstrapped_models}. The former
#' is a data frame with columns \code{lower_ci} and \code{upper_ci}. If
#'  \code{return_models = TRUE}, the latter is a list of
#'  fitted model objects from the bootstrap replicates; otherwise \code{NULL}.
#'
#' @import stats
#' @import parallel
#' @export
bootstrap_ci <- function(
  W,
  fitted_model,
  n_boot,
  m = NULL,
  alpha = 0.05,
  parallelize = FALSE,
  ncores = 5,
  seed = NULL,
  return_models = FALSE,
  verbose = FALSE,
  adjust = FALSE
) {
  n <- nrow(W)
  J <- ncol(W)

  if (is.null(m)) {
    m <- sqrt(n)
  }
  if (is.null(seed)) {
    seed <- 0
  }

  if (is.null(fitted_model$wts)) {
    fitted_model$wts <- rep(1, n * J)
  }
  wts <- fitted_model$wts
  set.seed(seed)
  boot_seeds <- sample(1:1e8, n_boot)

  boot_results <- vector(n_boot, mode = "list")

  if (!parallelize) {
    for (boot_iter in 1:n_boot) {
      if (verbose) {
        message(paste0("Bootstrap iteration ", boot_iter))
      }
      set.seed(boot_seeds[boot_iter])
      boot_weights <- rgamma(n, m / n)
      boot_weights <- boot_weights / sum(boot_weights)
      boot_weights <- rep(boot_weights, each = J)
      boot_weights <- boot_weights * wts
      boot_weights <- J * boot_weights / sum(boot_weights)

      boot_model <- estimate_parameters(
        W = W,
        X = fitted_model$X,
        Z = fitted_model$Z,
        Z_tilde = fitted_model$Z_tilde,
        Z_tilde_gamma_cols = fitted_model$Z_tilde_gamma_cols,
        gammas = fitted_model$gammas,
        gammas_fixed_indices = fitted_model$gammas_fixed_indices,
        P = fitted_model$P,
        P_fixed_indices = fitted_model$P_fixed_indices,
        B = fitted_model$B,
        B_fixed_indices = fitted_model$B_fixed_indices,
        X_tilde = fitted_model$X_tilde,
        P_tilde = fitted_model$P_tilde,
        P_tilde_fixed_indices = fitted_model$P_tilde_fixed_indices,
        gamma_tilde = fitted_model$gamma_tilde,
        gamma_tilde_fixed_indices = fitted_model$gamma_tilde_fixed_indices,
        alpha_tilde = fitted_model$alpha_tilde,
        Z_tilde_list = fitted_model$Z_tilde_list,
        barrier_t = 1, #starting value of reciprocal barrier penalty coef.
        barrier_scale = 10, #increments for value of barrier penalty
        max_barrier = 1e12, #maximum value of barrier_t
        initial_conv_tol = 1000,
        final_conv_tol = 0.1,
        constraint_tolerance = 1e-10,
        hessian_regularization = 0.01,
        criterion = "Poisson",
        profile_P = TRUE,
        verbose = verbose,
        wts = boot_weights,
        profiling_maxit = 25
      )

      boot_results[[boot_iter]] <- boot_model
    }
  }
  if (parallelize) {
    boot_weights <- lapply(1:n_boot, function(x) {
      # set.seed(x)
      bwts <- rgamma(n, m / n)
      bwts <- rep(bwts, each = J)
      bwts <- bwts * wts
      bwts <- J * bwts / sum(bwts)
      return(bwts)
    })

    boot_results <-
      parallel::mclapply(
        1:n_boot,
        function(k) {
          do_one_boot(
            W = W,
            fitted_model = fitted_model,
            m = m,
            seed = boot_seeds[k],
            boot_weights = boot_weights[[k]]
          )
        },
        mc.cores = ncores,
        mc.set.seed = TRUE
      )
  }

  boot_matrix <-
    do.call(
      cbind,
      lapply(1:n_boot, function(k) {
        matrix(
          sqrt(m) *
            (boot_results[[k]]$varying$value - fitted_model$varying$value),
          ncol = 1
        )
      })
    )

  if (adjust) {
    num_nonsillyparams <- sum(fitted_model$varying$param != "gamma")
    adjust_factor <- (n * J) / (n * J - num_nonsillyparams)
    boot_matrix <- boot_matrix * sqrt(adjust_factor)
  }

  lower_boot_quantiles <- apply(boot_matrix, 1, function(x) {
    quantile(x, alpha / 2)
  })

  upper_boot_quantiles <- apply(boot_matrix, 1, function(x) {
    quantile(x, 1 - alpha / 2)
  })

  summary_df <- fitted_model$varying

  summary_df$lower_ci <- summary_df$value - (1 / sqrt(n)) * upper_boot_quantiles
  summary_df$upper_ci <- summary_df$value - (1 / sqrt(n)) * lower_boot_quantiles

  summary_df$lower_ci[summary_df$param %in% c("P", "P_tilde")] <-
    pmax(summary_df$lower_ci[summary_df$param %in% c("P", "P_tilde")], 0)

  summary_df$upper_ci[summary_df$param %in% c("P", "P_tilde")] <-
    pmin(summary_df$upper_ci[summary_df$param %in% c("P", "P_tilde")], 1)

  if (return_models) {
    return(list("ci" = summary_df, "bootstrapped_models" = boot_results))
  } else {
    return(list("ci" = summary_df, "bootstrapped_models" = NULL))
  }
}
