#' Criterion Evaluation Function
#'
#' Evaluates the objective function used during model fitting under either a
#' Poisson likelihood criterion or a GMM criterion.
#'
#' @param W An \eqn{n \times J} matrix of numeric HTS output (e.g., read counts, coverages, etc.)
#' @param X The sample efficiency design -- an \eqn{n \times p} matrix
#' @param Z The sample-specimen design -- an \eqn{n \times K} matrix whose \eqn{ij}-th entry
#' indicates the proportional contribution of specimen \eqn{j} to sample \eqn{i}. Rows must
#' sum to 1 or be identically 0.
#' @param Z_tilde The spurious read design -- an \eqn{n \times \tilde{K}} matrix where
#' \eqn{\tilde{K}} is the number of spurious read sources modeled.
#' @param Z_tilde_gamma_cols A numeric vector containing the columns of Z_tilde which should be
#' multiplied by exp(gamma).
#' @param Z_tilde_list Optional list-form representation of `Z_tilde` or
#'   related auxiliary structures.
#' @param X_tilde Design matrix associated with spurious-read or auxiliary
#'   components of the mean model.
#' @param fixed_df A data frame of fixed model parameters.
#' @param varying_df A data frame of parameters currently being optimized,
#'   represented on the natural scale unless `lr_scale = TRUE`.
#' @param varying_lr_df Optional data frame of varying parameters on the
#'   log-ratio scale. Used when `lr_scale = TRUE`.
#' @param barrier_t Optional numeric tuning parameter for the log-ratio barrier
#'   penalty.
#' @param criterion Character string specifying the criterion to evaluate.
#'   Currently one of `"Poisson"` or `"GMM"`.
#' @param lr_scale Logical; if `TRUE`, convert log-ratio scale parameters to
#'   relative abundance scale before evaluation.
#' @param include_log_penalty Logical; if `TRUE`, include the log-barrier
#'   penalty when `lr_scale = TRUE`.
#' @param wts Optional numeric weights
#' @param gmm_inv_wts Optional inverse weighting matrix or vector for the GMM
#'   criterion. If `NULL`, it is estimated internally
#' @param return_gmm_inv_weights Logical; if `TRUE` and `criterion = "GMM"`,
#'   return both the GMM criterion value and the inverse weights.
#'
#' @return
#' If `criterion = "Poisson"`, a single numeric criterion value.
#'
#' If `criterion = "GMM"` and `return_gmm_inv_weights = FALSE`, a single numeric
#' criterion value.
#'
#' If `criterion = "GMM"` and `return_gmm_inv_weights = TRUE`, a list with
#' components:
#' \describe{
#'   \item{gmm_crit}{The numeric GMM criterion value.}
#'   \item{inv_wts}{The inverse weights used in the GMM calculation.}
#' }
#'
#' @export
evaluate_criterion_lr <- function(
  W,
  X,
  Z,
  Z_tilde,
  Z_tilde_gamma_cols,
  Z_tilde_list = NULL,
  X_tilde,
  fixed_df,
  varying_df,
  varying_lr_df = NULL,
  barrier_t = NULL,
  criterion = "Poisson",
  lr_scale = TRUE,
  include_log_penalty = TRUE,
  wts = NULL,
  gmm_inv_wts = NULL,
  return_gmm_inv_weights = FALSE
) {
  if (lr_scale) {
    varying_df <- lr_to_ra(fixed_df, varying_lr_df, varying_df)
  }

  params <- dataframes_to_parameters(fixed_df, varying_df)

  if (lr_scale) {
    if (include_log_penalty) {
      log_penalty <- calculate_log_penalty(varying_lr_df, fixed_df, barrier_t)
    } else {
      log_penalty <- 0
    }
  } else {
    log_penalty <- 0
  }

  means <- meaninate(
    gammas = params$gammas,
    B = params$B,
    X = X,
    Z = Z,
    P = params$P,
    X_tilde = X_tilde,
    Z_tilde = Z_tilde,
    Z_tilde_gamma_cols = Z_tilde_gamma_cols,
    Z_tilde_list = Z_tilde_list,
    P_tilde = params$P_tilde,
    gamma_tilde = params$gamma_tilde,
    alpha_tilde = params$alpha_tilde,
    return_separate = FALSE
  )

  if (criterion == "Poisson") {
    return(poisson_criterion(W = W, means = means, wts = wts) + log_penalty)
  }
  if (criterion == "GMM") {
    n <- nrow(W)
    W_long <- lapply(1:n, function(i) as.numeric(W[i, ]))
    W_long <- do.call(c, W_long)
    means_long <- lapply(1:n, function(i) as.numeric(means[i, ]))
    means_long <- do.call(c, means_long)

    if (is.null(gmm_inv_wts)) {
      stop("We lost get_gmm_inv_weights somewhere. Amy can find it if you need it.")
      # inv_wts <- get_gmm_inv_weights(W_long = W_long, means_long = means_long)
    } else {
      inv_wts <- gmm_inv_wts
    }

    if (!return_gmm_inv_weights) {
      return(
        gmm_criterion(
          W_long = W_long,
          means_long = means_long,
          inv_wts = inv_wts
        ) +
          log_penalty
      )
    } else {
      return(list(
        "gmm_crit" = gmm_criterion(
          W_long = W_long,
          means_long = means_long,
          inv_wts = inv_wts
        ) +
          log_penalty,
        "inv_wts" = inv_wts
      ))
    }
  }
}
