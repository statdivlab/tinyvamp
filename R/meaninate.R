#' Calculate the mean given parameters
#'
#' @param gammas A numeric vector of length n of starting values for read intensity parameter gamma
#' @param B A \eqn{p \times J} numeric matrix giving initial values for the sample efficiencies.
#' @param X The sample efficiency design -- an \eqn{n \times p} matrix
#' @param Z The sample-specimen design -- an \eqn{n \times K} matrix whose \eqn{ij}-th entry
#' indicates the proportional contribution of specimen \eqn{j} to sample \eqn{i}. Rows must
#' sum to 1 or be identically 0.
#' @param Z_tilde The spurious read design -- an \eqn{n \times \tilde{K}} matrix where
#' \eqn{\tilde{K}} is the number of spurious read sources modeled.
#' @param Z_tilde_gamma_cols A numeric vector containing the columns of Z_tilde which should be
#' multiplied by exp(gamma).
#' @param P A \eqn{K \times J} numeric matrix giving initial values for the relative abundance matrix.
#' @param X_tilde A \eqn{\tilde{K} \times p} matrix giving the spurious read source efficiency design matrix
#' @param P_tilde A \eqn{\tilde{K} \times J} numeric matrix giving initial values for the spurious read source relative abundances.
#' @param gamma_tilde A numeric vector of length \eqn{\tilde{K}} of starting values for spurious read intensity parameter gamma_tilde
#' @param alpha_tilde A numeric vector containing starting values of length \eqn{M}. If used, \code{Z_tilde_list} must be provided.
#' @param Z_tilde_list A list of length \eqn{M + 1} containing matrices \eqn{\tilde{Z}_1,\dots,\tilde{Z}_{M + 1}} to be linearly combined to
#' generate \code{Z_tilde}: \eqn{\tilde{Z} = \tilde{Z}_{(1)} + \sum_{m = 1}^M \exp(\tilde{\alpha}_m)\tilde{Z}_{(m + 1)}}. If used,
#' \code{alpha_tilde} must be provided.
#' @param return_separate Boolean. Return the summed mean, or separate the sample and contamination pieces. Defaults to FALSE.
#' @param exclude_gammas Boolean, defaults to FALSE. Should the gamma components be ignored?
#'
meaninate <- function(
  gammas,
  B,
  X,
  Z,
  P,
  X_tilde,
  Z_tilde = NULL,
  Z_tilde_gamma_cols,
  P_tilde,
  gamma_tilde,
  alpha_tilde = NULL,
  Z_tilde_list = NULL,
  return_separate = FALSE,
  exclude_gammas = FALSE
) {
  if (!is.null(alpha_tilde)) {
    Z_tilde <- construct_Z_tilde(Z_tilde_list, alpha_tilde)
  }

  J <- ncol(B)
  n <- nrow(X)

  if (!exclude_gammas) {
    #multiply appropriate columns of Z_tilde by exp(gamma)
    if (length(Z_tilde_gamma_cols > 0)) {
      for (colnum in Z_tilde_gamma_cols) {
        Z_tilde[, colnum] <- exp(gammas) * Z_tilde[, colnum]
      }
    }

    sample_part <-
      (Z %*% P) * (exp(gammas %*% matrix(1, nrow = 1, ncol = J) + X %*% B))
    spurious_part <- Z_tilde %*%
      (P_tilde *
        exp(gamma_tilde %*% matrix(1, nrow = 1, ncol = J) + X_tilde %*% B))
  } else {
    sample_part <-
      (Z %*% P) * (exp(X %*% B))
    spurious_part <- Z_tilde %*%
      (P_tilde *
        exp(gamma_tilde %*% matrix(1, nrow = 1, ncol = J) + X_tilde %*% B))
  }

  if (return_separate) {
    return(list("sample" = sample_part, "spurious" = spurious_part))
  } else {
    return(sample_part + spurious_part)
  }
}
