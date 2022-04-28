#' Calculate derivative of mu_ij with respect to an entry of P
#'
#' @param i The sample index (must be in 1, ..., n)
#' @param j The taxon index (must be in 1, ..., J)
#' @param m The row of P with respect to which to take derivative
#' @param gammas Numeric vector of read intensities
#' @param B Detection efficiency matrix
#' @param X The efficiency design matrix (n x p)
#' @param Z The sample design matrix (n x K)
#' @param P The sample relative abundance matrix (K x J)
#' @param X_tilde The spurious read efficiency design (K_tilde x p)
#' @param Z_tilde The spurious read design (n x K_tilde)
#' @param Z_tilde_gamma_cols Numeric vector containing indexes of columns of
#' Z_tilde to scale by exp(gamma); NULL if no columns to be scaled
#' @param P_tilde The spurious source relative abundance matrix (K_tilde x J)
#' @param gamma_tilde Spurious read intensity parameter
#'
#' @return A derivative d mu_ij / d P_kj
mu_d_P <- function(i,
                   j,
                   m,
                   gammas,
                   B,
                   X,
                   Z,
                   P){

  mu_deriv <-
    Z[i,m,drop = F]*exp(gammas[i] +
                          X[i,,drop = F]%*%B[,j,drop = F])

  return(mu_deriv)
}
