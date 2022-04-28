#' Calculate derivative of mu_ij with respect to ith entry of gamma
#'
#' @param i The sample index (must be in 1, ..., n)
#' @param j The taxon index (must be in 1, ..., J)
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
#' @return A derivative d mu_ij / d gamma_i

mu_d_gamma <- function(i,
                       j,
                       gammas,
                       B,
                       X,
                       Z,
                       P,
                       X_tilde,
                       Z_tilde,
                       Z_tilde_gamma_cols,
                       P_tilde,
                       gamma_tilde){

  mu_deriv <- Z[i,,drop = F]%*%P[,j,drop = F]*exp(gammas[i] +
                                                    X[i,,drop = F]%*%B[,j,drop = F])

  K_tilde <- dim(P_tilde)[1]
  for(k_tilde in 1:K_tilde){
    if(k_tilde %in% Z_tilde_gamma_cols){
      mu_deriv <- mu_deriv + exp(gammas[i])*
        (Z_tilde[i,k_tilde,drop = F])%*%
        (P_tilde[k_tilde,j,drop = F])*
        exp(gamma_tilde[k_tilde] +
              X_tilde[k_tilde,,drop = F]%*%B[,j,drop = F])
    }
  }
  return(mu_deriv)
}
