#' Calculate derivative of mu_ij with respect to a row of matrix-valued
#' parameter rho
#'
#' @param i The sample index (must be in 1, ..., n)
#' @param J The total number of taxa modeled
#' @param k Row index (which row of rho with respect to which to take derivative)
#' @param gammas Numeric vector of read intensities
#' @param B Detection efficiency matrix
#' @param X The efficiency design matrix (n x p)
#' @param Z The sample design matrix (n x K)
#' @param Ak_list List containing matrices that map back-transformed
#' rho to entries of P
#' @param rho_k Value of kth row of rho
#' @param fixed_P_multipliers Numeric vector of length K containing values in (0,1]
#' equal to 1 - sum(fixed relative abundances in row k of P)
#' Z_tilde to scale by exp(gamma); NULL if no columns to be scaled
#'
#' @return A derivative d mu_i / d rho_k
mu_d_rho_faster <- function(i,
                            J,
                            k,
                            gammas,
                            B,
                            X,
                            Z,
                            Ak_list,
                            rho_k,
                            fixed_P_multipliers,
                            proportion_scale = FALSE){
  dmu_i_dpk <- Matrix::Diagonal(x = as.numeric(Z[i,k]*exp(matrix(gammas[i],nrow = 1, ncol = J) + X[i,,drop=F]%*%B)))
  if(!proportion_scale){
  nu_k <- exp(rho_k)
  front_term <- 1/(1 + sum(nu_k))
  dpk_drhok <- Matrix::Matrix(fixed_P_multipliers[k]*(
    front_term*rbind(diag(nu_k), 0) -
      front_term^2*outer(c(nu_k,1),nu_k))
  )

  rho_jacob_i <- dmu_i_dpk%*%dpk_drhok

  return(rho_jacob_i)
  } else{
    return(dmu_i_dpk)
  }

}
