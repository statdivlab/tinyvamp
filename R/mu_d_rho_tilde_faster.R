#' Calculate derivative of mu_i with respect to a row of matrix-valued
#' parameter rho_tilde
#'
#' @param i The sample index (must be in 1, ..., n)
#' @param J The total number of taxa modeled
#' @param k_tilde Row index (which row of rho_tilde with respect to which to take derivative)
#' @param gammas Numeric vector of read intensities
#' @param B Detection efficiency matrix
#' @param rho_tilde_k Value of kth row of rho_tilde
#' @param A_tilde_k_list List containing matrices that map back-transformed
#' rho_tilde to entries of P_tilde
#' @param fixed_P_multipliers Numeric vector of length K containing values in (0,1]
#' equal to 1 - sum(fixed relative abundances in row k of P_tilde)
#' @param alpha_tilde A numeric vector containing starting values of length \eqn{M}. If used, \code{Z_tilde_list} must be provided.
#' @param Z_tilde_list A list of length \eqn{M + 1} containing matrices \eqn{\tilde{Z}_1,\dots,\tilde{Z}_{M + 1}} to be linearly combined to
#' generate \code{Z_tilde}: \eqn{\tilde{Z} = \tilde{Z}_{(1)} + \sum_{m = 1}^M \exp(\tilde{\alpha}_m)\tilde{Z}_{(m + 1)}}. If used,
#' \code{alpha_tilde} must be provided.
#' @param X_tilde The spurious read efficiency design (K_tilde x p)
#' @param Z_tilde The spurious read design (n x K_tilde)
#' @param Z_tilde_gamma_cols Numeric vector containing indexes of columns of
#' Z_tilde to scale by exp(gamma); NULL if no columns to be scaled
#'
#' @return A derivative d mu_i / d rho_tilde_k
mu_d_rho_tilde_faster <- function(i,
                                  J,
                                  k_tilde,
                                  gammas,
                                  B,
                                  rho_tilde_k,
                                  A_tilde_k_list,
                                  fixed_P_tilde_multipliers,
                                  X_tilde,
                                  Z_tilde,
                                  Z_tilde_gamma_cols,
                                  alpha_tilde = NULL,
                                  Z_tilde_list = NULL,
                                  gamma_tilde,
                                  proportion_scale = FALSE)
{

  if(!is.null(alpha_tilde)){
    Z_tilde <- construct_Z_tilde(Z_tilde_list,
                                 alpha_tilde)
  }

  for(zcol in Z_tilde_gamma_cols){
    Z_tilde[,zcol] <- exp(gammas)*Z_tilde[,zcol]
  }

  dmu_i_dptildektilde <- Matrix::Diagonal(x = as.numeric(Z_tilde[i,k_tilde]*
                                                           exp(matrix(gamma_tilde[k_tilde],nrow = 1, ncol = J) +
                                                                 X_tilde[k_tilde,,drop = F]%*%B)))

  if(!proportion_scale){
  nu_tilde_k <- exp(rho_tilde_k)

  front_term <- 1/(1 + sum(nu_tilde_k))
  dptildek_drhotildek <- Matrix::Matrix(fixed_P_tilde_multipliers[k_tilde]*(
    front_term*rbind(diag(nu_tilde_k), 0) -
      front_term^2*outer(c(nu_tilde_k,1),nu_tilde_k)
  ))

  rho_tilde_jacob_i <- dmu_i_dptildektilde%*%A_tilde_k_list[[k_tilde]]%*%dptildek_drhotildek

  return(rho_tilde_jacob_i)
  } else{
    return(dmu_i_dptildektilde)
  }

}
