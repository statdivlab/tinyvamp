
mu_d_gamma_faster <- function(i,
                              J,
                              gammas,
                              B,
                              X,
                              Z,
                              P,
                              X_tilde,
                              Z_tilde,
                              Z_tilde_gamma_cols,
                              alpha_tilde = NULL,
                              Z_tilde_list = NULL,
                              P_tilde,
                              gamma_tilde){

  if(!is.null(alpha_tilde)){
    Z_tilde <- construct_Z_tilde(Z_tilde_list,
                                 alpha_tilde)
  }

  mu_deriv <- Z[i,,drop = F]%*%P*exp(gammas[i] +
                                       X[i,,drop = F]%*%B)

  K_tilde <- dim(P_tilde)[1]
  for(k_tilde in 1:K_tilde){
    if(k_tilde %in% Z_tilde_gamma_cols){
      mu_deriv <- mu_deriv + exp(gammas[i])*
        (Z_tilde[i,k_tilde,drop = F])%*%
        (P_tilde[k_tilde,,drop = F])*
        exp(gamma_tilde[k_tilde] +
              X_tilde[k_tilde,,drop = F]%*%B)
    }
  }
  return(mu_deriv)
}
