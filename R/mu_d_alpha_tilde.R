

mu_d_alpha_tilde <- function(i,
                             J,
                             a_tilde,
                             gammas,
                             B,
                             X_tilde,
                             Z_tilde_gamma_cols,
                             alpha_tilde,
                             Z_tilde_list,
                             P_tilde,
                             gamma_tilde){

  Z_tilde_piece <- Z_tilde_list[[a_tilde + 1]]

  if(length(Z_tilde_gamma_cols >0)){
    for(colnum in Z_tilde_gamma_cols){
      Z_tilde_piece[,colnum] <- exp(gammas)* Z_tilde_piece[,colnum]
    }
  }

  return(as.numeric(
    exp(alpha_tilde[a_tilde])*Z_tilde_piece[i,]%*%(
    P_tilde*exp(matrix(gamma_tilde,ncol = 1)%*%matrix(1, nrow = 1, ncol = J) +
                  X_tilde%*%B))
  )
  )




}
