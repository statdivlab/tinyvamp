par_to_jacobian_row <- function(params,
                                param_status,
                                i,
                                j,
                                X,
                                Z,
                                X_tilde,
                                Z_tilde,
                                Z_tilde_gamma_cols){
  J <- ncol(params$P)

  ### jacobian row in P
  which_P_rows <- which(apply(param_status$P,1, max) ==1)
  P_jac_row <- numeric(J*length(which_P_rows))

  if(length(which_P_rows) >0){
    for(row_index in 1:length(which_P_rows)){
      P_jac_row[(row_index - 1)*J + j] <- mu_d_P(i,
                                                 j,
                                                 m = which_P_rows[row_index],
                                                 gammas = params$gammas,
                                                 B = params$B,
                                                 X = X,
                                                 Z = Z,
                                                 P = params$P)
    }
  }

  ### jacobian row in P_tilde
  which_P_tilde_rows <- which(apply(param_status$P_tilde,1, max) ==1)
  P_tilde_jac_row <- numeric(J*length(which_P_tilde_rows))

  if(length(which_P_tilde_rows)>0){
    for(row_index in 1:length(which_P_tilde_rows)){
      P_tilde_jac_row[(row_index - 1)*J + j] <- mu_d_P_tilde(i,
                                                             j,
                                                             k_tilde = which_P_tilde_rows[row_index],
                                                             gammas = params$gammas,
                                                             B = params$B,
                                                             X = X,
                                                             Z = Z,
                                                             P = params$P,
                                                             X_tilde = X_tilde,
                                                             Z_tilde = Z_tilde,
                                                             Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                                             P_tilde = params$P_tilde,
                                                             gamma_tilde = params$gamma_tilde)
    }
  }

  ### jacobian row in B
  which_B_rows <- which(apply(param_status$B,1, max) ==1)
  B_jac_row <- numeric((J - 1)*length(which_B_rows))

  if(length(which_B_rows)>0){
    if(j < J){ #row of jacobian is zero if j = J
      for(row_index in 1:length(which_B_rows)){
        B_jac_row[(row_index - 1)*(J - 1) + j] <- mu_d_beta(i,
                                                            j,
                                                            q = which_B_rows[row_index],
                                                            gammas = params$gammas,
                                                            B = params$B,
                                                            X = X,
                                                            Z = Z,
                                                            P = params$P,
                                                            X_tilde = X_tilde,
                                                            Z_tilde = Z_tilde,
                                                            Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                                            P_tilde = params$P_tilde,
                                                            gamma_tilde = params$gamma_tilde)
      }
    }
  }

  ### gamma jacobian row
  which_gammas <- which(apply(param_status$gammas,1, max) ==1)
  gammas_jac_row <- numeric(sum(param_status$gammas))

  if(length(which_gammas)>0){
    gammas_jac_row[i] <- mu_d_gamma(i,
                                    j,
                                    gammas = params$gammas,
                                    B = params$B,
                                    X = X,
                                    Z = Z,
                                    P = params$P,
                                    X_tilde = X_tilde,
                                    Z_tilde = Z_tilde,
                                    Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                    P_tilde = params$P_tilde,
                                    gamma_tilde = params$gamma_tilde)
  }
  gammas_jac_row <- gammas_jac_row[which_gammas] #???

  ### gamma tilde jacobian row
  which_gamma_tilde <- which(apply(param_status$gamma_tilde,1, max) ==1)
  gamma_tilde_jac_row <- numeric(sum(param_status$gamma_tilde))

  if(length(which_gamma_tilde)>0){
    for(row_index in 1:length(which_gamma_tilde)){
      gamma_tilde_jac_row[row_index] <- mu_d_gamma_tilde(i,
                                                         j,
                                                         k_tilde = which_gamma_tilde[row_index],
                                                         gammas = params$gammas,
                                                         B = params$B,
                                                         X = X,
                                                         Z = Z,
                                                         P = params$P,
                                                         X_tilde = X_tilde,
                                                         Z_tilde = Z_tilde,
                                                         Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                                         P_tilde = params$P_tilde,
                                                         gamma_tilde = params$gamma_tilde)
    }
  }

  return(c(P_jac_row,
           P_tilde_jac_row,
           B_jac_row,
           gammas_jac_row,
           gamma_tilde_jac_row))

}
