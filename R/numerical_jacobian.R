

numerical_jacobian <- function(varying_lr_df,
                               fixed_df,
                               varying_df,
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
                               Z_tilde_list = NULL
                               ){
  n <- nrow(Z)
  J <- ncol(P)
  npar <- nrow(varying_lr_df)
  mean_func <- function(x,index,i,j){
    temp_lr <- varying_lr_df
    temp_lr$value[index] <- x
    temp_params <- dataframes_to_parameters(fixed_df,
                                            lr_to_ra(fixed_df,
                                                     temp_lr,
                                                     varying_df))
    return(with(temp_params,meaninate(gammas,
                                      B,
                                      X,
                                      Z,
                                      P,
                                      X_tilde,
                                      Z_tilde = NULL,
                                      Z_tilde_gamma_cols,
                                      P_tilde,
                                      gamma_tilde,
                                      alpha_tilde = alpha_tilde,
                                      Z_tilde_list = Z_tilde_list,
                                      return_separate = FALSE,
                                      exclude_gammas = FALSE)[i,j]))
  }


num_jacob <- matrix(0,nrow = n*J, ncol = npar)

for(i in 1:n){
  print(paste("i = ", i, sep = "", collapse = ""))
  for(j in 1:J){
    print(paste("j = ", j, sep = "", collapse = ""))
    for(parindex in 1:npar){
      num_jacob[(i - 1)*J + j, parindex] <- numDeriv::grad(
        function(x) mean_func(x,index = parindex,i = i,j = j),
        varying_lr_df$value[parindex])

    }
  }
}

return(num_jacob)



}
