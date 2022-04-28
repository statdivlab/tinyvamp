mean_jac <- function(varying_df,
                     fixed_df,
                     X,
                     Z,
                     X_tilde,
                     Z_tilde,
                     Z_tilde_gamma_cols,
                     Z_tilde_list = NULL){



  params <- dataframes_to_parameters(fixed_df,
                                     varying_df)

  if(!is.null(params$alpha_tilde)){
    Z_tilde <- construct_Z_tilde(Z_tilde_list,
                                 params$alpha_tilde)
  }

  fixed_status <- fixed_df
  if(nrow(fixed_status)>0){
    fixed_status$value <- 0
  }
  varying_status <- varying_df
  varying_status$value <- 1

  param_status <- dataframes_to_parameters(fixed_status,
                                           varying_status)

  n <- nrow(X)
  J <- ncol(params$P)
  n_varying <- nrow(varying_df)

  jacobian <- matrix(nrow = n*J,
                     ncol = n_varying)
  counter <- 0
  for(i in 1:n){
    for(j in 1:J){
      counter <- counter + 1
      jacobian[counter,] <-
        par_to_jacobian_row(params,
                            param_status,
                            i,
                            j,
                            X,
                            Z,
                            X_tilde,
                            Z_tilde,
                            Z_tilde_gamma_cols)

    }
  }
  return(jacobian)
}
