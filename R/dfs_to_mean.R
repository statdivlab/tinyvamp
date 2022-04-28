

dfs_to_mean <- function(varying_df,
                        varying_lr_df = NULL,
                        fixed_df,
                        X,
                        Z,
                        X_tilde,
                        Z_tilde,
                        Z_tilde_gamma_cols){


  if(!is.null(varying_lr_df)){
    varying_df <- lr_to_ra(fixed_df,varying_lr_df,
                           varying_df)
  }

  params <- dataframes_to_parameters(fixed_df,
                                     varying_df)

  return(meaninate(gammas = params$gammas,
            B = B,
            X = X,
            Z = Z,
            P = params$P,
            X_tilde = X_tilde,
            Z_tilde = Z_tilde,
            Z_tilde_gamma_cols = Z_tilde_gamma_cols,
            P_tilde = params$P_tilde,
            gamma_tilde = params$gamma_tilde))
}


