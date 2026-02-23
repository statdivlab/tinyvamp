

calculate_U <- function(W,
                        full_model,
                        null_model,
                        training_indicator,
                        log_scale = TRUE){
  W_train <- W[training_indicator,,drop = FALSE]
  W_test <- W[!training_indicator,,drop = FALSE]
  
  ### may need to set some entries of P to be fixed if no observations in
  ### training set correspond to them
  
  if(is.null(full_model$Z_tilde)){
    full_model_Z_tilde <- NULL
    full_model_Z_tilde_test <- NULL
  } else{
    full_model_Z_tilde <- full_model$Z_tilde[training_indicator,,drop = FALSE]
    full_model_Z_tilde_test <-
      full_model$Z_tilde[!training_indicator,,drop = FALSE]
  }
  if(!is.null(full_model$Z_tilde_list)){
    full_model_Z_tilde_list <-
      lapply(1:length(full_model$Z_tilde_list),
             function(k)
               full_model$Z_tilde_list[[k]][training_indicator,,drop = FALSE])
    full_model_Z_tilde_list_test <-
      lapply(1:length(full_model$Z_tilde_list),
             function(k)
               full_model$Z_tilde_list[[k]][!training_indicator,,drop = FALSE])
  } else{
    full_model_Z_tilde_list <- full_model$Z_tilde_list
    full_model_Z_tilde_list_test <- full_model$Z_tilde_list
  }
  
  
  if(is.null(null_model$Z_tilde)){
    null_model_Z_tilde <- NULL
    null_model_Z_tilde_test <- NULL
  } else{
    null_model_Z_tilde <- null_model$Z_tilde[training_indicator,,drop = FALSE]
    null_model_Z_tilde_test <-
      null_model$Z_tilde[!training_indicator,,drop = FALSE]
  }
  if(!is.null(null_model$Z_tilde_list)){
    null_model_Z_tilde_list <-
      lapply(1:length(null_model$Z_tilde_list),
             function(k)
               null_model$Z_tilde_list[[k]][training_indicator,,drop = FALSE])
    null_model_Z_tilde_list_test <-
      lapply(1:length(null_model$Z_tilde_list),
             function(k)
               null_model$Z_tilde_list[[k]][!training_indicator,,drop = FALSE])
  } else{
    null_model_Z_tilde_list <- null_model$Z_tilde_list
    null_model_Z_tilde_list_test <- null_model$Z_tilde_list
  }
  
  full_training <-
    estimate_parameters(W = W_train,
                        X = full_model$X[training_indicator,,drop = FALSE],
                        Z = full_model$Z[training_indicator,,drop = FALSE],
                        Z_tilde = full_model_Z_tilde,
                        Z_tilde_gamma_cols = full_model$Z_tilde_gamma_cols,
                        gammas = full_model$gammas[training_indicator],
                        gammas_fixed_indices =
                          full_model$gammas_fixed_indices[training_indicator],
                        P = full_model$P,
                        P_fixed_indices = full_model$P_fixed_indices,
                        B = full_model$B,
                        B_fixed_indices = full_model$B_fixed_indices,
                        X_tilde = full_model$X_tilde,
                        P_tilde = full_model$P_tilde,
                        P_tilde_fixed_indices = full_model$P_tilde_fixed_indices,
                        gamma_tilde = full_model$gamma_tilde,
                        gamma_tilde_fixed_indices =
                          full_model$gamma_tilde_fixed_indices,
                        alpha_tilde = full_model$alpha_tilde,
                        Z_tilde_list = full_model_Z_tilde_list)
  
  null_training <-
    estimate_parameters(W = W_train,
                        X = null_model$X[training_indicator,,drop = FALSE],
                        Z = null_model$Z[training_indicator,,drop = FALSE],
                        Z_tilde = null_model_Z_tilde,
                        Z_tilde_gamma_cols = null_model$Z_tilde_gamma_cols,
                        gammas = null_model$gammas[training_indicator],
                        gammas_fixed_indices =
                          null_model$gammas_fixed_indices[training_indicator],
                        P = null_model$P,
                        P_fixed_indices = null_model$P_fixed_indices,
                        B = null_model$B,
                        B_fixed_indices = null_model$B_fixed_indices,
                        X_tilde = null_model$X_tilde,
                        P_tilde = null_model$P_tilde,
                        P_tilde_fixed_indices = null_model$P_tilde_fixed_indices,
                        gamma_tilde = null_model$gamma_tilde,
                        gamma_tilde_fixed_indices =
                          null_model$gamma_tilde_fixed_indices,
                        alpha_tilde = null_model$alpha_tilde,
                        Z_tilde_list = null_model_Z_tilde_list)
  
  full_starter_means <- meaninate(gammas = rep(0, sum(!training_indicator)),
                                  B = full_training$B,
                                  X =
                                    full_model$X[!training_indicator,,drop = FALSE],
                                  Z =
                                    full_model$Z[!training_indicator,,drop = FALSE],
                                  P = full_training$P,
                                  X_tilde = full_model$X_tilde,
                                  Z_tilde = full_model_Z_tilde_test,
                                  Z_tilde_gamma_cols =
                                    full_model$Z_tilde_gamma_cols,
                                  P_tilde = full_training$P_tilde,
                                  gamma_tilde = full_training$gamma_tilde,
                                  alpha_tilde = full_training$alpha_tilde,
                                  Z_tilde_list = full_model_Z_tilde_list_test)
  
  full_profile_gammas <- log(apply(W_test,1,sum)) -
    log(apply(full_starter_means,1,sum))
  
  null_starter_means <- meaninate(gammas = rep(0, sum(!training_indicator)),
                                  B = null_training$B,
                                  X =
                                    null_model$X[!training_indicator,,drop = FALSE],
                                  Z =
                                    null_model$Z[!training_indicator,,drop = FALSE],
                                  P = null_training$P,
                                  X_tilde = null_model$X_tilde,
                                  Z_tilde = null_model_Z_tilde_test,
                                  Z_tilde_gamma_cols =
                                    null_model$Z_tilde_gamma_cols,
                                  P_tilde = null_training$P_tilde,
                                  gamma_tilde = null_training$gamma_tilde,
                                  alpha_tilde = null_training$alpha_tilde,
                                  Z_tilde_list = null_model_Z_tilde_list_test)
  
  null_profile_gammas <- log(apply(W_test,1,sum)) -
    log(apply(null_starter_means,1,sum))
  
  
  full_means_test <- meaninate(gammas = full_profile_gammas,
                               B = full_training$B,
                               X =
                                 full_model$X[!training_indicator,,drop = FALSE],
                               Z =
                                 full_model$Z[!training_indicator,,drop = FALSE],
                               P = full_training$P,
                               X_tilde = full_model$X_tilde,
                               Z_tilde = full_model_Z_tilde_test,
                               Z_tilde_gamma_cols =
                                 full_model$Z_tilde_gamma_cols,
                               P_tilde = full_training$P_tilde,
                               gamma_tilde = full_training$gamma_tilde,
                               alpha_tilde = full_training$alpha_tilde,
                               Z_tilde_list = full_model_Z_tilde_list_test)
  
  null_means_test <- meaninate(gammas = null_profile_gammas,
                               B = null_training$B,
                               X =
                                 null_model$X[!training_indicator,,drop = FALSE],
                               Z =
                                 null_model$Z[!training_indicator,,drop = FALSE],
                               P = null_training$P,
                               X_tilde = null_model$X_tilde,
                               Z_tilde = null_model_Z_tilde_test,
                               Z_tilde_gamma_cols =
                                 null_model$Z_tilde_gamma_cols,
                               P_tilde = null_training$P_tilde,
                               gamma_tilde = null_training$gamma_tilde,
                               alpha_tilde = null_training$alpha_tilde,
                               Z_tilde_list = null_model_Z_tilde_list_test)
  
  l_HA <- poisson_criterion(W_test,full_means_test)
  l_H0 <- poisson_criterion(W_test,null_means_test)
  
  if(log_scale){
    return(l_HA - l_H0)
  } else{
    return(exp(l_HA - l_H0))
  }
  
}
