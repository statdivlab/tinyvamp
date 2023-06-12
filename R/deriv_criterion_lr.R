
deriv_criterion_lr <- function(W,
                               X,
                               Z,
                               Z_tilde = NULL,
                               Z_tilde_gamma_cols,
                               Z_tilde_list = NULL,
                               alpha_tilde = NULL,
                               X_tilde,
                               which_k_p,
                               which_k_p_tilde,
                               Ak_list,
                               A_tilde_k_list,
                               fixed_P_multipliers,
                               fixed_P_tilde_multipliers,
                               which_B_keep,
                               which_B_rows,
                               which_gammas,
                               which_gamma_tilde,
                               fixed_df,
                               varying_df,
                               varying_lr_df = NULL,
                               K,
                               K_tilde,
                               barrier_t = NULL,
                               criterion = "Poisson",
                               lr_scale = TRUE,
                               include_log_penalty_derivatives = TRUE,
                               return_info = FALSE,
                               wts = NULL,
                               gmm_inv_wts = NULL) {
  
  n <- nrow(W)
  J <- ncol(W)
  
  if (lr_scale) {
    varying_df <- lr_to_ra(fixed_df,
                           varying_lr_df,
                           varying_df)
  }
  
  # if (is.null(varying_lr_df)) {
  #   varying_lr_df <-   ra_to_lr(varying_df)
  # }
  
  params <- dataframes_to_parameters(fixed_df,
                                     varying_df)
  
  log_P <- lapply(1:K,
                  function(k) {
                    if (k %in% which_k_p) {
                      rho <- c(varying_lr_df$value[varying_lr_df$param == "rho" &
                                                     varying_lr_df$k == k],0);
                      rho_cent <- rho -
                        sum_of_logs(rho) +
                        log(fixed_P_multipliers[k])
                    } else{
                      rho_cent <- NULL
                    }
                    return(rho_cent)})
  
  log_P_tilde <- lapply(1:K_tilde,
                        function(k) {
                          if (k %in% which_k_p_tilde) {
                            rho <-
                              c(varying_lr_df$value[varying_lr_df$param == "rho_tilde" &
                                                      varying_lr_df$k == k],0);
                            rho_cent <- rho -
                              sum_of_logs(rho) +
                              log(fixed_P_tilde_multipliers[k])
                          } else{
                            rho_cent <- NULL
                          }
                          return(rho_cent)})
  
  params$log_P <- log_P
  params$log_P_tilde <- log_P_tilde
  
  jacobian <- mean_jac_lr_faster(fixed_df = fixed_df,
                                 varying_lr_df = varying_lr_df,
                                 varying_df = varying_df,
                                 which_k_p = which_k_p,
                                 which_k_p_tilde = which_k_p_tilde,
                                 which_B_rows = which_B_rows,
                                 which_B_keep = which_B_keep,
                                 which_gammas = which_gammas,
                                 which_gamma_tilde = which_gamma_tilde,
                                 Ak_list = Ak_list,
                                 A_tilde_k_list = A_tilde_k_list,
                                 fixed_P_multipliers = fixed_P_multipliers,
                                 fixed_P_tilde_multipliers = fixed_P_tilde_multipliers,
                                 params = params,
                                 X = X,
                                 Z = Z,
                                 K = K,
                                 K_tilde = K_tilde,
                                 X_tilde = X_tilde,
                                 Z_tilde = Z_tilde,
                                 Z_tilde_gamma_cols =Z_tilde_gamma_cols,
                                 Z_tilde_list = Z_tilde_list,
                                 sparse = TRUE)
  
  ############### inline test ################
  # num_jacob <- numerical_jacobian(varying_lr_df,
  #                                fixed_df,
  #                                varying_df,
  #                                gammas,
  #                                B,
  #                                X,
  #                                Z,
  #                                P,
  #                                X_tilde,
  #                                Z_tilde = Z_tilde,
  #                                Z_tilde_gamma_cols,
  #                                P_tilde,
  #                                gamma_tilde,
  #                                alpha_tilde = alpha_tilde,
  #                                Z_tilde_list = Z_tilde_list
  # )
  #
  #
  # plot(asinh(num_jacob), asinh(as.matrix(jacobian)),pch = ".")
  #
  # plot(apply(jacobian - num_jacob,2,function(x) asinh(max(abs(x)))))
  
  #################            ################
  
  
  
  means <- meaninate(gammas  = params$gammas,
                     B = params$B,
                     X = X,
                     Z = Z,
                     P = params$P,
                     X_tilde = X_tilde,
                     Z_tilde = Z_tilde,
                     Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                     P_tilde = params$P_tilde,
                     gamma_tilde = params$gamma_tilde,
                     Z_tilde_list = Z_tilde_list,
                     alpha_tilde = params$alpha_tilde,
                     return_separate = FALSE)
  
  W_long <- lapply(1:n,function(i) as.numeric(W[i,]))
  W_long <- do.call(c,W_long)
  means_long <- lapply(1:n, function(i) as.numeric(means[i,]))
  means_long <- do.call(c,means_long)
  
  ### only for Poisson:
  if (criterion == "Poisson") {
    V <- diag(1/means_long)
  } else{
    V <- diag(rep(1, length(means_long)))
  }
  
  if (!is.null(wts)) {
    diag(V) <- diag(V)*wts
  }
  if (!is.null(gmm_inv_wts)) {
    diag(V) <- diag(V)/gmm_inv_wts
  }
  # V <- as.matrix.csr(V,nrow = nrow(V))
  V <- as(V, "sparseMatrix")
  if (return_info) {
    V2 <- as(sqrt(abs(V)),"sparseMatrix")
  }
  # V2 <- as.matrix.csr(sqrt( diag(1/means_long)), nrow = nrow(V))
  # V_gm <- exp(mean(log(diag(V)[diag(V)!=0])))
  # V <- V/V_gm
  
  # jacobian_gm <- exp(mean(log(abs(jacobian)[jacobian!=0])))
  # jacobian <- jacobian/jacobian_gm
  
  
  
  crit_grad <- -Matrix::crossprod(jacobian, V) %*% (
    Matrix::Matrix(W_long - means_long,ncol = 1))
  
  # step_dir <- qr.solve(info,-crit_grad)
  
  # crit_grad <- crit_grad*V_gm*jacobian_gm
  if (return_info) {
    pre_info <- V2 %*% jacobian
    
    # crit_info <- t(pre_info) %*% pre_info
    if (is.null(wts)) {
      crit_info <- Matrix::crossprod(pre_info)
    } else{
      crit_info <- Matrix::crossprod(pre_info,diag(sign(wts)) %*% pre_info)
    }
  }
  
  ### log_penalties
  if (include_log_penalty_derivatives) {
    lp_grad <- rep(0, sum(!(varying_lr_df$param %in% c("rho","rho_tilde",
                                                      "alpha_tilde"))))
    if (return_info) {
      lp_hess <- matrix(0,  ncol = length(lp_grad),nrow = length(lp_grad))
    }
    which_rho_k <- unique(varying_lr_df$k[varying_lr_df$param == "rho"])
    which_rho_tilde_k <-
      unique(varying_lr_df$k[varying_lr_df$param == "rho_tilde"])
    
    for(k in which_rho_k) {
      temp_rho <- varying_lr_df$value[varying_lr_df$param == "rho" &
                                        varying_lr_df$k == k]
      
      #calculate and tack on grad of penalty
      lp_grad <- c(lp_grad, log_penalty_grad(temp_rho,barrier_t))
      
      if (return_info) {
        #calculate hessian of penalty
        new_hess_piece <- log_penalty_hess(temp_rho,barrier_t)
        
        #tack new block of penalty hessian onto existing penalty hessian
        hess_rows <- nrow(lp_hess)
        new_rows <- nrow(new_hess_piece)
        new_hess_piece <- cbind(matrix(0, nrow = new_rows,
                                       ncol = hess_rows),
                                new_hess_piece)
        lp_hess <- cbind(lp_hess,
                         matrix(0, nrow = hess_rows,
                                ncol = new_rows))
        lp_hess <- rbind(lp_hess,new_hess_piece)
      }
    }
    
    for(k in which_rho_tilde_k) {
      temp_rho_tilde <- varying_lr_df$value[varying_lr_df$param == "rho_tilde" &
                                              varying_lr_df$k == k]
      
      #calculate and tack on grad of penalty
      lp_grad <-  c(lp_grad, log_penalty_grad(temp_rho_tilde,barrier_t))
      
      if (return_info) {
        #calculate hessian of penalty
        new_hess_piece <- log_penalty_hess(temp_rho_tilde,barrier_t)
        
        #tack new block of penalty hessian onto existing penalty hessian
        hess_rows <- nrow(lp_hess)
        new_rows <- nrow(new_hess_piece)
        new_hess_piece <- cbind(matrix(0, nrow = new_rows,
                                       ncol = hess_rows),
                                new_hess_piece)
        lp_hess <- cbind(lp_hess,
                         matrix(0, nrow = hess_rows,
                                ncol = new_rows))
        lp_hess <- rbind(lp_hess,new_hess_piece)
      }
    }
    
  }
  
  #alpha_tilde stored *after* rho, rho_tilde in df -- have to tack on
  #zeroes for it
  n_alpha_tilde <- sum(varying_lr_df$param=="alpha_tilde")
  if (n_alpha_tilde>0) {
    if (include_log_penalty_derivatives) {
      lp_grad <- c(lp_grad,
                   rep(0, n_alpha_tilde))}
    
    if (return_info) {
      
      lp_hess <- cbind(lp_hess,matrix(0, ncol = n_alpha_tilde,
                                      nrow = nrow(lp_hess)))
      
      lp_hess <- rbind(lp_hess,matrix(0, nrow = n_alpha_tilde,
                                      ncol = ncol(lp_hess)))
    }
    
  }
  if (!return_info) {
    crit_info <- NULL
    lp_hess <- NULL
  }
  
  if (!include_log_penalty_derivatives) {
    return(list("grad" = crit_grad,
                "info" = crit_info))
  } else{
    return(list("grad" = crit_grad + lp_grad,
                "info" = crit_info + lp_hess))
  }
  
}
