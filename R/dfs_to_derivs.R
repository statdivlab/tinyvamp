

dfs_to_derivs <- function(varying_df,
              varying_lr_df = NULL,
              fixed_df,
              X,
              Z,
              X_tilde,
              Z_tilde,
              Z_tilde_gamma_cols,
              criterion = "Poisson",
              gmm_inv_wts = NULL){

  if(!is.null(varying_lr_df)){
    varying_df <- lr_to_ra(fixed_df,
                           varying_lr_df,
                           varying_df)
  }


  K <- max(c(varying_df$k[varying_df$param == "P"],
             fixed_df$k[fixed_df$param == "P"]))

  fixed_P_multipliers <- sapply(1:K, function(k)
    1 - sum(fixed_df$value[fixed_df$param == "P"&
                             fixed_df$k ==k]))

  K_tilde <- max(c(varying_df$k[varying_df$param == "P_tilde"],
                   fixed_df$k[fixed_df$param == "P_tilde"]))

  fixed_P_tilde_multipliers <- sapply(1:K_tilde, function(k)
    1 - sum(fixed_df$value[fixed_df$param == "P_tilde"&
                             fixed_df$k ==k]))

  Ak_list <- get_Ak_list(fixed_df,
                         varying_df,
                         varying_lr_df)

  A_tilde_k_list <- get_A_tilde_k_list(fixed_df,
                                       varying_df,
                                       varying_lr_df)

  # message("created Ak_list and A_tilde_k_list")

  #calculate at outset of optimization  pass as argument to mean_jac_lr, etc.
  which_k_p <- sapply(1:K, function(k) ifelse(is.null(Ak_list[[k]]),
                                              NA, k))

  which_k_p <- which_k_p[!is.na(which_k_p)]

  #calculate at outset of optimization  pass as argument to mean_jac_lr, etc.
  which_k_p_tilde <- sapply(1:K_tilde,
                            function(k) ifelse(
                              is.null(A_tilde_k_list[[k]]),
                              NA,k
                            ))

  which_k_p_tilde <- which_k_p_tilde[!is.na(which_k_p_tilde)]

  # message("saved which_k_p and which_k_p_tilde")

  #calculate at outset of optimization
  which_B_rows <- unique(varying_df$k[varying_df$param == "B"])
  which_B_rows <- which_B_rows[order(which_B_rows)]

  #calculate at outset of optimization
  which_B_keep <- lapply(which_B_rows,
                         function(k) sapply(1:(J - 1),
                                            function(j)
                                              j %in% varying_lr_df$j[
                                                varying_lr_df$param == "B" &
                                                  varying_lr_df$k == k]
                         ))
  which_B_keep <- do.call(rbind,which_B_keep)

  # message("saved which_B_keep")

  which_gammas <- unique(varying_df$k[varying_df$param == "gamma"])

  which_gamma_tilde <- unique(varying_df$k[varying_df$param == "gamma_tilde"])

  which_unconstrained <- varying_lr_df$param %in% c("B","gamma","gamma_tilde")
  which_rho <- varying_lr_df$param %in% c("rho")
  which_rho_tilde <- varying_lr_df$param %in% c("rho_tilde")
  npar <- nrow(varying_lr_df)

  return(deriv_criterion_lr(W = W,
                     X = X,
                     Z = Z,
                     which_k_p = which_k_p,
                     which_k_p_tilde = which_k_p_tilde,
                     fixed_P_multipliers = fixed_P_multipliers,
                     fixed_P_tilde_multipliers = fixed_P_tilde_multipliers,
                     which_B_keep = which_B_keep,
                     which_B_rows = which_B_rows,
                     which_gammas = which_gammas,
                     which_gamma_tilde = which_gamma_tilde,
                     Z_tilde = Z_tilde,
                     Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                     X_tilde = X_tilde,
                     Ak_list = Ak_list,
                     A_tilde_k_list = A_tilde_k_list,
                     fixed_df = fixed_df,
                     varying_df = varying_df,
                     varying_lr_df = varying_lr_df,
                     K = K,
                     K_tilde = K_tilde,
                     barrier_t = 1,
                     criterion = "Poisson",
                     lr_scale = TRUE,
                     include_log_penalty_derivatives = TRUE,
                     return_info = FALSE,
                     wts = NULL,
                     gmm_inv_wts = gmm_inv_wts))
}
