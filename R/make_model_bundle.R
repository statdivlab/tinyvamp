#
# make_model_bundle <- function(W,
#                               X,
#                               Z,
#                               Z_tilde,
#                               Z_tilde_gamma_cols,
#                               gammas,
#                               gammas_fixed_indices,
#                               P,
#                               P_fixed_indices,
#                               B,
#                               B_fixed_indices,
#                               X_tilde,
#                               P_tilde,
#                               P_tilde_fixed_indices,
#                               gamma_tilde,
#                               gamma_tilde_fixed_indices,
#                               criterion = "Poisson"){
#
#
#   model_bundle <- list()
#   model_bundle$W <- W
#   model_bundle$designs <- list(Z = Z,
#                                Z_tilde = Z_tilde,
#                                Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                                X_tilde = X_tilde)
#   model_bundle$parameters <-
#     list(gammas = gammas,
#          gammas_fixed_indices = gammas_fixed_indices,
#          P = P,
#          P_fixed_indices = P_fixed_indices,
#          B = B,
#          B_fixed_indices = B_fixed_indices,
#          P_tilde = P_tilde,
#          P_tilde_fixed_indices = P_tilde_fixed_indices,
#          gamma_tilde = gamma_tilde,
#          gamma_tilde_fixed_indices = gamma_tilde_fixed_indices)
#   model_bundle$criterion <- criterion
#
#
#   model_bundle$n <- nrow(W)
#   model_bundle$J <- ncol(W)
#   model_bundle$K <- NA
#   model_bundle$K_tilde <- NA
#
#   model_bundle$optimization_info <- list()
#   model_bundle$optimization_info$nsteps <-
#     ceiling(log(max_barrier/barrier_t)/log(barrier_scale)) + 1
#
#   model_bundle$optimization_info$tolerances <-
#     exp(seq(log(initial_conv_tol),log(final_conv_tol), length.out = nsteps))
#
#   model_bundle$parameter_dfs <-  parameters_to_dataframes(P,
#                               P_fixed_indices,
#                               P_tilde,
#                               P_tilde_fixed_indices,
#                               B,
#                               B_fixed_indices,
#                               gammas,
#                               gammas_fixed_indices,
#                               gamma_tilde,
#                               gamma_tilde_fixed_indices)
#
#
#
#
#   # message("created parameter dfs")
#
#   model_bundle$K <-
#     with(model_bundle$parameter_dfs,
#          max(c(varying_df$k[varying_df$param == "P"],
#              fixed_df$k[fixed_df$param == "P"])))
#
#   model_bundle$K_tilde <-
#     with(model_bundle$parameter_dfs,
#     max(c(varying_df$k[varying_df$param == "P_tilde"],
#           fixed_df$k[fixed_df$param == "P_tilde"]))
#     )
#
#   model_bundle$rho_P_conversion <-list()
#   model_bundle$rho_P_conversion$fixed_P_multipliers <-
#   sapply(1:model_bundle$K, function(k)
#     with(model_bundle$parameter_dfs,1 - sum(fixed_df$value[fixed_df$param == "P"&
#                              fixed_df$k ==k])))
#
#   model_bundle$rho_P_conversion$fixed_P_tilde_multipliers <-
#     sapply(1:model_bundle$K_tilde, function(k)
#     with(model_bundle$parameter_dfs, 1 - sum(fixed_df$value[fixed_df$param == "P_tilde"&
#                              fixed_df$k ==k])))
#
#
#   # create matrices to track rho-P and rho_tilde-P_tilde relationships
#   model_bundle$rho_P_conversion$Ak_list <-
#     with(model_bundle$parameter_dfs,
#          get_Ak_list(fixed_df,
#                          varying_df,
#                          varying_lr_df))
#
#   A_tilde_k_list <- get_A_tilde_k_list(fixed_df,
#                                        varying_df,
#                                        varying_lr_df)
#
#   # message("created Ak_list and A_tilde_k_list")
#
#   #calculate at outset of optimization  pass as argument to mean_jac_lr, etc.
#   which_k_p <- sapply(1:K, function(k) ifelse(is.null(Ak_list[[k]]),
#                                               NA, k))
#
#   which_k_p <- which_k_p[!is.na(which_k_p)]
#
#   #calculate at outset of optimization  pass as argument to mean_jac_lr, etc.
#   which_k_p_tilde <- sapply(1:K_tilde,
#                             function(k) ifelse(
#                               is.null(A_tilde_k_list[[k]]),
#                               NA,k
#                             ))
#
#   which_k_p_tilde <- which_k_p_tilde[!is.na(which_k_p_tilde)]
#
#   # message("saved which_k_p and which_k_p_tilde")
#
#   #calculate at outset of optimization
#   which_B_rows <- unique(varying_df$k[varying_df$param == "B"])
#   which_B_rows <- which_B_rows[order(which_B_rows)]
#
#   #calculate at outset of optimization
#   which_B_keep <- lapply(which_B_rows,
#                          function(k) sapply(1:(J - 1),
#                                             function(j)
#                                               j %in% varying_lr_df$j[
#                                                 varying_lr_df$param == "B" &
#                                                   varying_lr_df$k == k]
#                          ))
#   which_B_keep <- do.call(rbind,which_B_keep)
#
#   # message("saved which_B_keep")
#
#   which_gammas <- unique(varying_df$k[varying_df$param == "gamma"])
#
#   which_gamma_tilde <- unique(varying_df$k[varying_df$param == "gamma_tilde"])
#
#   which_unconstrained <- varying_lr_df$param %in% c("B","gamma","gamma_tilde")
#   which_rho <- varying_lr_df$param %in% c("rho")
#   which_rho_tilde <- varying_lr_df$param %in% c("rho_tilde")
#   npar <- nrow(varying_lr_df)
# }
