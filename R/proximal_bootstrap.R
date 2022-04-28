#
#
# proximal_bootstrap <- function(
#   n, #number of observations
#   alpha_n, #bootstrap scaling term
#   varying_df, #dataframe containing fitted parameter values
#   fixed_df, #dataframe containing fixed parameter values
#   W, #observations
#   X, #sample efficiency design
#   Z, #sample specimen design
#   Z_tilde = NULL,
#   Z_tilde_gamma_cols,
#   gammas_fixed_indices,
#   P_fixed_indices,
#   B_fixed_indices,
#   X_tilde,
#   P_tilde_fixed_indices,
#   gamma_tilde_fixed_indices,
#   Z_tilde_list = NULL,
#   constraint_tolerance = 1e-10,
#   hessian_regularization = 0.01,
#   criterion = "Poisson",
#   gmm_inv_wts = NULL,
#   nu = 1, #starting lagrangian penalty
#   mu = 1, #starting augmented lagrangian penalty
#   ){
#
#   params <- dataframes_to_parameters(fixed_df,
#                                      varying_df)
#
#   jacobian <- mean_jac_lr_faster(fixed_df = fixed_df,
#                                  varying_lr_df = varying_lr_df,
#                                  varying_df = varying_df,
#                                  which_k_p = which_k_p,
#                                  which_k_p_tilde = which_k_p_tilde,
#                                  which_B_rows = which_B_rows,
#                                  which_B_keep = which_B_keep,
#                                  which_gammas = which_gammas,
#                                  which_gamma_tilde = which_gamma_tilde,
#                                  Ak_list = Ak_list,
#                                  A_tilde_k_list = A_tilde_k_list,
#                                  fixed_P_multipliers = fixed_P_multipliers,
#                                  fixed_P_tilde_multipliers = fixed_P_tilde_multipliers,
#                                  params = params,
#                                  X = X,
#                                  Z = Z,
#                                  K = K,
#                                  K_tilde = K_tilde,
#                                  X_tilde = X_tilde,
#                                  Z_tilde = Z_tilde,
#                                  Z_tilde_gamma_cols =Z_tilde_gamma_cols,
#                                  Z_tilde_list = Z_tilde_list,
#                                  sparse = TRUE,
#                                  proportion_scale = TRUE,
#                                  P_fixed_indices = P_fixed_indices,
#                                  P_tilde_fixed_indices = P_tilde_fixed_indices)
#
#   means <- meaninate(gammas = params$gammas,
#                      B = params$B,
#                      X = X,
#                      Z = Z,
#                      P = params$P,
#                      X_tilde = X_tilde,
#                      Z_tilde = Z_tilde,
#                      Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                      P_tilde = params$P_tilde,
#                      gamma_tilde = params$gamma_tilde,
#                      alpha_tilde = params$alpha_tilde,
#                      Z_tilde_list = Z_tilde_list)
#
#   W_long <- lapply(1:n,function(i) as.numeric(W[i,]))
#   W_long <- do.call(c,W_long)
#   means_long <- lapply(1:n, function(i) as.numeric(means[i,]))
#   means_long <- do.call(c,means_long)
#
#   ### only for Poisson:
#   if(criterion == "Poisson"){
#     V <- diag(1/means_long)
#   } else{
#     V <- diag(rep(1, length(means_long)))
#   }
#
#   if(!is.null(gmm_inv_wts)){
#     diag(V) <- diag(V)/gmm_inv_wts
#   }
#
#   dli_dtheta <-
#     lapply(1:n, function(i)
#            -Matrix::crossprod(jacobian[1:J + (i - 1)*J,],
#                               V[1:J + (i - 1)*J,1:J + (i - 1)*J])%*%(
#     Matrix::Matrix(W_long[1:J + (i - 1)*J] - means_long[1:J + (i - 1)*J],
#                    ncol = 1)))
#
#   dli_dtheta <- do.call(cbind,dli_dtheta)
#   dli_dtheta <- dli_dtheta/n
#
#   lgrad <- apply(dli_dtheta,1,sum)
#
#   V2 <- as(sqrt(abs(V)),"sparseMatrix")
#   pre_info <- V2%*%jacobian
#   H_n <-  Matrix::crossprod(pre_info)
#
#
#   boot_weights <- rexp(n)
#   boot_weights <- boot_weights/sum(boot_weights)
#
#
#   lgrad_star <- apply(dli_dtheta%*%diag(boot_weights),1,sum)
#
#   diff_dls <- alpha_n*sqrt(n)*as.numeric(
#     (apply(dli_star_dtheta,1,sum) - apply(dli_dtheta,1,sum)))
#   diff_dls <- matrix(diff_dls,nrow = 1)
#
#   prox_crit <- function(x){
#     x_deviation <- matrix(x - varying_df$value,ncol = 1)
#
#     return(as.numeric(diff_dls%*%x_deviation +
#       0.5*t(x_deviation)%*%H_n%*%x_deviation))
#
#   }
#
#   K <- length(unique(varying_df$k[varying_df$param == "P"]))
#   K_tilde <- length(unique(varying_df$k[varying_df$param == "P_tilde"]))
#   aug_lag_params <- matrix(1,ncol = 2, nrow = n_simplex_constraints)
#
#   varying_p_k <- unique(varying_df$k[varying_df$param == "P"])
#   simplex_matrix_P <-
#     lapply(varying_p_k,
#            function(d) as.numeric(
#              (varying_df$k == d) & (varying_df$param == "P"))
#     )
#
#   simplex_matrix_P <- do.call(rbind,simplex_matrix_P)
#
#   varying_p_tilde_k <- unique(varying_df$k[varying_df$param == "P_tilde"])
#   simplex_matrix_P_tilde <-
#     lapply(varying_p_tilde_k,
#            function(d) as.numeric(
#              (varying_df$k == d) & (varying_df$param == "P_tilde"))
#     )
#
#   simplex_matrix_P_tilde <- do.call(rbind,simplex_matrix_P_tilde)
#
#   for(k in 1:maxit){
#     counter <- counter + 1
#     # print(counter)
#
#     new_x <- optim(varying_df$value,
#                    )
#
#     # print("checking constraints")
#
#     V <- abs(sum(x) - 1)
#
#     satisfied <- V< constraint_tolerance
#
#
#
#     if(V < constraint_tolerance){
#       return(x)
#     }
#
#     if(V < 0.25*previous_V){
#       nu <- nu + 2*mu*V
#     } else{
#       mu <- 2*mu
#     }
#     # print(mu)
#     # print(nu)
#     previous_V <- V
#
#   }
#
#
#
#
# }
