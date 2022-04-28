

# test_that("GMM derivatives in gamma and rho are correct in simple case", {
#
#
#   set.seed(0)
#   W <- matrix(sapply(1:10,function(x) rpois(1,1000)),
#               ncol = 5)
#   X <- matrix(1,ncol = 1, nrow = 2)
#   Z <- matrix(1,nrow = 2, ncol = 1)
#   Z_tilde <- matrix(0,nrow = 2, ncol = 1)
#   Z_tilde_gamma_cols <- 1
#   gammas <- apply(W,1,function(x) log(sum(x)))
#   gammas_fixed_indices <- rep(F,2)
#   P <- matrix(1/5, nrow = 1, ncol = 5)
#   P_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 5)
#   B <- matrix(0,ncol = 5, nrow = 1)
#   B_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
#   X_tilde <- matrix(0,ncol = 1, nrow = 1)
#   P_tilde <- matrix(1/5,ncol = 5, nrow = 1)
#   P_tilde_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
#   gamma_tilde <- matrix(0,nrow = 1, ncol = 1)
#   gamma_tilde_fixed_indices <- matrix(TRUE, nrow = 1, ncol = 1)
#
#
#   poisson_fit <- estimate_parameters(W = W,
#                                           X = X,
#                                           Z = Z,
#                                           Z_tilde = Z_tilde,
#                                           Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                                           gammas = gammas,
#                                           gammas_fixed_indices = gammas_fixed_indices,
#                                           P = P,
#                                           P_fixed_indices = P_fixed_indices,
#                                           B = B,
#                                           B_fixed_indices = B_fixed_indices,
#                                           X_tilde = X_tilde,
#                                           P_tilde = P_tilde,
#                                           P_tilde_fixed_indices = P_tilde_fixed_indices,
#                                           gamma_tilde = gamma_tilde,
#                                           gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
#                                           barrier_t = 1, #starting value of reciprocal barrier penalty coef.
#                                           barrier_scale = 10, #increments for value of barrier penalty
#                                           max_barrier = 1e12, #maximum value of barrier_t
#                                           initial_conv_tol = 1000,
#                                           final_conv_tol = 0.1,
#                                           final_f = 1e-6,
#                                           constraint_tolerance = 1e-15,
#                                           hessian_regularization = 0.01,
#                                           criterion = "Poisson",
#                                           subproblem_method = "Newton",
#                                           profile_P = FALSE,
#                                           profiling_maxit = 25
#   )
#
#
#   poisson_fit_params <- dataframes_to_parameters(fixed = poisson_fit$fixed,
#                                                  varying = poisson_fit$varying)
#
#   #update values of parameters
#   gammas <- poisson_fit_params$gammas
#   P <- poisson_fit_params$P
#   B <- poisson_fit_params$B
#   P_tilde <- poisson_fit_params$P_tilde
#   gamma_tilde <- poisson_fit_params$gamma_tilde
#
#
#   parameter_dfs <- parameters_to_dataframes(P,
#                                             P_fixed_indices,
#                                             P_tilde,
#                                             P_tilde_fixed_indices,
#                                             B,
#                                             B_fixed_indices,
#                                             gammas,
#                                             gammas_fixed_indices,
#                                             gamma_tilde,
#                                             gamma_tilde_fixed_indices)
#
#   varying_df <- parameter_dfs$varying
#   fixed_df <- parameter_dfs$fixed
#   varying_lr_df <- ra_to_lr(varying_df)
#
#
#   K <- max(c(varying_df$k[varying_df$param == "P"],
#              fixed_df$k[fixed_df$param == "P"]))
#
#   fixed_P_multipliers <- sapply(1:K, function(k)
#     1 - sum(fixed_df$value[fixed_df$param == "P"&
#                              fixed_df$k ==k]))
#
#   K_tilde <- max(c(varying_df$k[varying_df$param == "P_tilde"],
#                    fixed_df$k[fixed_df$param == "P_tilde"]))
#
#   fixed_P_tilde_multipliers <- sapply(1:K_tilde, function(k)
#     1 - sum(fixed_df$value[fixed_df$param == "P_tilde"&
#                              fixed_df$k ==k]))
#
#   Ak_list <- get_Ak_list(fixed_df,
#                          varying_df,
#                          varying_lr_df)
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
#
#
#   n <- nrow(W)
#   W_long <- lapply(1:n,function(i) as.numeric(W[i,]))
#   W_long <- do.call(c,W_long)
#   par_to_long_mean <- function(x){
#     temp_varying <- varying_lr_df
#     temp_varying$value <- x
#     temp_varying <- lr_to_ra(fixed_df, temp_varying, varying_df)
#     temp_params <- dataframes_to_parameters(fixed_df,temp_varying)
#
#     means <- meaninate(gammas = temp_params$gammas,
#                        B = temp_params$B,
#                        P = temp_params$P,
#                        X = X,
#                        Z = Z,
#                        X_tilde = X_tilde,
#                        Z_tilde = Z_tilde,
#                        Z_tilde_gamma_cols = 1,
#                        gamma_tilde = temp_params$gamma_tilde,
#                        P_tilde= temp_params$P_tilde)
#
#     means_long <- lapply(1:n, function(i) as.numeric(means[i,]))
#     means_long <- do.call(c,means_long)
#     return(means_long)
#
#   }
#   means_long <- par_to_long_mean(varying_lr_df$value)
#
#   gmm_inv_wts <- get_gmm_inv_weights(W_long, means_long)
#
#   analytical_deriv <-
#     deriv_criterion_lr(W = W,
#                        X = X,
#                        Z = Z,
#                        which_k_p = which_k_p,
#                        which_k_p_tilde = which_k_p_tilde,
#                        fixed_P_multipliers = fixed_P_multipliers,
#                        fixed_P_tilde_multipliers = fixed_P_tilde_multipliers,
#                        which_B_keep = which_B_keep,
#                        which_B_rows = which_B_rows,
#                        which_gammas = which_gammas,
#                        which_gamma_tilde = which_gamma_tilde,
#                        Z_tilde = Z_tilde,
#                        Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                        X_tilde = X_tilde,
#                        Ak_list = Ak_list,
#                        A_tilde_k_list = A_tilde_k_list,
#                        fixed_df = fixed_df,
#                        varying_df = varying_df,
#                        varying_lr_df = varying_lr_df,
#                        K = K,
#                        K_tilde = K_tilde,
#                        barrier_t = 1,
#                        criterion = "GMM",
#                        lr_scale = TRUE,
#                        include_log_penalty_derivatives = TRUE,
#                        return_info = FALSE,
#                        wts = NULL,
#                        gmm_inv_wts = gmm_inv_wts)
#
#
#
#
#
#   gmm_direct <- function(x){ gmm_criterion(W_long, par_to_long_mean(x),
#                                            gmm_inv_wts)}
#
#   numerical_deriv <- numDeriv::grad(gmm_direct,varying_lr_df$value)
#
#   expect_equal(numerical_deriv,as.numeric(analytical_deriv$grad),
#                tolerance = 1e-2)
#
#
#
#
#
# }
# )

test_that("Poisson derivatives in gamma and rho are correct in simple case", {



  set.seed(0)
  W <- matrix(sapply(1:10,function(x) rpois(1,1000)),
              ncol = 5)
  X <- matrix(1,ncol = 1, nrow = 2)
  Z <- matrix(1,nrow = 2, ncol = 1)
  Z_tilde <- matrix(0,nrow = 2, ncol = 1)
  Z_tilde_gamma_cols <- 1
  gammas <- apply(W,1,function(x) log(sum(x)))
  gammas_fixed_indices <- rep(F,2)
  P <- matrix(1/5, nrow = 1, ncol = 5)
  P_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 5)
  B <- matrix(0,ncol = 5, nrow = 1)
  B_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
  X_tilde <- matrix(0,ncol = 1, nrow = 1)
  P_tilde <- matrix(1/5,ncol = 5, nrow = 1)
  P_tilde_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
  gamma_tilde <- matrix(0,nrow = 1, ncol = 1)
  gamma_tilde_fixed_indices <- matrix(TRUE, nrow = 1, ncol = 1)


  parameter_dfs <- parameters_to_dataframes(P,
                                            P_fixed_indices,
                                            P_tilde,
                                            P_tilde_fixed_indices,
                                            B,
                                            B_fixed_indices,
                                            gammas,
                                            gammas_fixed_indices,
                                            gamma_tilde,
                                            gamma_tilde_fixed_indices)

  varying_df <- parameter_dfs$varying
  fixed_df <- parameter_dfs$fixed
  varying_lr_df <- ra_to_lr(varying_df)


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

  analytical_deriv <-
    deriv_criterion_lr(W = W,
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
                       gmm_inv_wts = NULL)


  par_to_mean <- function(x){
    temp_varying <- varying_lr_df
    temp_varying$value <- x
    temp_varying <- lr_to_ra(fixed_df, temp_varying, varying_df)
    temp_params <- dataframes_to_parameters(fixed_df,temp_varying)

    means <- meaninate(gammas = temp_params$gammas,
                       B = temp_params$B,
                       P = temp_params$P,
                       X = X,
                       Z = Z,
                       X_tilde = X_tilde,
                       Z_tilde = Z_tilde,
                       Z_tilde_gamma_cols = 1,
                       gamma_tilde = temp_params$gamma_tilde,
                       P_tilde= temp_params$P_tilde)

    return(means)

  }

  pois_direct <- function(x){ poisson_criterion(W, par_to_mean(x))}

  numerical_deriv <- numDeriv::grad(pois_direct,varying_lr_df$value)

  expect_equal(numerical_deriv,as.numeric(analytical_deriv$grad),
               tolerance = 1e-2)
}
)


# test_that("GMM derivatives in gamma, rho, B, gamma_tilde,
# and rho_tilde are correct in simple case", {
#
#
#   set.seed(0)
#   W <- matrix(sapply(1:10,function(x) rpois(1,1000)),
#               ncol = 5)
#   X <- matrix(1,ncol = 1, nrow = 2)
#   Z <- matrix(1,nrow = 2, ncol = 1)
#   Z_tilde <- matrix(1,nrow = 2, ncol = 1)
#   Z_tilde_gamma_cols <- 1
#   gammas <- apply(W,1,function(x) log(sum(x)))
#   gammas_fixed_indices <- rep(F,2)
#   P <- matrix(1/5, nrow = 1, ncol = 5)
#   P_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 5)
#   B <- matrix(0,ncol = 5, nrow = 1)
#   B_fixed_indices <- matrix(FALSE, ncol = 5, nrow = 1)
#   B_fixed_indices[,5] <- TRUE
#   X_tilde <- matrix(0,ncol = 1, nrow = 1)
#   P_tilde <- matrix(1/5,ncol = 5, nrow = 1)
#   P_tilde_fixed_indices <- matrix(FALSE, ncol = 5, nrow = 1)
#   gamma_tilde <- matrix(0,nrow = 1, ncol = 1)
#   gamma_tilde_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 1)
#
#
#   parameter_dfs <- parameters_to_dataframes(P,
#                                             P_fixed_indices,
#                                             P_tilde,
#                                             P_tilde_fixed_indices,
#                                             B,
#                                             B_fixed_indices,
#                                             gammas,
#                                             gammas_fixed_indices,
#                                             gamma_tilde,
#                                             gamma_tilde_fixed_indices)
#
#   varying_df <- parameter_dfs$varying
#   fixed_df <- parameter_dfs$fixed
#   varying_lr_df <- ra_to_lr(varying_df)
#
#
#   K <- max(c(varying_df$k[varying_df$param == "P"],
#              fixed_df$k[fixed_df$param == "P"]))
#
#   fixed_P_multipliers <- sapply(1:K, function(k)
#     1 - sum(fixed_df$value[fixed_df$param == "P"&
#                              fixed_df$k ==k]))
#
#   K_tilde <- max(c(varying_df$k[varying_df$param == "P_tilde"],
#                    fixed_df$k[fixed_df$param == "P_tilde"]))
#
#   fixed_P_tilde_multipliers <- sapply(1:K_tilde, function(k)
#     1 - sum(fixed_df$value[fixed_df$param == "P_tilde"&
#                              fixed_df$k ==k]))
#
#   Ak_list <- get_Ak_list(fixed_df,
#                          varying_df,
#                          varying_lr_df)
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
#
#   n <- nrow(W)
#   W_long <- lapply(1:n,function(i) as.numeric(W[i,]))
#   W_long <- do.call(c,W_long)
#   par_to_long_mean <- function(x){
#     temp_varying <- varying_lr_df
#     temp_varying$value <- x
#     temp_varying <- lr_to_ra(fixed_df, temp_varying, varying_df)
#     temp_params <- dataframes_to_parameters(fixed_df,temp_varying)
#
#     means <- meaninate(gammas = temp_params$gammas,
#                        B = temp_params$B,
#                        P = temp_params$P,
#                        X = X,
#                        Z = Z,
#                        X_tilde = X_tilde,
#                        Z_tilde = Z_tilde,
#                        Z_tilde_gamma_cols = 1,
#                        gamma_tilde = temp_params$gamma_tilde,
#                        P_tilde= temp_params$P_tilde)
#
#     means_long <- lapply(1:n, function(i) as.numeric(means[i,]))
#     means_long <- do.call(c,means_long)
#     return(means_long)
#
#   }
#   means_long <- par_to_long_mean(varying_lr_df$value)
#
#   gmm_inv_wts <- get_gmm_inv_weights(W_long, means_long)
#
#   analytical_deriv <-
#     deriv_criterion_lr(W = W,
#                        X = X,
#                        Z = Z,
#                        which_k_p = which_k_p,
#                        which_k_p_tilde = which_k_p_tilde,
#                        fixed_P_multipliers = fixed_P_multipliers,
#                        fixed_P_tilde_multipliers = fixed_P_tilde_multipliers,
#                        which_B_keep = which_B_keep,
#                        which_B_rows = which_B_rows,
#                        which_gammas = which_gammas,
#                        which_gamma_tilde = which_gamma_tilde,
#                        Z_tilde = Z_tilde,
#                        Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                        X_tilde = X_tilde,
#                        Ak_list = Ak_list,
#                        A_tilde_k_list = A_tilde_k_list,
#                        fixed_df = fixed_df,
#                        varying_df = varying_df,
#                        varying_lr_df = varying_lr_df,
#                        K = K,
#                        K_tilde = K_tilde,
#                        barrier_t = 1,
#                        criterion = "GMM",
#                        lr_scale = TRUE,
#                        include_log_penalty_derivatives = TRUE,
#                        return_info = FALSE,
#                        wts = NULL,
#                        gmm_inv_wts = gmm_inv_wts)
#
#
#
#   gmm_direct <- function(x){ gmm_criterion(W_long, par_to_long_mean(x),
#                                            gmm_inv_wts)}
#
#   numerical_deriv <- numDeriv::grad(gmm_direct,varying_lr_df$value)
#
#
#   analytical_deriv <- as.numeric(analytical_deriv$grad)
#
#   expect_equal(numerical_deriv,analytical_deriv)
#
#
#
# }
# )

test_that("Poisson derivatives in gamma, rho, and alpha_tilde are correct
          when alpha_tilde present", {



  set.seed(0)
  W <- matrix(sapply(1:10,function(x) rpois(1,1000)),
              ncol = 5)
  X <- matrix(1,ncol = 1, nrow = 2)
  Z <- matrix(1,nrow = 2, ncol = 1)
  # Z_tilde <- matrix(1,nrow = 2, ncol = 1)
  Z_tilde_gamma_cols <- 1
  Z_tilde_list <- list(matrix(0,nrow = 2,ncol = 1),
                       matrix(1,nrow = 2,ncol = 1))
  alpha_tilde <- 3
  gammas <- apply(W,1,function(x) log(sum(x)))
  gammas_fixed_indices <- rep(F,2)
  P <- matrix(1/5, nrow = 1, ncol = 5)
  P_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 5)
  B <- matrix(0,ncol = 5, nrow = 1)
  B_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
  X_tilde <- matrix(0,ncol = 1, nrow = 1)
  P_tilde <- matrix(1/5,ncol = 5, nrow = 1)
  P_tilde_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
  gamma_tilde <- matrix(0,nrow = 1, ncol = 1)
  gamma_tilde_fixed_indices <- matrix(TRUE, nrow = 1, ncol = 1)


  parameter_dfs <- parameters_to_dataframes(P,
                                            P_fixed_indices,
                                            P_tilde,
                                            P_tilde_fixed_indices,
                                            B,
                                            B_fixed_indices,
                                            gammas,
                                            gammas_fixed_indices,
                                            gamma_tilde,
                                            gamma_tilde_fixed_indices,
                                            alpha_tilde)

  varying_df <- parameter_dfs$varying
  fixed_df <- parameter_dfs$fixed
  varying_lr_df <- ra_to_lr(varying_df)


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

  analytical_deriv <-
    deriv_criterion_lr(W = W,
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
                       # Z_tilde = Z_tilde,
                       Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                       Z_tilde_list = Z_tilde_list,
                       alpha_tilde = alpha_tilde,
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
                       gmm_inv_wts = NULL)


  par_to_mean <- function(x){
    temp_varying <- varying_lr_df
    temp_varying$value <- x
    temp_varying <- lr_to_ra(fixed_df, temp_varying, varying_df)
    temp_params <- dataframes_to_parameters(fixed_df,temp_varying)

    means <- meaninate(gammas = temp_params$gammas,
                       B = temp_params$B,
                       P = temp_params$P,
                       X = X,
                       Z = Z,
                       X_tilde = X_tilde,
                       # Z_tilde = Z_tilde,
                       Z_tilde_list = Z_tilde_list,
                       Z_tilde_gamma_cols = 1,
                       gamma_tilde = temp_params$gamma_tilde,
                       P_tilde= temp_params$P_tilde,
                       alpha_tilde = temp_params$alpha_tilde)

    return(means)

  }

  pois_direct <- function(x){ poisson_criterion(W, par_to_mean(x))}

  numerical_deriv <- numDeriv::grad(pois_direct,varying_lr_df$value)

  expect_equal(numerical_deriv,as.numeric(analytical_deriv$grad),
               tolerance = 1e-2)

          }
)

# test_that("GMM derivatives in gamma, rho, B, gamma_tilde, rho_tilde,
# and alpha_tilde are correct in simple case", {
#
#
#   set.seed(0)
#   W <- matrix(sapply(1:10,function(x) rpois(1,1000)),
#               ncol = 5)
#   X <- matrix(1,ncol = 1, nrow = 2)
#   Z <- matrix(1,nrow = 2, ncol = 1)
#   Z_tilde <- NULL
#   Z_tilde_list <- list(matrix(c(1,0),nrow = 2, ncol = 1),
#                        matrix(c(0,1),nrow = 2, ncol = 1))
#   alpha_tilde <- 3
#   Z_tilde_gamma_cols <- 1
#   gammas <- apply(W,1,function(x) log(sum(x)))
#   gammas_fixed_indices <- rep(F,2)
#   P <- matrix(1/5, nrow = 1, ncol = 5)
#   P_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 5)
#   B <- matrix(0,ncol = 5, nrow = 1)
#   B_fixed_indices <- matrix(FALSE, ncol = 5, nrow = 1)
#   B_fixed_indices[,5] <- TRUE
#   X_tilde <- matrix(0,ncol = 1, nrow = 1)
#   P_tilde <- matrix(1/5,ncol = 5, nrow = 1)
#   P_tilde_fixed_indices <- matrix(FALSE, ncol = 5, nrow = 1)
#   gamma_tilde <- matrix(0,nrow = 1, ncol = 1)
#   gamma_tilde_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 1)
#
#
#   parameter_dfs <- parameters_to_dataframes(P,
#                                             P_fixed_indices,
#                                             P_tilde,
#                                             P_tilde_fixed_indices,
#                                             B,
#                                             B_fixed_indices,
#                                             gammas,
#                                             gammas_fixed_indices,
#                                             gamma_tilde,
#                                             gamma_tilde_fixed_indices,
#                                             alpha_tilde = alpha_tilde)
#
#   varying_df <- parameter_dfs$varying
#   fixed_df <- parameter_dfs$fixed
#   varying_lr_df <- ra_to_lr(varying_df)
#
#
#   K <- max(c(varying_df$k[varying_df$param == "P"],
#              fixed_df$k[fixed_df$param == "P"]))
#
#   fixed_P_multipliers <- sapply(1:K, function(k)
#     1 - sum(fixed_df$value[fixed_df$param == "P"&
#                              fixed_df$k ==k]))
#
#   K_tilde <- max(c(varying_df$k[varying_df$param == "P_tilde"],
#                    fixed_df$k[fixed_df$param == "P_tilde"]))
#
#   fixed_P_tilde_multipliers <- sapply(1:K_tilde, function(k)
#     1 - sum(fixed_df$value[fixed_df$param == "P_tilde"&
#                              fixed_df$k ==k]))
#
#   Ak_list <- get_Ak_list(fixed_df,
#                          varying_df,
#                          varying_lr_df)
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
#
#   n <- nrow(W)
#   W_long <- lapply(1:n,function(i) as.numeric(W[i,]))
#   W_long <- do.call(c,W_long)
#   par_to_long_mean <- function(x){
#     temp_varying <- varying_lr_df
#     temp_varying$value <- x
#     temp_varying <- lr_to_ra(fixed_df, temp_varying, varying_df)
#     temp_params <- dataframes_to_parameters(fixed_df,temp_varying)
#
#     means <- meaninate(gammas = temp_params$gammas,
#                        B = temp_params$B,
#                        P = temp_params$P,
#                        X = X,
#                        Z = Z,
#                        X_tilde = X_tilde,
#                        Z_tilde = Z_tilde,
#                        Z_tilde_gamma_cols = 1,
#                        Z_tilde_list = Z_tilde_list,
#                        gamma_tilde = temp_params$gamma_tilde,
#                        P_tilde= temp_params$P_tilde,
#                        alpha_tilde = temp_params$alpha_tilde)
#
#     means_long <- lapply(1:n, function(i) as.numeric(means[i,]))
#     means_long <- do.call(c,means_long)
#     return(means_long)
#
#   }
#   means_long <- par_to_long_mean(varying_lr_df$value)
#
#   gmm_inv_wts <- get_gmm_inv_weights(W_long, means_long)
#
#   analytical_deriv <-
#     deriv_criterion_lr(W = W,
#                        X = X,
#                        Z = Z,
#                        which_k_p = which_k_p,
#                        which_k_p_tilde = which_k_p_tilde,
#                        fixed_P_multipliers = fixed_P_multipliers,
#                        fixed_P_tilde_multipliers = fixed_P_tilde_multipliers,
#                        which_B_keep = which_B_keep,
#                        which_B_rows = which_B_rows,
#                        which_gammas = which_gammas,
#                        which_gamma_tilde = which_gamma_tilde,
#                        Z_tilde = Z_tilde,
#                        Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                        Z_tilde_list = Z_tilde_list,
#                        alpha_tilde = alpha_tilde,
#                        X_tilde = X_tilde,
#                        Ak_list = Ak_list,
#                        A_tilde_k_list = A_tilde_k_list,
#                        fixed_df = fixed_df,
#                        varying_df = varying_df,
#                        varying_lr_df = varying_lr_df,
#                        K = K,
#                        K_tilde = K_tilde,
#                        barrier_t = 1,
#                        criterion = "GMM",
#                        lr_scale = TRUE,
#                        include_log_penalty_derivatives = TRUE,
#                        return_info = FALSE,
#                        wts = NULL,
#                        gmm_inv_wts = gmm_inv_wts)
#
#
#
#   gmm_direct <- function(x){ gmm_criterion(W_long, par_to_long_mean(x),
#                                            gmm_inv_wts)}
#
#   numerical_deriv <- numDeriv::grad(gmm_direct,varying_lr_df$value)
#
#
#   analytical_deriv <- as.numeric(analytical_deriv$grad)
#
#   expect_equal(numerical_deriv,analytical_deriv,
#                tolerance = 1e-2)
#
#
# }
# )

test_that("Poisson derivatives correct in realistic setting with alpha_tilde", {


  set.seed(0)
  p_mock1 <- c(rep(0,5),rep(1/5,5))
  p_mock2 <- c(rep(1/5,5),rep(0,5))
  p_contam <- c(rexp(5),rep(0,5))
  p_contam <- p_contam/sum(p_contam)
  p_true <- rep(1/10,10)
  dilutions <- rep(3^(1:5),3)
  W <- matrix(NA,nrow = 15, ncol = 10)
  for(i in 1:15){
    if(i<6){
      W[i,] <- round(rexp(1,3^((i - 1)%%5)*1/10000)*(p_mock1 + dilutions[i]*p_contam),0)
    } else{
      if(i<11){
        W[i,] <- round(rexp(1,3^((i- 1)%%5)*1/10000)*(p_mock2 + dilutions[i]*p_contam),0)

      } else{
        W[i,] <- round(rexp(1,3^((i-1)%%5)*1/10000)*(p_true + dilutions[i]*p_contam),0)

      }
    }
  }
  X <- matrix(0,ncol = 1, nrow = 15)
  Z <- cbind(c(rep(1,5),rep(0,10)),
             c(rep(0,5),rep(1,5), rep(0,5)),
             c(rep(0,10),rep(1,5)))
  Z_tilde <- matrix(dilutions/exp(mean(log(dilutions))), ncol = 1)
  Z_tilde_list <- list(Z_tilde*matrix(c(rep(1,5),rep(0,10))),
                       Z_tilde*matrix(c(rep(0,5),rep(1,5),rep(0,5))),
                       Z_tilde*matrix(c(rep(0,10),rep(1,5))))

  Z_tilde_gamma_cols <- 1
  gammas <- apply(W,1,function(x) log(sum(x)))
  gammas_fixed_indices <- rep(F,length(gammas))
  P <- rbind(p_mock1,
             p_mock2,
             rep(.1,10))
  P_fixed_indices <- matrix(FALSE, nrow = 3, ncol = 10)
  P_fixed_indices[1:2,] <- TRUE
  B <- matrix(0,ncol = 10, nrow = 1)
  B_fixed_indices <- matrix(TRUE, ncol = 10, nrow = 1)
  X_tilde <- matrix(0,ncol = 1, nrow = 1)
  P_tilde <- matrix(1/10,ncol = 10, nrow = 1)
  P_tilde_fixed_indices <- matrix(FALSE, ncol = 10, nrow = 1)
  gamma_tilde <- matrix(0,nrow = 1, ncol = 1)
  gamma_tilde_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 1)
  alpha_tilde <- c(0,0)

  parameter_dfs <- parameters_to_dataframes(P,
                                            P_fixed_indices,
                                            P_tilde,
                                            P_tilde_fixed_indices,
                                            B,
                                            B_fixed_indices,
                                            gammas,
                                            gammas_fixed_indices,
                                            gamma_tilde,
                                            gamma_tilde_fixed_indices,
                                            alpha_tilde = alpha_tilde)

  varying_df <- parameter_dfs$varying
  fixed_df <- parameter_dfs$fixed
  varying_lr_df <- ra_to_lr(varying_df)



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

            analytical_deriv <-
              deriv_criterion_lr(W = W,
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
                                 # Z_tilde = Z_tilde,
                                 Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                 Z_tilde_list = Z_tilde_list,
                                 alpha_tilde = alpha_tilde,
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
                                 include_log_penalty_derivatives = FALSE,
                                 return_info = FALSE,
                                 wts = NULL,
                                 gmm_inv_wts = NULL)


            par_to_mean <- function(x){
              temp_varying <- varying_lr_df
              temp_varying$value <- x
              temp_varying <- lr_to_ra(fixed_df, temp_varying, varying_df)
              temp_params <- dataframes_to_parameters(fixed_df,temp_varying)

              means <- meaninate(gammas = temp_params$gammas,
                                 B = temp_params$B,
                                 P = temp_params$P,
                                 X = X,
                                 Z = Z,
                                 X_tilde = X_tilde,
                                 # Z_tilde = Z_tilde,
                                 Z_tilde_list = Z_tilde_list,
                                 Z_tilde_gamma_cols = 1,
                                 gamma_tilde = temp_params$gamma_tilde,
                                 P_tilde= temp_params$P_tilde,
                                 alpha_tilde = temp_params$alpha_tilde)

              return(means)

            }

            pois_direct <- function(x){ poisson_criterion(W, par_to_mean(x))}

            numerical_deriv <- numDeriv::grad(pois_direct,varying_lr_df$value)

            # plot(1:nrow(varying_lr_df),analytical_deriv$grad - numerical_deriv)
            #
            # data.frame(param = varying_lr_df$param,
            #            nderiv = numerical_deriv,
            #            aderiv = as.numeric(analytical_deriv$grad),
            #            index = 1:nrow(varying_lr_df)) %>%
            #   ggplot() +
            #   geom_point(aes(x = index, y= asinh(aderiv-nderiv),color = param))+
            #   theme_bw()

            expect_equal(numerical_deriv,as.numeric(analytical_deriv$grad))

          }
)



