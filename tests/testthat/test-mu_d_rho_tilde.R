test_that("rho derivative is correct when all Z_tilde rows are multiplied by exp(gamma)", {
  require(Matrix)
  gammas <-  4.53
  B <-  matrix(c(rnorm(4),0),nrow = 1)
  X_tilde <-  matrix(1,nrow = 1, ncol = 1)
  Z_tilde <-  matrix(1, nrow = 1, ncol = 1)
  P_tilde <-  matrix((5:1)/15, nrow = 1, ncol = 5)
  rho_tilde_k = log(P_tilde[1,1:4]/P_tilde[1,5])
  gamma_tilde <- 1
  Z_tilde_gamma_cols <- 1

  function_value <- mu_d_rho_tilde_faster(i = 1,
                                    J = 5,
                                    k_tilde = 1,
                                    gammas = gammas,
                                    B = B,
                                    rho_tilde_k = rho_tilde_k,
                                    A_tilde_k_list = list(diag(5)),
                                    fixed_P_tilde_multipliers = 1,
                                    X_tilde = X_tilde,
                                    Z_tilde = Z_tilde,
                                    Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                    gamma_tilde = gamma_tilde)

  #convert to matrix from dgeMatrix
  function_value <- as.matrix(function_value)

  #remove empty dimnames
  dimnames(function_value) <- NULL


  direct_d_mu_dP_tilde_k <- diag(as.numeric((exp(gammas)*Z_tilde)%*%exp(X_tilde%*%B + gamma_tilde)))
  direct_dP_tilde_k_drho_tilde_k <- cbind(diag(exp(rho_tilde_k)/(sum(c(1,exp(rho_tilde_k))))),0) -
    outer(c(exp(rho_tilde_k))/(sum(c(1,exp(rho_tilde_k)))),c(exp(rho_tilde_k),1)/(sum(c(1,exp(rho_tilde_k)))))
  direct_calculation <- t(direct_dP_tilde_k_drho_tilde_k%*%direct_d_mu_dP_tilde_k)

  expect_equal(function_value, direct_calculation)
})

test_that("rho_tilde derivative correct in realistic setting with alpha_tilde", {


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
  Z_tilde <- NULL
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

  rho_tilde_k = log(P_tilde[1,1:9]/P_tilde[1,10])
  i <- 1

  function_value <- mu_d_rho_tilde_faster(i = i,
                                          J = 10,
                                          k_tilde = 1,
                                          gammas = gammas,
                                          B = B,
                                          rho_tilde_k = rho_tilde_k,
                                          A_tilde_k_list = A_tilde_k_list,
                                          fixed_P_tilde_multipliers = 1,
                                          X_tilde = X_tilde,
                                          Z_tilde = Z_tilde,
                                          Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                                          gamma_tilde = gamma_tilde,
                                          alpha_tilde = alpha_tilde,
                                          Z_tilde_list = Z_tilde_list)

  #convert to matrix from dgeMatrix
  function_value <- as.matrix(function_value)

  #remove empty dimnames
  dimnames(function_value) <- NULL

  Z_tilde <- construct_Z_tilde(Z_tilde_list,
                               alpha_tilde)


  direct_d_mu_dP_tilde_k <- diag(as.numeric((exp(gammas[i])*Z_tilde[i,,drop = FALSE])%*%exp(X_tilde%*%B + gamma_tilde%*%matrix(1,ncol = 10,
                                                                                                           nrow = 1))))
  direct_dP_tilde_k_drho_tilde_k <- cbind(diag(exp(rho_tilde_k)/(sum(c(1,exp(rho_tilde_k))))),0) -
    outer(c(exp(rho_tilde_k))/(sum(c(1,exp(rho_tilde_k)))),c(exp(rho_tilde_k),1)/(sum(c(1,exp(rho_tilde_k)))))
  direct_calculation <- t(direct_dP_tilde_k_drho_tilde_k%*%direct_d_mu_dP_tilde_k)

  expect_equal(function_value, direct_calculation)

  get_mean <- function(rho_tilde_k){
    temp_lr_df <- varying_lr_df
    temp_lr_df$value[temp_lr_df$param == "rho_tilde"] <- rho_tilde_k
    temp_params <- dataframes_to_parameters(fixed_df,
                                       lr_to_ra(fixed_df,
                                                temp_lr_df,
                                                varying_df))

   return( meaninate(gammas = temp_params$gammas,
              B = temp_params$B,
              Z = Z,
              P = P,
              X = X,
              X_tilde = X_tilde,
              Z_tilde = NULL,
              Z_tilde_gamma_cols = 1,
              P_tilde = temp_params$P_tilde,
              gamma_tilde = temp_params$gamma_tilde,
              alpha_tilde = temp_params$alpha_tilde,
              Z_tilde_list = Z_tilde_list)[1,1])


  }

  num_grad <-
    numDeriv::grad(get_mean,
                   varying_lr_df$value[varying_lr_df$param == "rho_tilde"])

  expect_equal(function_value[1,],num_grad)

}
)
