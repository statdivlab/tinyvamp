test_that("jacobian correct for B and gammas in simple example", {

  X <- matrix(1)
  X_tilde <- matrix(1)
  Z <- matrix(1)
  Z_tilde <- matrix(1)
  Z_tilde_gamma_cols <- 1

  P <- matrix(rep(1/8, 8), nrow = 1)
  P_fixed_indices <- matrix(rep(TRUE,8), nrow = 1)
  P_tilde <- matrix((1:8)/sum(1:8), nrow = 1)
  P_tilde_fixed_indices <- matrix(c(rep(TRUE,8)), nrow = 1)
  B <- matrix(c(-3:3,0),nrow = 1)
  B_fixed_indices <- matrix(c(rep(FALSE,7),TRUE), nrow = 1)
  gammas <- 8
  gammas_fixed_indices <- FALSE
  gamma_tilde <- matrix(log(100))
  gamma_tilde_fixed_indices <- matrix(TRUE)

param_dfs <- parameters_to_dataframes(P = P,
                         P_tilde = P_tilde,
                         P_fixed_indices = P_fixed_indices,
                         P_tilde_fixed_indices = P_tilde_fixed_indices,
                         B = B,
                         B_fixed_indices = B_fixed_indices,
                         gammas = gammas,
                         gammas_fixed_indices = gammas_fixed_indices,
                         gamma_tilde = gamma_tilde,
                         gamma_tilde_fixed_indices = gamma_tilde_fixed_indices)

params <- dataframes_to_parameters(param_dfs$fixed,
                                   param_dfs$varying)




  function_output <-
    mean_jac_lr_faster(fixed_df = param_dfs$fixed,
                       varying_lr_df = ra_to_lr(param_dfs$varying),
                       varying_df = param_dfs$varying,
                       which_k_p = NULL,
                       which_k_p_tilde = NULL,
                       which_B_rows = 1,
                       which_B_keep = matrix(rep(TRUE,7),nrow= 1),
                       which_gammas = 1,
                       which_gamma_tilde = 1,
                       params = params,
                       Ak_list = list(NULL),
                       A_tilde_k_list = list(NULL),
                       fixed_P_multipliers = 1,
                       fixed_P_tilde_multipliers = 1,
                       K = 1,
                       K_tilde = 1,
                       X = X,
                       Z = Z,
                       X_tilde = X_tilde,
                       Z_tilde = Z_tilde,
                       Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                       sparse = TRUE)


  means <-
    P*exp(gammas + B) + P_tilde*exp(as.numeric(gamma_tilde) + gammas + B)

  means <- as.numeric(means)

  direct_comp <- cbind(diag(means)[,1:7],means)
  colnames(direct_comp) <- NULL

  expect_equal(max(abs(as.matrix(function_output) -  direct_comp)), 0) 

})

test_that("jacobian correct for rho in simple example", {

  X <- matrix(1)
  X_tilde <- matrix(1)
  Z <- matrix(1)
  Z_tilde <- matrix(1)
  Z_tilde_gamma_cols <- 1

  P <- matrix(rep(1/8, 8), nrow = 1)
  P_fixed_indices <- matrix(rep(FALSE,8), nrow = 1)
  P_tilde <- matrix((1:8)/sum(1:8), nrow = 1)
  P_tilde_fixed_indices <- matrix(c(rep(TRUE,8)), nrow = 1)
  B <- matrix(c(-3:3,0),nrow = 1)
  B_fixed_indices <- matrix(c(rep(TRUE,7),TRUE), nrow = 1)
  gammas <- 8
  gammas_fixed_indices <- TRUE
  gamma_tilde <- matrix(log(100))
  gamma_tilde_fixed_indices <- matrix(TRUE)

  param_dfs <- parameters_to_dataframes(P = P,
                                        P_tilde = P_tilde,
                                        P_fixed_indices = P_fixed_indices,
                                        P_tilde_fixed_indices = P_tilde_fixed_indices,
                                        B = B,
                                        B_fixed_indices = B_fixed_indices,
                                        gammas = gammas,
                                        gammas_fixed_indices = gammas_fixed_indices,
                                        gamma_tilde = gamma_tilde,
                                        gamma_tilde_fixed_indices = gamma_tilde_fixed_indices)

  params <- dataframes_to_parameters(param_dfs$fixed,
                                     param_dfs$varying)




  function_output <-
    mean_jac_lr_faster(fixed_df = param_dfs$fixed,
                       varying_lr_df = ra_to_lr(param_dfs$varying),
                       varying_df = param_dfs$varying,
                       which_k_p = 1,
                       which_k_p_tilde = NULL,
                       which_B_rows = 1,
                       which_B_keep = matrix(rep(FALSE,7),nrow= 1),
                       which_gammas = NULL,
                       which_gamma_tilde = NULL,
                       params = params,
                       Ak_list = list(diag(8)),
                       A_tilde_k_list = list(NULL),
                       fixed_P_multipliers = 1,
                       fixed_P_tilde_multipliers = 0,
                       K = 1,
                       K_tilde = 1,
                       X = X,
                       Z = Z,
                       X_tilde = X_tilde,
                       Z_tilde = Z_tilde,
                       Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                       sparse = TRUE)


 premultiplier <- (1/8)*diag(rep(1,8))[1:7,] -
   (1/64)*matrix(1,nrow = 7, ncol = 1) %*% matrix(1, nrow = 1, ncol = 8)


function_output <- as.matrix(function_output)
dimnames(function_output) <- NULL


  expect_equal(function_output,
               t(premultiplier%*%diag(as.numeric(exp(gammas + B)))))

})

test_that("jacobian correct for rho_tilde in simple example", {

  X <- matrix(1)
  X_tilde <- matrix(1)
  Z <- matrix(1)
  Z_tilde <- matrix(1)
  Z_tilde_gamma_cols <- 1

  P <- matrix(rep(1/8, 8), nrow = 1)
  P_fixed_indices <- matrix(rep(TRUE,8), nrow = 1)
  P_tilde <- matrix(rep(1/8, 8), nrow = 1)
  P_tilde_fixed_indices <- matrix(c(rep(FALSE,8)), nrow = 1)
  B <- matrix(c(-3:3,0),nrow = 1)
  B_fixed_indices <- matrix(c(rep(TRUE,7),TRUE), nrow = 1)
  gammas <- 8
  gammas_fixed_indices <- TRUE
  gamma_tilde <- matrix(log(100))
  gamma_tilde_fixed_indices <- matrix(TRUE)

  param_dfs <- parameters_to_dataframes(P = P,
                                        P_tilde = P_tilde,
                                        P_fixed_indices = P_fixed_indices,
                                        P_tilde_fixed_indices = P_tilde_fixed_indices,
                                        B = B,
                                        B_fixed_indices = B_fixed_indices,
                                        gammas = gammas,
                                        gammas_fixed_indices = gammas_fixed_indices,
                                        gamma_tilde = gamma_tilde,
                                        gamma_tilde_fixed_indices = gamma_tilde_fixed_indices)

  params <- dataframes_to_parameters(param_dfs$fixed,
                                     param_dfs$varying)




  function_output <-
    mean_jac_lr_faster(fixed_df = param_dfs$fixed,
                       varying_lr_df = ra_to_lr(param_dfs$varying),
                       varying_df = param_dfs$varying,
                       which_k_p = NULL,
                       which_k_p_tilde = 1,
                       which_B_rows = 1,
                       which_B_keep = matrix(rep(FALSE,7),nrow= 1),
                       which_gammas = NULL,
                       which_gamma_tilde = NULL,
                       params = params,
                       Ak_list = list(NULL),
                       A_tilde_k_list = list(diag(8)),
                       fixed_P_multipliers = 0,
                       fixed_P_tilde_multipliers = 1,
                       K = 1,
                       K_tilde = 1,
                       X = X,
                       Z = Z,
                       X_tilde = X_tilde,
                       Z_tilde = Z_tilde,
                       Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                       sparse = TRUE)


  premultiplier <- (1/8)*diag(rep(1,8))[1:7,] -
    (1/64)*matrix(1,nrow = 7, ncol = 1) %*% matrix(1, nrow = 1, ncol = 8)


  function_output <- as.matrix(function_output)
  dimnames(function_output) <- NULL


  expect_equal(function_output,
               t(premultiplier%*%diag(as.numeric(exp(gammas + as.numeric(gamma_tilde) + B)))))

})



test_that("Jacobian correct in realistic setting with alpha_tilde", {

  skip("This test take a moment (and is redundant if rho_tilde derivative tests are passing)")

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

  params <- dataframes_to_parameters(fixed_df,varying_df)



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

  analytical_jacobian <-
  mean_jac_lr_faster(fixed_df = fixed_df,
                     varying_lr_df = varying_lr_df,
                     varying_df = varying_df,
                     which_k_p = which_k_p,
                     which_k_p_tilde = which_k_p_tilde,
                     which_B_rows = which_B_rows,
                     which_B_keep = which_B_keep,
                     which_gammas = which_gammas,
                     which_gamma_tilde = which_gamma_tilde,
                     params = params,
                     Ak_list = Ak_list,
                     A_tilde_k_list = A_tilde_k_list,
                     fixed_P_multipliers = fixed_P_multipliers,
                     fixed_P_tilde_multipliers = fixed_P_tilde_multipliers,
                     K = K,
                     K_tilde = K_tilde,
                     X = X,
                     Z = Z,
                     X_tilde = X_tilde,
                     Z_tilde = NULL,
                     Z_tilde_gamma_cols = 1,
                     Z_tilde_list = Z_tilde_list,
                     sparse = TRUE)


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



  numerical_deriv <-
    lapply(1:nrow(W),
           function(i) do.call(rbind,lapply(1:ncol(W),
                                        function(j)
                              matrix(numDeriv::grad(function(x)
    par_to_mean(x)[i,j],varying_lr_df$value),
    nrow = 1))))

  numerical_deriv <- do.call(rbind, numerical_deriv)

  analytical_jacobian <- as.matrix(analytical_jacobian)
  dimnames(analytical_jacobian) <- NULL
  expect_equal(numerical_deriv,   analytical_jacobian)

}
)
