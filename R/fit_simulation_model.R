

fit_simulation_model <- function(W,
                                 estimator,
                                 return_variance = FALSE){

  n <- nrow(W)/16
  J <- ncol(W)


  Z_tilde <- do.call(rbind,lapply(1:n,
                                  function(x) matrix(rep(c(1,9,81,729),
                                                         4),
                                                     ncol = 1)))

  Z_tilde_list <- list(Z_tilde*rbind(matrix(1,nrow = 4*n,ncol = 1),
                                     matrix(0, nrow = 3*4*n,ncol = 1)),
                       Z_tilde*rbind(matrix(0,nrow = 4*n,ncol = 1),
                                     matrix(1, nrow = 4*n,ncol = 1),
                                     matrix(0,nrow = 2*4*n,ncol = 1)),
                       Z_tilde*rbind(matrix(0,nrow = 2*4*n,ncol = 1),
                                     matrix(1, nrow = 4*n,ncol = 1),
                                     matrix(0,nrow = 4*n,ncol = 1)),
                       Z_tilde*rbind(matrix(0,nrow = 3*4*n,ncol = 1),
                                     matrix(1, nrow = 4*n,ncol = 1)))

  Z_tilde <- NULL

  Z_tilde_gamma_cols <- 1

  alpha_tilde <- c(0,0,0)
  gamma_tilde <- matrix(-5,ncol = 1, nrow = 1)

  ### generate Z
  Z <- do.call(rbind,lapply(1:4,
                            function(x) do.call(rbind,
                                                lapply(1:(4*n),function(k) matrix(
                                                  as.numeric(x == 1:4),nrow = 1
                                                )))))

  X <- matrix(1,nrow = n*16,ncol = 1)

  B <- matrix(0,ncol = J, nrow = 1)
  B_fixed_indices <- matrix(c(rep(FALSE,J - 1), TRUE), nrow = 1)

  P <- matrix(1/J, nrow = 4, ncol = J)
  P[1,] <- 2^(seq(0,4,length.out = J))
  P[1,] <- P[1,]/sum(P[1,])
  P[2,] <- P[1,J:1]

  P_fixed_indices <- rbind(matrix(TRUE, nrow = 1, ncol = J),
                           matrix(TRUE,nrow = 1, ncol = J),
                           matrix(FALSE, nrow = 1, ncol = J),
                           matrix(FALSE, nrow = 1, ncol = J))

  P_tilde <- matrix(1/J,ncol = J, nrow = 1)
  P_tilde_fixed_indices <- matrix(FALSE, nrow = 1, ncol = J)
  X_tilde <- matrix(1,ncol = 1, nrow = 1)
  gamma_tilde <- matrix(-3.7,ncol = 1,nrow = 1)
  gamma_tilde_fixed_indices <- matrix(FALSE, ncol = 1, nrow =1)

  gammas <- apply(W, 1, function(x) log(sum(x)))
  gammas_fixed_indices <- rep(FALSE, 16*n)

  fitted_model <- estimate_parameters(W = W,
                      X = X,
                      Z = Z,
                      Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                      gammas = gammas,
                      gammas_fixed_indices = gammas_fixed_indices,
                      P = P,
                      P_fixed_indices = P_fixed_indices,
                      B = B,
                      B_fixed_indices = B_fixed_indices,
                      X_tilde = X_tilde,
                      P_tilde = P_tilde,
                      P_tilde_fixed_indices = P_tilde_fixed_indices,
                      gamma_tilde = gamma_tilde,
                      gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
                      alpha_tilde = alpha_tilde,
                      Z_tilde_list = Z_tilde_list,
                                  barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                                  barrier_scale = 10, #increments for value of barrier penalty
                                  max_barrier = 1e12, #maximum value of barrier_t
                                  initial_conv_tol = 1000,
                                  final_conv_tol = 0.1,
                                  final_f = 1e-6,
                                  constraint_tolerance = 1e-10,
                                  hessian_regularization = 0.01,
                                  criterion = estimator,
                                  subproblem_method = "Newton",
                                  profile_P = TRUE,
                                  profiling_maxit = 25,
                                  wts = NULL,
                                  verbose = TRUE,
                                  bootstrap_failure_cutoff = NULL,
                      return_variance = return_variance
  )



return(fitted_model)

}
