test_that("gmm criterion behaves predictably in the inv_wts argument", {
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

  param_dfs <- parameters_to_dataframes(P,
                                        P_fixed_indices,
                                        P_tilde,
                                        P_tilde_fixed_indices,
                                        B,
                                        B_fixed_indices,
                                        gammas,
                                        gammas_fixed_indices,
                                        gamma_tilde,
                                        gamma_tilde_fixed_indices)

  means <- meaninate(gammas = gammas,
                     B = B,
                     X = X,
                     Z = Z,
                     P = P,
                     X_tilde = X_tilde,
                     Z_tilde = Z_tilde,
                     Z_tilde_gamma_cols = 1,
                     P_tilde = P_tilde,
                     gamma_tilde = gamma_tilde)
  n <- nrow(W)
  W_long <- lapply(1:n,function(i) as.numeric(W[i,]))
  W_long <- do.call(c,W_long)
  means_long <- lapply(1:n, function(i) as.numeric(means[i,]))
  means_long <- do.call(c,means_long)

  sse <- sum((W_long - means_long)^2)

  expect_equal(0.5*sse,gmm_criterion(W_long,means_long,
                                 rep(1,length(W_long))))

  expect_true(0.5*sse > gmm_criterion(W_long,means_long,
                                  rep(2,length(W_long))))


  W_long[1] <- means_long[1]

  expect_equal(gmm_criterion(W_long,means_long,c(0,rep(1,length(W_long) - 1))),
               gmm_criterion(W_long,means_long,c(1,rep(1,length(W_long) - 1))))

})
