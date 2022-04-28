test_that("rho derivative is correct when all Z_tilde rows are multiplied by exp(gamma)", {
  gammas <-  4.53
  B <-  matrix(c(rnorm(4),0),nrow = 1)
  X <-  matrix(1, nrow = 1, ncol = 1)
  Z <-  matrix(1, nrow = 1, ncol = 1)
  P <-  matrix((1:5)/15, nrow = 1, ncol = 5)
  X_tilde <-  matrix(1,nrow = 1, ncol = 1)
  Z_tilde <-  matrix(1, nrow = 1, ncol = 1)
  P_tilde <-  matrix((5:1)/15, nrow = 1, ncol = 5)
  rho_k = log(P[1,1:4]/P[1,5])
  gamma_tilde <- 1
  function_value <- mu_d_rho_faster(i = 1,
                           J = 5,
                           k = 1,
                           gammas = gammas,
                           B = B,
                           X = X,
                           Z = Z,
                           rho_k = rho_k,
                           Ak_list = list(diag(5)),
                           fixed_P_multipliers = 1
                          )

  #convert to matrix from dgeMatrix
  function_value <- as.matrix(function_value)

  #remove empty dimnames
  dimnames(function_value) <- NULL


  direct_d_mu_dPk <- diag(as.numeric((Z%*%exp(X%*%B + gammas))))
  direct_dPk_drho_k <- cbind(diag(exp(rho_k)/(sum(c(1,exp(rho_k))))),0) -
    outer(c(exp(rho_k))/(sum(c(1,exp(rho_k)))),c(exp(rho_k),1)/(sum(c(1,exp(rho_k)))))
  direct_calculation <- t(direct_dPk_drho_k%*%direct_d_mu_dPk)

  expect_equal(function_value, direct_calculation)
})
