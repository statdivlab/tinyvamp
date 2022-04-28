test_that("gamma derivative is correct when all Z_tilde rows are multiplied by exp(gamma)", {
  gammas <-  4.53
  B <-  matrix(c(rnorm(4),0),nrow = 1)
  X <-  matrix(1, nrow = 1, ncol = 1)
  Z <-  matrix(1, nrow = 1, ncol = 1)
  P <-  matrix((1:5)/15, nrow = 1, ncol = 5)
  X_tilde <-  matrix(1,nrow = 1, ncol = 1)
  Z_tilde <-  matrix(1, nrow = 1, ncol = 1)
  P_tilde <-  matrix((5:1)/15, nrow = 1, ncol = 5)
  gamma_tilde <- 1
function_value <- mu_d_gamma(i = 1,
           j = 5,
           gammas = gammas,
           B = B,
           X = X,
           Z = Z,
           P = P,
           X_tilde = X_tilde,
           Z_tilde = Z_tilde,
           Z_tilde_gamma_cols = 1,
           P_tilde = P_tilde,
           gamma_tilde = gamma_tilde)

direct_calculation <- ((Z%*%P)*exp(X%*%B + gammas) +
  Z_tilde%*%(P_tilde*exp(X_tilde%*%B + gamma_tilde + gammas)))[5]

expect_equal(as.numeric(function_value), direct_calculation)
})
