test_that("P derivative is correct in simple setting", {
  gammas <-  4.53
  B <-  matrix(c(rnorm(4),0),nrow = 1)
  X <-  matrix(1, nrow = 1, ncol = 1)
  Z <-  matrix(1, nrow = 1, ncol = 1)
  P <-  matrix((1:5)/15, nrow = 1, ncol = 5)
  X_tilde <-  matrix(1,nrow = 1, ncol = 1)
  Z_tilde <-  matrix(1, nrow = 1, ncol = 1)
  P_tilde <-  matrix((5:1)/15, nrow = 1, ncol = 5)
  gamma_tilde <- 1
  function_value <- mu_d_P(i = 1,
                                     j = 5,
                                     m = 1,
                                     gammas = gammas,
                                     B = B,
                                     X = X,
                                     Z = Z,
                                     P = P)

  direct_calculation <- (Z%*%exp(X%*%B + gammas))[5]

  expect_equal(as.numeric(function_value), direct_calculation)
})
