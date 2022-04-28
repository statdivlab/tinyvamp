test_that("Beta derivative is correct in no-contamination case", {

  set.seed(9382)
  B <- matrix(c(rnorm(6),0), nrow = 1)
  X <-  matrix(1,ncol = 1, nrow = 1)
  Z <-  matrix(1, ncol = 1, nrow = 1)
  P <-  matrix(1/7, nrow =1 , ncol = 7)
  X_tilde <-  matrix(1, ncol = 1, nrow = 1)
  Z_tilde <-  matrix(1, ncol = 1, nrow = 1)
  P_tilde <-  matrix(1/7, ncol = 7, nrow = 1)
  deriv_from_fn <- mu_d_beta(i = 1,
            j = 2,
            q = 1,
            gammas = 8,
            B = B,
            X = X,
            Z = Z,
            P = P,
            X_tilde = X_tilde,
            Z_tilde = Z_tilde,
            Z_tilde_gamma_cols = 1,
            P_tilde = P_tilde,
            gamma_tilde = matrix(-100, ncol = 1, nrow = 1))

 direct_deriv <- ((Z%*%P)*(exp(X%*%B + 8)))[1,2]
 expect_equal(as.numeric(deriv_from_fn), direct_deriv)
})
