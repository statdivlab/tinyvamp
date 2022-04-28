test_that("mean makes sense in sample-read-only setting", {
  gammas <- log(3*(1:5))
  B <- matrix(rep(0,3),nrow = 1)
  X <- matrix(1,nrow = 5,ncol = 1)
  Z <- matrix(1,nrow = 5, ncol = 1)
  P <- matrix((1:3)/6,nrow = 1)
  X_tilde <- matrix(0,ncol = 1, nrow = 1)
  P_tilde <- P
  Z_tilde <- Z*0
  gamma_tilde <- matrix(0,ncol = 1, nrow = 1)
  means <- meaninate(gammas = gammas,
            B = B,
            Z = Z,
            X = Z,
            P = P,
            Z_tilde = Z_tilde,
            X_tilde = X_tilde,
            Z_tilde_gamma_cols = 1,
            P_tilde = P_tilde,
            gamma_tilde = gamma_tilde)

  #make sure meaninate returns a matrix
  expect_true(is.matrix(means))
  #make sure its dimensions are correct
  expect_equal(dim(means),c(5,3))
  #test that P is being used correctly here
  expect_equal(means[1,3]/means[1,1],3)
  #test that gammas are being used correctly here
  expect_equal

})
