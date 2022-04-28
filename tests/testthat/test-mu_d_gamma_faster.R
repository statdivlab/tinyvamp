
test_that("mu_d_gamma_faster output equal to mean when gamma_tilde negligible",{
set.seed(75943)
B <- matrix(c(rnorm(6),0),nrow = 1)
dgamma <- mu_d_gamma_faster(i = 1,
                  J = 7,
                  gammas = 1,
                  B = B,
                  X = matrix(1, ncol = 1),
                  Z = matrix(1, ncol = 1),
                  P = matrix((1:7)/(28),nrow = 1),
                  X_tilde = matrix(1, ncol = 1),
                  Z_tilde = matrix(1, ncol = 1),
                  Z_tilde_gamma_cols = NULL,
                  P_tilde = matrix((7:1)/28, nrow = 1),
                  gamma_tilde = -100)


mean_term <- meaninate(gammas = 1,
          B = B,
          X = matrix(1, ncol = 1),
          Z = matrix(1, ncol = 1),
          P = matrix((1:7)/(28),nrow = 1),
          X_tilde = matrix(1, ncol = 1),
          Z_tilde = matrix(1, ncol = 1),
          Z_tilde_gamma_cols = NULL,
          P_tilde = matrix((7:1)/28, nrow = 1),
          gamma_tilde = -100)

expect_equal(mean_term, dgamma)
})


test_that("mu_d_gamma_faster output equal to mean when gamma_tilde non-negligible
          and all columns of Z_tilde in Z_tilde_gamma_tilde_cols",{
  set.seed(742343)
  B <- matrix(c(rnorm(6),0),nrow = 1)
  dgamma <- mu_d_gamma_faster(i = 1,
                              J = 7,
                              gammas = 1,
                              B = B,
                              X = matrix(1, ncol = 1),
                              Z = matrix(1, ncol = 1),
                              P = matrix((1:7)/(28),nrow = 1),
                              X_tilde = matrix(1, ncol = 1),
                              Z_tilde = matrix(1, ncol = 1),
                              Z_tilde_gamma_cols = 1,
                              P_tilde = matrix((7:1)/28, nrow = 1),
                              gamma_tilde = 2)


  mean_term <- meaninate(gammas = 1,
                         B = B,
                         X = matrix(1, ncol = 1),
                         Z = matrix(1, ncol = 1),
                         P = matrix((1:7)/(28),nrow = 1),
                         X_tilde = matrix(1, ncol = 1),
                         Z_tilde = matrix(1, ncol = 1),
                         Z_tilde_gamma_cols = 1,
                         P_tilde = matrix((7:1)/28, nrow = 1),
                         gamma_tilde = 2)

  expect_equal(mean_term, dgamma)
})
