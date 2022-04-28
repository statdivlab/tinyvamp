# test_that("same optimum achieved for no contamination case", {

  # set.seed(0)
  # B <- matrix(c(-3:3,0),nrow = 1)
  # X <- matrix(1,ncol = 1,
  #             nrow = 10)
  # P <- matrix(rep(1/8,8),nrow = 1)
  # Z <- matrix(1,nrow = 10, ncol = 1)
  # Z_tilde <- matrix(1,nrow = 1,
  #                   ncol =1)
  # X_tilde <- matrix(0,ncol = 1, nrow = 1)
  # gammas <- matrix(rnorm(10,12),ncol = 1)
  #
  # means_W <- (Z%*%P)*exp(gammas%*%matrix(1, nrow = 1, ncol = 8) + X%*%B)
  #
  # W <- apply(means_W, c(1,2), function(x) rpois(1, x))
  #
  #
  # X <- matrix(rep(1),nrow = nrow(W))
  # X_tilde <- matrix(0,nrow = 1, ncol = 1)
  # Z_tilde <- matrix(0, nrow = nrow(W), ncol = 1)
  # P_tilde <-matrix( rep(1/ncol(W),ncol(W)),
  #                   nrow = 1)
  #
  # gamma_tilde <- matrix(0,nrow = 1)
  #
  # gammas_fixed_indices <- matrix(rep(F, length(gammas)),ncol = 1)
  #
  # gamma_tilde_fixed_indices <- rep(T, length(gamma_tilde))
  #
  # J <- ncol(W)
  #
  # P_fixed_indices <- matrix(ncol = ncol(P),
  #                           nrow = nrow(P),
  #                           TRUE)
  #
  # P_tilde_fixed_indices <- matrix(ncol = ncol(P_tilde),
  #                                 nrow = nrow(P_tilde),
  #                                 TRUE)
  #
  #
  # B_fixed_indices <- matrix(ncol = ncol(B),
  #                           nrow = nrow(B),
  #                           TRUE)
  # B_fixed_indices[,(J - 7):(J -1)] <- FALSE
  # Z_tilde_gamma_cols <- 1
  #
  # estimate_newton <- estimate_parameters(W = W,
  #                                          X = X,
  #                                          Z = Z,
  #                                          Z_tilde = Z_tilde,
  #                                          Z_tilde_gamma_cols = Z_tilde_gamma_cols,
  #                                          gammas = gammas,
  #                                          gammas_fixed_indices = gammas_fixed_indices,
  #                                          P = P,
  #                                          P_fixed_indices = P_fixed_indices,
  #                                          B = B,
  #                                          B_fixed_indices = B_fixed_indices,
  #                                          X_tilde = X_tilde,
  #                                          P_tilde = P_tilde,
  #                                          P_tilde_fixed_indices = P_tilde_fixed_indices,
  #                                          gamma_tilde = gamma_tilde,
  #                                          gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
  #                                          barrier_t = 1, #starting value of reciprocal barrier penalty coef.
  #                                          barrier_scale = 10, #increments for value of barrier penalty
  #                                          max_barrier = 1e12, #maximum value of rbpc
  #                                          final_epsilon = .01,
  #                                          final_f = 1e-9,
  #                                          subproblem_method = "Newton")
#
#   estimate_trust <- estimate_parameters(W = W,
#                                          X = X,
#                                          Z = Z,
#                                          Z_tilde = Z_tilde,
#                                          Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                                          gammas = gammas,
#                                          gammas_fixed_indices = gammas_fixed_indices,
#                                          P = P,
#                                          P_fixed_indices = P_fixed_indices,
#                                          B = B,
#                                          B_fixed_indices = B_fixed_indices,
#                                          X_tilde = X_tilde,
#                                          P_tilde = P_tilde,
#                                          P_tilde_fixed_indices = P_tilde_fixed_indices,
#                                          gamma_tilde = gamma_tilde,
#                                          gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
#                                          barrier_t = 1, #starting value of reciprocal barrier penalty coef.
#                                          barrier_scale = 10, #increments for value of barrier penalty
#                                          max_barrier = 1e12, #maximum value of rbpc
#                                          final_epsilon = .01,
#                                          final_f = 1e-8,
#                                          subproblem_method = "trust")
#
#   estimate_bfgs <- estimate_parameters(W = W,
#                                         X = X,
#                                         Z = Z,
#                                         Z_tilde = Z_tilde,
#                                         Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                                         gammas = gammas,
#                                         gammas_fixed_indices = gammas_fixed_indices,
#                                         P = P,
#                                         P_fixed_indices = P_fixed_indices,
#                                         B = B,
#                                         B_fixed_indices = B_fixed_indices,
#                                         X_tilde = X_tilde,
#                                         P_tilde = P_tilde,
#                                         P_tilde_fixed_indices = P_tilde_fixed_indices,
#                                         gamma_tilde = gamma_tilde,
#                                         gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
#                                         barrier_t = 1, #starting value of reciprocal barrier penalty coef.
#                                         barrier_scale = 10, #increments for value of barrier penalty
#                                         max_barrier = 1e12, #maximum value of rbpc
#                                         final_epsilon = .01,
#                                         final_f = 1e-8,
#                                         subproblem_method = "BFGS")
#
#
#
#
# expect_equal(estimate_newton$varying$value,
#              estimate_trust$varying$value)
# expect_equal(estimate_trust$varying$value,
#              estimate_bfgs$varying$value,
#              tolerance = 1e-5)
# }
# )
#
# test_that("same optimum achieved, contamination case", {
#
#   set.seed(0)
#   B <- matrix(c(-3:3,0),nrow = 1)
#   X <- matrix(1,ncol = 1,
#               nrow = 10)
#   P <- rbind(matrix(c(0.05,0.2,0.05,0.2,0.05,0.2,0.05,.2),nrow = 1),
#              matrix(c(0.2,0.05,0.2,0.05,0.2,0.05,.2,0.05),nrow = 1))
#
#   Z <- cbind(matrix(rep(c(0,1),5),nrow = 10, ncol = 1),
#              matrix(rep(c(1,0),5),nrow = 10, ncol = 1))
#   Z_tilde <- matrix(1,ncol = 1, nrow = 10)
#   P_tilde <- matrix((1:8)/(sum(1:8)),nrow = 1)
#   X_tilde <- matrix(0,ncol = 1, nrow = 1)
#   gammas <- matrix(rnorm(10,12),ncol = 1)
#   gamma_tilde <- matrix(-1,ncol = 1, nrow = 1)
#
#   means_W <- (Z%*%P)*exp(gammas%*%matrix(1, nrow = 1, ncol = 8) + X%*%B) +
#     (exp(gammas)*Z_tilde)%*%(P_tilde*exp(gamma_tilde%*%matrix(1,ncol = 8,nrow = 1)))
#
#   W <- apply(means_W, c(1,2), function(x) rpois(1, x))
#
#
#   gammas_fixed_indices <- matrix(rep(F, length(gammas)),ncol = 1)
#
#   gamma_tilde_fixed_indices <- rep(F, length(gamma_tilde))
#
#   J <- ncol(W)
#
#   P_fixed_indices <- matrix(ncol = ncol(P),
#                             nrow = nrow(P),
#                             TRUE)
#
#   P_tilde_fixed_indices <- matrix(ncol = ncol(P_tilde),
#                                   nrow = nrow(P_tilde),
#                                   FALSE)
#
#
#   B_fixed_indices <- matrix(ncol = ncol(B),
#                             nrow = nrow(B),
#                             TRUE)
#   B_fixed_indices[,(J - 7):(J -1)] <- FALSE
#   Z_tilde_gamma_cols <- 1
#
#   ### set starting values for parameters to be estimated
#   P_tilde_start <- P_tilde
#   P_tilde_start[] <- rexp(8) %>% (function(x) x/sum(x))
#
#   gammas_start <- gammas
#   gammas_start[] <- apply(W,1, function(x) log(sum(x)))
#
#   B_start <- B
#   B_start[] <- 0
#
#   gamma_tilde_start <- gamma_tilde
#   gamma_tilde_start[] <- 0
#
#
#   estimate_newton <-  estimate_parameters(W = W,
#                                           X = X,
#                                           Z = Z,
#                                           Z_tilde = Z_tilde,
#                                           Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                                           gammas = gammas_start,
#                                           gammas_fixed_indices = gammas_fixed_indices,
#                                           P = P,
#                                           P_fixed_indices = P_fixed_indices,
#                                           B = B_start,
#                                           B_fixed_indices = B_fixed_indices,
#                                           X_tilde = X_tilde,
#                                           P_tilde = P_tilde_start,
#                                           P_tilde_fixed_indices = P_tilde_fixed_indices,
#                                           gamma_tilde = gamma_tilde,
#                                           gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
#                                           barrier_t = 1, #starting value of reciprocal barrier penalty coef.
#                                           barrier_scale = 10, #increments for value of barrier penalty
#                                           max_barrier = 1e20, #maximum value of rbpc
#                                           final_epsilon = .01,
#                                           final_f = 1e-9,
#                                           subproblem_method = "Newton")
#
#   estimate_trust <- estimate_parameters(W = W,
#                                         X = X,
#                                         Z = Z,
#                                         Z_tilde = Z_tilde,
#                                         Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                                         gammas = gammas_start,
#                                         gammas_fixed_indices = gammas_fixed_indices,
#                                         P = P,
#                                         P_fixed_indices = P_fixed_indices,
#                                         B = B_start,
#                                         B_fixed_indices = B_fixed_indices,
#                                         X_tilde = X_tilde,
#                                         P_tilde = P_tilde_start,
#                                         P_tilde_fixed_indices = P_tilde_fixed_indices,
#                                         gamma_tilde = gamma_tilde,
#                                         gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
#                                         barrier_t = 1, #starting value of reciprocal barrier penalty coef.
#                                         barrier_scale = 10, #increments for value of barrier penalty
#                                         max_barrier = 1e20, #maximum value of rbpc
#                                         final_epsilon = .01,
#                                         final_f = 1e-9,
#                                         subproblem_method = "trust")
#
#   estimate_bfgs <- estimate_parameters(W = W,
#                                        X = X,
#                                        Z = Z,
#                                        Z_tilde = Z_tilde,
#                                        Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                                        gammas = gammas_start,
#                                        gammas_fixed_indices = gammas_fixed_indices,
#                                        P = P,
#                                        P_fixed_indices = P_fixed_indices,
#                                        B = B_start,
#                                        B_fixed_indices = B_fixed_indices,
#                                        X_tilde = X_tilde,
#                                        P_tilde = P_tilde_start,
#                                        P_tilde_fixed_indices = P_tilde_fixed_indices,
#                                        gamma_tilde = gamma_tilde,
#                                        gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
#                                        barrier_t = 1, #starting value of reciprocal barrier penalty coef.
#                                        barrier_scale = 10, #increments for value of barrier penalty
#                                        max_barrier = 1e20, #maximum value of rbpc
#                                        final_epsilon = .01,
#                                        final_f = 1e-9,
#                                        subproblem_method = "BFGS")
#
#
#
#
#   expect_equal(estimate_newton$varying$value,
#                estimate_trust$varying$value,
#                tolerance = 1e-1)
#   expect_equal(estimate_trust$varying$value,
#                estimate_bfgs$varying$value,
#                tolerance = 1e-1)
#   expect_equal(estimate_newton$varying$value,
#                estimate_bfgs$varying$value,
#                tolerance = 1e-2)
# }
# )
# ### add tests for mean specification -- esp Z!
