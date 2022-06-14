test_that("When MLE has closed form, estimate_parameters finds it", {

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


  numerical_result <- estimate_parameters(W = W,
                      X = X,
                      Z = Z,
                      Z_tilde = Z_tilde,
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
                      barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                      barrier_scale = 10, #increments for value of barrier penalty
                      max_barrier = 1e12, #maximum value of barrier_t
                      initial_conv_tol = 1000,
                      final_conv_tol = 0.1,
                      constraint_tolerance = 1e-15,
                      hessian_regularization = 0.01,
                      criterion = "Poisson",
                      profile_P = FALSE,
                      profiling_maxit = 25
  )
  P_analytical <- apply(W,2,sum) %>% (function(x) x/sum(x))

  P_numerical <- numerical_result$varying %>% with(value[param == "P"])

  expect_equal(P_analytical,P_numerical)

})

test_that("When MLE has closed form,reweighted estimator is close", {

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


  numerical_result <- estimate_parameters(W = W,
                                          X = X,
                                          Z = Z,
                                          Z_tilde = Z_tilde,
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
                                          barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                                          barrier_scale = 10, #increments for value of barrier penalty
                                          max_barrier = 1e12, #maximum value of barrier_t
                                          initial_conv_tol = 1000,
                                          final_conv_tol = 0.1,
                                          
                                          constraint_tolerance = 1e-15,
                                          hessian_regularization = 0.01,
                                          criterion = "reweighted_Poisson",
                                          
                                          profile_P = FALSE,
                                          profiling_maxit = 25,
                                          return_variance = TRUE
  )
  P_analytical <- apply(W,2,sum) %>% (function(x) x/sum(x))

  P_numerical <- numerical_result$varying %>% with(value[param == "P"])

  expect_equal(P_analytical,P_numerical,tolerance = .01)

  # numerical_result$variance_function %>%
  #   ggplot() +
  #   geom_point(aes(x = mean, y = squerror)) +
  #   geom_line(aes(x = mean, y = estd_var))+
  #   theme_bw()

})

test_that("When MLE has closed form, estimate_parameters finds it and profiling step doesn't fail", {

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


  numerical_result <- estimate_parameters(W = W,
                                          X = X,
                                          Z = Z,
                                          Z_tilde = Z_tilde,
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
                                          barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                                          barrier_scale = 10, #increments for value of barrier penalty
                                          max_barrier = 1e12, #maximum value of barrier_t
                                          initial_conv_tol = 1000,
                                          final_conv_tol = 0.1,
                                          constraint_tolerance = 1e-15,
                                          hessian_regularization = 0.01,
                                          
                                          profile_P = TRUE,
                                          profiling_maxit = 25
  )
  P_analytical <- apply(W,2,sum) %>% (function(x) x/sum(x))

  P_numerical <- numerical_result$varying %>% with(value[param == "P"])

  expect_equal(P_analytical,P_numerical)

})


test_that("When MLE has closed form at boundary, estimate_parameters finds it and profiling step doesn't fail", {

  set.seed(0)
  W <- matrix(sapply(1:10,function(x) rpois(1,1000)),
              ncol = 5)
  W[,3] <- 0 #taxon 3 not observed
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


  numerical_result <- estimate_parameters(W = W,
                                          X = X,
                                          Z = Z,
                                          Z_tilde = Z_tilde,
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
                                          barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                                          barrier_scale = 10, #increments for value of barrier penalty
                                          max_barrier = 1e12, #maximum value of barrier_t
                                          initial_conv_tol = 1000,
                                          final_conv_tol = 0.1,
                                          constraint_tolerance = 1e-10,
                                          hessian_regularization = 0.01,
                                          
                                          profile_P = TRUE,
                                          profiling_maxit = 25
  )
  P_analytical <- apply(W,2,sum) %>% (function(x) x/sum(x))

  P_numerical <- numerical_result$varying %>% with(value[param == "P"])

  expect_equal(P_analytical,P_numerical)

})


test_that("Setting all weights equal to 1 does not affect ability to
          find MLE (in case where it has closed form)", {

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


  numerical_result <- estimate_parameters(W = W,
                                          X = X,
                                          Z = Z,
                                          Z_tilde = Z_tilde,
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
                                          barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                                          barrier_scale = 10, #increments for value of barrier penalty
                                          max_barrier = 1e12, #maximum value of barrier_t
                                          initial_conv_tol = 1000,
                                          final_conv_tol = 0.1,
                                          constraint_tolerance = 1e-15,
                                          hessian_regularization = 0.01,
                                          
                                          profile_P = FALSE,
                                          profiling_maxit = 25,
                                          wts = rep(1,10)
  )
  P_analytical <- apply(W,2,sum) %>% (function(x) x/sum(x))

  P_numerical <- numerical_result$varying %>% with(value[param == "P"])

  expect_equal(P_analytical,P_numerical)

})

# test_that("GMM returns something reasonable in simple setting", {
#
# set.seed(0)
# W <- matrix(sapply(1:10,function(x) rpois(1,1000)),
#             ncol = 5)
# X <- matrix(1,ncol = 1, nrow = 2)
# Z <- matrix(1,nrow = 2, ncol = 1)
# Z_tilde <- matrix(0,nrow = 2, ncol = 1)
# Z_tilde_gamma_cols <- 1
# gammas <- apply(W,1,function(x) log(sum(x)))
# gammas_fixed_indices <- rep(F,2)
# P <- matrix(1/5, nrow = 1, ncol = 5)
# P_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 5)
# B <- matrix(0,ncol = 5, nrow = 1)
# B_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
# X_tilde <- matrix(0,ncol = 1, nrow = 1)
# P_tilde <- matrix(1/5,ncol = 5, nrow = 1)
# P_tilde_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
# gamma_tilde <- matrix(0,nrow = 1, ncol = 1)
# gamma_tilde_fixed_indices <- matrix(TRUE, nrow = 1, ncol = 1)
#
#
# numerical_result <- estimate_parameters(W = W,
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
#                                         max_barrier = 1e12, #maximum value of barrier_t
#                                         initial_conv_tol = 1000,
#                                         final_conv_tol = 0.1,
#                                         
#                                         constraint_tolerance = 1e-15,
#                                         hessian_regularization = 0.01,
#                                         criterion = "GMM",
#                                         
#                                         profile_P = FALSE,
#                                         profiling_maxit = 25
# )
# P_analytical <- apply(W,2,sum) %>% (function(x) x/sum(x))
#
# P_numerical <- numerical_result$varying %>% with(value[param == "P"])
#
# expect_equal(P_analytical,P_numerical,tolerance = .1)
#
#
# })
#
# test_that("GMM returns something reasonable in simple setting", {
#
#   set.seed(0)
#   n_samples <- 20
#   W <- lapply(1:n_samples, function(x)
#     matrix(rnbinom(5,mu = 300*x,size = 10), nrow = 1))
#   W <- do.call(rbind,W)
#
#   X <- matrix(1,ncol = 1, nrow = n_samples)
#   Z <- matrix(1,nrow = n_samples, ncol = 1)
#   Z_tilde <- matrix(0,nrow = n_samples, ncol = 1)
#   Z_tilde_gamma_cols <- 1
#   gammas <- apply(W,1,function(x) log(sum(x)))
#   gammas_fixed_indices <- rep(F,n_samples)
#   P <- matrix(1/5, nrow = 1, ncol = 5)
#   P_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 5)
#   B <- matrix(0,ncol = 5, nrow = 1)
#   B_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
#   X_tilde <- matrix(0,ncol = 1, nrow = 1)
#   P_tilde <- matrix(1/5,ncol = 5, nrow = 1)
#   P_tilde_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
#   gamma_tilde <- matrix(0,nrow = 1, ncol = 1)
#   gamma_tilde_fixed_indices <- matrix(TRUE, nrow = 1, ncol = 1)
#
#
#   numerical_result <- estimate_parameters(W = W,
#                                           X = X,
#                                           Z = Z,
#                                           Z_tilde = Z_tilde,
#                                           Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                                           gammas = gammas,
#                                           gammas_fixed_indices = gammas_fixed_indices,
#                                           P = P,
#                                           P_fixed_indices = P_fixed_indices,
#                                           B = B,
#                                           B_fixed_indices = B_fixed_indices,
#                                           X_tilde = X_tilde,
#                                           P_tilde = P_tilde,
#                                           P_tilde_fixed_indices = P_tilde_fixed_indices,
#                                           gamma_tilde = gamma_tilde,
#                                           gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
#                                           barrier_t = 1, #starting value of reciprocal barrier penalty coef.
#                                           barrier_scale = 10, #increments for value of barrier penalty
#                                           max_barrier = 1e12, #maximum value of barrier_t
#                                           initial_conv_tol = 1000,
#                                           final_conv_tol = 0.1,
#                                           
#                                           constraint_tolerance = 1e-15,
#                                           hessian_regularization = 0.01,
#                                           criterion = "GMM",
#                                           
#                                           profile_P = FALSE,
#                                           profiling_maxit = 25
#   )
#   P_pois <- apply(W,2,sum) %>% (function(x) x/sum(x))
#
#   P_gmm <- numerical_result$varying %>% with(value[param == "P"])
#
#   expect_true(sd(P_pois)>sd(P_gmm))
#
#
# })

#
# test_that("When MLE has closed form, GMM also works", {
#
#   set.seed(0)
#   W <- matrix(sapply(1:10,function(x) rpois(1,1000)),
#               ncol = 5)
#   X <- matrix(1,ncol = 1, nrow = 2)
#   Z <- matrix(1,nrow = 2, ncol = 1)
#   Z_tilde <- matrix(0,nrow = 2, ncol = 1)
#   Z_tilde_gamma_cols <- 1
#   gammas <- apply(W,1,function(x) log(sum(x)))
#   gammas_fixed_indices <- rep(F,2)
#   P <- matrix(1/5, nrow = 1, ncol = 5)
#   P_fixed_indices <- matrix(FALSE, nrow = 1, ncol = 5)
#   B <- matrix(0,ncol = 5, nrow = 1)
#   B_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
#   X_tilde <- matrix(0,ncol = 1, nrow = 1)
#   P_tilde <- matrix(1/5,ncol = 5, nrow = 1)
#   P_tilde_fixed_indices <- matrix(TRUE, ncol = 5, nrow = 1)
#   gamma_tilde <- matrix(0,nrow = 1, ncol = 1)
#   gamma_tilde_fixed_indices <- matrix(TRUE, nrow = 1, ncol = 1)
#
#
#   numerical_result <- estimate_parameters(W = W,
#                                           X = X,
#                                           Z = Z,
#                                           Z_tilde = Z_tilde,
#                                           Z_tilde_gamma_cols = Z_tilde_gamma_cols,
#                                           gammas = gammas,
#                                           gammas_fixed_indices = gammas_fixed_indices,
#                                           P = P,
#                                           P_fixed_indices = P_fixed_indices,
#                                           B = B,
#                                           B_fixed_indices = B_fixed_indices,
#                                           X_tilde = X_tilde,
#                                           P_tilde = P_tilde,
#                                           P_tilde_fixed_indices = P_tilde_fixed_indices,
#                                           gamma_tilde = gamma_tilde,
#                                           gamma_tilde_fixed_indices = gamma_tilde_fixed_indices,
#                                           barrier_t = 1, #starting value of reciprocal barrier penalty coef.
#                                           barrier_scale = 10, #increments for value of barrier penalty
#                                           max_barrier = 1e12, #maximum value of barrier_t
#                                           initial_conv_tol = 1000,
#                                           final_conv_tol = 0.1,
#                                           
#                                           constraint_tolerance = 1e-15,
#                                           hessian_regularization = 0.01,
#                                           
#                                           profile_P = TRUE,
#                                           profiling_maxit = 25,
#                                           criterion = "GMM"
#   )
#   P_analytical <- apply(W,2,sum) %>% (function(x) x/sum(x))
#
#   P_numerical <- numerical_result$varying %>% with(value[param == "P"])
#
#   expect_equal(P_analytical,P_numerical,tolerance = .1)
#
# })


test_that("Contamination setting with alpha_tilde works", {

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


  numerical_result <-
    estimate_parameters(W = W,
                        X = X,
                        Z = Z,
                        Z_tilde = NULL,
                        Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                        Z_tilde_list = Z_tilde_list,
                        alpha_tilde = alpha_tilde,
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
                        barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                        barrier_scale = 10, #increments for value of barrier penalty
                        max_barrier = 1e12, #maximum value of barrier_t
                        initial_conv_tol = 1000,
                        final_conv_tol = 0.1,
                        constraint_tolerance = 1e-15,
                        hessian_regularization = 0.01,
                        criterion = "Poisson",
                        
                        profile_P = FALSE,
                        profiling_maxit = 25
  )

  P_numerical <- numerical_result$varying %>% with(value[param == "P"])

  expect_equal(P_numerical, p_true, tolerance = 1e-1)

  P_tilde_numerical <- numerical_result$varying %>%
    with(value[param == "P_tilde"])

  expect_equal(P_tilde_numerical, p_contam, tolerance = .01)

})


test_that("Contamination setting with alpha_tilde works when alpha_tilde is not zero", {

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
        W[i,] <- round(rexp(1,3^((i- 1)%%5)*1/10000)*(p_mock2 + exp(0.25)*dilutions[i]*p_contam),0)

      } else{
        W[i,] <- round(rexp(1,3^((i-1)%%5)*1/10000)*(p_true + exp(-1.5)*dilutions[i]*p_contam),0)

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


  numerical_result <-
    estimate_parameters(W = W,
                        X = X,
                        Z = Z,
                        Z_tilde = NULL,
                        Z_tilde_gamma_cols = Z_tilde_gamma_cols,
                        Z_tilde_list = Z_tilde_list,
                        alpha_tilde = alpha_tilde,
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
                        barrier_t = 1, #starting value of reciprocal barrier penalty coef.
                        barrier_scale = 10, #increments for value of barrier penalty
                        max_barrier = 1e12, #maximum value of barrier_t
                        initial_conv_tol = 1000,
                        final_conv_tol = 0.1,
                        
                        constraint_tolerance = 1e-15,
                        hessian_regularization = 0.01,
                        criterion = "Poisson",
                        
                        profile_P = FALSE,
                        profiling_maxit = 25
    )

  P_numerical <- numerical_result$varying %>% with(value[param == "P"])

  expect_equal(P_numerical, p_true, tolerance = 1e-1)

  P_tilde_numerical <- numerical_result$varying %>%
    with(value[param == "P_tilde"])

  expect_equal(P_tilde_numerical, p_contam, tolerance = .01)

  alpha_tildes <-
    numerical_result$varying %>% with(value[param == "alpha_tilde"])

  expect_equal(alpha_tildes, c(0.25,-1.5), tolerance = .01)

})
