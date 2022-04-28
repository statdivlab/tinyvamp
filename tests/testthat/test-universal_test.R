# test_that("multiplication works", {
#   n <- 10
#   sim_p_10 <- numeric(100)
#   set.seed(4939323)
#   for(sim in 1:100){
#     print(paste("Simulation ", sim, sep = "", collapse = ""))
#     W <- simulate_simple_data(matrix(0, nrow = 1, ncol = 2),
#                               distrib = "Poisson",
#                               n = n,
#                               gamma_mean = 11)
#
#     full_model <- fit_simple_model(W,
#                                      B_fixed_at_zero = FALSE)
#
#     null_model <- fit_simple_model(W,
#                                    B_fixed_at_zero = TRUE)
#
#
#
#
#     sim_p_10[sim] <- universal_test(W,
#                                     full_model,
#                                     null_model)
#
#
#   }
#
#   qqplot(sim_p_10,runif(10000),type = "s")
#
#   n <- 10
#   sim_p_10_nb <- numeric(100)
#   set.seed(4939323)
#   for(sim in 1:100){
#     print(paste("Simulation ", sim, sep = "", collapse = ""))
#     W <- simulate_simple_data(matrix(0, nrow = 1, ncol = 2),
#                               distrib = "nb10",
#                               n = n,
#                               gamma_mean = 11)
#
#     full_model <- fit_simple_model(W,
#                                    B_fixed_at_zero = FALSE)
#
#     null_model <- fit_simple_model(W,
#                                    B_fixed_at_zero = TRUE)
#
#
#
#
#     sim_p_10_nb[sim] <- universal_test(W,
#                                     full_model,
#                                     null_model,
#                                     parallelize = TRUE)
#
#     print(sim_p_10_nb[sim])
#
#
#   }
#
#   qqplot(sim_p_10,runif(10000),type = "s")
# })
