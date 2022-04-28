# # # # # test_that("multiplication works", {
# # # # #
# # # # n <- 10
# # # # sim_p_10 <- numeric(100)
# # # # set.seed(4939323)
# # # # for(sim in 1:100){
# # # #   print(paste("Simulation ", sim, sep = "", collapse = ""))
# # # #   W <- simulate_simple_data(matrix(0, nrow = 1, ncol = 2),
# # # #                             distrib = "nb10",
# # # #                             n = n,
# # # #                             gamma_mean = 11)
# # # #
# # # #   fitted_model <- fit_simple_model(W,
# # # #                                    B_fixed_at_zero = FALSE,
# # # #                                    reweight = TRUE)
# # # #
# # # #   null_param <- fitted_model
# # # #   null_param$B[] <- 0
# # # #   null_param$B_fixed_indices[] <- TRUE
# # # #
# # # #   boot_fit <- bootstrap_lrt(W = W,
# # # #                             fitted_model = fitted_model,
# # # #                             null_param = null_param,
# # # #                             n_boot = 100,
# # # #                             m = n^(3/4),
# # # #                             recalculate_W0 = FALSE,
# # # #                             parallelize = TRUE,
# # # #                             ncores = 7,
# # # #                             save_models = FALSE)
# # # #
# # # #   sim_p_10[sim] <- boot_fit$boot_pval
# # # #
# # # #   print(sim_p_10[sim])
# # # #   hist(boot_fit$boot_lr_stats,breaks = 7)
# # # #   abline(v = boot_fit$observed_lr_stat, col = "red")
# # # #
# # # # }
# # # # #
# # # # qs <- seq(0.01,.99,by = .01)
# # # # #
# # # # plot(qs, sapply(qs,function(k) quantile(sim_p_10,k)), type = "s",
# # # #      ylim = c(0,1),
# # # #      xlim = c(0,1))
# # # # points(qs, sapply(qs,function(k) quantileCI(sim_p_10,k,
# # # #                                            method = "asymptotic")$conf.int[1]),
# # # #      type = "s",lty = 2)
# # # # points(qs, sapply(qs,function(k) quantileCI(sim_p_10,k,
# # # #                                             method = "asymptotic")$conf.int[2]),
# # # #        type = "s",lty = 2)
# # # # abline(a = 0, b = 1, lty = 2)
# # # ################################ error in how weights are assigned?
# n <- 25
# sim_p_25 <- numeric(1000)
# sim_p_25_weighted <- numeric(1000)
# qs <- seq(0,1,by = .01)
# set.seed(4939323)
# boot_fits <- vector(1000,mode = "list")
# boot_fits_weighted <- boot_fits
# for(sim in 1:1000){
#   print(paste("Simulation ", sim, sep = "", collapse = ""))
#   W <- simulate_simple_data(matrix(c(0,0), nrow = 1, ncol = 2),
#                             distrib = "nb10",
#                             n = n,
#                             gamma_mean = 11)
#
#   fitted_model <- fit_simple_model(W = W,
#                                    B_fixed_at_zero = FALSE,
#                                    reweight = FALSE)
#
#   null_param <- fitted_model
#   null_param$B[] <- 0
#   null_param$B_fixed_indices[] <- TRUE
#
#   boot_fit <- bootstrap_lrt(W = W,
#                             fitted_model = fitted_model,
#                             null_param = null_param,
#                             n_boot = 1000,
#                             m = sqrt(m),
#                             recalculate_W0 = FALSE,
#                             parallelize = TRUE,
#                             ncores = 5,
#                             save_models = FALSE)
#
#   boot_fits[[sim]] <- boot_fit
#
#   sim_p_25[sim] <- boot_fit$boot_pval
#
#   fitted_model_weighted <- fit_simple_model(W = W,
#                                    B_fixed_at_zero = FALSE,
#                                    reweight = TRUE)
#
#   null_param <- fitted_model_weighted
#   null_param$B[] <- 0
#   null_param$B_fixed_indices[] <- TRUE
#
#   boot_fit_weighted <- bootstrap_lrt(W = W,
#                             fitted_model = fitted_model_weighted,
#                             null_param = null_param,
#                             n_boot = 1000,
#                             m = sqrt(m),
#                             recalculate_W0 = FALSE,
#                             parallelize = TRUE,
#                             ncores = 5,
#                             save_models = FALSE)
#
#   boot_fits_weighted[[sim]] <- boot_fit_weighted
#
#   sim_p_25_weighted[sim] <- boot_fit_weighted$boot_pval
#
#   print(sim_p_25[sim])
#   print(sim_p_25_weighted[sim])
#
#
#   plot(qs,sapply(qs, function(k) quantile(sim_p_25[1:sim],k)),
#                  type = "s",
#        xlim = c(0,1),
#        ylim = c(0,1))
#
#   lines(qs,sapply(qs, function(k) quantile(sim_p_25_weighted[1:sim],k)),
#        type = "s",
#        xlim = c(0,1),
#        ylim = c(0,1),
#        col = "red")
#
#   abline(a = 0, b = 1, lty = 2)
#
#
#
# }
#
# qs <- seq(0,1,by = .01)
# #
# plot(qs, sapply(qs,function(k) quantile(sim_p_25,k)), type = "s")
# abline(a = 0, b = 1, lty = 2)
# #
# # boot_lrs <- lapply(1:100, function(k) boot_fits[[k]]$boot_lr_stats)
# # obs_lrs <- sapply(1:100, function(k) boot_fits[[k]]$observed_lr_stat)
# # plot(1:100, log(sapply(1:100, function(k) median(boot_lrs[[k]]))/obs_lrs))
# # abline(a = 0, b= 1, lty = 2)
# # })
