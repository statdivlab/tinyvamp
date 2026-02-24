# test_that("We are able to fit a simulation model (Poisson likelihood) to
# Poisson simulation data", {
#   W <- simulate_paper_data(n = 1,
#                                   J = 5,
#                                   B_multiplier = 1,
#                                   distrib = "Poisson",
#                                   seed = 0)
# 
# 
#   fitted.model <- try(fit_simulation_model(W,"Poisson"))
#   expect_type(fitted.model,"list")
# })
# 
# test_that("We are able to fit a simulation model (Poisson likelihood) to
# negative binomial simulation data", {
#   W <- simulate_paper_data(n = 1,
#                            J = 5,
#                            B_multiplier = 1,
#                            distrib = "NB",
#                            seed = 0)
# 
# 
#   fitted.model <- try(fit_simulation_model(W,"Poisson"))
#   expect_type(fitted.model,"list")
# })
# 
# test_that("We are able to fit a simulation model (Poisson likelihood)
# to Poisson simulation data", {
#   W <- simulate_paper_data(n = 1,
#                            J = 5,
#                            B_multiplier = 1,
#                            distrib = "Poisson",
#                            seed = 0)
# 
# 
#   fitted.model <- try(fit_simulation_model(W,"reweighted_Poisson"))
#   expect_type(fitted.model,"list")
# })
# 
# test_that("We are able to fit a simulation model (reweighted Poisson likelihood) to
# negative binomial simulation data", {
#   W <- simulate_paper_data(n = 1,
#                            J = 5,
#                            B_multiplier = 1,
#                            distrib = "NB",
#                            seed = 0)
# 
# 
#   fitted.model <- try(fit_simulation_model(W,"reweighted_Poisson",
#                                            return_variance= TRUE))
#   expect_type(fitted.model,"list")
# 
#   # fitted.model$variance_function %>%
#   #   ggplot() +
#   #   geom_point(aes(x= mean, y = squerror),
#   #              size = 0.5) +
#   #   geom_line(aes(x = mean, y = estd_var),
#   #             color = "red") +
#   #   theme_bw()+
#   #   scale_y_sqrt() +
#   #   scale_x_sqrt()
# })
# 
# test_that("We get different estimates from Poisson and
# reweighted Poisson estimators fit to negative binomial data", {
#   W <- simulate_paper_data(n = 3,
#                            J = 20,
#                            B_multiplier = 1,
#                            distrib = "NB",
#                            seed = 0)
# 
# 
#   poisson_fit <- try(fit_simulation_model(W,"Poisson"))
#   reweighted_fit <- try(fit_simulation_model(W,"reweighted_Poisson",
#                                              return_variance = TRUE))
# 
#   expect_true(
#     mean(abs(poisson_fit$varying$value - reweighted_fit$varying$value))
#     >0.01)
# #
# #
# #   poisson_cis <- bootstrap_ci(W,
# #                               fitted_model = poisson_fit,
# #                               n_boot = 100,
# #                               verbose= TRUE,
# #                               parallelize = TRUE)
# #
# #   reweighted_cis <- bootstrap_ci(W,
# #                                  fitted_model = reweighted_fit,
# #                                  n_boot = 100,
# #                                  verbose = TRUE,
# #                                  parallelize = TRUE
# #                                  )
# #
# # poisson_cis$ci$method <- "Poisson"
# # reweighted_cis$ci$method <- "Reweighted"
# #
# # rbind(poisson_cis$ci,
# #       reweighted_cis$ci) %>%
# #   dplyr::filter(param == "B") %>%
# #   ggplot() +
# #   geom_errorbar(aes(x = j, ymin = lower_ci, ymax= upper_ci, color = method),
# #                 position = position_dodge(0.5)) +
# #   facet_wrap(~k)+
# #   theme_bw()
# #
# # rbind(poisson_cis$ci,
# #       reweighted_cis$ci) %>%
# #   dplyr::group_by(param, method) %>%
# #   dplyr::summarize(mean_width = mean(abs(upper_ci - lower_ci))) %>%
# #   ggplot() +
# #   geom_point(aes(x= param, y = mean_width, color = method),
# #              position = position_dodge(0.5)) +
# #   scale_y_log10()+
# #   theme_bw()
# #
# # rbind(poisson_cis$ci,
# #       reweighted_cis$ci) %>%
# #   mutate(width = upper_ci - lower_ci) %>%
# #   filter(param == "P_tilde") %>%
# #   select(param, j, method, width) %>%
# #   pivot_wider(values_from = width, names_from = method) %>%
# #   ggplot(aes(x = j, y = Poisson/Reweighted)) +
# #   geom_point() +
# #   theme_bw() +
# #   scale_y_log10()
# #
# #
# #
# # reweighted_fit$variance_function %>%
# #   ggplot() +
# #   geom_point(aes(x= mean,y = squerror)) +
# #   geom_line(aes(x = mean,y = estd_var),color="red") +
# #   scale_y_log10() +
# #   scale_x_log10() +
# #   theme_bw()
# 
# 
# })
