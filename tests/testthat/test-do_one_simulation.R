test_that("do_one_simulation returns list of lists", {

  tinysim <- do_one_simulation(n = 1,
                                J = 5,
                                distrib = "NB",
                                B_multiplier = 0,
                                seed = 1,
                                label = "test",
                                n_boot = 5,
                                verbose =FALSE,
                                # load_tinyvamp = FALSE,
                                folder_name = "test",
                                return_dont_save = TRUE)

  expect_type(tinysim$poisson_lrt,"list")
  expect_type(tinysim$poisson_ci,"list")
  expect_type(tinysim$reweighted_lrt,"list")
  expect_type(tinysim$reweighted_ci,"list")

})



# test_that("Simulation that was stalling on Bayes
# is not stalling",{
#
#   please_dont_stall <-
#     do_one_simulation(n = 1,
#                       J = 20,
#                       distrib = "NB",
#                       B_multiplier = 0,
#                       seed = 1,
#                       label = "testrun",
#                       n_boot = 5,
#                       load_tinyvamp = FALSE,
#                       folder_name,
#                       return_dont_save = TRUE,
#                       parallelize = TRUE,
#                       verbose= TRUE)
#
#   expect_type(please_dont_stall$poisson_lrt,"list")
#   expect_type(please_dont_stall$poisson_ci,"list")
#   expect_type(please_dont_stall$reweighted_lrt,"list")
#   expect_type(please_dont_stall$reweighted_ci,"list")
#
# }
#           )

test_that("Simulation that failed bc of log penalty
gradient off-by-one-type error now succeeds.",{

  failing_on_bayes <- do_one_simulation(n = 1,
                    J = 5,
                    distrib = "Poisson",
                    B_multiplier = 1,
                    seed = 1,
                    label = "test",
                    n_boot = 5,
                    parallelize= FALSE,
                    # load_tinyvamp = FALSE,
                    folder_name = "test",
                    verbose = FALSE,
                    return_dont_save = TRUE)

  expect_type(failing_on_bayes$poisson_lrt,"list")
  expect_type(failing_on_bayes$poisson_ci,"list")
  expect_type(failing_on_bayes$reweighted_lrt,"list")
  expect_type(failing_on_bayes$reweighted_ci,"list")

})


# test_that("Simulation that might be stalled on Bayes
# actually runs.",{
#
#   possibly_failing_on_bayes <- do_one_simulation(n = 3,
#                                         J = 5,
#                                         distrib = "NB",
#                                         B_multiplier = 1,
#                                         seed = 1,
#                                         label = "test",
#                                         n_boot = 5,
#                                         parallelize= FALSE,
#                                         load_tinyvamp = FALSE,
#                                         folder_name = "test",
#                                         verbose = TRUE,
#                                         return_dont_save = TRUE)
#
#   expect_type(possibly_failing_on_bayes$poisson_lrt,"list")
#   expect_type(possibly_failing_on_bayes$poisson_ci,"list")
#   expect_type(possibly_failing_on_bayes$reweighted_lrt,"list")
#   expect_type(possibly_failing_on_bayes$reweighted_ci,"list")
#
# })

# test_that("Yet another thing that stalled on Bayes doesn't stall anymore.", {
#
#   stalling_on_bayes <- do_one_simulation(n = 1,
#                                         J = 5,
#                                         distrib = "Poisson",
#                                         B_multiplier = 1,
#                                         seed = 1,
#                                         label = "test",
#                                         n_boot = 5,
#                                         parallelize= FALSE,
#                                         load_tinyvamp = FALSE,
#                                         folder_name = "test",
#                                         return_dont_save = TRUE,
#                                         verbose = TRUE)
#
#
# })

test_that("Another simulation stalled on Bayes can run and
          gives different estimates depending on estimator used.",{

  stalling_out_on_bayes <-
    do_one_simulation(n = 3,
                      J = 20,
                      distrib = "Poisson",
                      B_multiplier = 0,
                      seed = 67791629,
                      label = "trying_for_hundo",
                      n_boot = 10,
                      parallelize= FALSE,
                      # load_tinyvamp = FALSE,
                      folder_name = "test",
                      verbose = TRUE,
                      return_dont_save = TRUE)

  expect_true(mean(stalling_out_on_bayes$poisson_lrt$boot_lr_stats) !=
                mean(stalling_out_on_bayes$reweighted_lrt$boot_lr_stats))

  expect_true(sd(stalling_out_on_bayes$poisson_ci$ci$lower_ci -
    stalling_out_on_bayes$reweighted_ci$ci$lower_ci) >1e-4)

  expect_true(sd(stalling_out_on_bayes$poisson_ci$ci$upper_ci -
                   stalling_out_on_bayes$reweighted_ci$ci$upper_ci) >1e-4)

  expect_true(sd(stalling_out_on_bayes$poisson_ci$ci$value -
    stalling_out_on_bayes$reweighted_ci$ci$value)>1e-4)




})
