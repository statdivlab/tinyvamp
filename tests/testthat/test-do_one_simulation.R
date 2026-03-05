test_that("do_one_simulation returns list of lists", {
  
  the_n_boot <- 3
  tinysim <- do_one_simulation(n = 1,
                               J = 5,
                               distrib = "NB",
                               B_multiplier = 0,
                               seed = 1,
                               label = "test",
                               n_boot = the_n_boot,
                               parallelize=FALSE,
                               verbose =FALSE,
                               folder_name = "test",
                               return_dont_save = TRUE)
  expect_equal(length(tinysim$poisson_lrt$boot_lr_stats), the_n_boot)
  expect_equal(length(tinysim$reweighted_lrt$boot_lr_stats), the_n_boot)
  
  expect_true("data.frame" %in% class(tinysim$poisson_ci$ci))
  expect_true("data.frame" %in% class(tinysim$reweighted_ci$ci))
  
})

test_that("Simulations that previously failed now succeed.",{
  
  skip("This test is very slow")
  
  failing_on_bayes <- do_one_simulation(n = 1,
                                        J = 5,
                                        distrib = "Poisson",
                                        B_multiplier = 1,
                                        seed = 1,
                                        label = "test",
                                        n_boot = 5,
                                        parallelize= FALSE,
                                        folder_name = "test",
                                        verbose = FALSE,
                                        return_dont_save = TRUE)
  
  expect_type(failing_on_bayes$poisson_lrt,"list")
  expect_type(failing_on_bayes$poisson_ci,"list")
  expect_type(failing_on_bayes$reweighted_lrt,"list")
  expect_type(failing_on_bayes$reweighted_ci,"list")
  
  stalling_out_on_bayes <- do_one_simulation(n = 3,
                                             J = 20,
                                             distrib = "Poisson",
                                             B_multiplier = 0,
                                             seed = 67791629,
                                             label = "trying_for_hundo",
                                             n_boot = 10,
                                             parallelize= FALSE,
                                             folder_name = "test",
                                             verbose = FALSE,
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
