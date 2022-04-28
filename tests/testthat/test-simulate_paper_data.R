test_that("Simulated data is reproducible", {
  W1 <- simulate_paper_data(n = 1,
                           J = 5,
                           B_multiplier = 1,
                           distrib = "Poisson",
                           seed = 0)

  W2 <- simulate_paper_data(n = 1,
                           J = 5,
                           B_multiplier = 1,
                           distrib = "Poisson",
                           seed = 0)

  expect_equal(W1,W2)

  W1 <- simulate_paper_data(n = 1,
                            J = 5,
                            B_multiplier = 1,
                            distrib = "NB",
                            seed = 0)

  W2 <- simulate_paper_data(n = 1,
                            J = 5,
                            B_multiplier = 1,
                            distrib = "NB",
                            seed = 0)

  expect_equal(W1,W2)
})

test_that("Simulated data is reasonably unique", {
  W1 <- simulate_paper_data(n = 1,
                            J = 5,
                            B_multiplier = 1,
                            distrib = "Poisson",
                            seed = 0)

  W2 <- simulate_paper_data(n = 1,
                            J = 5,
                            B_multiplier = 1,
                            distrib = "Poisson",
                            seed = 1)

  expect_true(all(W1 != W2))

  W1 <- simulate_paper_data(n = 1,
                            J = 5,
                            B_multiplier = 1,
                            distrib = "NB",
                            seed = 0)

  W2 <- simulate_paper_data(n = 1,
                            J = 5,
                            B_multiplier = 1,
                            distrib = "NB",
                            seed = 1)

  expect_true(all(W1 != W2))
})


