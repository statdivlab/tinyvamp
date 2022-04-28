test_that("safe_divide works when numerator = denominator = 0", {

  expect_equal(safe_divide(0,0),1)
})

test_that("safe_divide works when numerator and denominator are positive", {
  expect_equal(safe_divide(1,2),0.5)
})

test_that("safe_divide returns correct penalty when numerator nonzero and denominator 0", {
  expect_equal(safe_divide(1,0,100),100)
})
