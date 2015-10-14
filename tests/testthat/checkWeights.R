# Tests for checkWeights

test_that("checkWeights passes when the input is a proper set of weights",
{
  expect_true(checkWeights(c(0.213, 0.467)))
  expect_true(checkWeights(1:2))
})

test_that("checkWeights fails when the input is not correct", 
{
  expect_error(checkWeights("cake"))
  expect_error(checkWeights(42))
  expect_error(checkWeights(c(0, 1)))
  expect_error(checkWeights(c(-10, 1)))
  expect_error(checkMaf(NaN))
  expect_error(checkWeights(NULL))
})