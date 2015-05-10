# Tests for checkMiss

test_that("checkMiss passes when the input is a proper missingness rate",
{
  expect_true(checkMiss(0.352))
  expect_warning(checkMiss(1), TRUE)
  expect_message(checkMiss(0), TRUE)
  expect_identical(checkMiss(NaN), TRUE)
})

test_that("checkMiss fails when the input is not correct", 
{
  expect_error(checkMiss("cake"))
  expect_error(checkMiss(c(0.1, 0.2, 0.3)))
  expect_error(checkMiss(42))
  expect_error(checkMiss(-5))
  expect_error(checkMiss(NA))
  expect_error(checkMiss(NULL))
})