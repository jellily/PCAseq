# Tests for check.maf

test_that("check.miss passes when the input is a proper missingness rate",
{
  expect_true(check.miss(0.352))
  expect_warning(check.miss(1), TRUE)
  expect_message(check.miss(0), TRUE)
  expect_identical(check.miss(NaN), TRUE)
})

test_that("check.miss fails when the input is not correct", 
{
  expect_error(check.miss("cake"))
  expect_error(check.miss(c(0.1, 0.2, 0.3)))
  expect_error(check.miss(42))
  expect_error(check.miss(-5))
  expect_error(check.miss(NA))
  expect_error(check.miss(NULL))
})