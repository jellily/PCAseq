# Tests for check.bool

test_that("check.bool passes when the input is a logical variable",
{
  expect_true(check.bool(TRUE))
  expect_true(check.bool(FALSE))
})

test_that("check.bool fails when the input is a logical variable", 
{
  expect_error(check.bool("TRUE"), "Input is not logical. See help for details.")
  expect_error(check.bool(42),"Input is not logical. See help for details.")
  expect_error(check.bool(c("cake", "cookie")), "Input is not logical. See help for details.")
  expect_error(check.bool(NA), "Input set to NA. See help for details.")
  expect_error(check.bool(NaN), "Input is not logical. See help for details.")
  expect_error(check.bool(NULL), "Input is not logical. See help for details.")
})