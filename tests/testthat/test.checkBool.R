# Tests for checkBool

test_that("checkBool passes when the input is a logical variable",
{
  expect_true(checkBool(TRUE))
  expect_true(checkBool(FALSE))
})

test_that("checkBool fails when the input is a logical variable", 
{
  expect_error(checkBool("TRUE"), "Input is not logical. See help for details.")
  expect_error(checkBool(42),"Input is not logical. See help for details.")
  expect_error(checkBool(c("cake", "cookie")), "Input is not logical. See help for details.")
  expect_error(checkBool(NA), "Input set to NA. See help for details.")
  expect_error(checkBool(NaN), "Input is not logical. See help for details.")
  expect_error(checkBool(NULL), "Input is not logical. See help for details.")
})