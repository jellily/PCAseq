# Tests for checkEcnt

test_that("checkEcnt passes when the input is a proper missingness rate",
{
  expect_true(checkEcnt(1))
  expect_true(checkEcnt(57))
})

test_that("checkEcnt fails when the input is not correct", 
{
  expect_error(checkEcnt("cake"))
  expect_error(checkEcnt(1:3))
  expect_error(checkEcnt(0))
  expect_error(checkEcnt(-5))
  expect_error(checkEcnt(NA))
  expect_error(checkEcnt(NULL))
})