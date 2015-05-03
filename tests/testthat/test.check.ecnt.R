# Tests for check.maf

test_that("check.miss passes when the input is a proper missingness rate",
{
  expect_true(check.ecnt(1))
  expect_true(check.ecnt(57))
})

test_that("check.miss fails when the input is not correct", 
{
  expect_error(check.ecnt("cake"))
  expect_error(check.ecnt(1:3))
  expect_error(check.ecnt(0))
  expect_error(check.ecnt(-5))
  expect_error(check.ecnt(NA))
  expect_error(check.ecnt(NULL))
})