# Tests for check.maf

test_that("check.maf passes when the input is a proper MAF",
{
  expect_true(check.maf(0.352))
  expect_true(check.maf(c(0.213, 0.467)))
  expect_true(check.maf(c(0, 0.5)))
  expect_true(check.maf(NaN))
})

test_that("check.maf fails when the input is not correct", 
{
  expect_error(check.maf("cake"))
  expect_error(check.maf(c(0.1, 0.2, 0.3)))
  expect_error(check.maf(42))
  expect_error(check.maf(c(-1, 0.4)))
  expect_error(check.maf(c(0, 27)))
  expect_error(check.maf(c(0.2, 0.2)))
  expect_error(check.maf(c(0.2, 0.19)))
  expect_error(check.maf(NA))
  expect_error(check.maf(NULL))
})