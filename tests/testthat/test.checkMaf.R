# Tests for checkMaf

test_that("checkMaf passes when the input is a proper MAF",
{
  expect_true(checkMaf(0.352))
  expect_true(checkMaf(c(0.213, 0.467)))
  expect_true(checkMaf(c(0, 0.5)))
  expect_true(checkMaf(NaN))
})

test_that("checkMaf fails when the input is not correct", 
{
  expect_error(checkMaf("cake"))
  expect_error(checkMaf(c(0.1, 0.2, 0.3)))
  expect_error(checkMaf(42))
  expect_error(checkMaf(c(-1, 0.4)))
  expect_error(checkMaf(c(0, 27)))
  expect_error(checkMaf(c(0.2, 0.2)))
  expect_error(checkMaf(c(0.2, 0.19)))
  expect_error(checkMaf(NA))
  expect_error(checkMaf(NULL))
})