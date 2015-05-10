# Tests for check.bool

test_that("getSnps returns the appropriate list of vectors",
{
  expect_equal(getSnps(1, 3, 12), 1:3)
  expect_equal(getSnps(2, 3, 12), 4:6)
  expect_equal(getSnps(2, 3, 5), 4:5)
})
