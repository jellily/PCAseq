# Tests for getIndex

test_that("getIndex returns the appropriate list of vectors",
{
  expect_equal(getIndex(1, 3, 12), 1:3)
  expect_equal(getIndex(2, 3, 12), 4:6)
  expect_equal(getIndex(2, 3, 5), 4:5)
})
