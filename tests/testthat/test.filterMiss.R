# Tests for filterMiss

test <- matrix(c(1, 2, 1, 0 , 3, 1, 0, 0 ,0, 0, 3, 0), ncol = 3)
test_that("filterMiss removes snps with more than a specified missing rate",
{
  expect_identical(filterMiss(test, 0.4), 1:4)
  expect_identical(filterMiss(test, 0.1), 2:4)
})

