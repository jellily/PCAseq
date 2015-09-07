# Tests for getMiss


test <- matrix(c(1, 2, 1, 0 , NA, 1, 0, 0 , 0, NA, 0, 0), ncol = 3)
test_that("getMiss calculates the correct missingness rate",
{
  expect_identical(getMissRate(test), c(1/3, 1/3, 0, 0))
})

test <- matrix(c(1, 2, 1, 0 , 2, 1, 0, 0 , 0, 1, 0, 0), ncol = 3)
test_that("filterMaf removes snps by maf",
{
  expect_identical(getMissRate(test), c(0, 0, 0, 0))
})

test <- matrix(rep(NA, 12), ncol = 3)
test_that("filterMaf removes snps by maf",
{
  expect_identical(getMissRate(test), c(1, 1, 1, 1))
})



