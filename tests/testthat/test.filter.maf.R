# Tests for filter.maf


test <- matrix(c(1, 2, 1, 0 , 1, 1, 0, 0 ,0, 0, 0, 0), ncol = 3)
test_that("filter.maf removes snps by maf",
{
  expect_identical(filter.maf(test, 0.2), test[,1:2])
})

test <- matrix(c(1, 2, 1, 0 , 1, 1, 0, 0 ,0, 1, 1, 0), ncol = 3)
test_that("filter.maf removes snps by maf",
{
  expect_identical(filter.maf(test, 0), test)
})

test <- matrix(c(2, 2, 2, 2 , 1, 1, 0, 0 ,2, 2, 2, 1), ncol = 3)
test_that("filter.maf removes snps by maf",
{
  expect_identical(filter.maf(test, c(0.1, 0.4)), test[,2:3]) 
})


test <- matrix(c(2, 2, 2, 2 , 1, 1, 0, 0 ,2, 2, 2, 1), ncol = 3)
test_that("filter.maf returns an error when all data is removed",
{
  expect_error(filter.maf(test, c(0.1, 0.2)))
  expect_error(filter.maf(test, 1))
})
