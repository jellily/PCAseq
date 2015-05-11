# Tests for filterMaf


test <- matrix(c(1, 2, 1, 0 , 1, 1, 0, 0 ,0, 0, 0, 0), ncol = 3)
test_that("filterMaf removes snps by maf",
{
  expect_identical(filterMaf(test, 0.2), test[1:2,])
})

test <- matrix(c(1, 2, 1, 0 , 1, 1, 0, 0 ,0, 1, 1, 0), ncol = 3)
test_that("filterMaf removes snps by maf",
{
  expect_identical(filterMaf(test, 0), test[1:3,])
})

test <- matrix(c(2, 0, 2, 2 , 1, 1, 0, 0 ,2, 2, 2, 1), ncol = 3)
test_that("filterMaf removes snps by maf",
{
  expect_identical(filterMaf(test, c(0.1, 0.4)), test[c(1,3),]) 
})


test <- matrix(c(1, 0, 0, 1 , 0, 0, 0, 0 ,1, 0, 0, 1), ncol = 3)
test_that("filterMaf returns an error when all data is removed",
{
  expect_error(filterMaf(test, c(0.1, 0.3)))
  expect_error(filterMaf(test, 1))
})
