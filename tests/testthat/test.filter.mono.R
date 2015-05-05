# Tests for filter.mono


test <- matrix(c(1, 2, 1, 0 , 1, 1, 0, 0 ,0, 0, 0, 0), ncol = 3)
test_that("filter.mono removes monomorphic snps",
{
  expect_identical(filter.mono(test), test[,1:2])
})


test <- matrix(c(1, 2, 1, 0 , 1, 1, 0, 0 ,0, 1, 1, 0), ncol = 3)
test_that("filter.mono removes monomorphic snps",
{
  expect_identical(filter.mono(test), test)
})


test <- matrix(c(2, 2, 2, 2 , 1, 1, 0, 0 ,2, 2, 2, 1), ncol = 3)
test_that("filter.mono removes monomorphic snps",
{
  expect_identical(filter.mono(test), test[,2:3])
})

test <- matrix(c(2, 0, 2, 2 , 1, 1, 1, 1 ,2, 2, 2, 1), ncol = 3)
test_that("filter.mono removes monomorphic snps",
{
  expect_identical(filter.mono(test), test)
})
