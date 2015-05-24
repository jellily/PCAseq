# Tests for filterMono
test <- matrix(c(1, 2, 1, 0 , 1, 1, 0, 0 ,0, 0, 0, 0), ncol = 3)
test_that("filterMono removes monomorphic snps",
{
  expect_identical(filterMono(test), test[1:3,])
})


test <- matrix(c(1, 2, 1, 0 , 1, 1, 0, 1 ,0, 1, 1, 0), ncol = 3)
test_that("filterMono removes monomorphic snps",
{
  expect_identical(filterMono(test), test)
})


test <- matrix(c(2, 1, 2, 2 , 1, 1, 0, 2 ,2, 2, 2, 2), ncol = 3)
test_that("filterMono removes monomorphic snps",
{
  expect_identical(filterMono(test), test[1:3,])
})

test <- matrix(c(2, 0, 2, 1 , 1, 1, 1, 1 ,2, 2, 2, 1), ncol = 3)
test_that("filterMono removes monomorphic snps",
{
  expect_identical(filterMono(test), test)
})
