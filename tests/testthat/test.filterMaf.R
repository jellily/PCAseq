# Tests for filterMaf


test <- matrix(c(1, 2, 1, 0 , 1, 1, 0, 0 ,0, 0, 0, 0), ncol = 3)
test_that("filterMaf removes snps by maf",
{
  expect_equal(filterMaf(test, "(0, 0.2]"), 3)
})


test_that("filterMaf removes snps by maf",
{
  expect_equal(filterMaf(test, "[0, 0.2]"), 3:4)
})


test_that("filterMaf removes snps by maf",
{
  expect_equal(filterMaf(test, "(0, 0.5)"), c(1,3))
})


test_that("filterMaf removes snps by maf",
{
  expect_equal(filterMaf(test, "[0, 0.5]"), 1:4)
})


