# Tests for calcMaf

test <- c(0, 0.1, 0.75, 1)
test_that("filterMaf removes snps by maf",
{
  expect_equal(calcMaf(test), c(0, 0.1, 0.25, 0))
})