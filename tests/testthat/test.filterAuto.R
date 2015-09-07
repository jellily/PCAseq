# Tests for filterAuto
test <- matrix(c(1, 2, 1, 0 , 1, 1, 0, 0 ,0, 0, 0, 0), ncol = 3)
test_that("filterAuto removes non-autosomal snps",
{
  expect_identical(filterAuto(snpChromosome = 1:4), 1:4)
  expect_identical(filterAuto(snpChromosome = c(1:3, 23)), 1:3)
  expect_identical(filterAuto(snpChromosome = c(1:3, "M")), 1:3)
})
