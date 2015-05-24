# Tests for filterAuto


test <- matrix(c(1, 2, 1, 0 , 1, 1, 0, 0 ,0, 0, 0, 0), ncol = 3)
test_that("filterAuto removes non-autosomal snps",
{
  expect_identical(filterAuto(test, snpChromosome = 1:4), test)
  expect_identical(filterAuto(test, snpChromosome = c(1:3, 23)), test[1:3, ])
  expect_identical(filterAuto(test, snpChromosome = c(1:3, "M")), test[1:3, ])
})
