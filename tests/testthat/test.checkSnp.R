# Tests for checkSnp

test_that("checkSnp passes when the input is a valid SNP id vector",
{
  expect_true(checkSnp(1:3, 1:4))
  expect_true(checkSnp(1:3, 1:3))
  expect_true(checkSnp(c("snp1", "snp2"), c("snp1", "snp2", "snp3")))
})

test_that("checkSnp fails when the input is not a valid SNP id vector", 
{
  expect_error(checkSnp(c("snp1", "snp5"),c("snp1", "snp2", "snp3")) , "Snp.id vector specified has SNP IDs not in the data set.")
  expect_error(checkSnp(c(),c("snp1", "snp2", "snp3")) , "Snp.id vector specified has length of 0.")
  expect_error(checkSnp(c("snp1", "snp2", "snp3"),c("snp1", "snp2")) , "More SNP IDs given than are in the genotype data set.")
  
})