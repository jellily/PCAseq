# Tests for checkSamp

test_that("checkSamp passes when the input is a valid Samp id vector",
{
  expect_true(checkSamp(1:3, 1:4))
  expect_true(checkSamp(1:3, 1:3))
  expect_true(checkSamp(c("Samp1", "Samp2"), c("Samp1", "Samp2", "Samp3")))
})

test_that("checkSamp fails when the input is not a valid Sample id vector", 
{
  expect_error(checkSamp(c("Samp1", "Samp5"),c("Samp1", "Samp2", "Samp3")) , "Sample.id vector specified has sample IDs not in the data set.")
  expect_error(checkSamp(c(),c("Samp1", "Samp2", "Samp3")) , "Sample.id vector specified has length of 0.")
  expect_error(checkSamp(c("Samp1", "Samp2", "Samp3"),c("Samp1", "Samp2")) , "More sample IDs given than are in the genotype data set.")
  
})