# Tests for seqPCAClass

grm <- matrix(1:9, ncol = 3)
method <- "pcaseq"
eigen.res <- eigen(grm)
sample.id <- c("A", "C", "C")
snp.id <- c("s1", "s2", "s3")

res <- seqPCAClass(grm, method, eigen.res, sample.id, snp.id, eigen.cnt = 2, FALSE)
test_that("seqPCAClass passes when the input is correct",
{
  expect_is(res, "seqPCAClass")
  expect_identical(res$genmat, NULL)
  expect_equal(length(res$eigenval), 2)
  expect_equal(dim(res$eigenvect)[2], 2)
})

res <- seqPCAClass(grm, method, eigen.res, sample.id, snp.id, eigen.cnt = 2, TRUE)
test_that("seqPCAClass passes when the input is correct",
{
  expect_is(res, "seqPCAClass")
  expect_identical(res$genmat, grm)
  expect_equal(length(res$eigenval), 2)
  expect_equal(dim(res$eigenvect)[2], 2)
})

res <- seqPCAClass(grm, method, eigen.res, sample.id, snp.id, eigen.cnt = 5, FALSE)
test_that("seqPCAClass passes when the input is correct",
{
  expect_is(res, "seqPCAClass")
  expect_identical(res$genmat, NULL)
  expect_equal(length(res$eigenval), length(sample.id))
  expect_equal(dim(res$eigenvect)[2], length(sample.id))
})


res <- seqPCAClass(grm, method, eigen.res, sample.id, snp.id, eigen.cnt = 0, FALSE)
test_that("seqPCAClass passes when the input is correct",
{
  expect_is(res, "seqPCAClass")
  expect_identical(res$genmat, NULL)
  expect_equal(length(res$eigenval), length(sample.id))
  expect_equal(dim(res$eigenvect)[2], length(sample.id))
})