# Tests for seqPCAClass

grm <- matrix(1:9, ncol = 3)
weights <- c(1, 1)
eigen.res <- eigen(grm)
sample.id <- c("A", "C", "C")
snp.id <- c("s1", "s2", "s3")

res <- seqPCAClass(grm, weights, maf = NaN, eigen.res, sample.id, snp.id,
                   eigenCnt = 2, FALSE)
test_that("seqPCAClass passes when the input is correct",
{
  expect_is(res, "seqPCAClass")
  expect_identical(res$genmat, NULL)
  expect_equal(length(res$eigenval), 2)
  expect_equal(dim(res$eigenvect)[2], 2)
})

res <- seqPCAClass(grm, method, maf = NaN, eigen.res, sample.id, snp.id,
                   eigenCnt = 2, TRUE)
test_that("seqPCAClass passes when the input is correct",
{
  expect_is(res, "seqPCAClass")
  expect_identical(res$genmat, grm)
  expect_equal(length(res$eigenval), 2)
  expect_equal(dim(res$eigenvect)[2], 2)
})

test_that("seqPCAClass passes when the input is correct",
{
  expect_warning(seqPCAClass(grm, method, maf = NaN, eigen.res, sample.id,
                             snp.id, eigenCnt = 5, FALSE),
                 "Number of eigenvectors and values to return is more than the
            dimensions of the GRM. All eigenvalues and vectors will be
            returned.")
})


res <- seqPCAClass(grm, method, maf = NaN, eigen.res, sample.id, snp.id, 
                   eigenCnt = 0, FALSE)
test_that("seqPCAClass passes when the input is correct",
{
  expect_is(res, "seqPCAClass")
  expect_identical(res$genmat, NULL)
  expect_equal(length(res$eigenval), length(sample.id))
  expect_equal(dim(res$eigenvect)[2], length(sample.id))
})


res <- seqPCAClass(grm, method, maf = NaN, eigen.res, sample.id, snp.id, 
                   eigenCnt = 0, FALSE)
test_that("seqPCAClass passes when the input is correct",
{
  expect_is(res, "seqPCAClass")
  expect_identical(res$genmat, NULL)
  expect_equal(length(res$eigenval), length(sample.id))
  expect_equal(dim(res$eigenvect)[2], length(sample.id))
})