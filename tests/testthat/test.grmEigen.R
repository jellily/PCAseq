# Test for grmEigen

# Exmple 1
# read in the solutions
snprelate.res <- readRDS("snprelate-example-1.txt")
my.res <- readRDS("eigen-example-1.txt")

# open the test file
test.file <- snpgdsOpen("ex1.gds")
samp.id <- read.gdsn(index.gdsn(test.file, "sample.id"))
snp.id <- read.gdsn(index.gdsn(test.file, "snp.id"))

test_that("grmEigen returns the appropriate GRM matrix and that it is equivalent to SNPRelate",
{
  expect_equal(grmEigen(test.file, sampleId = samp.id, snpId = snp.id, autosomeOnly = FALSE, 
                             removeMonosnp = FALSE, maf = NaN, missingRate = NaN, transpose = FALSE), 
                   my.res)
  expect_equal(grmEigen(test.file, sampleId = samp.id, snpId = snp.id, autosomeOnly = FALSE, 
                         removeMonosnp = FALSE, maf = NaN, missingRate = NaN, transpose = FALSE), 
               snprelate.res)
})
snpgdsClose(test.file)

# 
# 
# # Exmple 1 Transposed
# # read in the solutions
# snprelate.res <- readRDS("snprelate-example-1t.txt")
# my.res <- readRDS("eigen-example-1t.txt")
# 
# # open the test file
# test.file <- snpgdsOpen("ex1t.gds")
# samp.id <- read.gdsn(index.gdsn(test.file, "sample.id"))
# snp.id <- read.gdsn(index.gdsn(test.file, "snp.id"))
# 
# test_that("grmEigen returns the appropriate GRM matrix and that it is equivalent to SNPRelate",
# {
#   expect_equal(grmEigen(test.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
#                          remove.monosnp = FALSE, maf = NaN, missing.rate = NaN, transpose = TRUE), 
#                my.res)
#   expect_equal(grmEigen(test.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
#                          remove.monosnp = FALSE, maf = NaN, missing.rate = NaN, transpose = TRUE), 
#                snprelate.res)
# })
# snpgdsClose(test.file)
# 
# 
# 
# # Example 2
# # read in the solutions
# snprelate.res <- readRDS("snprelate-example-2.txt")
# my.res <- readRDS("eigen-example-2.txt")
# 
# # open the test file
# test.file <- snpgdsOpen("ex2.gds")
# samp.id <- read.gdsn(index.gdsn(test.file, "sample.id"))
# snp.id <- read.gdsn(index.gdsn(test.file, "snp.id"))
# 
# test_that("grmEigen returns the appropriate GRM matrix",
# {
#   expect_equal(grmEigen(test.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
#                          remove.monosnp = FALSE, maf = NaN, missing.rate = NaN, transpose = FALSE), 
#                my.res)
#   expect_equal(grmEigen(test.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
#                          remove.monosnp = FALSE, maf = NaN, missing.rate = NaN, transpose = FALSE), 
#                snprelate.res)
# })
# snpgdsClose(test.file)



# Example 2 Transposed
# read in the solutions
# snprelate.res <- readRDS("snprelate-example-2t.txt")
# my.res <- readRDS("eigen-example-2t.txt")
# 
# # open the test file
# test.file <- snpgdsOpen("ex2t.gds")
# samp.id <- read.gdsn(index.gdsn(test.file, "sample.id"))
# snp.id <- read.gdsn(index.gdsn(test.file, "snp.id"))
# 
# test_that("grmEigen returns the appropriate GRM matrix",
# {
#   expect_equal(grmEigen(test.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
#                          remove.monosnp = FALSE, maf = NaN, missing.rate = NaN, transpose = TRUE), 
#                my.res)
#   expect_equal(grmEigen(test.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
#                          remove.monosnp = FALSE, maf = NaN, missing.rate = NaN, transpose = TRUE), 
#                snprelate.res)
# })
# snpgdsClose(test.file)