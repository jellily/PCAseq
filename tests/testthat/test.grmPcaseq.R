# Test for grmPcaseq

# Exmple 1
# read in the solutions
my.res <- readRDS("pcaseq-example-1.txt")

# open the test file
test.file <- snpgdsOpen("ex1.gds")
samp.id <- read.gdsn(index.gdsn(test.file, "sample.id"))
snp.id <- read.gdsn(index.gdsn(test.file, "snp.id"))

test_that("grmPcaseq returns the appropriate GRM matrix",
{
  expect_equal(grmPcaseq(test.file, sampleId = samp.id, snpId = snp.id, autosomeOnly = FALSE, 
                         removeMonosnp = FALSE, maf = NaN, missingRate = NaN, transpose = FALSE)[[1]], 
              my.res)
})
snpgdsClose(test.file)

# 
# 
# # Example 1 Transposed
# my.res <- readRDS("pcaseq-example-1t.txt")
# 
# # open the test file
# test.file <- snpgdsOpen("ex1t.gds")
# samp.id <- read.gdsn(index.gdsn(test.file, "sample.id"))
# snp.id <- read.gdsn(index.gdsn(test.file, "snp.id"))
# 
# test_that("grmPcaseq returns the appropriate GRM matrix",
# {
#   expect_equal(grmPcaseq(test.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
#                           remove.monosnp = FALSE, maf = NaN, missing.rate = NaN, transpose = TRUE), 
#                my.res)
# })
# snpgdsClose(test.file)
# 
# 
# 
# # Example 2
# # read in the solutions
# my.res <- readRDS("pcaseq-example-2.txt")
# 
# # open the test file
# test.file <- snpgdsOpen("ex2.gds")
# samp.id <- read.gdsn(index.gdsn(test.file, "sample.id"))
# snp.id <- read.gdsn(index.gdsn(test.file, "snp.id"))
# 
# test_that("grm.eigen returns the appropriate GRM matrix",
# {
#   expect_equal(grmPcaseq(test.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
#                          remove.monosnp = FALSE, maf = NaN, missing.rate = NaN, transpose = FALSE), 
#                my.res)
# })
# snpgdsClose(test.file)
# 
# 
# # Example 2 Transposed
# my.res <- readRDS("pcaseq-example-2t.txt")
# 
# # open the test file
# test.file <- snpgdsOpen("ex2t.gds")
# samp.id <- read.gdsn(index.gdsn(test.file, "sample.id"))
# snp.id <- read.gdsn(index.gdsn(test.file, "snp.id"))
# 
# test_that("grmPcaseq returns the appropriate GRM matrix",
# {
#   expect_equal(grmPcaseq(test.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
#                           remove.monosnp = FALSE, maf = NaN, missing.rate = NaN, transpose = TRUE), 
#                my.res)
# })
# snpgdsClose(test.file)

# Example 4
# test.file <- snpgdsOpen("ex4.gds", allow.duplicate = TRUE)
# samp.id <- read.gdsn(index.gdsn(test.file, "sample.id"))
# snp.id <- read.gdsn(index.gdsn(test.file, "snp.id"))
# 
# test_that("grmPcaseq returns an error when the GRM is empty (matrix of zeros)",
# {
#   expect_error(grmPcaseq(test.file, sampleId = samp.id, snpId = snp.id, autosomeOnly = FALSE,
#                removeMonosnp = FALSE, maf = NaN, missingRate = NaN, transpose = FALSE))
# })
# 
# snpgdsClose(test.file)