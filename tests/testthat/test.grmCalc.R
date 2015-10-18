# Test for grmCalc

# Exmple 1
# read in the solutions
#snprelate.res <- readRDS("snprelate-example-1.txt")
res1 <- read.table("eigen-example-1.txt")
res2 <- read.table("pcaseq-example-1.txt")

# open the test file
test.file <- snpgdsOpen("ex1.gds")
samp.id <- read.gdsn(index.gdsn(test.file, "sample.id"))
snp.id <- read.gdsn(index.gdsn(test.file, "snp.id"))

test_that("grmEigen returns the appropriate GRM matrix and that it is equivalent to SNPRelate",
{
  expect_equivalent(grmCalc(test.file, c(0.5, 0.5), sampleId = samp.id,
                            snpId = snp.id, autosomeOnly = FALSE,
                            removeMonosnp = FALSE, maf = NaN,
                            missingRate = NaN, transpose = FALSE)[[1]],
                    as.matrix(my.res1))

  expect_equivalent(grmCalc(test.file, c(1, 1), sampleId = samp.id,
                            snpId = snp.id, autosomeOnly = FALSE,
                            removeMonosnp = FALSE, maf = NaN,
                            missingRate = NaN, transpose = FALSE)[[1]],
                    as.matrix(res2))
})
snpgdsClose(test.file)
