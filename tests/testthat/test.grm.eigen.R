# Test for grm.eigen

# Exmple 1
test <- t(matrix(c(0,0,1, 0,1,0, 1,1,1), ncol = 3, byrow = TRUE))
# create a test gds file
snpgdsCreateGeno("test.gds", genmat = test,
                 sample.id = 1:3, snp.id = 1:3,
                 snp.chromosome = rep(1,3),
                 snp.position = 1:3,
                 snp.allele = rep(0, 3) )
test.file <- snpgdsOpen("test.gds")

samp.id <- read.gdsn(index.gdsn(test.file, "sample.id"))
snp.id <- read.gdsn(index.gdsn(test.file, "snp.id"))

# SNPRelate reports a GRM that is multiplied by (n-1)/trace(GRM)
my.eigen <- grm.eigen(test.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
          remove.monosnp = FALSE, maf = NaN, missing.rate = NaN)
scale <- (ncol(my.eigen)-1)/matrix.trace(my.eigen)

snprelate.res <- snpgdsPCA(test.file, genmat.only = TRUE)

test_that("grm.eigen returns the appropriate GRM matrix",
{
  expect_equal(grm.eigen(test.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
                             remove.monosnp = FALSE, maf = NaN, missing.rate = NaN), 
                   snprelate.res$genmat/scale)
})
unlink("test.gds")


# Example 2
data(hapmap_geno)
snpgdsCreateGeno("test2.gds", genmat = hapmap_geno$genotype[1:5, 1:3],
                 sample.id = hapmap_geno$sample.id[1:3], snp.id = hapmap_geno$snp.id[1:5],
                 snp.chromosome = hapmap_geno$snp.chromosome[1:5],
                 snp.position = hapmap_geno$snp.position[1:5],
                 snp.allele = hapmap_geno$snp.allele[1:5])
test2.file <- snpgdsOpen("test2.gds")

samp.id <- read.gdsn(index.gdsn(test2.file, "sample.id"))
snp.id <- read.gdsn(index.gdsn(test2.file, "snp.id"))

# SNPRelate reports a GRM that is multiplied by (n-1)/trace(GRM)
my.res <- grm.eigen(test2.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
                    remove.monosnp = FALSE, maf = NaN, missing.rate = NaN) 
scale <- (ncol(my.res)-1)/matrix.trace(my.res)

snprelate.res <- snpgdsPCA(test2.file, genmat.only = TRUE)

test_that("grm.eigen returns the appropriate GRM matrix",
{
  expect_equal(grm.eigen(test2.file, sample.id = samp.id, snp.id = snp.id, autosome.only = FALSE, 
                         remove.monosnp = FALSE, maf = NaN, missing.rate = NaN), 
               snprelate.res$genmat/scale)
})
unlink("test2.gds")