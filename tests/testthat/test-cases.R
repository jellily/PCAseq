# This code generates the test gds files and 
# manually calculates the solutions for testing
# the GRM calculations
library(SNPRelate)
library(matrixcalc)

# Example 1
test <- t(matrix(c(0,0,1, 0,1,0, 1,1,1), ncol = 3, byrow = TRUE))

# create a test gds file
snpgdsCreateGeno("ex1.gds", genmat = test,
                 sample.id = 1:3, snp.id = 1:3,
                 snp.chromosome = rep(1,3),
                 snp.position = 1:3,
                 snp.allele = rep(0, 3) )

# EIGENSTRAT GRM Calculation by hand
solution <- sweep(test, MARGIN = 1,  STATS = rowMeans(test))
sigma <- sqrt(0.5*rowMeans(test)*(1-0.5*rowMeans(test)))
solution <- sweep(solution, MARGIN = 1, STATS = sigma, FUN = "/")
solution <- crossprod(solution)/dim(solution)[1]
saveRDS(solution, "eigen-example-1.txt")

# SNPRelate reports a GRM that is multiplied by (n-1)/trace(GRM)
scale <- (ncol(solution)-1)/matrix.trace(solution)
test.file <- snpgdsOpen("ex1.gds")
snprelate.res <- snpgdsPCA(test.file, genmat.only = TRUE)
saveRDS(snprelate.res$genmat/scale, "snprelate-example-1.txt")
snpgdsClose(test.file)

# PCA-seq by Hand
solution <- sweep(test, MARGIN = 1,  STATS = rowMeans(test))
solution <- crossprod(solution)
sigma <- sum(rowMeans(test)*(1-0.5*rowMeans(test)))
solution <- solution/sigma
saveRDS(solution, "pcaseq-example-1.txt")


# Example 1 Transposed
test <- matrix(c(0,0,1, 0,1,0, 1,1,1), ncol = 3, byrow = TRUE)

# create a test gds file
snpgdsCreateGeno("ex1t.gds", genmat = test,
                 sample.id = 1:3, snp.id = 1:3,
                 snp.chromosome = rep(1,3),
                 snp.position = 1:3,
                 snp.allele = rep(0, 3), snpfirstdim = FALSE)


# EIGENSTRAT GRM Calculation by hand
solution <- sweep(test, MARGIN = 2,  STATS = colMeans(test))
sigma <- sqrt(0.5*colMeans(test)*(1-0.5*colMeans(test)))
solution <- sweep(solution, MARGIN = 2, STATS = sigma, FUN = "/")
solution <- tcrossprod(solution)/dim(solution)[2]
saveRDS(solution, "eigen-example-1t.txt")

# SNPRelate reports a GRM that is multiplied by (n-1)/trace(GRM)
scale <- (ncol(solution)-1)/matrix.trace(solution)
test.file <- snpgdsOpen("ex1t.gds")
snprelate.res <- snpgdsPCA(test.file, genmat.only = TRUE)
saveRDS(snprelate.res$genmat/scale, "snprelate-example-1t.txt")

# PCA-seq by Hand
solution <- sweep(test, MARGIN = 2,  STATS = colMeans(test))
solution <- tcrossprod(solution)
sigma <- sum(colMeans(test)*(1-0.5*colMeans(test)))
solution <- solution/sigma
saveRDS(solution, "pcaseq-example-1t.txt")



# Example 2
data(hapmap_geno)
test <- hapmap_geno$genotype[1:5, 1:3]
snpgdsCreateGeno("ex2.gds", genmat = test,
                 sample.id = hapmap_geno$sample.id[1:3], snp.id = hapmap_geno$snp.id[1:5],
                 snp.chromosome = hapmap_geno$snp.chromosome[1:5],
                 snp.position = hapmap_geno$snp.position[1:5],
                 snp.allele = hapmap_geno$snp.allele[1:5])

# EIGENSTRAT GRM Calculation by hand
solution <- sweep(test, MARGIN = 1,  STATS = rowMeans(test))
sigma <- sqrt(0.5*rowMeans(test)*(1-0.5*rowMeans(test)))
solution <- sweep(solution, MARGIN = 1, STATS = sigma, FUN = "/")
solution <- crossprod(solution)/dim(solution)[1]
saveRDS(solution, "eigen-example-2.txt")

# SNPRelate reports a GRM that is multiplied by (n-1)/trace(GRM)
scale <- (ncol(solution)-1)/matrix.trace(solution)
test.file <- snpgdsOpen("ex2.gds")
snprelate.res <- snpgdsPCA(test.file, genmat.only = TRUE)
saveRDS(snprelate.res$genmat/scale, "snprelate-example-2.txt")
snpgdsClose(test.file)

# PCA-seq by Hand
solution <- sweep(test, MARGIN = 1,  STATS = rowMeans(test))
solution <- crossprod(solution)
sigma <- sum(rowMeans(test)*(1-0.5*rowMeans(test)))
solution <- solution/sigma
saveRDS(solution, "pcaseq-example-2.txt")


# Example 2 Transposed
test <- t(hapmap_geno$genotype[1:5, 1:3])
snpgdsCreateGeno("ex2t.gds", genmat = test,
                 sample.id = hapmap_geno$sample.id[1:3], snp.id = hapmap_geno$snp.id[1:5],
                 snp.chromosome = hapmap_geno$snp.chromosome[1:5],
                 snp.position = hapmap_geno$snp.position[1:5],
                 snp.allele = hapmap_geno$snp.allele[1:5], snpfirstdim = FALSE)

# EIGENSTRAT GRM Calculation by hand
solution <- sweep(test, MARGIN = 2,  STATS = colMeans(test))
sigma <- sqrt(0.5*colMeans(test)*(1-0.5*colMeans(test)))
solution <- sweep(solution, MARGIN = 2, STATS = sigma, FUN = "/")
solution <- tcrossprod(solution)/dim(solution)[1]
saveRDS(solution, "eigen-example-2t.txt")

# SNPRelate reports a GRM that is multiplied by (n-1)/trace(GRM)
scale <- (ncol(solution)-1)/matrix.trace(solution)
test.file <- snpgdsOpen("ex2t.gds")
snprelate.res <- snpgdsPCA(test.file, genmat.only = TRUE)
saveRDS(snprelate.res$genmat/scale, "snprelate-example-2t.txt")
snpgdsClose(test.file)

# PCA-seq by Hand
solution <- sweep(test, MARGIN = 2,  STATS = colMeans(test))
solution <- tcrossprod(solution)
sigma <- sum(colMeans(test)*(1-0.5*colMeans(test)))
solution <- solution/sigma
saveRDS(solution, "pcaseq-example-2t.txt")
