setwd("~/Documents/School/Research/PCAseq/R")

library(testthat)
library(SNPRelate)
library(matrixcalc)

source("grm.pcaseq.R")
source("grm.eigen.R")
source("get.snps.R")
source("filter.snps.R")
source("filter.mono.R")
source("filter.miss.R")
source("filter.maf.R")
source("filter.auto.R")
source("check.bool.R")
source("check.maf.R")
source("check.ecnt.R")
source("check.miss.R")
source("make.snpgdsPCAClass.R")
source("run.grm.R")
source("seqGRM.R")
source("seqPCA.R")

# Example 1
test <- t(matrix(c(0,0,1, 0,1,0, 1,1,1), ncol = 3, byrow = TRUE))

# create a test gds file
snpgdsCreateGeno("test.gds", genmat = test,
                 sample.id = 1:3, snp.id = 1:3,
                 snp.chromosome = rep(1,3),
                 snp.position = 1:3,
                 snp.allele = rep(0, 3) )

my.pca <- seqPCA("test.gds", method = "eigen", need.genmat = TRUE)
