setwd("~/Documents/School/Research/PCAseq/R")

library(testthat)
library(SNPRelate)
library(matrixcalc)

source("grmPcaseq.R")
source("grmEigen.R")
source("getSnps.R")
source("filterSnps.R")
source("filterMono.R")
source("filterMiss.R")
source("filterMaf.R")
source("filterAuto.R")
source("checkBool.R")
source("checkMaf.R")
source("checkEcnt.R")
source("checkMiss.R")
source("seqPCAClass.R")
source("runGrm.R")
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

#test.file <- snpgdsOpen("test.gds")
#rez <- snpgdsPCA(test.file, genmat.only = TRUE)

my.pca <- seqPCA("test.gds", method = "eigen", need.genmat = TRUE)
my.grm <- seqGRM("test.gds", method = "eigen")

# demonstrate compatibility with SNPRelate
snpgdsPCASNPLoading(my.pca, snpgdsOpen("test.gds"))
snpgdsPCACorr(my.pca, snpgdsOpen("test.gds"))
