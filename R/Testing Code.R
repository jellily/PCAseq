library(SNPRelate)
setwd("~/Documents/School/Research/Sequence Data")

# Test the GRM functions on the GAW 19 Data
snpgdsBED2GDS(bed.fn = "chr1-clean.bed",
              fam.fn = "chr1-clean.fam",
              bim.fn = "chr1-clean.bim",
              out.gdsfn = "chr1-clean.gds")

geno.dat <- openfn.gds("chr1-clean.gds")
  


eigen.test <- grm.eigen(geno.dat)
pcaseq.test <- grm.pcaseq(geno.dat)

# Need to test the data checks by passing in things
# that do not meet the requirements


# test the popStruct function
eigen.test <- popStruct(geno.dat, method = "eigen")
pcaseq.test <- popStruct(geno.dat, method = "pcaseq")