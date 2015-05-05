# seqGRM
# this funciton needs work beuase it doesn't interact properly with make.snpgdsPCAClass
# there are several slots from that class which are empty if we are only making the GRM

seqGRM <- function(gdsobj, sample.id = NULL, snp.id = NULL, autosome.only = TRUE,
                   remove.monosnp = TRUE, maf = NaN, missing.rate = NaN,
                   verbose = TRUE, method)
{
  # Check the inputs for the appropriate classes and values:
  check.bool(autosome.only)
  check.bool(remove.monosnp)
  check.bool(verbose)
  
  check.maf(maf)
  
  check.miss(missing.rate)
  
  
  # Find the GRM
  grm <- run.grm(gdsobj, sample.id, snp.id, autosome.only, remove.monosnp, maf, missing.rate, method)
  
  # Return the appropriate object
  # eigen.cnt and need.genmat
  make.snpgdsPCAClass(grm, eigen.res = NULL, sample.id, snp.id, eigen.cnt = 32, need.genmat = TRUE)
}
