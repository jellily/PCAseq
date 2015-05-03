# seqGRM

seqGRM <- function(gdsobj, sample.id = NULL, snp.id = NULL, autosome.only = TRUE,
                   remove.monosnp = TRUE, maf = NaN, missing.rate = NaN,
                   verbose = TRUE)
{
  # Check the inputs for the appropriate classes and values:
  check.gdsobj(gdsobj)
  
  check.bool(autosome.only)
  check.bool(remove.monosnp)
  check.bool(verbose)
  
  check.maf(maf)
  
  check.miss(missing.rate)
  
  
  # Find the GRM
  run.grm(gdsobj, sample.id, snp.id, autosome.only, remove.monosnp, maf, missing.rate)
  
  # Return the appropriate object
  make.snpgdsPCAClass(grm, eigen.res, sample.id, snp.id)
}
