# seqPCA

seqPCA <- function(gdsobj, method, sample.id = NULL, snp.id = NULL, autosome.only = TRUE,
          remove.monosnp = TRUE, maf = NaN, missing.rate = NaN, eigen.cnt = 32,
          need.genmat = FALSE, verbose = TRUE)
{
  # Check the inputs for the appropriate classes and values
  check.bool(autosome.only)
  check.bool(remove.monosnp)
  check.bool(need.genmat)
  check.bool(verbose)
  
  check.ecnt(eigen.cnt)
  
  check.maf(maf)
  
  check.miss(missing.rate)
  
  
  # Find the GRM
  grm <- run.grm(gdsobj, sample.id, snp.id, autosome.only, remove.monosnp, maf, missing.rate, method)
  
  # Check if the GRM only has one entry
  if (dim(grm)[1] == 1 | dim(grm)[2] == 1 | class(grm) != "matrix")
  {
    warning("GRM has only one entry.")
  }
  
  # Find the eigendecomposition
  eigen.res <- eigen(grm)
  
  # Return the appropriate object
  make.snpgdsPCAClass(grm, eigen.res, sample.id, snp.id, eigen.cnt, need.genmat)
}
