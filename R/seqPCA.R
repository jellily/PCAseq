# seqPCA

seqPCA <- function(gdsobj, method, sample.id = NULL, snp.id = NULL, 
                   autosome.only = TRUE, remove.monosnp = TRUE, maf = NaN, 
                   missing.rate = NaN, eigen.cnt = 32, need.genmat = FALSE, 
                   verbose = TRUE){
  
  # Check the inputs for the appropriate classes and values
  checkBool(autosome.only)
  checkBool(remove.monosnp)
  checkBool(need.genmat)
  checkBool(verbose)
  
  checkEcnt(eigen.cnt)
  
  checkMaf(maf)
  
  checkMiss(missing.rate)
  
  
  # Find the GRM
  grm <- runGRM(gdsobj, method, sample.id, snp.id, autosome.only, 
                remove.monosnp, maf, missing.rate)
  
  # Check if the GRM only has one entry
  if (dim(grm)[1] == 1 | dim(grm)[2] == 1 | class(grm) != "matrix")
  {
    warning("GRM has only one entry.")
  }
  
  # Find the eigendecomposition
  eigenRes <- eigen(grm)
  
  # Return the appropriate object
  seqPCAClass(grm, method, maf, eigenRes, sample.id, snp.id, eigen.cnt, 
              need.genmat)
}
