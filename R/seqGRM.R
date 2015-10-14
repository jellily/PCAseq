seqGRM <- function(gdsobj, weights, sample.id = NULL, snp.id = NULL,
                   autosome.only = TRUE, remove.monosnp = TRUE, maf = NaN,
                   missing.rate = NaN, verbose = TRUE){

  # Check the inputs for the appropriate classes and values:
  checkWeights(weights)
  
  checkBool(autosome.only)
  checkBool(remove.monosnp)
  checkBool(verbose)

  checkMaf(maf)

  checkMiss(missing.rate)


  # Find the GRM
  grm <- runGRM(gdsobj, weights, sample.id, snp.id, autosome.only, 
                remove.monosnp, maf, missing.rate)

  # check if the GRM only has one entry
  if (dim(grm)[1] == 1 | dim(grm)[2] == 1 | class(grm) != "matrix") {
    warning("GRM has only one entry.")
  }

  # Return the appropriate object
  # eigen.cnt and need.genmat
  seqPCAClass(grm, weights, maf, sample.id, snp.id, eigenRes = NULL,
              eigenCnt = 32, needGenmat = TRUE)

}
