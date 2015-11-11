# seqPCAClass ------------------------------------------------------------------
seqPCAClass <- function(grm, weights, maf, eigenRes, sampleId, snpId, eigenCnt,
                        needGenmat){

  # find the trace of the GRM
  if (requireNamespace("matrixcalc", quietly = TRUE)) {
    matTrace <- matrixcalc::matrix.trace(grm)
  } else {
    matTrace <- NA
    warning("Package matrixcalc needed to calculate the trace of the GRM.")
  }

  # if needGenmat is FALSE,
  # set genmat to NULL
  if (!needGenmat) {
    grm <- NULL
  }

  # if eigen.cnt is greater than the total number
  # use the total number & issue a warning
  if (eigenCnt > length(eigenRes$values)) {
    warning("Number of eigenvectors and values to return is more than the
            dimensions of the GRM. All eigenvalues and vectors will be
            returned.")
    eigenCnt <- length(eigenRes$values)
  } else if(eigenCnt == 0) {
    eigenCnt <- length(eigenRes$values)
  }

  if (is.null(eigenRes)) {
    eigenval <- NULL
    eigenvect <- NULL
    varprop <- NULL
  } else {
    eigenval <- eigenRes$values[1:eigenCnt]
    eigenvect <- eigenRes$vectors[ , 1:eigenCnt]
    varprop <- eigenval / sum(eigenRes$values)
  }

  # create the object of class snpgdsPCAClass
  pcaRes <- list("weights" = weights,
                 "maf" = maf,
                 "sample.id" = sampleId,
                 "snp.id" = snpId,
                 "eigenval" = eigenval,
                 "eigenvect" = eigenvect,
                 "varprop" = varprop,
                 "TraceXTX" = matTrace,
                 "Bayesian" = FALSE,
                 "genmat" = grm)
  class(pcaRes) <- list("seqPCAClass", "snpgdsPCAClass")

  return(pcaRes)
}