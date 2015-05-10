# make.snpgdsPCAClass

seqPCAClass <- function(grm, method, eigen.res, sample.id, snp.id, eigen.cnt, need.genmat)
{
  # find the trace of the GRM
  mat.trace <- matrix.trace(grm)
  
  # if need.genmat is FALSE,
  # set genmat to NULL
  if (!need.genmat){
    grm <- NULL
  }
  
  # if eigen.cnt is greater than the total number
  # use the total number & issue a warning
  if (eigen.cnt > length(eigen.res$values)){
    warning("Number of eigenvectors and values to return is more than the dimensions of the GRM. All eigenvalues and vectors will be returned.")
    eigen.cnt <- length(eigen.res$values)
  }
  
  if (is.null(eigen.res)){
    eigenval <- NULL
    eigenvect <- NULL
    varprop <- NULL
  }else{
    eigenval <- eigen.res$values[1:eigen.cnt]
    eigenvect <- eigen.res$vectors[ , 1:eigen.cnt]
    varprop <- eigenval/sum(eigen.res$values)
  }
  # create the object of class snpgdsPCAClass
  pca.res <- list("method" = method,
                  "sample.id" = sample.id,
                  "snp.id" = snp.id,
                  "eigenval" = eigenval,
                  "eigenvect" = eigenvect,
                  "varprop" = varprop,
                  "TraceXTX" = mat.trace,
                  "Bayesian" = FALSE,
                  "genmat" = grm)
  class(pca.res) <- list("seqPCAClass", "snpgdsPCAClass")
  
  return(pca.res)  
}