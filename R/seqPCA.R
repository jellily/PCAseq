#' Run PCA-seq on a GDS file
#'
#' This function calculates the Genetic Relatedness Matrix (GRM) on a GDS file
#' using the PCA-seq method for sequence data.
#'
#' @aliases seqGRM
#'
#' @param gdsobj an object of the class SNPGDSFileClass, a SNP GDS file.
#' @param method a string, either "eigen" or "pcaseq", indicating which method
#' to use to calculate the GRM; see Details.
#' @param sample.id a vector of sample ids specifying the samples to use for
#' analysis; if NULL, all samples are used.
#' @param snp.id a vector of SNP ids specifying the SNPs to use for analysis;
#' if NULL, all SNPs are used.
#' @param autosome.only if TRUE, use autosomal SNPs only; if it is a numeric or
#' character vector, keep SNPs according to the specified chromosomes.
#' @param remove.monosnp if TURE, remove monomorphic SNPs.
#' @param maf if one number is specified, use SNPs with MAF greater than or
#' equal to this value; if a numeric vector of length two is specified, only
#' SNPs with MAFs in (min, max) are taken.
#' @param missing.rate to use the SNPs with missing rates less than or equal to
#' missing.rate; if NaN, no misisng threshold.
#' @param eigen.cnt the number of eigen vectors and values to return; if zero,
#' return all eigenvalues and vectors.
#' @param need.genmat if TRUE, return the genetic relatedness matrix.
#' @param verbose Not supported.
#'
#' @details If \code{method} is "eigen", the GRM is calculated using the
#' EIGENSTRAT method as given in Patterson et al 2006. If \code{method} is
#' "pcaseq", the GRM is calculated using the PCA-seq method.
#'
#' @section Note:
#' If you need to run the EIGENSTRAT method on a very large data set and do not
#' need to subset by both a minimum and maximum MAF, the
#' \code{\link[SNPRelate]{snpgdsPCA}} function will be faster.
#'
#'
#' @return Return a \code{snpPCAClass} object, a list with the follow slots:
#' \describe{
#' \item{\code{method}}{the method used to calculate the GRM}
#' \item{\code{maf}}{the MAF cutoffs used}
#' \item{\code{sample.id}}{the sample ids used in the analysis}
#' \item{\code{snp.id}}{the SNP ids used in the analysis}
#' \item{\code{eigenval}}{eigenvalues}
#' \item{\code{eigenvect}}{a matrix of eigenvectors of dimensions # of samples by
#' \code{eigen.cnt}}
#' \item{\code{varprop}}{the proportion of the variance explained by each principal
#' component}
#' \item{\code{TraceXTX}}{the trace of the genetic relateness matrix}
#' \item{\code{Bayesian}}{indicates Bayes normalization; set to FALSE, as this is not
#' currently supported}
#' \item{\code{genmat}}{the genetic relateness matrix}}

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
  grmRes <- runGRM(gdsobj, method, sample.id, snp.id, autosome.only,
                remove.monosnp, maf, missing.rate)

  grm <- grmRes[[1]]
  sampleId <- grmRes[[2]]
  snpId <- grmRes[[3]] # fix this to actully return the list of SNPs use
  
  # Check if the GRM only has one entry
  if (dim(grm)[1] == 1 | dim(grm)[2] == 1 | class(grm) != "matrix"){
    warning("GRM has only one entry.")
  }

  # Find the eigendecomposition
  eigenRes <- eigen(grm)
  
  # Return the appropriate object
  seqPCAClass(grm, method, maf, eigenRes, sampleId, snpId, eigen.cnt,
              need.genmat)
}
