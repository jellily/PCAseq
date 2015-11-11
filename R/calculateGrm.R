# Functions to calculate the GRMs

# runGRM -----------------------------------------------------------------------
# Open the GDS file and check the orientation of the genotype data; call the
# appropriate GRM method
runGRM <- function(gdsobj, weights, sampleId, snpId, autosomeOnly,
                   removeMonosnp, maf, missingRate) {
  # Open the file
  genoDat <- snpgdsOpen(gdsobj)

  # If no sample ids are given, take all of them
  # otherwise check that the sample.id vector
  # has no duplicates and contains sample ids in the
  # file sample.id vector
  samples <- read.gdsn(index.gdsn(genoDat, "sample.id"))

  if (is.null(sampleId)) {
    sampleId <- samples
  } else {
    checkSamp(sampleId, samples)
  }

  # if no SNP IDs are given, take all of them
  # (note SNP filtering happens at a later time point)
  snps <- read.gdsn(index.gdsn(genoDat, "snp.id"))

  if (is.null(snpId)) {
    snpId <- snps
  } else {
    checkSnp(snpId, snps)
  }

  # Determine the orientation of the file -- the code to calculate
  # the GRM is written for SNP x SAMPLE data, if the data is not
  # in that order, it is transposed as it is read in
  transpose <- identical(names(get.attr.gdsn(index.gdsn(genoDat, "genotype"))),
                 "sample.order")

  # call the appropraite function based on the method
  grm <-  grmCalc(genoDat, weights, sampleId, snpId, autosomeOnly,
                                 removeMonosnp, maf, missingRate, transpose)

  # Close the file
  snpgdsClose(genoDat)

  return(list(grm[[1]], samples, snpId[[2]]))
}


# grmCalc ---------------------------------------------------------------------
# Calculate the GRM
grmCalc <- function(genoDat, weights, sampleId, snpId, autosomeOnly,
                    removeMonosnp, maf, missingRate, transpose){

  # constants
  nBlocks <- 5000
  byRows <- 1
  nCopies <- 2

  nSubj <- length(sampleId)
  nSnps <- length(snpId)

  alpha <- weights[1]
  beta <- weights[2]

  emptyMat <- matrix(0, nrow = nSubj, ncol = nSubj)
  max <- ceiling(nSnps / nBlocks)  # maximum number of blocks to loop over

  # create empty grm & vector to count the number of snps used
  grm <- matrix(0, nrow = nSubj, ncol = nSubj)

  # Loop through the SNPs in blocks of size nblock
  for(i in 1:max) {

    message(paste("Computing GRM: Block", i, "of", max))

    # Get the SNPs to choose
    snps <- getIndex(i, nBlocks, nSnps)

    # Read in the relevant SNP data, subsetting by subject
    snpDat <- snpgdsGetGeno(genoDat, snp.id = snpId[snps], sample.id = sampleId)
    snpChrom <- read.gdsn(index.gdsn(genoDat,"snp.chromosome"))[snps]

    # Transpose data into SNPs x Samples
    if ( isTRUE(transpose) ) {
      snpDat <- t(snpDat)
    }

    # Filter the data
    snpInfo <- filterSnps(snps, snpDat, autosomeOnly, removeMonosnp,
                          missingRate, maf, snpChrom)
    snps <- snpInfo[[1]]
    snpDat <- snpInfo[[2]]

    # check to make sure there are still SNPs in the data set
    if ( !(identical(class(snpDat), "matrix")) | (dim(snpDat)[1] == 0) ) {
      message("No data remains in this block after filtering. Going to next
                block.")
     next
    } else {
      alleleFreq <- (1 / nCopies) * rowMeans(snpDat, na.rm = TRUE)
      
      # Estimate the variance at each SNP
      genoCent <- sweep(snpDat, byRows, STATS = nCopies * alleleFreq)
      
      weights <- betaWeights(alleleFreq, alpha, beta)
      
      # Find the empirical correlation matrix
      zee <- sweep(genoCent, byRows, STATS = weights, FUN = "*")
      grm <- grm + crossprod(zee)
    }
  }

  if (identical(grm, emptyMat)) {
    stop("GRM is the zero matrix. Perhaps all of the SNPs were removed when
         filtering or there is no variability in the genotype data.")
  } else {
    return(list(grm, snps))
  }
}


# getIndex ---------------------------------------------------------------------
# function to calculate the index for the current block of SNPs
getIndex <- function(i, nBlock, nSnps){

  index <- (1:nBlock) + (i - 1) * nBlock

  if (index[nBlock] > nSnps) {
    index <- index[1]:nSnps
  }

  return(index)
}


# betaWeights ------------------------------------------------------------------
# function to find the weights for each SNP by MAF
betaWeights <- function(alleleFreq, alpha, beta) {

  # find the minor allele frequency
  minorFreq <- calcMaf(alleleFreq)

  return(dbeta(minorFreq, alpha, beta))
}