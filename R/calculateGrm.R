# Functions to calculate the GRMs

# runGRM -----------------------------------------------------------------------
# Open the GDS file and check the orientation of the genotype data; call the 
# appropriate GRM method

runGRM <- function(gdsobj, method, sampleId, snpId, autosomeOnly, removeMonosnp, 
                    maf, missingRate){
  # Open the file
  genoDat <- snpgdsOpen(gdsobj)
  
  # If no sample ids are given, take all of them
  # otherwise check that the sample.id vector
  # has no duplicates and contains sample ids in the
  # file sample.id vector
  if (is.null(sampleId)){
    sampleId <- read.gdsn(index.gdsn(genoDat, "sample.id"))
  } else {
    checkSamp(sampleId, read.gdsn(index.gdsn(genoDat, "sample.id")))
  }
  
  # if no SNP ids are given, take all of them
  # (note SNP filtering happens at a later time point)
  if (is.null(snpId)){
    snpId <- read.gdsn(index.gdsn(genoDat, "snp.id"))
  } else {
    checkSnp(snpId, read.gdsn(index.gdsn(genoDat, "snp.id")))
  }
  
  # Determine the orientation of the file -- the code to calculate
  # the GRM is written for SNP x SAMPLE data, if the data is not
  # in that order, it is transposed as it is read in
  if (names(get.attr.gdsn(index.gdsn(genoDat, "genotype"))) == "sample.order"){
    transpose <- TRUE
  } else {
    transpose <- FALSE
  }
  
  # call the appropraite function based on the method
  grm <- switch(method, 
                eigen = grmEigen(genoDat, sampleId, snpId, autosomeOnly,
                                 removeMonosnp, maf, missingRate, transpose),
                pcaseq = grmPcaseq(genoDat, sampleId, snpId, autosomeOnly,
                                   removeMonosnp, maf, missingRate, transpose))
  
  # Close the file
  snpgdsClose(genoDat)
  
  return(grm)
}


# grmEigen ---------------------------------------------------------------------
# Calculate the GRM using the method in Price et al 2006

grmEigen <- function(genoDat, sampleId, snpId, autosomeOnly, removeMonosnp, 
                     maf, missingRate, transpose){
  # constants
  nBlocks <- 5000
  byRows <- 1
  nCopies <- 2
  nSubj <- length(sampleId)
  nSnps <- length(snpId)
  max <- floor(nSnps / nBlocks) + 1  # maximum number of blocks to loop over
  
  # create empty grm & vector to count the number of snps used
  grm <- matrix(0, nrow = nSubj, ncol = nSubj)
  totalSnps <- 0
  
  # Loop through the SNPs in blocks of size nblock
  for(i in 1:max)
  {
    message(paste("Computing GRM: Block", i, "of", max))
    
    # Get the SNPs to choose
    snps <- getIndex(i, nBlocks, nSnps)
    
    # Read in the relevant SNP data, subsetting by subject
    snpDat <- snpgdsGetGeno(genoDat, snp.id = snpId[snps], sample.id = sampleId)  
    
    # Transpose data into SNPs x Samples
    if (transpose == TRUE)
    {
      snpDat <- t(snpDat)
    }
    
    # Filter the data
    snpDat <- filterSnps(snpDat, autosomeOnly, removeMonosnp, missingRate, 
                          maf, read.gdsn(index.gdsn(genoDat, 
                                                    "snp.chromosome")))
    
    totalSnps <- dim(snpDat)[byRows] # total # of snps
    alleleFreq <- (1 / nCopies) * rowMeans(snpDat, na.rm = TRUE)
    
    # Estimate the variance at each SNP
    genoCent <- sweep(snpDat, byRows, STATS = nCopies * alleleFreq)
    sigmaEigen <- sqrt(alleleFreq * (1 - alleleFreq)) # standard deviation
    
    # Find the empirical correlation matrix
    zee <- sweep(genoCent, byRows, STATS = sigmaEigen, FUN = "/")
    grm <- grm + crossprod(zee)
  }
  grm <- (1 / totalSnps) * grm
  return(grm)
}

# grmPcaseq --------------------------------------------------------------------

grmPcaseq <- function(genoDat, sampleId, snpId, autosomeOnly, removeMonosnp, 
                      maf, missingRate, transpose)
{
  # constants
  nBlocks <- 5000
  byRows <- 1
  nCopies <- 2
  nSubj <- length(sampleId)
  nSnps <- length(snpId)
  max <- floor(nSnps / nBlocks) + 1  #  maximum number of blocks to loop over
  
  # create empty grm & vector to store the allele frequencies
  grm <- matrix(0, nrow = nSubj, ncol = nSubj)
  allelePs <- c()
  
  # Loop through the SNPs in blocks of size nblock
  for(i in 1:max)
  {
    message(paste("Computing GRM: Block", i, "of", max))
    
    # Get the SNPs to choose
    snps <- getIndex(i, nBlocks, nSnps)
    
    # Read in the relevant SNP data, subsetting by sample if necessary
    snpDat <- snpgdsGetGeno(genoDat, snp.id = snpId[snps], sample.id = sampleId)  
    
    # Transpose data into SNPs x Samples
    if (transpose)
    {
      snpDat <- t(snpDat)
    }
    
    # Filter the data
    snpDat <- filterSnps(snpDat, autosomeOnly, removeMonosnp, missingRate, 
                          maf, read.gdsn(index.gdsn(genoDat, "snp.chromosome")))
    
    # Calculate the SNP allele frequencies
    alleleFreq <- (1 / nCopies) * rowMeans(snpDat, na.rm = TRUE) 
    
    # Estimate the variance at each SNP
    genoCent <- sweep(snpDat, byRows, STATS = 2 * alleleFreq)
    allelePs <- c(allelePs, alleleFreq)
    
    # Find the empirical correlation matrix
    grm <- grm + crossprod(genoCent) 
  }
  
  sigmaPcaseq <- sum(nCopies * allelePs * (1 - allelePs)) # standard deviation
  grm <- grm / sigmaPcaseq
  return(grm)
}

# getIndex ---------------------------------------------------------------------
# function to calculate the index for the current block of SNPs

getIndex <- function(i, nBlock, nSnps)
{
  index <- (1:nBlock) + (i - 1) * nBlock
  
  if (index[nBlock] > nSnps){
    index <- index[1]:nSnps
  }
  
  return(index)
}