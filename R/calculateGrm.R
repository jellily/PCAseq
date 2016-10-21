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

  return(list(grm[[1]], sampleId, grm[[2]]))
}



# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                         decreasing=FALSE, head=FALSE, n=5) {
  napply <- function(names, fn) sapply(names, function(x)
    fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
    capture.output(format(utils::object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
    as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
  if (!missing(order.by))
    out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
    out <- head(out, n)
  out
}

# shorthand
lsos <- function(..., n=10) {
  .ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n)
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
  
  alpha <- weights[1]
  beta <- weights[2]
  
  # get the index of the subjects
  subj <- read.gdsn(index.gdsn(genoDat, "sample.id"))
  subj <- which(subj %in% sampleId)
  snps <- read.gdsn(index.gdsn(genoDat, "snp.id")) # the letter SNP codes
  
  keepSnps <- c()
  nSnps <- length(snps)
  #maxBlocks <- ceiling(nSnps / nBlocks)  # maximum number of blocks to loop over
  maxBlocks <- 20
  # create empty grm & vector to count the number of snps used
  grm <- matrix(0, nrow = nSubj, ncol = nSubj)
  
  print(mem_used())
  # Loop through the SNPs in blocks of size nblock
  for(i in 1:maxBlocks) {
    
    message(paste("Computing GRM: Block", i, "of", maxBlocks))
    
    # Get the SNPs to choose
    first <- 1 + (i - 1) * nBlocks
    count <- nBlocks
    
    if((nBlocks + first - 1) > nSnps){
      count <- nSnps - first + 1  
    }
    
    
    # Transpose data into SNPs x Samples
    if ( isTRUE(transpose) ) {
      snpDat <- read.gdsn(index.gdsn(genoDat, "genotype"), start = c(1, first), count = c(-1, count)) 
      snpDat <- t(snpDat)
    } else {
      snpDat <- read.gdsn(index.gdsn(genoDat, "genotype"), start = c(first,1), count = c(count,-1)) 
    }

    # Read in the chromosome ids: might be faster to move this into the autosome only function,
    # so that it's only done when it needs to be done...
    snpChrom <- read.gdsn(index.gdsn(genoDat, "snp.chromosome"), start = first, count = count)

    # Filter the data
    # by subject
    snpDat <- snpDat[ , subj] # subset by subject ID

    # by SNP
    snpIndex <- filterSnps(snpDat, autosomeOnly, removeMonosnp,
                           missingRate, maf, snpChrom)
    snpIndex <- snpIndex & snps[first:(first + count - 1)] %in% snpId
    snpDat <- snpDat[snpIndex, ] # subset by SNP ID
    keepSnps <- c(keepSnps, snps[which(snpIndex) + (i-1) * nBlocks])

    # check to make sure there are still SNPs in the data set
    if ( !(identical(class(snpDat), "matrix")) | (dim(snpDat)[1] == 0) ) {
      message("No data remains in this block after filtering. Going to next
              block.")
      next
    } else {

      alleleFreq <- (1 / nCopies) * rowMeans(snpDat, na.rm = TRUE)
      weights <- betaWeights(alleleFreq, alpha, beta)

      # Estimate the variance at each SNP
      genoCent <- sweep(snpDat, byRows, STATS = nCopies * alleleFreq, check.margin = FALSE)
      
      # Find the empirical correlation matrix
      zee <- sweep(genoCent, byRows, STATS = weights, FUN = "*", check.margin = FALSE)

      grm <- grm + crossprod(zee)
      print(mem_used())
      print(lsos())
    }
  }
  
  print(mem_used())
  
  if (identical(grm, matrix(0, nrow = nSubj, ncol = nSubj))) {
    stop("GRM is the zero matrix. Perhaps all of the SNPs were removed when
         filtering or there is no variability in the genotype data.")
  } else {
    return(list(grm, keepSnps))
  }
}

# betaWeights ------------------------------------------------------------------
# function to find the weights for each SNP by MAF
betaWeights <- function(alleleFreq, alpha, beta) {

  # find the minor allele frequency
  return(dbeta(calcMaf(alleleFreq), alpha, beta))
}