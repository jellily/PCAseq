# filterSnps -------------------------------------------------------------------
# Subset the genotype data set by removing SNPs based on parameters passed in
filterSnps <- function(snpDat, autosomeOnly, removeMonosnp, missingRate,
                       maf, snpChromosome){ 
  
  # set up empty vectors
  aoCheck <- rep(TRUE, nrow(snpDat))
  msCheck <- rep(TRUE, nrow(snpDat))
  mrCheck <- rep(TRUE, nrow(snpDat))
  maCheck <- rep(TRUE, nrow(snpDat))
  
  # remove sex chromosome snps
  if (autosomeOnly) {
    aoCheck <- 1:length(snpChromosome) %in% filterAuto(snpChromosome) 
  }
  
  # remove monomorphic snps
  if (removeMonosnp) {
    msCheck <- 1:nrow(snpDat) %in% filterMono(snpDat)
  }
  
  # remove snps with too much missingness
  if (!is.nan(missingRate)) {
    mrCheck <- 1:nrow(snpDat) %in% filterMiss(snpDat, missingRate) 
  }
  
  # filter based on MAF
  if (!is.na(maf)) {
    maCheck <- 1:nrow(snpDat) %in% filterMaf(snpDat, maf) 
  }
  
  return(aoCheck & msCheck & mrCheck & maCheck)
}


# filterMono -------------------------------------------------------------------
# Remove monomorphic snps
filterMono <-function(snpDat){

  # find the allele frequencies
  alleleFreq <- 0.5*rowMeans(snpDat, na.rm = TRUE)
  
  # remove monomorphic SNPs
  snps <- which(alleleFreq > 0 & alleleFreq < 1)

  return(snps) # returns the indices for the current block index
}


# filterAuto -------------------------------------------------------------------
# Subset SNPs to only those on the autosomal chromosomes
filterAuto <-function(snpChromosome){

  autosomeCodes <- as.character(1:22)

  # Convert the vector of SNP choromosome labels to numeric
  # Any character label will become NA
  snpChromosome <- as.character(snpChromosome)

  # Select only those SNPs with chromosome labels 1-22 (autosomes)
  snps <- which(snpChromosome %in% autosomeCodes)

  return(snps) # returns the indices for the current block indexes
}


# filterMiss -------------------------------------------------------------------
# Subset the SNPs to those with missingness rates less than missingRate
filterMiss <-function(snpDat, missingRate){

  # replace the missing code 3 with NA
  snpDat[snpDat == 3] <- NA

  # find the proportion missing for each SNP
  missing <- getMissRate(snpDat)
  snps <- which(missing <= missingRate)

  return(snps) # returns the indices for the current block index
}

# filterMaf --------------------------------------------------------------------
# Subset the SNPs to those with MAFs either greater than maf (if a single value)
# or in the range specified (if two values)
filterMaf <- function(snpDat, maf) {

  # find the loci MAFs
  MAFs <- calcMaf(0.5*rowMeans(snpDat, na.rm = TRUE))

  mafMin <- getBound(maf, "min")
  mafMax <- getBound(maf, "max")

  lowerBound <- substr(maf, 1, 1)
  upperBound <- substr(maf, nchar(maf), nchar(maf))

  # subset based on the interval specified
  # four options ( ), ( ], [ ), and [ ]
  if(lowerBound  == "(" & upperBound == ")") {
    snps <- which(MAFs > mafMin & MAFs < mafMax)
  } else if (lowerBound == "(" & upperBound == "]") {
    snps <- which(MAFs > mafMin & MAFs <= mafMax)
  } else if (lowerBound == "[" & upperBound == ")") {
    snps <- which(MAFs >= mafMin & MAFs < mafMax)
  } else if(lowerBound == "[" & upperBound == "]") {
    snps <- which(MAFs >= mafMin & MAFs <= mafMax)
  } else {
    stop("There was an error with the MAF boundaries please check that they are 
         properly specified.")
  }

  return(snps) # returns the indices for the current block index
}

# calcMAF ----------------------------------------------------------------------
# calculate the minor allele frequency
calcMaf <- function(alleleFreq) {
  return(pmin(alleleFreq, 1-alleleFreq))
}

# getMissRate ------------------------------------------------------------------
# Get the missingness rate for the SNPs
getMissRate <- function(snps){
  byRows <- 1

  # calculate the proportion missing for one snp
  return(apply(snps, byRows, FUN = mean(ifelse(is.na(snp), 1, 0))))
}
