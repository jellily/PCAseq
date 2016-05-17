# filterSnps -------------------------------------------------------------------
# Subset the genotype data set by removing SNPs based on parameters passed in
filterSnps <- function(snps, snpDat, autosomeOnly, removeMonosnp, missingRate,
                       maf, snpChromosome){  
  print(removeMonosnp)
  # remove sex chromosome snps
  if (autosomeOnly) {
    temp <- filterAuto(snpChromosome)
    snps <- snps %in% temp
  }
  
  # remove monomorphic snps
  if (removeMonosnp) {
    print("cake!")
    temp <- filterMono(snpDat)
    snps <- snps %in% temp
  }
  
  # remove snps with too much missingness
  if (!is.nan(missingRate)) {
    temp <- filterMiss(snpDat, missingRate)
    snps <- snps %in% temp
  }
  
  # filter based on MAF
  if (!is.na(maf)) {
    temp <- filterMaf(snpDat, maf)
    snps <- snps %in% temp
  }
  
  return(snps)
}


# filterMono -------------------------------------------------------------------
# Remove monomorphic snps
filterMono <-function(snpDat){

  # find the allele frequencies
  alleleFreq <- 0.5*rowMeans(snpDat, na.rm = TRUE)

  # remove monomorphic SNPs
  snps <- which(alleleFreq > 0 & alleleFreq < 1)

  return(snps)
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

  return(snps)
}


# filterMiss -------------------------------------------------------------------
# Subset the SNPs to those with missingness rates less than missingRate
filterMiss <-function(snpDat, missingRate){

  # replace the missing code 3 with NA
  snpDat[snpDat == 3] <- NA

  # find the proportion missing for each SNP
  missing <- getMissRate(snpDat)
  snps <- which(missing <= missingRate)

  return(snps)
}

# filterMaf --------------------------------------------------------------------
# Subset the SNPs to those with MAFs either greater than maf (if a single value)
# or in the range specified (if two values)
filterMaf <- function(snpDat, maf) {

  # find the allele frequencies
  alleleFreq <- 0.5*rowMeans(snpDat, na.rm = TRUE)
  alleleMAFs <- calcMaf(alleleFreq)

  mafMin <- getBound(maf, "min")
  mafMax <- getBound(maf, "max")

  lowerBound <- substr(maf, 1, 1)
  upperBound <- substr(maf, nchar(maf), nchar(maf))

  # subset based on the interval specified
  # four options ( ), ( ], [ ), and [ ]
  if(lowerBound  == "(" & upperBound == ")") {
    snps <- which(alleleMAFs > mafMin & alleleMAFs < mafMax)
  } else if (lowerBound == "(" & upperBound == "]") {
    snps <- which(alleleMAFs > mafMin & alleleMAFs <= mafMax)
  } else if (lowerBound == "[" & upperBound == ")") {
    snps <- which(alleleFreq >= mafMin & alleleFreq < mafMax)
  } else if(lowerBound == "[" & upperBound == "]") {
    snps <- which(alleleFreq >= mafMin & alleleFreq <= mafMax)
  } else {
    stop("There was an error with the MAF boundaries please check that they are 
         properly specified.")
  }

  return(snps)
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

  # function to calculate the proportion missing for one snp
  propMiss <- function(snp){
    mean(ifelse(is.na(snp), 1, 0))
  }

  missing <- apply(snps, byRows, FUN = propMiss)

  return(missing)
}
