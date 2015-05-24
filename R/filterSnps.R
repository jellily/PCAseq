# filterSnps -------------------------------------------------------------------
# Subset the genotype data set by removing SNPs based on parameters passed in

filterSnps <- function(snps, autosomeOnly, removeMonosnp, missingRate, maf,
                       snpChromosome){

  # remove sex chromosome snps
  if (autosomeOnly){
    snps <- filterAuto(snps, snpChromosome)
  }

  # remove monomorphic snps
  if (removeMonosnp){
    snps <- filterMono(snps)
  }

  # remove snps with too much missingness
  if (!is.nan(missingRate)){
    snps <- filterMiss(snps, missing.rate)
  }

  # filter based on MAF
  if (length(maf) != 1){
    snps <- filterMaf(snps, maf)
  }

  return(snps)
}


# filterMono -------------------------------------------------------------------
# Remove monomorphic snps

filterMono <-function(snps){

  # find the allele frequencies
  alleleFreq <- 0.5*rowMeans(snps, na.rm = TRUE)

  # remove monomorphic SNPs
  snps <- snps[alleleFreq > 0 & alleleFreq < 1, ]

  return(snps)
}

# filterAuto -------------------------------------------------------------------
# Subset SNPs to only those on the autosomal chromosomes

filterAuto <-function(snps, snpChromosome){

  autosomeCodes <- as.character(1:22)

  # Convert the vector of SNP choromosome labels to numeric
  # Any character label will become NA
  snpChromosome <- as.character(snpChromosome)

  # Select only those SNPs with chromosome labels 1-22 (autosomes)
  snps <- snps[which(snpChromosome %in% autosomeCodes), ]

  return(snps)
}


# filterMiss -------------------------------------------------------------------
# Subset the SNPs to those with missingness rates less than missingRate

filterMiss <-function(snps, missingRate){

  # replace the missing code 3 with NA
  snps[snps == 3] <- NA

  # find the proportion missing for each SNP
  missing <- getMissRate(snps)
  snps <- snps[missing <= missingRate, ]

  return(snps)
}

# filterMaf --------------------------------------------------------------------
# Subset the SNPs to those with MAFs either greater than maf (if a single value)
# or in the range specified (if two values)

filterMaf <-function(snps, maf){

  # find the allele frequencies
  alleleFreq <- 0.5*rowMeans(snps, na.rm = TRUE)

  # different filtering based on what value is specified
  if (length(maf) == 1){
    snps <- snps[alleleFreq > maf & alleleFreq < (1 - maf), ]
  } else { #  if length(maf) == 2
    min <- maf[1]
    max <- maf[2]

    snps <- snps[alleleFreq >= min & alleleFreq <= max |
                   alleleFreq >= (1 - max) & alleleFreq <= (1 - min), ]
  }

  return(snps)
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
