#### grm.calc  ####
#
# Function to calculate the empirical genotype relation matrix
# 
# Arguments:
# geno -- matrix of genotype data; subjects in columns
#
# Returns: matrix

# the pca-seq calculations are off here
# need to go back and fix those

grm.calc <- function(genodat, freq.min, freq.max)
{
  # Get SNP & Sample ids
  snpid <- read.gdsn(index.gdsn(genodat, "snp.id"))
  sampid <- read.gdsn(index.gdsn(genodat, "sample.id"))
  
  # Number of SNPs & Samples
  nloci <- length(snpid)
  nsamp <- length(sampid)
  
  # Loop over the SNPs, reading in 5,000 at a time
  grm.pcaseq <- matrix(0, nrow = nsamp, ncol = nsamp)
  total.snps <- 0
  sigma.hat <- 0
  
  max <- floor(length(snpid)/5000)
  
  for(i in 0:max)
  {
    # Get the SNP indexes
    snps <- (1:5000) + i*5000
    if(snps[5000] > nloci)
    {
      snps <- snps[1]:nloci
    }
    
    # Read in the relevant SNP data
    snpdat <- snpgdsGetGeno(genodat, snp.id = snpid[snps])
    
    # Check for missing genotypes
    if(any(is.na(snpdat)))
    {
      stop("Missing genotype data. See help(popStruct) for details.")
    }
    
    # check for genotype data in the wrong format
    # includes letters, non-integers, etc.
    if( any( !(snpdat %in% 0:2) ) )
    {
      stop("Genotype data is not 0, 1, or 2. See help(popStruct) for details.")
    }
    
    # check that there are at least 2 subjects in the data set
    if( dim(snpdat)[1] < 2)
    {
      stop("Dataset includes only 1 subject.")
    }  
    
    # Get the allele frequencies
    allele.freq <- 0.5*colMeans(snpdat)
    
    # Remove monomorphic & singleton SNPs
    single.freq <- 2/nsamp
    snpdat <- snpdat[ , allele.freq > single.freq 
                     & allele.freq < (1 - single.freq)]
    
    
    
    # Check that there are still SNPs in the data set
    # if there are none, go to the next block of SNPs
    if( dim(snpdat)[2] == 0 )
    {
      next
    }
    
    # Recalculate the allele frequencies
    allele.freq <- 0.5*colMeans(snpdat)
    
    # Calculate the number of SNPs actually used
    total.snps <- total.snps + ncol(snpdat)
    
    # Estimated variance at each SNP
    if(method == "eigen")
    {
      geno.cent <- (snpdat - 2*allele.freq)/(4*allele.freq*(1-allele.freq))
    }else if(method == "pcaseq")
    {
      geno.cent <- (snpdat - 2*allele.freq)
      sigma.hat <- sigma.hat + sum(allele.freq*(1-allele.freq))
    }else
    {
      stop("Methods other than EIGENSTRAT (eigen) and PCA-seq (pcaseq) are not supported.")
    }
    

    grm <- grm + tcrossprod(geno.cent)
  }
  
  if( method == "eigen")
  {
    grm <- grm/nloci
  }else if(method = "pcaseq")
  {
    grm <- grm/4*sigma.hat
  }else
  {
    stop("Methods other than EIGENSTRAT (eigen) and PCA-seq (pcaseq) are not supported.")
  }
  
  return(grm) 
}

