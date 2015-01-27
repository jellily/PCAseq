#### grm.eigen ####
#
# Function to calculate the empirical genetic relatedness matrix
# using the usual estimator of relateness, an average of the covariance
#
# Arguments:
# geno -- matrix of genotype data; subjects in columns
#
# Returns: matrix of kinship estimates 


grm.eigen <- function(genodat)
{
  #   # check for missing genotypes
  #   if(any(is.na(genodat)))
  #   {
  #     stop("Missing genotype data. See help(popStruct) for details.")
  #   }
  #   
  #   # check for genotype data in the wrong format
  #   # includes letters, non-integers, etc.
  #   if( any( !(genodat %in% 0:2) ) )
  #   {
  #     stop("Genotype data is not 0, 1, or 2. See help(popStruct) for details.")
  #   }
  #   
  
  # Get SNP & Sample ids
  snpid <- read.gdsn(index.gdsn(genodat, "snp.id"))
  sampid <- read.gdsn(index.gdsn(genodat, "sample.id"))
  
  # Number of SNPs & Samples
  nloci <- length(snpid)
  nsamp <- length(sampid)
  
  #   # check that there are at least 2 subjects in the data set
  #   if( dim(genodat)[2] < 2)
  #   {
  #     stop("Data set includes only 1 subject.")
  #   }
  
  # Loop over the SNPs, reading in 5,000 at a time
  grm.eigen <- matrix(0, nrow = nsamp, ncol = nsamp)
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
    snp.dat <- snpgdsGetGeno(geno.dat, snp.id = snpid[snps])
    
    # Get the allele frequencies
    allele.freq <- 0.5*colMeans(snp.dat, na.rm = TRUE)
    
    # Remove monomorphic & singleton SNPs
    single.freq <- 2/nsamp
    snp.dat <- snp.dat[ , allele.freq > single.freq 
                       & allele.freq < (1 - single.freq)]
    
    # Recalculate the allele frequencies
    allele.freq <- 0.5*colMeans(snp.dat, na.rm = TRUE)
    
    # Calculate the number of SNPs actually used
    total.snps <- total.snps + ncol(snp.dat)
    
    # Estimated variance at each SNP
    geno.cent <- (snp.dat - 2*allele.freq)
    
    # find the standard deviation
    sigma.hat <- sqrt( 4*allele.freq*(1-allele.freq) )
    
    # For each block matrix of snps,
    # Find the empirical correlation matrix
    grm.eigen <- grm.eigen + tcrossprod(geno.cent/sigma.hat)
  }
  
  grm.eigen <- (1/nloci) * grm.eigen
  
  return(grm.eigen) 
}

