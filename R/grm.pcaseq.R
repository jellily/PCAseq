#### grm.pcaseq  ####
#
# Function to calculate the empirical genotype relation matrix
# 
# Arguments:
# geno -- matrix of genotype data; subjects in columns
#
# Returns: matrix

# Method of Moments Estimator
grm.pcaseq <- function(genodat)
{
  # Get SNP & Sample ids
  snpid <- read.gdsn(index.gdsn(genodat, "snp.id"))
  sampid <- read.gdsn(index.gdsn(genodat, "sample.id"))
  
  # Number of SNPs & Samples
  nloci <- length(snpid)
  nsamp <- length(sampid)
  
  # Loop over the SNPs, reading in 5,000 at a time
  grm.eigen <- matrix(0, nrow = nsamp, ncol = nsamp)
  total.snps <- 0
  sigma.hat <- 0
  
  for(i in 0:max)
  {
    # Get the number of SNPs to choose
    snps <- (1:5000) + i*5000
    if(snps[5000] > length(snpid))
    {
      snps <- snps[1]:length(snpid)
    }
    
    # Remove monomorphic & singleton SNPs
    single.freq <- 2/nsamp
    snp.dat <- snp.dat[ , allele.freq > single.freq 
                       & allele.freq < (1 - single.freq)]
    
    # Get the allele frequencies
    allele.freq[snps] <- 0.5*colMeans(snp.dat, na.rm = TRUE)
    
    # Calculate the number of SNPs actually used
    total.snps <- total.snps + ncol(snp.dat)
    
    # Estimated variance at each SNP
    geno.cent <- snp.dat - 2* allele.freq[snps]

      
    # For each block matrix of snps,
    # find the empirical correlation matrix
    grm.pcaseq <- grm.pcaseq + tcrossprod(geno.cent)
  }
    
  sigma.hat <- sigma.hat + sum(allele.freq*(1-allele.freq))
    
  # close the file
  closefn.gds(geno.dat)
}
  
  # take the ratio of the two estimators
  grm.pcaseq <- grm.pcaseq/(4*sigma.hat)
  
  return(cov) 
}
