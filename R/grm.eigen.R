#### grm.eigen ####
#
# Function to calculate the empirical genetic relatedness matrix
# using the usual estimator of relateness, an average of the covariance
#
# Arguments:
# geno -- matrix of genotype data; subjects in columns
#
# Returns: matrix of kinship estimates 


grm.eigen <- function(geno)
{
  ns <- dim(geno)[2]
  
  stopifnot( ns > 0 )
  
  # allele freq est for each SNP
  pA <- 0.5*rowMeans(geno, na.rm = TRUE)
  
  # remove monomorphic SNPs
  geno <- subset(geno, pA > 0)
  pA <- 0.5*rowMeans(geno, na.rm = TRUE)
  
  # number of SNPs
  nsnps <- dim(geno)[1]
  
  # estimated variance at each SNP
  sigma.hat <- sqrt(2*pA*(1-pA))
  Zu <- (geno-2*pA)/sigma.hat
  #Zu[which(is.na(Zu))] <- 0
  
  # empirical correlation matrix
  Psiu <- (1/nsnps)*crossprod(Zu)
  
  return(Psiu) 
}