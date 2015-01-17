#### grm.pcaseq  ####
#
# Function to calculate the empirical genotype relation matrix
# 
# Arguments:
# geno -- matrix of genotype data; subjects in columns
#
# Returns: matrix

# Method of Moments Estimator
grm.pcaseq <- function(geno)
{
  ns <- dim(geno)[2]
  
  stopifnot( ns > 0 )
  
  # allele freq est for each SNP
  pA <- 0.5*rowMeans(geno, na.rm = TRUE)
  
  nsnps <- length(pA)
  
  # estimated variance at each SNP
  Zu <- (geno-2*pA)
  
  # empirical correlation matrix
  Psiu <- (1/nsnps)*crossprod(Zu)
  
  # MoM estimator of the variance
  sigma.hat <- mean( 4*pA*(1-pA) )
  
  Psiu<- Psiu/sigma.hat
  
  return(Psiu) 
}

