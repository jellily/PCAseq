# grm.pcaseq

grmPcaseq <- function(geno.dat, sample.id, snp.id, autosome.only,remove.monosnp, 
                      maf, missing.rate, transpose)
{
  # constants
  nblock <- 5000
  byrows <- 1
  ncopies <- 2
  nsubj <- length(sample.id)
  nsnps <- length(snp.id)
  max <- floor(nsnps/nblock) + 1  # maximum number of blocks to loop over
  
  # create empty grm & vector to count the number of snps used
  grm <- matrix(0, nrow = nsubj, ncol = nsubj)
  allele.ps <- c()
  
  # Loop through the SNPs in blocks of size nblock
  for(i in 1:max)
  {
    message(paste("Computing GRM: Block", i, "of", max))
    
    # Get the SNPs to choose
    snps <- getSnps(i, nblock, nsnps)
    
    # Read in the relevant SNP data, subsetting by sample if necessary
    snp.dat <- snpgdsGetGeno(geno.dat, snp.id = snp.id[snps], sample.id = sample.id)  

    # Transpose data into SNPs x Samples
    if (transpose)
    {
      snp.dat <- t(snp.dat)
    }
    
    # Filter the data
    snp.dat <- filterSnps(snp.dat, autosome.only, remove.monosnp, missing.rate, 
                           maf, read.gdsn(index.gdsn(geno.dat, "snp.chromosome")))
    
    # Calculate the SNP allele frequencies
    allele.freq <- (1/ncopies)*rowMeans(snp.dat, na.rm = TRUE)  # snp allele frequencies
    
    # Estimate the variance at each SNP
    geno.cent <- sweep(snp.dat, byrows, STATS = 2*allele.freq)
    allele.ps <- c(allele.ps, allele.freq)
    
    # Find the empirical correlation matrix
    grm <- grm + crossprod(geno.cent) 
  }
  
  sigma.pcaseq <- sum(ncopies*allele.ps*(1-allele.ps)) # standard deviation
  grm <- grm/sigma.pcaseq
  return(grm)
}
