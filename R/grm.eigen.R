# grm.eigen

grm.eigen <- function(geno.dat, sample.id, snp.id, autosome.only, remove.monosnp, 
                      maf, missing.rate)
{
  # constants
  nblock <- 5000
  bycolumns <- 2
  ncopies <- 2
  nsubj <- length(sample.id)
  nsnps <- length(snp.id)
  max <- floor(nsnps/nblock) + 1  # maximum number of blocks to loop over

  # create empty grm & vector to count the number of snps used
  grm <- matrix(0, nrow = nsubj, ncol = nsubj)
  total.snps <- 0
  
  # Loop through the SNPs in blocks of size nblock
  for(i in 1:max)
  {
    message(paste("Computing GRM: Block", i, "of", max))
    
    # Get the number of SNPs to choose
    snps <- snp.id[get.snps(i, nblock, nsnps)]
    
    # Read in the relevant SNP data, subsetting by subject
    snp.dat <- snpgdsGetGeno(geno.dat, snp.id = snp.id[snps], sample.id = sample.id)  
    
    # Filter the data
    snp.dat <- filter.snps(snp.dat, autsome.only, remove.monosnp, missing.rate, maf)
    
    total.snps <- dim(snp.dat)[bycolumns] # total # of snps
    allele.freq <- (1/ncopies)*colMeans(snp.dat)  # snp allele frequencies
    
    # Estimate the variance at each SNP
    geno.cent <- sweep(snp.dat, bycolumns, STATS = ncopies*allele.freq)   
    sigma.eigen <- sqrt(allele.freq[snps]*(1-allele.freq)) # standard deviation
    
    # Find the empirical correlation matrix
    zee <- sweep(geno.cent, bycolumns, STATS = sigma.eigen, FUN = "/")
    grm <- grm + tcrossprod(zee) 
  }
    
  grm <- (1/total.snps)*grm
  return(grm)
}
