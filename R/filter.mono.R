# filter.mono

filter.mono <-function(snps){
  
  # find the allele frequencies
  allele.freq <- 0.5*colMeans(snps)
  
  # remove monomorphic SNPs
  snps <- snps[ ,allele.freq > 0 & allele.freq < 1]
  
  # check to make sure there are still SNPs in the data set
  if (class(snps) != "matrix"){
    stop("All SNPs are monomorphic. No data remains removing monomorphic SNPs.")
  }else if(ncol(snps) == 0){
    stop("All SNPs are monomorphic. No data remains removing monomorphic SNPs.")
  }else if (ncol(snps) < 50){
    message("Fewer than 50 SNPs remain after removing monomorphic SNPs.")
    return(snps)
  }else{
    return(snps)
  }
}