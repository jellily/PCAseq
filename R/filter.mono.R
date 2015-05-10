# filter.mono

filter.mono <-function(snps){
  
  # find the allele frequencies
  allele.freq <- 0.5*rowMeans(snps, na.rm = TRUE)
  
  # remove monomorphic SNPs
  snps <- snps[allele.freq > 0 & allele.freq < 1, ]
  
  # check to make sure there are still SNPs in the data set
  if (class(snps) != "matrix" | nrow(snps) == 0){
    stop("All SNPs are monomorphic. No data remains after removing monomorphic SNPs.")
  }else if (nrow(snps) < 50){
    message("Fewer than 50 SNPs remain after removing monomorphic SNPs.")
    return(snps)
  }else{
    return(snps)
  }
}