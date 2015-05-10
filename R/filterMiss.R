# filter.miss

filterMiss <-function(snps, missing.rate){
  
  # replace the missing code 3 with
  # NA
  #snp[snp == 3] <- NA
  
  # find the proportion missing for each SNP
  #missing <- get.missing(snps)
  
  #snps <- snps[misisng <= missing.rate, ]
  
  # check to make sure there are still SNPs in the data set
  #if (class(snps) != "matrix" | nrow(snps) == 0){
  #  stop("All SNPs have missing rates above specified threshold. No data remains after missingness filtering.")
  #}else if (nrow(snps) < 50){
  #  message("Fewer than 50 SNPs remain after filtering by missing rate.")
  #  return(snps)
  #}else{
  #  return(snps)
  #}
  
  return(snps)
}