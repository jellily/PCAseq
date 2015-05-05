# filter.maf

filter.maf <-function(snps, maf){
  
  # find the allele frequencies
  allele.freq <- 0.5*colMeans(snps)
  
  # different filtering based on what value is specified
  if (length(maf) == 1){
    snps <- snps[, allele.freq > maf & allele.freq < (1-maf)]
  }else{ #if length(maf) == 2
    min <- maf[1]
    max <- maf[2]

    snps <- snps[, allele.freq >= min & allele.freq <= max |
                   allele.freq >= (1-max) & allele.freq <= (1-min)]
  }

  # check to make sure there are still SNPs in the data set
  if (class(snps) != "matrix"){
    stop("All SNPs have MAF below specified threshold. No data remains after MAF filtering.")
  }else if(ncol(snps) == 0){
    stop("All SNPs have MAF below specified threshold. No data remains after MAF filtering.")
  }else if (ncol(snps) < 50){
    message("Fewer than 50 SNPs remain after filtering by MAF.")
    return(snps)
  }else{
    return(snps)
  }
}