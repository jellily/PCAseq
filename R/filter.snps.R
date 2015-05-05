# filter.snps

filter.snps <- function(snps, autosome.only, remove.monosnp, missing.rate, maf){
  # remove monomorphic snps
  if (remove.monosnp){
    snps <- filter.mono(snps)
  }
  
  # remove sex chromosome snps
  if (autosome.only){
    snps <- filter.auto(snps)
  }
  
  # remove snps with too much missingness
  if (!is.nan(missing.rate)){
    snps <- filter.miss(snps, missing.rate)
  }
  
  # filter based on MAF
  if (!is.nan(maf)){
    snps <- filter.maf(snps, maf)
  }
  
  return(snps)
}