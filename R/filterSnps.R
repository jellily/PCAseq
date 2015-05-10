# filter.snps

filterSnps <- function(snps, autosome.only, remove.monosnp, missing.rate, maf, snp.chromosome){
  # remove monomorphic snps
  if (remove.monosnp){
    snps <- filterMono(snps)
  }
  
  # remove sex chromosome snps
  if (autosome.only){
    snps <- filterAuto(snps, snp.chromosome)
  }
  
  # remove snps with too much missingness
  if (!is.nan(missing.rate)){
    snps <- filterMiss(snps, missing.rate)
  }
  
  # filter based on MAF
  if (!is.nan(maf)){
    snps <- filterMaf(snps, maf)
  }
  
  return(snps)
}