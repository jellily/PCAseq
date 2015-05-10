# filter.auto

filter.auto <-function(snps, snp.chromosome){
  snp.chromsome <- as.numeric(snp.chromosome)
  snps <- snps[as.numeric(snp.chromosome) %in% 1:22, ]
  return(snps)
}