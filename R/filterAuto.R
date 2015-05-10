# filter.auto

filterAuto <-function(snps, snp.chromosome){
  
  # Convert the vector of SNP choromosome labels to numeric
  # Any character label will become NA
  snp.chromsome <- as.numeric(snp.chromosome)
  
  # Select only those SNPs with chromosome labels 1-22 (autosomes)
  snps <- snps[as.numeric(snp.chromosome) %in% 1:22, ]
  
  return(snps)
}