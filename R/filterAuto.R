# filter.auto

filterAuto <-function(snps, snp.chromosome){
  
  autosome.codes <- as.character(1:22)
  
  # Convert the vector of SNP choromosome labels to numeric
  # Any character label will become NA
  snp.chromsome <- as.character(snp.chromosome)
  
  # Select only those SNPs with chromosome labels 1-22 (autosomes)
  snps <- snps[snp.chromosome %in% autosome.codes, ]
  
  return(snps)
}