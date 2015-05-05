# run.grm

run.grm <- function(gdsobj, sample.id, snp.id, autosome.only, remove.monosnp, 
                    maf, missing.rate, method)
{
  # Open the file
  geno.dat <- snpgdsOpen(gdsobj)
  
  # If no sample ids are given, take all of them
  if (is.null(sample.id)){
    sample.id <- read.gdsn(index.gdsn(geno.dat, "sample.id"))
  }
  
  # if no SNP ids are given, take all of them
  # (note SNP filtering happens at a later time point)
  if (is.null(snp.id)){
    snp.id <- read.gdsn(index.gdsn(geno.dat, "snp.id"))
  }
  
  # call the appropraite function based on the method
  grm <- switch(method, 
                eigen = grm.eigen(geno.dat, sample.id, snp.id, autosome.only,remove.monosnp, 
                                  maf, missing.rate),
                pcaseq = grm.pcaseq(geno.dat, sample.id, snp.id, autosome.only,remove.monosnp, 
                                     maf, missing.rate))
  
  # Close the file
  snpgdsClose(geno.dat)
  return(grm)
}