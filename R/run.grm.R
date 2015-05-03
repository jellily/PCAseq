# run.grm

run.grm <- function(gdsobj, sample.id, snp.id, autosome.only, remove.monosnp, 
                    maf, missing.rate, method)
{
  # Open the file
  geno.dat <- snpgdsOpen(gdsobj)
  
  # Find the number of subjects and use it to create
  # an empty GRM matrix
  if(is.null(sample.id))
  {
    sample.id <- read.gdsn(index.gdsn(geno.dat, "sample.id"))
  }
  
  # call the appropraite function based on the method
  grm <- switch(method, 
                eigen = grm.eigen(geno.dat, sample.id, snp.id, autosome.only,remove.monosnp, 
                                  maf, missing.rate),
                pcaseq = grm.pcaseq(geno.dat, sample.id, snp.id, autosome.only,remove.monosnp, 
                                     maf, missing.rate))

  return(grm)
}