# run.grm

runGrm <- function(gdsobj, sample.id, snp.id, autosome.only, remove.monosnp, 
                    maf, missing.rate, method)
{
  # Open the file
  geno.dat <- snpgdsOpen(gdsobj)
  
  # If no sample ids are given, take all of them
  # otherwise check that the sample.id vector
  # has no duplicates and contains sample ids in the
  # file sample.id vector
  if (is.null(sample.id)){
    sample.id <- read.gdsn(index.gdsn(geno.dat, "sample.id"))
  }else{
    checkSamp(sample.id, read.gdsn(index.gdsn(geno.dat, "sample.id")))
  }
  
  # if no SNP ids are given, take all of them
  # (note SNP filtering happens at a later time point)
  if (is.null(snp.id)){
    snp.id <- read.gdsn(index.gdsn(geno.dat, "snp.id"))
  }else{
    checkSnp(snp.id, read.gdsn(index.gdsn(geno.dat, "snp.id")))
  }
  
  # Determine the orientation of the file -- the code to calculate
  # the GRM is written for SNP x SAMPLE data, if the data is not
  # in that order, it is transposed as it is read in
  if (names(get.attr.gdsn(index.gdsn(geno.dat, "genotype"))) == "sample.order"){
    transpose <- TRUE
  }else{
    transpose <- FALSE
  }
  
  # call the appropraite function based on the method
  grm <- switch(method, 
                eigen = grmEigen(geno.dat, sample.id, snp.id, autosome.only,remove.monosnp, 
                                  maf, missing.rate, transpose),
                pcaseq = grmPcaseq(geno.dat, sample.id, snp.id, autosome.only,remove.monosnp, 
                                     maf, missing.rate, transpose))
  
  # Close the file
  snpgdsClose(geno.dat)
  return(grm)
}