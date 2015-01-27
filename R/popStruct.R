# function to find the population structure

# genodat: genotype data should be a matrix with 0,1,2 for the genotype
# method: method to calculate the genetic relatedness matrix string: pcaseq, eigenstrat
# npcs: number of Principle Componenets (eigen vectors & eigenvalues) 
#      to return (integer between 1 and # of subjects)
# not sure what else you would want to see here


popStruct <- function(genodat, method, npcs)
{ 
  # validity checks
  # check if the data set is in GDS format
  if( class(genodat) != "gds.class" )
  {
    stop("Genotype data is not GDS format. See help(popStruct) for details.")
  }
  
  # based on method, call the appropriate function
  grm <- switch(method,
         eigen = grm.eigen(genodat, allele.freq),
         pcaseq = grm.pcaseq(genodat, allele.freq))
  
  return(grm)
}