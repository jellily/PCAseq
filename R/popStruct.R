# function to find the population structure

# genodat: genotype data should be a matrix with 0,1,2 for the genotype
# method: method to calculate the genetic relatedness matrix string: pcaseq, eigenstrat
# npcs: number of Principle Componenets (eigen vectors & eigenvalues) 
#      to return (integer between 1 and # of subjects)
# freq.min: mimum minor allele frequency
# freq.max: maximum minor allele frequency


popStruct <- function(genodat, method, npcs, freq.min = NA, freq.max = NA)
{ 
  # validity checks
  # check if the data set is in GDS format
  if( class(genodat) != "gds.class" )
  {
    stop("Genotype data is not GDS format. See help(popStruct) for details.")
  }
  
  # based on method, call the appropriate function
  
  # find the principle components
  
  # return the grm and the specified number of principle components
  
}


