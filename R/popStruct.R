# function to find the population structure

# genodat: genotype data should be a matrix with 0,1,2 for the genotype
# method: method to calculate the genetic relatedness matrix string: pcaseq, eigenstrat
# npcs: number of Principle Componenets (eigen vectors & eigenvalues) 
#      to return (integer between 1 and # of subjects)
# not sure what else you would want to see here


popStruct <- function(genodat, method, npcs)
{
  
  # validity checks
  # check if the data set is a vector
  if( class(genodat) == "numeric" | class(genodat) == "data.frame" )
  {
    stop("Genotype data is not a matrix. See help(popStruct) for details.")
  }
  
  # check for missing genotypes
  if(any(is.na(genodat)))
  {
    stop("Missing genotype data.")
  }
  
  # check for genotype data in the wrong format
  # includes letters, non-integers, etc.
  if( any( !(genodat %in% 0:2) ) )
  {
    stop("Genotype data is not 0, 1, or 2. See help(popStruct) for details.")
  }
  
  # check that there are at least 2 subjects in the data set
  if( dim(genodat)[2] < 2)
  {
    stop("Data set includes only 1 subject.")
  }
  
  # based on method, call the appropriate function
  grm <- switch(method,
         eigen = grm.eigen(genodat),
         pcaseq = grm.pcaseq(genodat))
  
  return(grm)
}
