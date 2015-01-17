# function to find the population structure

# genodat: genotype data should be a matrix with 0,1,2 for the genotype
# method: method to calculate the genetic relatedness matrix string: pcaseq, eigenstrat
# npcs: number of Principle Componenets (eigen vectors & eigenvalues) 
#      to return (integer between 1 and # of subjects)
# not sure what else you would want to see here

popStruct <- function(genodat, method, npcs)
{
  # check a bunch of stuff about the genotype 
  
}