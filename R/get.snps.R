get.snps <- function(i, nblock, nsnps)
{
  index <- (1:nblock) + (i-1)*nblock
  if(index[nblock] > nsnps)
  {
    index <- index[1]:nsnps
  }
  
  return(index)
}