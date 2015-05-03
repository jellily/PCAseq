get.snps <- function( i, nblock, nsnps)
{
  snps <- (1:nblock) + (i-1)*nblock
  if(snps[nblock] > nsnps)
  {
    snps <- snps[1]:nsnps
  }
  
  return(snps)
}