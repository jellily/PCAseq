# getMiss

getMiss <- function(snps)
{
  byrows <- 1
  
  # function to calculate the proportion missing for one snp
  propMiss <- function(snp)
  {
    mean(ifelse(is.na(snp), 1, 0))
  }
  
  missing <- apply(snps, byrows, FUN = propMiss)
  
  return(missing)
}