# get.miss

get.miss(snps)
{
  byrows <- 1
  missing <- apply(snps, byrows, FUN = mean(ifelse(is.na(x), 1, 0)) )
  return(missing)
}