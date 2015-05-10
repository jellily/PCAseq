# check.vector function

checkSnp <- function(user.snp, data.snp)
{
  if (length(user.snp) > length(data.snp)){
    stop("More SNP IDs given than are in the genotype data set.")
  }else if (length(user.snp) <= 0){
    stop("Snp.id vector specified has length of 0.")
  }else if (any(!(user.snp %in% data.snp)))
  {
    stop("Snp.id vector specified has SNP IDs not in the data set.")
  }else{
    return(TRUE)
  }
  
}