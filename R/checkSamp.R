# check.vector function

checkSamp <- function(user.samp, data.samp)
{
  if (length(user.samp) > length(data.samp)){
    stop("More sample IDs given than are in the genotype data set.")
  }else if (length(user.samp) <= 0){
    stop("Sample.id vector specified has length of 0.")
  }else if (any(!(user.samp %in% data.samp)))
  {
    stop("Sample.id vector specified has sample IDs not in the data set.")
  }else{
    return(TRUE)
  }
  
}