# check.maf

checkMaf <- function(maf)
{
  if(is.numeric(maf)) # if maf is a numeric value (vector or singleton)
  {
    if(length(maf) == 1) # has length one
    {
      if(is.nan(maf))
      {
        return(TRUE)
      }else if(maf >= 0 & maf <= 0.5) # and is in [0, 0.5]
      {
        return(TRUE)
      }else  # otherwise, an error
      {
        stop("MAF needs to be between 0 and 0.5. See help(seqPCA) for more details.")
      }
    }else if(length(maf) == 2) # OR has length two
    {
      if(maf[1] >= 0 & maf[2] > maf[1] & maf[2] <= 0.5) # and 0 <= maf[1] < maf[2] <= 0.5
      {
        return(TRUE)
      }else # otherwise return errors
      {
        stop("If MAF is a vector of length 2, please specify (min, max) in [0, 0.5]. See help(seqPCA) for more details.") 
      }
    }else
    {
      stop("MAF needs to be of length 1 or 2. See help(seqPCA) for more details.")
    }
  }else
  {
    stop("MAF should be numeric. See help(seqPCA) for more details.")
  }

}