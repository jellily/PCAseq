# check.bool

check.bool <- function(bool)
{
  if(is.logical(bool)) # if bool is a logical
  {
    if(!is.na(bool)) # and is not NA
    {
      return(TRUE)
    }else  # otherwise, errors
    {
      stop("Input set to NA. See help for details.")
    }
    
  }else
  {
    stop("Input is not logical. See help for details.")
  }
}