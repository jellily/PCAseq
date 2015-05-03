# check.ecnt

check.ecnt <- function(ecnt)
{
  if(is.numeric(ecnt)) # if ecnt is a number
  {
    if(length(ecnt) == 1) # and there is only one value
    {
     if(as.integer(ecnt) == ecnt) # and it is an integer
     {
       if(ecnt > 0) # and that value is greater than zero
       {
         return(TRUE)
       }else # otherwise, errors
       {
         stop("Number of eigenvectors and eigenvalues should be a counting number.")
       }
     }else
     {
       stop("Number of eigenvectors to return should be a counting number.")
     }
    }else
    {
      stop("Can only specify one value for the number of eigenvectors and eigenvalues.")
    }
  }else
  {
    stop("Number of eigenvalues and eigenvectors to return should be a counting number." )
  }
}