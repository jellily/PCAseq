# check.miss

checkMiss <- function(miss){
  if (is.numeric(miss)){ # if the missing rate is a number
    if (length(miss) == 1){ # and only one is given
      if (is.nan(miss)){ # and it is NaN
        return(TRUE)
      }else 
      {
        if (miss > 0 & miss < 1){ # or if the missing rate is in (0,1)
          return(TRUE)
        }else if(miss == 0){ # or is 0
          message("Missing rate set to 0, all SNPs with missing data will be removed.")
        }else if (miss == 1){ # or is 1
          warning("Missing rate set to 1. To turn off missingness removal, set missing rate to NaN.")
        } else{  # otherwise, errors
          stop("Missing rate should be in (0, 1)")
        }
      }
    }else{
      stop("Only a single missing rate can be specified.")
    }
  }else{
    stop("Missing rate should be numeric or NaN. See help for details.")
  }
}