# Functions to check the validity of parameters passed in by the user


# Check logicals ---------------------------------------------------------------
# Checks for a valid logical that is not NA.

checkBool <- function(bool){
  if (is.logical(bool)){
    if (!is.na(bool)){
      return(TRUE)
    }else{
      stop("Input set to NA. See help for details.")
    }   
  } else {
    stop("Input is not logical. See help for details.")
  }
} 


# Check eigen.cnt --------------------------------------------------------------
# Check that the eigen.cnt parameter for a single integer number that is 0 or 
# greater.

checkEcnt <- function(ecnt){
  if (is.numeric(ecnt)){
    if (length(ecnt) == 1){
      if(as.integer(ecnt) == ecnt){ 
        return(TRUE)
      } else {
        stop("Number of eigenvectors to return should be a counting number.")
      }
    } else {
      stop("Can only specify one value for the number of eigenvectors 
           and eigenvalues.")
    }
  } else {
    stop("Number of eigenvalues and eigenvectors to return should be a 
         counting number." )
  }
}


# Check maf --------------------------------------------------------------------
# Checks the maf paramter for either a single number in [0,0.5] or two values
# in [0, 0.5] that are not equal with the first value strictly less than the
# second value.

checkMaf <- function(maf){
  if (is.numeric(maf)){
    if (length(maf) == 1){
      if (is.nan(maf)){
        return(TRUE)
      } else if(maf >= 0 & maf <= 0.5) {
        return(TRUE)
      } else {
        stop("MAF needs to be between 0 and 0.5. See help(seqPCA) for more 
             details.")
      }
    } else if(length(maf) == 2) {
      if(maf[1] >= 0 & maf[2] > maf[1] & maf[2] <= 0.5) {
        return(TRUE)
      } else {
        stop("If MAF is a vector of length 2, please specify (min, max) in 
             [0, 0.5]. See help(seqPCA) for more details.") 
      }
    } else {
      stop("MAF needs to be of length 1 or 2. See help(seqPCA) for more 
           details.")
    }
  } else {
    stop("MAF should be numeric. See help(seqPCA) for more details.")
  }
  
}


# Check missing.rate -----------------------------------------------------------
# Check that miissing.rate is a single number in [0, 1) or is NaN.

checkMiss <- function(miss){
  if (is.numeric(miss)){ 
    if (length(miss) == 1){ 
      if (is.nan(miss)){ 
        return(TRUE)
      } else {
        if (miss > 0 & miss < 1){ 
          return(TRUE)
        }else if(miss == 0){ 
          message("Missing rate set to 0, all SNPs with missing data will be
                  removed.")
        }else if (miss == 1){ # or is 1
          warning("Missing rate set to 1. To turn off missingness removal, set 
                  missing rate to NaN.")
        } else {  
          stop("Missing rate should be in [0, 1)")
        }
      }
    }else {
      stop("Only a single missing rate can be specified.")
    }
  }else  {
    stop("Missing rate should be numeric or NaN. See help for details.")
  }
}


# Check sample.id --------------------------------------------------------------
# check that the sample.id vector is no longer than the number of sample IDs in 
# the GDS file, has at least one entry, and does not contain sample IDs not in 
# the GDS file.

checkSamp <- function(userSamp, dataSamp)
{
  if (length(userSamp) > length(dataSamp)){
    stop("More sample IDs given than are in the genotype data set.")
  }else if (length(userSamp) <= 0){
    stop("Sample.id vector specified has length of 0.")
  }else if (any(!(userSamp %in% dataSamp)))
  {
    stop("Sample.id vector specified has sample IDs not in the data set.")
  }else{
    return(TRUE)
  }
  
}


# Check snp.id -----------------------------------------------------------------
# Check that the the snp.id vector is no longer than the number of SNP IDs in
# the GDS file, has at least one entry, and does not contain SNP IDs not in the
# GDS file.
checkSnp <- function(userSnp, dataSnp)
{
  if (length(userSnp) > length(dataSnp)){
    stop("More SNP IDs given than are in the genotype data set.")
  }else if (length(userSnp) <= 0){
    stop("Snp.id vector specified has length of 0.")
  }else if (any(!(userSnp %in% dataSnp)))
  {
    stop("Snp.id vector specified has SNP IDs not in the data set.")
  }else{
    return(TRUE)
  }
  
}