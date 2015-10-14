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
      if(as.integer(ecnt) == ecnt & ecnt >= 0){
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
# Checks the maf paramter for either NA or a string that has the form of an
# interval with two numbers between 0 and 0.5 that are not equal and with the
# first strictly less than the second

checkMaf <- function(maf){
  if (!is.na(maf) & !is.character(maf) | is.nan(maf)){
    stop("MAF should be NA or a character. See help(seqPCA) for more details.")
  } else  if(is.character(maf)){
    mafMin <- as.numeric(mafBound(maf, 1))
    mafMax <- as.numeric(mafBound(maf, 2))
    
    if(mafMin < 0 | mafMax > 0.5 | mafMin >= mafMax)
    {
      stop("MAF bounds should be between 0 and 0.5 and given in min, max order.")
    } else {
      return(TRUE)
    } 
  } else if(is.na(maf)){
    return(TRUE)
  }
}


# function to get the MAF endpoints and check for consistency
mafBound <- function(maf, num){
  mafs <- unlist(strsplit(maf, split = ","))
  mafVal <- mafs[num]
  
  mafVal <- gsub("[", "", mafVal, fixed = TRUE)
  mafVal <- gsub("]", "", mafVal, fixed = TRUE)
  mafVal <- gsub("(", "", mafVal, fixed = TRUE)
  mafVal <- gsub(")", "", mafVal, fixed = TRUE)
  
  mafVal <- as.numeric(mafVal)
  
  return(mafVal)
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
checkSnp <- function(userSnp, dataSnp) {
  if (length(userSnp) > length(dataSnp)) {
    stop("More SNP IDs given than are in the genotype data set.")
  } else if (length(userSnp) <= 0){
    stop("Snp.id vector specified has length of 0.")
  } else if (any(!(userSnp %in% dataSnp)))
  {
    stop("Snp.id vector specified has SNP IDs not in the data set.")
  } else {
    return(TRUE)
  }

}


# Check weights ----------------------------------------------------------------
checkWeights <- function(weights) {
  if ( !(class(weights) %in% c("integer", "numeric")) ) {
    stop("Weights should be a vector of two numbers.")
  } else {
    if(length(weights) != 2) {
      stop("Weights should be a vector of length two.")
    } else {
      if(weights[1] > 0 & !is.nan(weights[1]) &
           weights[2] > 0 & !is.nan(weights[2])){
        return(TRUE)
      } else {
        stop("Weights should be positive numbers.")
      }
    }
  }
}