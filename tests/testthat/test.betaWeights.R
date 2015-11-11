betaWeights <- function(alleleFreq, alpha, beta) {

  # find the minor allele frequency
  minorFreq <- calcMaf(alleleFreq)

  return(sqrt(dbeta(minorFreq, alpha, beta)))
}


alleleFreq <- seq(0, 1, length.out = 0.001)
alpha <- seq(0.01, 2, length.out = 5)
beta <- seq(1, 50, length.out = 10)

dat <- expand.grid(alleleFreq, alpha, beta)
names(dat) <- c("alleleFreq", "alpha", "beta")

weights <- mapply(FUN = betaWeights, dat$alleleFreq, dat$alpha, dat$beta)