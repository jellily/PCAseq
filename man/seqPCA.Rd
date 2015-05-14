\name{seqPCA}
\alias{seqPCA}
\alias{seqGRM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Principal Component Analysis (PCA) with Sequence Data}
\description{
Calculate the eigenvectors and eigenvalues for principal component analysis of sequence data. This function is designed to extend the functionality of the SNPRelate package by implementing the PCA-seq method. 
}
\usage{
seqPCA(gdsobj, method, sample.id = NULL, snp.id = NULL, autosome.only = TRUE,
          remove.monosnp = TRUE, maf = NaN, missing.rate = NaN, eigen.cnt = 32,
          need.genmat = FALSE, verbose = TRUE)}
\usage{
seqGRM(gdsobj, method, sample.id = NULL, snp.id = NULL, autosome.only = TRUE,
                   remove.monosnp = TRUE, maf = NaN, missing.rate = NaN,
                   verbose = TRUE)
}

\arguments{
  \item{gdsobj}{an object of the class SNPGDSFileClass, a SNP GDS file}
  \item{method}{a string indicating which method to use to calculate the GRM; see Details}
  \item{sample.id}{a vector of sample ids specifying the samples to use for analysis; if NULL, all samples are used}
  \item{snp.id}{a vector of SNP ids specifying the SNps to use for analysis; if NULL, all SNPs are used}
  \item{autosome.only}{if TRUE, use autosomal SNPs only; if it is a numeric or character vector, keep SNPs according to the specified chromosomes}
  \item{remove.monosnp}{if TURE, remove monomorphic SNPs}
  \item{maf}{if one number is specified, use SNPs with MAF greater than or equal to this value; if a numeric vector of length two is specified, only SNPs with MAFs in (min, max) are taken}
  \item{missing.rate}{to use the SNPs with missing rates less than or equal to missing.rate; if NaN, no misisng threshold}
  \item{eigen.cnt}{the number of eigen vectors and values to return; if zero, return all eigenvalues and vectors}
  \item{need.genmat}{if TRUE, return the genetic relatedness matrix}
  \item{verbose}{not supported}
}
\details{
Some more details about the package here...}
\value{
Return a snpPCAClass object, a list with the follow slots:
\item{sample.id}{the sample ids used in the analysis}
\item{snp.id}{the SNP ids used in the analysis}
\item{eigenval}{eigenvalues}
\item{eigenvect}{a matrix of eigenvectors of dimensions # of samples by eigen.cnt}
\item{varprop}{the proportion of the variance explained by each principal component}
\item{TraceXTX}{the trace of the genetic relateness matrix}
\item{Bayesian}{indicates Bayes normalization; set to FALSE, as this is not currently supported}
\item{genmat}{the genetic relateness matrix}
}

\author{
Jennifer L Kirk
}


\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
