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
          need.genmat = FALSE, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gdsobj}{an object of the class SNPGDSFileClass, a SNP GDS file}
  \item{method}{}
  \item{sample.id}{a vector of sample ids specifying the samples to use for analysis; if NULL, all samples are used}
  \item{snp.id}{a vector of SNP ids specifying the SNps to use for analysis; if NULL, all SNPs are used}
  \item{autosome.only}{not supported}
  \item{remove.monosnp}{if TURE, remove monomorphic SNPs}
  \item{maf}{if one number is specified, use SNPs with MAF greater than or equal to this value; if a numeric vector of length two is specified, only SNPs with MAFs from (min, max) are taken}
  \item{missing.rate}{not supported}
  \item{need.genmat}{if TRUE, return the genetic relatedness matrix}
  \item{verbose}{not supported}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
Return a snpgdsPCAClass object, and it is a list:
\item{sample.id}{the sample ids used in the analysis}
\item{snp.id}{the SNP ids used in the analysis}
\item{eigenval}{eigenvalues}
\item{eigenvect}{a matrix of eigenvectors of dimensions # of samples by eigen.cnt}
\item{varprop}{the proportion of the variance explained by each principal component}
\item{TraceXTX}{the trace of the genetic relateness matrix}
\item(Bayesian){indicates Bayes normalization; set to FALSE, as this is not currently supported}
\item{genmat}{the genetic relateness matrix}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (x) 
{
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
