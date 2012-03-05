\name{make.sigma.B.full}
\encoding{latin1}
\Rdversion{1.1}
\alias{make.sigma.B.full}

\title{
  The Covariance Matrix for Temporal Trends and beta Fields.
}
\description{
  Computes the full (cross)-covariance matrix for the beta fields
  multiplied by the temporal trends; i.e. F*sigmaB*F' where
  \code{sigmaB} is a block diagonal matrix, usually the output from
  \cr \code{\link{make.sigma.B}}.
}
\usage{
make.sigma.B.full(sill, range, loc.ind1, loc.ind2 = loc.ind1,
                  F1, F2 = F1, dists)
}
\arguments{
  \item{sill, range}{
    Vectors with the sill and range parameters for each of the block
    matrices.
  }
  \item{loc.ind1,loc.ind2}{
    Vectors indicating which location the rows in \code{F1} and
    \code{F2} corresponds to, usually \code{mesa.data.model$obs$idx}.
  }
  \item{F1,F2}{
    (number of obs.) - by - (number of temporal trends) matrices
    containing the temporal trends. Usually \code{mesa.data.model$F},
    where \code{\link{mesa.data.model}} is obtained from
    \code{\link{create.data.model}}.
  }
  \item{dists}{
    Distance matrix giving the distances between the (observed)
    locations.
  }
}
\value{
  Returns the covariance or cross-covariance matrix F'*sigmaB*F.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
\code{\link{make.sigma.B}}, \code{\link{make.sigma.nu}} 

This function is called by \code{\link{loglike}}.
}
\examples{
##First create some random locations
x <- rnorm(3)
y <- rnorm(3)

##compute distance matrix
D <- as.matrix(dist( cbind(x,y) ))

##create a trend
trend <- cbind(rep(1,5),sin(1:5))
##an index of locations
idx <- c(rep(1:3,3),1:2,2:3)
##expand the F matrix to match the locations/times in idx/T.
F <- trend[c(rep(1:3,each=3),4,4,5,5),]

##compute F*sigmaB*F'
FsigmaF <- make.sigma.B.full(c(2,1), c(.3,1), loc.ind1=idx,
                             F1=F, dists=D)

##create the full F matrix
Falt <- matrix(0,dim(F)[1],dim(F)[2]*max(idx))
for(i in 1:dim(F)[1]){
  Falt[i,c(idx[i],idx[i]+max(idx))] <- F[i,]
}
##and compute F*sigmaB*F' using full matrices
sigma.B <- make.sigma.B(c(2,1), c(.3,1), D)
FsigmaF.alt <- Falt \%*\% sigma.B \%*\% t(Falt)

##compare the results
range(FsigmaF - FsigmaF.alt)
\dontshow{
if( max(abs(FsigmaF - FsigmaF.alt)) > 1e-10 ){
  stop("make.sigma.B.full: Results not equal.")
}
}
##use the function to compute a cross covariance
FsigmaF.cross <- make.sigma.B.full(c(2,1), c(.3,1), loc.ind1=idx[1:6],
      loc.ind2=idx[7:13], F1=F[1:6,], F2=F[7:13,], dists=D)

##compare to the relevant part of the full matrix
range(FsigmaF.cross-FsigmaF[1:6,7:13])
\dontshow{
if( max(abs(FsigmaF.cross-FsigmaF[1:6,7:13])) > 1e-10 ){
  stop("make.sigma.B.full: Results not equal - cross-covariance")
}
}
}
