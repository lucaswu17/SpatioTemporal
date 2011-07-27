\name{make.sigma.nu}
\encoding{latin1}
\Rdversion{1.1}
\alias{make.sigma.nu}
\alias{make.sigma.nu.cross.cov}

\title{
  Crete (Cross)-Covariance Matrices for the Residual Field
}
\description{
  Creates covariance and cross-covariance matrices for the space-time
  residual nu field.
}
\usage{
make.sigma.nu(sill, nugget, range, block.sizes, loc.index, dists)

make.sigma.nu.cross.cov(sill, nugget, range, loc.ind1, loc.ind2 =
   loc.ind1, loc.ind2.to.1 = 1:max(loc.ind2), T1, T2 = T1, dists)
}
\arguments{
  \item{sill, nugget, range}{
    Sill, nugget, and range parameters for the exponential covariance
    function.
  }
  \item{block.sizes}{
    A vector of containing the size(s) of each of the diagonal blocks,
    usually \cr \code{mesa.data.model$nt}.
  }
  \item{loc.index, loc.ind1, loc.ind2}{
    Vectors indicating the location of each element in the covariance
    matrix, used to index the \code{dist}-matrix to determine the
    distance between locations, usually \code{mesa.data.model$obs$idx}.
  }
  \item{dists}{
    Distance matrix giving the distances between the (observed)
    locations.
  }
  \item{loc.ind2.to.1}{
    A vector that translates locations in \code{loc.ind2} to locations
    in \code{loc.ind1}; used to determine if points are collocated.
  }
  \item{T1, T2}{
    Vectors with the observation times of each space-time point; used to
    determine if observtations are coincidental, usually
    \code{mesa.data.model$obs$date}. 
  }
}
\value{
  Returns the block-diagonal covariance or cross-covariance matrix for
  the spatio-temporal residual nu field.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom} and Adam Szpiro
}
\seealso{
\code{\link{make.sigma.B.full}}, \code{\link{make.sigma.B}} 

This function is called by \code{\link{loglike}}.
}
\examples{
##First create some random locations
x <- rnorm(5)
y <- rnorm(5)

##compute distance matrix
D <- as.matrix(dist( cbind(x,y) ))

#a vector of locations
I <- c(1,2,3,1,4,4,3,2,1,1)
T <- c(1,1,1,2,2,3,3,3,3,4)

##create a block diagonal matrix consisting of four parts with
##exponential covariance using same sill, range, and nugget, but
##different locations in each block.
sigma.nu <- make.sigma.nu(sill=2, nugget=.1, range=.4,
    block.sizes=c(3,2,4,1), dists=D, loc.index=I)
print(sigma.nu)

sigma.nu.c <- make.sigma.nu.cross.cov(sill=2, nugget=.1, range=.4,
     loc.ind1=I, loc.ind2=I[7:10], T1=T, T2=T[7:10], dists=D)
print(sigma.nu.c)
##compare the cross-covariance with the relevant part of sigma.nu
range(sigma.nu.c-sigma.nu[,7:10])
\dontshow{
if( abs(max(sigma.nu.c-sigma.nu[,7:10])) > 1e-10 ){
  stop("make.sigma.nu.cross.cov 1: Results not equal")
}
}
##an alternative showing the use of loc.ind2.to.1
sigma.nu.c <- make.sigma.nu.cross.cov(sill=2, nugget=.1, range=.4,
     loc.ind1=I, loc.ind2=1:2, loc.ind2.to.1=4:3,
     T1=T, T2=c(3,3), dists=D[,4:3])
print(sigma.nu.c)
##compare the cross-covariance with the relevant part of sigma.nu
range(sigma.nu.c-sigma.nu[,6:7])
\dontshow{
if( abs(max(sigma.nu.c-sigma.nu[,6:7])) > 1e-10 ){
  stop("make.sigma.nu.cross.cov 1: Results not equal")
}
}
}
