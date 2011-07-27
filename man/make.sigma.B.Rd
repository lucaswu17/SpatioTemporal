\name{make.sigma.B}
\encoding{latin1}
\Rdversion{1.1}
\alias{make.sigma.B}

\title{
  Create Covariance Matrix for the beta-Fields
}
\description{
  Creates the block diagonal covariance matrix for the beta fields.
}
\usage{
make.sigma.B(sill, range, dists)
}
\arguments{
  \item{sill, range}{
    Vectors with the sill and range parameters for each of the block
    matrices.
  }
  \item{dists}{
    Distance matrix giving the distances between the (observed)
    locations.
  }
}
\value{
  Returns a block diagonal covariance matrix with one block for each
  of the elements in the vectors \code{sill} and \code{range}. Each
  block contains the exponential covariance matrix for the elements in
  the \code{dists} matrix.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom} and Adam Szpiro
}
\seealso{
\code{\link{make.sigma.B.full}}, \code{\link{make.sigma.nu}}.

This function is called by \code{\link{loglike}}.
}
\examples{
##First create some random locations
x <- rnorm(5)
y <- rnorm(5)

##compute distance matrix
D <- as.matrix(dist( cbind(x,y) ))

##create a block diagonal matrix consisting of three parts with
##different sill and range (nugget=0) and exponential covariance
##function
sill <- c(2,1,3)
range <- c(.3,1,.5)
make.sigma.B(sill, range, D)
\dontshow{
B <- make.sigma.B(sill, range, D)
B.alt <- matrix(0,length(sill)*dim(D)[1],length(sill)*dim(D)[1])
for(i in 1:length(sill)){
  B.alt[((i-1)*dim(D)[1]+1):(i*dim(D)[1]),((i-1)*dim(D)[1]+1):(i*dim(D)[1])] <- sill[i]*exp(-D/range[i])
}
if( abs(max(B-B.alt)) > 1e-10){
  stop("make.sigma.B: Results not equal")
}
}%\dontshow
}%\examples