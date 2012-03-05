\name{dot.prod}
\encoding{latin1}
\Rdversion{1.1}
\alias{dot.prod}
\alias{norm2}
\title{
  Computes Inner Product and Squared 2-norm
}
\description{
  \code{dot.prod} computes the inner (or dot/scalar) product between two
  vectors.

  \code{norm2} computes the squared 2-norm of all the elements in a
  matrix or vector.
}
\usage{
dot.prod(v1, v2)

norm2(v1)
}
\arguments{
  \item{v1,v2}{Two vectors}
}
\details{
  If the vectors are of unequal length \code{dot.prod} will give a
  warning and then truncates the longer vector, discarding any excess
  elements before the computations.
}
\value{
  \code{dot.prod} returns the inner product of v1 and v2. \code{norm2}
  returns the squared 2-norm of all elements in v1.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
   Block matrices are created by \code{\link{make.sigma.B}}, 
   \code{\link{make.sigma.B.full}} and \cr
   \code{\link{make.sigma.nu}}.
  
  Other block matrix functions \code{\link{makeCholBlock}}, 
  \code{\link{block.mult}}, \code{\link{calc.tF.times.mat}}, \cr
  \code{\link{calc.iS.X}}, and \code{\link{sumLogDiag}}.
  
  This function is called by \code{\link{loglike}}.
}
\examples{
##Create two vectors of equal length
v1 <- rnorm(10)
v2 <- rnorm(10)

##compute the inner product between the vectors
dot.prod(v1,v2)
##or
sum(v1*v2)

##compute the square 2-norm of v1
norm2(v1)
##or
dot.prod(v1,v1)
##or
sum(v1*v1)

##If the vectors are of unequal length the longer vector
##gets truncated (with a warning). 
dot.prod(v1,c(v2,2))
\dontshow{
if( abs(dot.prod(v1,v2)-sum(v1*v2)) > 1e-10 ){
  stop("dot.prod: Results not equal")
}
if( (abs(norm2(v1) - dot.prod(v1,v1)) > 1e-10) ||
    (abs(norm2(v1) - sum(v1*v1)) > 1e-10) ){
  stop("norm2: Results not equal")
}
}%\dontshow
}%\examples
