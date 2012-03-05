\name{sumLogDiag}
\encoding{latin1}
\Rdversion{1.1}
\alias{sumLogDiag}
\title{
  Sum the Logarithm of Diagonal Elements
}
\description{
  Computes the sum of the logarithm of the diagonal elements in a
  matrix. This corresponds to the logarithm of the determinant for a
  Cholesky factor. Behaviour is undefined for any elements that are
  <=0.
}
\usage{
sumLogDiag(mat)
}
\arguments{
  \item{mat}{
    A square matrix (preferably a Cholesky factor).
  }
}
\value{
  Sum of the logarithm of the diagonal elements.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  \code{\link{chol}} or \code{\link{makeCholBlock}} for computation of
  Cholesky factors.
  
  
  Other block matrix functions \code{\link{dot.prod}}, 
  \code{\link{block.mult}}, \code{\link{calc.tF.times.mat}}, 
  and \cr \code{\link{calc.iS.X}}.
  
  This function is used by \code{\link{loglike}}.
}
\examples{
#Create a covariance matrix
S <- cbind(c(2,1),c(1,2))
#compute Cholesky factor
R <- chol(S)
#compute determinant
log(det(R))
#compare with sum of the logarithm of diagonal elements
sumLogDiag(R)
\dontshow{
if( abs(log(det(R)) - sumLogDiag(R)) > 1e-10 ){
  stop("sumLogDiag: Results not equal")
}
}%\dontshow
}%\examples