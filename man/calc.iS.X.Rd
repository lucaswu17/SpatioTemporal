\name{calc.iS.X}
\encoding{latin1}
\Rdversion{1.1}
\alias{calc.iS.X}
\alias{calc.X.iS.X}
\alias{calc.mu.B}

\title{
  Matrix Multiplication with Block Matrices
}
\description{
  Computes either the product between a block diagonal, square matrix
  \code{iS} and a block matrix \code{X}; the quadratic form of a block
  diagonal, square matrix, \code{t(X)*iS*X}; or a block matrix
  multiplied by a vector, \code{X*alpha}.
}
\usage{
calc.iS.X(X, iS, dim)

calc.X.iS.X(X, iS.X, dim)

calc.mu.B(dim, X, alpha)
}
\arguments{
  \item{X}{
    A list of matrices with which to form the block matrix. Usually \cr
    \code{mesa.data.model$X}, where \code{\link{mesa.data.model}} is
    obtained from \cr \code{\link{create.data.model}}. The list should
    contain \code{dim$m} matrices, with the i:th matrix being of size
    \code{dim$n} - by - \code{dim$p[i]}.
  }
  \item{iS}{
    A block diagonal, square matrix, with \code{dim$m} blocks
    each of size \cr \code{dim$n} - by - \code{dim$n}.
  }
  \item{iS.X}{
    Matrix containing the product of \code{iS} and \code{X}. Output from
    \code{calc.iS.X}.
  }
  \item{dim}{
    Information regarding parameter dimension, output from
    \code{\link{loglike.dim}}.
  }
  \item{alpha}{
    A list of \code{dim$m} vectors, with the i:th vector being
    of length \code{dim$p[i]}.
  }
}
\value{
  Returns iS*X, X'*iS*X, or X*alpha.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  Block matrices can be created by \code{\link{make.sigma.B}}, 
  \code{\link{make.sigma.B.full}} and \cr
  \code{\link{make.sigma.nu}}.
  
  Other block matrix functions \code{\link{makeCholBlock}}, 
  \code{\link{block.mult}}, \code{\link{calc.tF.times.mat}}, \cr
  \code{\link{dot.prod}}, \code{\link{sumLogDiag}}.
  
  This function is called by \code{\link{loglike}}.
}
\examples{
## Create a block diagonal matrix, ...
iS <- rbind(c(2,1,0,0),c(1,3,0,0),
            c(0,0,3,2),c(0,0,2,4))
## ... a block matrix ...
X <- list(matrix(c(1,2)),matrix(c(2,2,3,4),2,2))
## ... with alternative form, ...
Xt <- rbind(cbind(X[[1]],matrix(0,2,2)),
            cbind(matrix(0,2,1),X[[2]]))
## ... and a vector alpha.
alpha <- list(c(1),c(-2,1))
## Create matching dimension structure.
dim <- list(p=c(1,2), m=2, n=2)

## Compute iS * X
iS.X <- calc.iS.X(X,iS,dim)
## or
iS \%*\% Xt
\dontshow{
if( max(abs(iS.X - (iS \%*\% Xt))) > 1e-13 ){
  stop("calc.iS.X: Results not equal")
}
}
## Compute X'* iS * X
calc.X.iS.X(X, iS.X, dim)
## or
t(Xt) \%*\% iS \%*\% Xt
\dontshow{
if( max(abs(calc.X.iS.X(X, iS.X, dim) -
            (t(Xt) \%*\% iS \%*\% Xt))) > 1e-13 ){
  stop("calc.X.iS.X: Results not equal")
}
}
## Compute X* alpha
calc.mu.B(dim, X, alpha)
## or
cbind(X[[1]] \%*\% alpha[[1]], X[[2]] \%*\% alpha[[2]])
\dontshow{
if( max(abs(cbind(X[[1]] \%*\% alpha[[1]], X[[2]] \%*\% alpha[[2]]) -
            calc.mu.B(dim, X, alpha))) > 1e-13 ){
  stop("calc.mu.B: Results not equal")
}
}
}
