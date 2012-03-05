\name{block.mult}
\encoding{latin1}
\Rdversion{1.1}
\alias{block.mult}
\title{
  Multiplication of Block Diagonal Matrix and Vector
}
\description{
  Computes the matrix product between a block diagonal square matrix and
  a column vector (or matrix).
}
\usage{
block.mult(mat, X, n.x = dim(X)[2], n.blocks = 1, block.sizes =
           rep(dim(mat)[1]/n.blocks, n.blocks)) 
}
\arguments{
  \item{mat}{
    A block diagonal, square matrix.
  }
  \item{X}{
    Vector or matrix to multiply by \code{mat}.
  }
  \item{n.x}{
    Number of columns in \code{X}, defaults correctly if \code{X} is a
    matrix.
  }
  \item{n.blocks}{
    Number of diagonal blocks in \code{mat}. Defaults to 1
    (i.e. a full matrix) if not given.
  }
  \item{block.sizes}{
    A vector of length \code{n.blocks} with the size of each of the
    diagonal blocks. If not given it will assume equal-sized blocks.
  }
}
\value{
  Returns mat * X.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
Block matrices van be created by \code{\link{make.sigma.B}}, 
\code{\link{make.sigma.B.full}} and \cr
\code{\link{make.sigma.nu}}.
  
  Other block matrix functions \code{\link{makeCholBlock}}, 
  \code{\link{calc.iS.X}}, \code{\link{calc.tF.times.mat}}, \cr
  \code{\link{dot.prod}}, \code{\link{sumLogDiag}}.
  
  This function is called by \code{\link{loglike}}.
}
\examples{
#create a matrix
mat <- cbind(c(1,0,0),c(0,2,1),c(0,1,2))
#define the number of blocks and block sizes
block.sizes <- c(1,2)
n.blocks <- length(block.sizes)

#define a X vector
X <- c(1,2,3)

#compute mat \%*\% X
block.mult(mat, X, n.x = 1, n.blocks, block.sizes)
#or the old fashioned way
mat \%*\% X
\dontshow{
if( max(abs(block.mult(mat, X, n.x = 1, n.blocks, block.sizes) -
    (mat \%*\% X))) > 1e-13 ){
stop("block.mult: Results not equal")
}
}
}
