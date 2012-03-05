\name{makeCholBlock}
\encoding{latin1}
\Rdversion{1.1}
\alias{makeCholBlock}
\alias{invCholBlock}
\alias{solveTriBlock}
\title{
  Computations for Block Diagonal Matrices
}
\description{
  Provides block diagonal version of the base package functions
  \code{\link{chol}}, \code{\link{chol2inv}}, and \cr
  \code{\link{backsolve}}.
  
  Computes the Cholesky factor, the matrix inverse and solves matrix
  equation systems for block diagonal matrices.
}
\usage{
makeCholBlock(mat, n.blocks = 1, block.sizes =
                rep(dim(mat)[1]/n.blocks, n.blocks))

invCholBlock(R, n.blocks = 1, block.sizes =
               rep(dim(R)[1]/n.blocks, n.blocks))

solveTriBlock(R, B, n.x = dim(B)[2], n.blocks = 1,
                block.sizes = rep(dim(R)[1]/n.blocks, n.blocks),
                transpose = FALSE)
}
\arguments{
  \item{mat}{
    A block diagonal, square, positive definite matrix.
  }
  \item{R}{
    Upper right block diagonal Cholesky factor. The output from
    \code{\link{chol}} or \cr \code{makeCholBlock}.
  }
  \item{n.blocks}{
    Number of diagonal blocks in \code{mat} (or \code{R}). Defaults to 1
    (i.e. a full matrix) if not given.
  }
  \item{block.sizes}{
    A vector of length \code{n.blocks} with the size of each of the
    diagonal blocks. If not given it will assume equal size blocks.
  }
  \item{B}{
    Vector or matrix containg the right hand side of the equations
    system to be solved.
  }
  \item{n.x}{
    Number of columns in \code{B}, defaults correctly if \code{B} is a
    matrix.
  }
  \item{transpose}{
    Transpose \code{R} before solving the equation system. Controls if
    we solve the equations system given by R*x = B or R'*x=B.
  }
}
\details{
  \code{makeCholBlock} computes the Cholesky factor of a block
  diagonal matrix using the block diagonal structure to speed up
  computations.

  \code{invCholBlock} uses the Cholesky factor from
  \code{makeCholBlock} to compute the inverse of \code{mat}.

  \code{solveTriBlock} solves equation systems based on the Cholesky
  factor, using the block diagonal structure to speed up computations
  (c.f. \code{\link{backsolve}}). The function solves equations of
  the form R*x = B, and R'*x = B with respect to x, where the transpose
  is controlled by the parameter \code{transpose}. Aplying the function
  twice solves mat*x=B, see the examples.

  For all three functions the block diagonal structure of the matrix is
  defined by two input variables, the number of blocks \code{n.blocks},
  and the size of each block \code{block.sizes}. The size of the
  matrices must match the total number of blocks,
  i.e. \code{sum(block.sizes)} \emph{must} equal \code{dim(mat)}.

  The functions can be  used for full matrices by setting the number of
  blocks to 1. 
}
\value{
  \code{makeCholBlock} gives the Cholesky factor and
  \code{invCholBlock} gives the inverse of the matrix \code{mat}.
  \code{solveTriBlock} gives to answer to the equation system.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom} and Adam Szpiro
}
\seealso{
  Other block matrix functions \code{\link{dot.prod}}, 
  \code{\link{block.mult}}, \code{\link{calc.tF.times.mat}}, \cr
  \code{\link{calc.iS.X}}, and \code{\link{sumLogDiag}}.
  
  This function is used by \code{\link{loglike}}.
}
\examples{
#create a matrix
mat <- cbind(c(1,0,0),c(0,2,1),c(0,1,2))
#define the number of blocks and block sizes
block.sizes <- c(1,2)
n.blocks <- length(block.sizes)

#Compute the Cholesky factor
R <- makeCholBlock(mat, n.blocks, block.sizes)
#and the matrix inverse
i.mat <- invCholBlock(R, n.blocks, block.sizes)
#compare to the alternative
i.mat-solve(mat)
\dontshow{
if( max(abs(R-chol(mat))) > 1e-10 ){
  stop("makeCholBlock: Results not equal")
}
if( max(abs(i.mat-solve(mat))) > 1e-10 ){
  stop("invCholBlock: Results not equal")
}
}
#define a B vector
B <- c(1,2,3)
#solve the equation system (we need n.x since B is not a matrix)
x1 <- solveTriBlock(R, B, n.x=1, n.blocks, block.sizes, tr=TRUE)
x2 <- solveTriBlock(R, x1, n.x=1, n.blocks, block.sizes, tr=FALSE)
print(x2)
#compare to the alternative
print(solve(mat,B))
range(x2-solve(mat,B))
\dontshow{
if( max(abs(x2-solve(mat,B))) > 1e-10 ){
  stop("solveTriBlock: Results not equal")
}
}
#compute the quadratic form B'*i.mat*B
norm2(x1)
#compare to the alternative
t(B) \%*\% i.mat \%*\% B
\dontshow{
if( abs(norm2(x1) - t(B) \%*\% i.mat \%*\% B) > 1e-10 ){
  stop("solveTriBlock: Results not equal for norm")
}
}%\dontshow
}%\examples
