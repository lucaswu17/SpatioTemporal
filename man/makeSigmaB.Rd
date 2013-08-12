\name{makeSigmaB}
\alias{makeSigmaB}
\title{Create Block Covariance Matrix (Equal Block Sizes)}
\usage{
  makeSigmaB(pars, dist, type = "exp", nugget = 0,
    symmetry = dim(dist)[1] == dim(dist)[2],
    ind2.to.1 = 1:dim(dist)[2], sparse = FALSE)
}
\arguments{
  \item{pars}{List of parameters for each block; if not a
  list a single block matrix is assumed. Should match
  parameters suggested by \code{\link{parsCovFuns}}.}

  \item{dist}{Distance matrix.}

  \item{type}{Name(s) of covariance functions, see
  \code{\link{namesCovFuns}}.}

  \item{nugget}{Vector of nugget(s) to add to the diagonal
  of each matrix.}

  \item{symmetry}{\code{TRUE}/\code{FALSE} flag if the
  \code{dist} is symmetric, resulting in a symmetric
  covariance matrix.}

  \item{ind2.to.1}{Vectors, that for each index along the
  second dimension gives a first dimension index, used only
  if \code{symmetry=FALSE} to determine which covariances
  should have an added nugget (collocated sites).}

  \item{sparse}{If \code{TRUE}, return a block diagonal
  sparse matrix, see \code{\link[Matrix:bdiag]{bdiag}}.}
}
\value{
  Block covariance matrix of size
  \code{dim(dist)*n.blocks}.
}
\description{
  Function that creates a block covariance matrix with
  equal sized blocks. Used to construct the Sigma_B matrix.
}
\details{
  Any parameters given as scalars will be \code{rep}-ed to
  match \code{length(pars)}.
}
\examples{
##First create some random locations
x <- rnorm(5)
y <- rnorm(5)

##compute distance matrix
D <- crossDist( cbind(x,y) )

##create a block diagonal matrix exponential covariance matrix
##with different range, sill, and nugget
pars <- list(c(.3,2), c(2,1), c(1,3))
nugget <- c(.5,0,1)

Sigma1 <- makeSigmaB(pars, D, type="exp", nugget=nugget)
##or using different covariance functions for each block
Sigma2 <- makeSigmaB(pars, D, type=c("exp","exp2","cubic"),
                     nugget=nugget)

##make a cross-covariance matrix
Dcross <- D[1:3,c(1,1,2,2)]
Sigma.cross <- makeSigmaB(pars, Dcross, type="exp", nugget=nugget,
                          ind2.to.1=c(1,1,2,2))

\dontshow{
  Sigma.alt <- matrix(0, length(pars)*dim(D)[1], length(pars)*dim(D)[1])
  for(i in 1:length(pars)){
    Ind <- (1:dim(D)[1]) + (i-1)*dim(D)[1]
    Sigma.alt[Ind, Ind] <- pars[[i]][2]*exp(-D/pars[[i]][1])
    diag(Sigma.alt[Ind, Ind]) <- diag(Sigma.alt[Ind, Ind])+nugget[i]
  }
  if( abs(max(Sigma1-Sigma.alt)) > 1e-13){
    stop("makeSigmaB: Results not equal, covariance")
  }
  Ind <- c(1,1,2,2)
  Sigma.alt.cross <- Sigma.alt[c(1:3,6:8,11:13),c(Ind, 5+Ind, 10+Ind)]
  if( abs(max(Sigma.cross-Sigma.alt.cross)) > 1e-13){
    stop("makeSigmaB: Results not equal, cross-covariance")
  }

  Sigma2.s <- makeSigmaB(pars, D, type=c("exp","exp2","cubic"),
                         nugget=nugget, sparse=TRUE)
  if( max(abs(Sigma2.s-Sigma2)) > 1e-14 ){
    stop("makeSigmaB: Sparse matrix not equal, covariance")
  }
  Sigma.cross.s <- makeSigmaB(pars, Dcross, type="exp", nugget=nugget,
                              ind2.to.1=c(1,1,2,2), sparse=TRUE)
  if( max(abs(Sigma.cross.s-Sigma.cross)) > 1e-14 ){
    stop("makeSigmaB: Sparse matrix not equal, cross-covariance")
  }
}

}
\author{
  Johan Lindström
}
\seealso{
  Other block matrix functions: \code{\link{blockMult}},
  \code{\link{calc.FX}}, \code{\link{calc.FXtF2}},
  \code{\link{calc.iS.X}}, \code{\link{calc.mu.B}},
  \code{\link{calc.tFX}}, \code{\link{calc.tFXF}},
  \code{\link{calc.X.iS.X}}, \code{\link{invCholBlock}},
  \code{\link{makeCholBlock}}, \code{\link{makeSigmaNu}},
  \code{\link{solveTriBlock}}

  Other covariance functions: \code{\link{crossDist}},
  \code{\link{evalCovFuns}}, \code{\link{makeSigmaNu}},
  \code{\link{namesCovFuns}}, \code{\link{parsCovFuns}},
  \code{\link{updateCovf}}
}

