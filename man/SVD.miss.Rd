\name{SVD.miss}
\encoding{latin1}
\Rdversion{1.1}
\alias{SVD.miss}
\title{
  Replace Missing Values in a Data Matrix Using Iterative svd.
}
\description{
  Function that completes a data matrix using iterative svd as described
  in Fuentes et. al. (2006). The
  function iterates between computing the svd for the matrix and
  replacing the missing values by linear regression of the columns onto
  the first \code{ncomp} svd components. As initial replacement for the 
  missing values regression on the column averages are used. 
  The function \emph{will fail} if entire rows and/or columns are missing 
  from the data matrix.
}
\usage{
SVD.miss(X, niter = 25, ncomp = 4, conv.reldiff = 0.001)
}
\arguments{
  \item{X}{
    Data matrix, with missing values marked by \code{NA}.
  }
  \item{niter}{
    Maximum number of iterations to run before exiting, \code{Inf} will
    run until the \code{conv.reldiff} criteria is fulfilled.
  }
  \item{ncomp}{
    Number of svd components to use in the reconstruction.
  }
  \item{conv.reldiff}{
    Assume the iterative procedure has converged when the relative
    difference between two consecutive iterations is less than
    \code{conv.reldiff}.
  }
}
\value{
  Returns a list with the following components
\item{Xfill}{The completed data matrix with missing values replaced
  by fitting the data to the \code{ncomp} most important svd components}
\item{svd}{The result of svd on the completed data matrix,
  i.e. \code{svd(Xfill)}}
\item{status}{A vector of status variables: \code{diff} contains the absolute
  difference between the two last iterations, \code{rel.diff} contains the
  relative difference, and the number of iterations
  run are contained in \code{n.iter}}
}
\references{
  M. Fuentes, P. Guttorp, and P. D. Sampson. (2006) Using Transforms to
  Analyze Space-Time Processes in Statistical methods for
  spatio-temporal systems (B. \enc{Finkenstädt}{Finkenstadt}, L. Held,
  V. Isham eds.) 77-150
}
\author{
  Paul D. Sampson and \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  See also \code{\link{SVD.smooth}}, \code{\link{SVD.smooth.cv}}, and
  \code{\link{svd}}.

  \code{\link{create.data.matrix}} can be used to create a data matrix
  from a \code{\link{mesa.data}} object.

  Used by \code{\link{calc.smooth.trends}} to compute smooth trends for
  \code{\link{mesa.data}}.
}
\examples{
#create a data matrix
t <- seq(0,4*pi,len=50)
X.org <- matrix(cbind(cos(t),sin(2*t)),length(t),2) \%*\%
  matrix(rnorm(20),2,10)
#add some normal errors
X <- X.org + .25*rnorm(length(X.org))
#and mark some data as missing
X[runif(length(X))<.25] <- NA

#run the missing data svd
res <- SVD.miss(X, niter=100, ncomp=2)
#look at the status
res$status
#plot the first column of the data matrix
plot(t,X[,1])
lines(t,res$Xfill[,1])
lines(t,X.org[,1],col=2)
\dontshow{
if( max(abs(res$Xfill-X),na.rm=TRUE) > 1e-10 ){
  stop("SVD.miss: Points not interpolated")
}
}
}
