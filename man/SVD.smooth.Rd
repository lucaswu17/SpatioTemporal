\name{SVD.smooth}
\encoding{latin1}
\Rdversion{1.1}
\alias{SVD.smooth}
\title{
  Smooth Basis Functions for Data Matrix with Missing Values
}
\description{
  Function that computes smooth functions for a data matrix with missing
  values, as described in Fuentes et. al. (2006). The function uses
  \code{\link{SVD.miss}} to complete the matrix and then computes smooth
  basis functions by applying \code{\link{smooth.spline}} to the svd of
  the completed data matrix. 
}
\usage{
SVD.smooth(data, n.basis, date.ind=NA, scale.data=TRUE,
           niter=100, conv.reldiff=0.001, df=NULL, spar=NULL)
}
\arguments{
  \item{data}{
    Data matrix, with missing values marked by \code{NA}.
  }
  \item{n.basis}{
    Number of smooth basis functions to compute, will be passed as
    \code{ncomp} to \code{\link{SVD.miss}}.
  }
  \item{date.ind}{
    Vector giving the observation time of each row in \code{data}, used
    as \code{x} in \cr \code{\link{smooth.spline}} when computing the
    smooth basis functions. If set to \code{NA} it will try to coerce
    the rownames of \code{data} to dates \cr
    \code{date.ind = as.Date(rownames(data))}. \cr
    If \code{data} has no rownames it will assume equidistant data
    points, \cr
    \code{date.ind = 1:dim(data)[1]}.
  }
  \item{scale.data}{
    If \code{TRUE}, will use \code{\link{scale}} to scale the 
    data matrix before calling \code{\link{SVD.miss}}.
  }
  \item{niter, conv.reldiff}{
    Controls convergence for \code{\link{SVD.miss}}.
  }
  \item{df, spar}{
    The desired degrees of freedom/smoothing parameter for the spline,
    \cr see \code{\link{smooth.spline}}
  }
}
\value{
  The function returns a matrix where each column is a smooth basis
  function based on the svd of the completed data matrix. The left most
  column contains the smooth of the most important svd.
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
  See also \code{\link{SVD.smooth}}, \code{\link{SVD.miss}},
  \code{\link{svd}}, and \code{\link{smooth.spline}}.

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

#compute two smooth basis functions
res <- SVD.smooth(X, n.basis=2, niter=100)

#plot the two smooth basis functions
par(mfcol=c(3,2),mar=c(4,4,.5,.5))
plot(t,res[,1],ylim=range(res),type="l")
lines(t,res[,2],col=2)
#and some of the data fitted to the smooths
for(i in 1:5)\{
  plot(t,X[,i])
  lines(t,predict.lm(lm(X[,i]~res),data.frame(res)))
  lines(t,X.org[,i],col=2)
\}
}
