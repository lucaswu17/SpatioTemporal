\name{SVD.smooth.cv}
\encoding{latin1}
\Rdversion{1.1}
\alias{SVD.smooth.cv}
\title{
  Cross-validation for Smooth Basis Functions
}
\description{
  Function that uses cross-validation to evaluate the number of smooth
  functions needed to describe a data matrix with missing values
  (c.f. Fuentes et. al. (2006). Can be used to determine reasonable
  values for \code{n.basis} in \code{\link{SVD.smooth}}.

  The function uses leave-one-column-out cross-validation; holding one
  column out from \code{data}, calling \code{\link{SVD.smooth}}, and
  then regressing the held out column on the resulting smooth
  functions. Cross-validation statistics computed include RMSE,
  R-squared and BIC.
}
\usage{
SVD.smooth.cv(data, n.basis, date.ind=NA, scale.data=TRUE,
              niter=100, conv.reldiff=0.001, df=NULL, spar=NULL)
}
\arguments{
  \item{data}{
    Data matrix, with missing values marked by \code{NA}.
  }
  \item{n.basis}{
    A vector with the number of basis functions for which to run the
    \cr cross-validation.
  }
  \item{date.ind}{
    Vector giving the observation time of each row in \code{data}.
    Same as \code{date.ind} in \code{\link{SVD.smooth}}.
  }
  \item{scale.data}{
    Uses \code{\link{scale}} to scale the data matrix before calling
    \code{\link{SVD.miss}}.
  }
  \item{niter, conv.reldiff}{
    Controls convergence for \code{\link{SVD.miss}}.
  }
  \item{df, spar}{
    The desired degrees of freedom/smoothing parameter for the spline,
    \cr see \code{\link{smooth.spline}}.
  }
}
\value{
  Returns a list with the following components:
\item{CV.stat}{A data.frame with cross-validation statistics for each of
  the number of basis functions evaluated. Contains \code{RMSE},
  \code{R2} (R-squared) and \code{BIC} (Bayesian information
  criterion).}
\item{BIC.all}{A data.frame with the individual BIC values for each
  column in the data matrix and for each number of basis functions
  evaluated. Can be used to evaluate how many of the columns would
  benefit from using more/fewer basis functions.}
\item{smooth.SVD}{A list with \code{length(n.basis)} components. Each
  component contains an array where \code{smooth.SVD[[j]][,,i]} is the
  result of \code{\link{SVD.smooth}} applied to \code{data[,-i]} with
  \code{n.basis[j]} smooth functions. Can be used to investigate how
  different the smooth functions vary when different columns of the data
  matrix are dropped.}
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
X <- matrix(cbind(cos(t),sin(2*t)),length(t),2) \%*\%
  matrix(rnorm(20),2,10)
#add some normal errors
X <- X + .25*rnorm(length(X))
#and mark some data as missing
X[runif(length(X))<.25] <- NA

#compute cross-validation for 1 to 4 basis functions
res.cv <- SVD.smooth.cv(X, n.basis=1:4, niter=100)

#plot cross-validation statistics
par(mfcol=c(2,2),mar=c(4,4,.5,.5))
plot(res.cv$CV.stat$RMSE,type="l",ylab="RMSE")
plot(res.cv$CV.stat$R2,type="l",ylab="R2")
plot(res.cv$CV.stat$BIC,type="l",ylab="BIC")

#plot the BIC for each column, illustrates how many columns that
#improve and how many that worsen when increasing the number of
#basis functions.
pairs(res.cv$BIC.all,panel=function(x,y)\{points(x,y);abline(0,1)\})
}
