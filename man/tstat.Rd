\name{tstat}
\encoding{latin1}
\Rdversion{1.1}
\alias{tstat}
\alias{CVbasics}

\title{
  Basic diagnostic statistics and summaries for SpatioTemporal model output.
}
\description{
  Calculation of t-statistics and Hessian eigenvalues of model fits. 'tstat' prints out a standard R-style summary of point estimates, their SEs and t-statistics.
}
\usage{
tstat(pardat,alphas=TRUE)

CVbasics(cvout)


}
\arguments{
  \item{pardat}{List which must include a \code{'par.all'} vector and a
    square matrix \code{'hessian.all'}, of compatible sizes. Typically,
    this list would be the \code{'res.best'} component of the output of
    \code{\link{fit.mesa.model}}, when the latter is run with the option
    \code{hessian.all=TRUE}.}

  \item{alphas}{logical: should we add the point estimates as well in
    the return? Defaults to \code{TRUE}.} 
  \item{cvout}{List, the output of \code{\link{estimateCV}}, run with
    the option \code{hessian.all=TRUE}.} 
}

\details{
These functions enable (until further package upgrade) a quick, standard summary of model performance. The \code{tstat} function returns, if \code{alphas=TRUE} (default), a 3-column summary commonly used in \code{R}, with parameter names in rows, and point estimates, SEs and t-statistics in columns. Shown are estimates for both LUR and kriging parameters.

The \code{CVbasics} function is intended to diagnose cross-validation estimation. Besides point-estimates and t-statistics, the eigenvalues of the Hessian matrix (inverse of the parameter variance estimates) are also calculated and returned.}

\seealso{
 \code{\link{fit.mesa.model}}, \code{\link{estimateCV}}.
 }
 
 \author{Assaf P. Oron}

 \value{
 
 \code{tstat} returns as a default a 3-column numerical data frame, with parameter names in rows, and point estimates, SEs and t-statistics in columns. If \code{alphas=FALSE}, it only returns a vector with the t-statistics.
 
 \code{CVbasics} returns a list with 3 matrices, each of dimension (# of parameters) X (# of CV groups). The matrix names are 'alphas' with point estimates, 'tees' with t-statistics, and 'eigens' with Hessian eigenvalues. The latter is sorted largest to smallest, rather than to match any particular parameter.
 }
 
 \examples{
 
 ##load a model object
data(mesa.data.model)

##Set up initial parameter values for optimization
dimm <- loglike.dim(mesa.data.model)
x.init <- as.matrix(cbind(rep(2,dimm$nparam.cov),c(rep(c(1,-3),dimm$m+1),-3)))

################ tstat example
				
 \dontrun{
##estimate parameters
##This may take a while...
par.est <- fit.mesa.model(x.init, mesa.data.model, type="p",
      hessian.all=TRUE, control=list(trace=3,maxit=1000))
}
##Let's load precomputed results instead.
data(mesa.data.res)
par.est <- mesa.data.res$par.est

#### summary

tstat(par.est$res.best)
}
