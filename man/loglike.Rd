\name{loglike}
\encoding{latin1}
\Rdversion{1.1}
\alias{loglike}
\alias{loglike.naive}
\title{
  Compute the Log-likelihood for the Spatio-Temporal Model
}
\description{
  Computes the log-likelihood for the spatio-temporal model.
  \code{loglike} uses an optimised version of the log-likelihood, while
  \code{loglike.naive} uses the naive (slow) version and is included
  mainly for testing and speed checks.
}
\usage{
loglike(x = NA, mesa.data.model, type = "p")

loglike.naive(x = NA, mesa.data.model, type = "p")
}
\arguments{
  \item{x}{
    Point at which to compute the log-likelihood, should be only
    \emph{log}-covariance parameters if \code{type=c("p","r")} and
    regression parameters followed by \emph{log}-covariance parameters
    if \code{type="f"}. If \code{x=NA} the function acts as an alias for
    \code{\link{loglike.var.names}} returning the expected names of the
    input parameters.
  }
  \item{mesa.data.model}{
    Data structure holding observations, and information regarding which
    \cr geographic and spatio-temporal covariates to use when fitting
    the model. See \cr \code{\link{create.data.model}} and
    \code{\link{mesa.data.model}}.
  }
  \item{type}{
    A single character indicating the type of log-likelihood to
    compute. Valid options are "f", "p", and "r", for \emph{full},
    \emph{profile} or \emph{restricted maximum likelihood}
    (REML). Computation of the \emph{full} log-likelihood requires both
    regression and \emph{log}-covariance parameters. \emph{Profile} and
    \emph{REML} require only \emph{log}-covariance parameters.
  }
}
\value{
  Returns the log-likelihood of the spatio temporal model.

  If \code{x = NA} it returns reasonable parameter names, by calling \cr
  \code{loglike.var.names(mesa.data.model, all=(type=="f"))}.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom} and Adam Szpiro
}
\section{Warning}{
  \code{loglike.naive} may take long to run. However for some
  problems with many locations and short time series \code{loglike.naive}
  can be faster than \code{loglike}.
}
\seealso{
  Expected names for \code{x} are given by
  \code{\link{loglike.var.names}}, with \code{\link{loglike.dim}} giving
  number of parameters.

  Gradients and Hessian of the log-likelihood are computed by
  \code{\link{loglike.grad}} and \cr \code{\link{loglike.hessian}}.
  
  For optimization functions see \code{\link{create.data.model}},
  \code{\link{run.MCMC}}, \code{\link{fit.mesa.model}}, and 
  \code{\link{cond.expectation}}.
}
\examples{
##load the data
data(mesa.data.model)

##Compute dimensions for the data structure
dim <- loglike.dim(mesa.data.model)

##Find out in which order parameters should be given
loglike(NA,mesa.data.model)

##Let's create random vectors of values
x <- runif(dim$nparam.cov)
x.all <- runif(dim$nparam)

##Evaluate the log-likelihood for these values
loglike(x.all, mesa.data.model, "f")
loglike(x, mesa.data.model, "p")
loglike(x, mesa.data.model, "r")
\dontshow{
##check that profile and full give the same results
data(mesa.data.res)
x.all <- mesa.data.res$par.est$res.best$par.all
x <- mesa.data.res$par.est$res.best$par
if( abs(loglike(x.all, mesa.data.model, "f") -
        loglike(x, mesa.data.model, "p")) > 1e-8 ){
  stop("loglike: full and profile not equal")
}
}%\dontshow
}%\examples
