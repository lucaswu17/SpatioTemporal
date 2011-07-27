\name{loglike.grad}
\encoding{latin1}
\Rdversion{1.1}
\alias{loglike.grad}
\alias{loglike.hessian}
\alias{loglike.naive.grad}
\alias{loglike.naive.hessian}
\title{
  Compute Gradient and Hessian for the Log-likelihood
}
\description{
  Computes finite difference gradients and hessians for the
  log-likelihood functions \code{\link{loglike}} and
  \code{\link{loglike.naive}}.
}
\usage{
loglike.grad(x, mesa.data.model, type = "p", h = 0.001,
             diff.type = 0)

loglike.naive.grad(x, mesa.data.model, type = "p", h = 0.001,
                   diff.type = 0) 

loglike.hessian(x, mesa.data.model, type = "p", h = 0.001)

loglike.naive.hessian(x, mesa.data.model, type = "p", h = 0.001)
}
\arguments{
  \item{x}{
    Point at which to compute the gradient or hessian. See
    \code{\link{loglike}}.
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
  \item{h, diff.type}{
    Step length and type of finite difference to use when computing
    gradients; positive values of \code{diff.type} gives forward
    differences, \code{0} gives central differences, and negative values
    gives backward differences. See \code{\link{gen.gradient}} for
    details.
  }
}
\details{
  Uses \code{\link{gen.gradient}} and \code{\link{gen.hessian}} to
  compute finite difference derivatives of the log-likelihood function
  in \code{\link{loglike}} and \code{\link{loglike.naive}}.
  
  Used by the model fitting function \code{\link{fit.mesa.model}} and
  provided for users who want to implement their own model fitting.
}
\value{
  Returns the gradient (as a vector) or hessian matrix for the
  \code{\link{loglike}} and \code{\link{loglike.naive}} functions.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\section{Warning}{\code{loglike.naive.grad} and
  \code{loglike.naive.hessian} may take \strong{very} long time to run,
  use with extreme care.
}
\seealso{
  Computes gradients and hessians for the log-likelihood function
  \code{\link{loglike}}.

  Expected names for \code{x} are given by
  \code{\link{loglike.var.names}}.
  
  Used by the estimation functions \code{\link{fit.mesa.model}} and
  \code{\link{run.MCMC}}.

  For general computation of gradient and hessian see
  \code{\link{gen.gradient}} and \code{\link{gen.hessian}}.
  
  For further log-likelihood computations see \code{\link{loglike}}, 
  \code{\link{loglike.dim}}, and \cr
  \code{\link{loglike.var.names}}.
}
\examples{
\dontrun{
##load the data
data(mesa.data.model)

##Compute dimensions for the data structure
dim <- loglike.dim(mesa.data.model)

##Let's create random vectors of values
x <- runif(dim$nparam.cov)
x.all <- runif(dim$nparam)

##Compute the gradients
Gf <- loglike.grad(x.all, mesa.data.model, "f")
Gp <- loglike.grad(x, mesa.data.model, "p")
Gr <- loglike.grad(x, mesa.data.model, "r")

##And the Hessian, this may take some time...
Hf <- loglike.hessian(x.all, mesa.data.model, "f")
Hp <- loglike.hessian(x, mesa.data.model, "p")
Hr <- loglike.hessian(x, mesa.data.model, "r")
}%\dontrun
}%\examples
