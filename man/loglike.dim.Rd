\name{loglike.dim}
\encoding{latin1}
\Rdversion{1.1}
\alias{loglike.dim}

\title{
  Dimensions of the Model Data Structure
}
\description{
  Function that computes the dimension of several objects in the
  \code{\link{mesa.data.model}} data structure.
}
\usage{
loglike.dim(mesa.data.model)
}
\arguments{
  \item{mesa.data.model}{
    Data structure holding observations, and information regarding which
    \cr geographic and spatio-temporal covariates to use when fitting
    the model. See \cr \code{\link{create.data.model}} and
    \code{\link{mesa.data.model}}.
  }
}
\value{
  Returns a list with the following components:
  \item{T}{Number of distinct time points in the data}
  \item{m}{Number of temporal basis functions including the intercept.}
  \item{n}{Number of distinct locations in the data}
  \item{p}{A vector with \code{m} elements, containg the number of
    geographic covariates \cr (LUR:s) for each of the temporal basis
    functions.}
  \item{L}{Number of spatio-temporal covariates}
  \item{nparam}{Total number of parameters, including covariance and
    regression parameters, (\code{L+sum(p)+2*m+3})}
  \item{nparam.cov}{Number of covariance parameters (\code{2*m+3})}
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  See \code{\link{mesa.data.model}} for details on the data
  structure. Used internally by several functions in the package.
  
  For further log-likelihood computations see \code{\link{loglike}}, 
  \code{\link{loglike.grad}}, and \cr \code{\link{loglike.var.names}}.
}
\examples{
##load the data
data(mesa.data.model)

##compute dimensions for the data structure
loglike.dim(mesa.data.model)
}
