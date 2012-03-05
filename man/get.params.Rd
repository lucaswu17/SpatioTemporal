\name{get.params}
\encoding{latin1}
\Rdversion{1.1}
\alias{get.params}

\title{
  Extract Parameters from a Vector
}
\description{
  Extracts parameters for the log-likelihood from a parameter vector and
  separates regression parameters and \emph{log}-covariance parameters.
}
\usage{
get.params(x, dimensions)
}
\arguments{
  \item{x}{
    A vector containing regression and \emph{log}-covariance parameters.
    The ordering of the elements should be \strong{exactly} that
    indicated by \code{\link{loglike.var.names}}.
  }
  \item{dimensions}{
    Information regarding parameter dimension, output from
    \code{\link{loglike.dim}}.
  }
}
\value{
  Returns a list containing:
  \item{gamma}{A column matrix of regression coefficients for the
    spatio-temporal covariate(s).
  }
  \item{alpha}{A list of column matrices. Each column matrix contains
    the regression coefficients for geographic covariates in the
    corresponding beta field.
  }
  \item{range, sill}{Vectors with the range and sill for each of the
    beta fields. These are back transformed assuming that \code{x}
    contained log range and sill.
  }
  \item{phi.nu}{Range, sill and nugget for the residual fields. Back
    transformed assuming that \code{x} contained log parameters.
  }
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  Used internally by \code{\link{loglike}} and
  \code{\link{simulateMesaData}}.
}
\examples{
##load the data
data(mesa.data.model)

##Compute dimensions for the data structure
dim <- loglike.dim(mesa.data.model)

##Let's create random parameter vectors ...
x <- runif(dim$nparam.cov)
x.all <- runif(dim$nparam)

##... and pick them apart
get.params(x, dim)
get.params(x.all, dim)
}
