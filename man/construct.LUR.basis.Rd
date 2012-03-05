\name{construct.LUR.basis}
\encoding{latin1}
\Rdversion{1.1}
\alias{construct.LUR.basis}
\alias{construct.ST.basis}

\title{
  Extract and Create Covariate Matrices
}
\description{
  Extracts the requested geographic and spatio-temporal covariates from
  \code{\link{mesa.data}} and formats them into suitable matrices.
}
\usage{
construct.LUR.basis(mesa.data, LUR)

construct.ST.basis(mesa.data, ST.Ind)
}
\arguments{
  \item{mesa.data}{
    Data structure holding observations, observation locations, smooth
    temporal trends, geographic covariates and spatio-temporal
    covariates. See \code{\link{mesa.data}}.
  }
  \item{LUR}{
    A vector or list of vectors indicating which geographic covariates
    to use, see \code{\link{default.LUR.list}}.
  }
  \item{ST.Ind}{
    A vector indicating which spatio-temporal covariates to use, see \cr
    \code{\link{default.ST.list}}.
  }
}
\value{
  \code{construct.LUR.basis} returns a list of matrices, where each
  matrix contains geographic covariates for the corresponding temporal
  trend, starting with the intercept.

  \code{construct.ST.basis} returns a matrix with one column for each
  spatio-temporal covariate. Each row contains the covariate(s) for the
  corresponding observation in \code{mesa.data$obs$obs}.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  Used internally by: \code{\link{create.data.model}}.
  
  See \code{\link{default.LUR.list}} and \code{\link{default.ST.list}}
  for details on the variable specification.

For other data functions, see \code{\link{mesa.data.raw}},
\code{\link{mesa.data.model}}, \code{\link{default.LUR.list}}, \cr
and \code{\link{create.data.matrix}}.
}
\examples{
##load the data
data(mesa.data)

##create a simple set of covariates
LUR.1 <- construct.LUR.basis(mesa.data, list(c(1:3),1,2))
str(LUR.1)

##or a structure with the same covariates for all
##temporal trends
LUR.2 <- construct.LUR.basis(mesa.data, c(1:2))
str(LUR.2)

##or a structure with only intercept for the temporal trends
LUR.3 <- construct.LUR.basis(mesa.data, list(c(1:3),NULL,NULL))
str(LUR.3)

##Ask for covariates by name
LUR.4 <- construct.LUR.basis(mesa.data, "log10.m.to.a1")
str(LUR.4)

##Ask for non-existent covariate by name (produces a warning)
LUR.5 <- construct.LUR.basis(mesa.data, "log10.m.to.a4")
str(LUR.5)

##create a simple set of spatio-temporal covariates
ST.1 <- construct.ST.basis(mesa.data, 1)
str(ST.1)

##or create a empty set of spatio-temporal covariates
ST.2 <- construct.ST.basis(mesa.data, NULL)
str(ST.2)
}
