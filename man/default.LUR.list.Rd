\name{default.LUR.list}
\encoding{latin1}
\Rdversion{1.1}
\alias{default.LUR.list}
\alias{default.ST.list}

\title{
  Error Checking and Default Values for Covariate Selection
}
\description{
  Functions that compare the covariates requested (both geographic and
  spatio-temporal) with those available in \code{\link{mesa.data}} and
  modifies the covariate request to ensure that only available
  covariates are used.
}
\usage{
default.LUR.list(mesa.data, LUR = NA)

default.ST.list(mesa.data, ST.Ind = NA)
}
\arguments{
  \item{mesa.data}{
    Data structure holding observations, observation locations, smooth
    temporal trends, geographic covariates and spatio-temporal
    covariates. See \code{\link{mesa.data}}.
  }
  \item{LUR}{
    A vector or list of vectors indicating which geographic covariates
    to use.

    If \code{LUR = NA} all covariates are used for all the temporal
    trends. If a vector of integers or characters (names) is used to
    specify covariates, these covariates will be used for all the
    temporal trends.

    A list of vectors allows for different covariates for each temporal
    trend. In this case the list needs to contain one vector for each of
    the temporal trends, starting with coefficients for the intercept.
    
    Setting any element(s) of the list to \code{NULL} implies \emph{only
      an intercept} for the corresponding temporal trend(s).
    
    If covariates are specified using names these should match \cr
    \code{names(mesa.data$LUR)}, unmatched covariates are dropped.
  }

  \item{ST.Ind}{
    A vector of integers or characters (names) indicating which of the
    \cr spatio-temporal covariates in \code{\link{mesa.data}} to
    include. \code{ST.Ind = NA} will include all covariates,
    \code{ST.Ind = NULL} gives \emph{no} covariates.

    If covariates are specified using names these should match \cr
    \code{dimnames(mesa.data$SpatioTemp)}, unmatched elements are
    dropped.
  }
}
\value{
  Returns either a \code{LUR}-list or \code{ST.Ind}-vector that refers to
  the covariate elements of \code{\link{mesa.data}} by \emph{name},
  ensuring that all elements referred to actually exist in
  \code{\link{mesa.data}}.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  Used internally by: \code{\link{create.data.model}},
  \code{\link{construct.LUR.basis}}, and \cr
  \code{\link{construct.ST.basis}}.
  
  See also \code{\link{mesa.data}} and \code{\link{create.data.model}}.

  For other data functions, see \code{\link{mesa.data.raw}},
  \code{\link{mesa.data.model}} and \cr
  \code{\link{create.data.matrix}}.
}
\examples{
##load the data
data(mesa.data)

##create a simple set of covariates
default.LUR.list(mesa.data, list(c(1:3),1,2))

##or a structure with the same covariates for all
##temporal trends
default.LUR.list(mesa.data, c(1,2))

##or a structure with only intercept for the temporal trends
default.LUR.list(mesa.data, list(c(1:3),NULL,NULL))

##Ask for covariates by name
default.LUR.list(mesa.data, list(c("log10.m.to.a1","log10.m.to.a2"),
                                 "log10.m.to.a1","log10.m.to.a1"))

##Ask for non-existent covariate by name (produces three warnings, one
##for each temporal trend)
default.LUR.list(mesa.data, "log10.m.to.a4")

##create a simple set of spatio-temporal covariates
default.ST.list(mesa.data, 1)

##or create a empty set of spatio-temporal covariates
default.ST.list(mesa.data, NULL)
}
