\name{print.SVDcv}
\alias{print.SVDcv}
\title{Print details for \code{SVDcv} object}
\usage{
  \method{print}{SVDcv} (x, ...)
}
\arguments{
  \item{x}{\code{SVDcv} object to print information for.}

  \item{...}{ignored additional arguments.}
}
\value{
  Nothing
}
\description{
  \code{\link[base:print]{print}} method for class
  \code{SVDcv}, prints cross-validation statistics.
}
\examples{
## end of SVDsmooth example
}
\author{
  Johan Lindström
}
\seealso{
  Other SVDcv methods: \code{\link{plot.SVDcv}},
  \code{\link{SVDsmooth}}, \code{\link{SVDsmoothCV}}

  Other SVD for missing data:
  \code{\link{calcSmoothTrends}}, \code{\link{plot.SVDcv}},
  \code{\link{SVDmiss}}, \code{\link{SVDsmooth}},
  \code{\link{SVDsmoothCV}},
  \code{\link{updateSTdataTrend}}
}

