\name{print.SVDcv}
\alias{print.SVDcv}
\alias{summary.SVDcv}
\title{Print details for \code{SVDcv} object}
\usage{
  \method{print}{SVDcv} (x, ...)

  \method{summary}{SVDcv} (object, ...)
}
\arguments{
  \item{x}{\code{SVDcv} object to print information for.}

  \item{...}{ignored additional arguments.}

  \item{object}{\code{SVDcv} object to compute summary
  for.}
}
\value{
  Nothing
}
\description{
  \code{\link[base:print]{print}} and
  \code{\link[base:summary]{summary}} methods for class
  \code{SVDcv}, prints cross-validation statistics.
}
\examples{
##See SVDsmooth example
}
\author{
  Johan Lindström
}
\seealso{
  Other SVD for missing data: \code{\link{boxplot.SVDcv}},
  \code{\link{calcSmoothTrends}}, \code{\link{plot.SVDcv}},
  \code{\link{SVDmiss}}, \code{\link{SVDsmooth}},
  \code{\link{SVDsmoothCV}},
  \code{\link{updateSTdataTrend}},
  \code{\link{updateTrend}},
  \code{\link{updateTrend.STdata}},
  \code{\link{updateTrend.STmodel}}

  Other SVDcv methods: \code{\link{boxplot.SVDcv}},
  \code{\link{plot.SVDcv}}, \code{\link{SVDsmooth}},
  \code{\link{SVDsmoothCV}}
}

