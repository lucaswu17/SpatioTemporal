\name{updateSTdataTrend}
\alias{updateSTdataTrend}
\title{Update Trend in \code{STdata} Object}
\usage{
  updateSTdataTrend(STdata, n.basis = 0,
    extra.dates = NULL, ...)
}
\arguments{
  \item{STdata}{A \code{STdata} object, see
  \code{\link{mesa.data}}.}

  \item{n.basis}{number of basis functions for the temporal
  trend}

  \item{extra.dates}{Additional dates for which smooth
  trends should be computed.}

  \item{...}{Additional parameters passed to
  \code{\link{calcSmoothTrends}}.}
}
\value{
  Returns a modfied version of the input, with and
  added/altered smooth trend.
}
\description{
  Updates/sets the temporal trend for \code{STdata}
  objects. It also checks that the spatio-temporal
  covariate exists for all dates in the trend (mainly an
  issue if \code{extra.dates!=NULL} adding additional times
  at which to do predictions.
}
\examples{
##load data
data(mesa.data)

##default data and time trend for one location
par(mfrow=c(4,1),mar=c(2.5,2.5,3,1))
plot(mesa.data)

##let's try with no trend
mesa.data <- updateSTdataTrend(mesa.data, n.basis=0)
plot(mesa.data)

##...and just one basis function, based on only AQS sites
subset <- mesa.data$covars$ID[mesa.data$covars$type=="AQS"]
mesa.data <- updateSTdataTrend(mesa.data, n.basis=1, subset=subset)
plot(mesa.data)

##Five basis functions, based on only AQS sites and much less smooth
mesa.data <- updateSTdataTrend(mesa.data, n.basis=5, subset=subset, df=100)
plot(mesa.data)
}
\author{
  Johan Lindström
}
\seealso{
  Other STdata functions: \code{\link{createDataMatrix}},
  \code{\link{createSTdata}}, \code{\link{createSTmodel}},
  \code{\link{c.STmodel}}, \code{\link{detrendSTdata}},
  \code{\link{removeSTcovarMean}}

  Other SVD for missing data:
  \code{\link{calcSmoothTrends}}, \code{\link{plot.SVDcv}},
  \code{\link{print.SVDcv}}, \code{\link{SVDmiss}},
  \code{\link{SVDsmooth}}, \code{\link{SVDsmoothCV}}
}

