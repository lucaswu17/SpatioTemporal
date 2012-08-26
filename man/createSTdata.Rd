\name{createSTdata}
\alias{createSTdata}
\title{Construct STdata Object}
\usage{
  createSTdata(obs, covars, SpatioTemporal = NULL,
    transform.obs = function(x) {     return(x) },
    mean.0.ST = FALSE, n.basis = NULL, extra.dates = NULL,
    ..., detrend = FALSE, region = NULL, method = NULL)
}
\arguments{
  \item{obs}{Either a data.frame with fields \code{obs},
  \code{date}, \code{ID} giving obsevations, time-points
  and location names; or a matrix, e.g. output from
  \code{\link{createDataMatrix}}.}

  \item{covars}{matrix/data.frame of covariates; should
  include both geographic covariates and coordinates of all
  locations, see \code{\link{stCheckCovars}}.}

  \item{SpatioTemporal}{possible spatio-temporal covariate,
  see \code{\link{stCheckSTcovars}}.}

  \item{transform.obs}{function to apply to the
  observations, defaults to an identity transform. Possible
  options are \code{\link[base:log]{log}},
  \code{\link[base:sqrt]{sqrt}}, and
  \code{\link[base:exp]{exp}}.}

  \item{mean.0.ST}{Call \code{\link{removeSTcovarMean}} to
  produce a mean-zero spatio-temporal covariate?}

  \item{n.basis}{Number of temporal components in the
  smooth trends computed by
  \code{\link{updateSTdataTrend}}, if \code{NULL} no trend
  is computed (implies only a constant).}

  \item{extra.dates}{Additional dates for which smooth
  trends should be computed, used by
  \code{\link{updateSTdataTrend}}.}

  \item{...}{Additional parameters passed to
  \code{\link{updateSTdataTrend}}.}

  \item{detrend}{Use \code{\link{detrendSTdata}} to remove
  a termporal trend from the observations; requires
  \code{n.basis!=NULL}.}

  \item{region,method}{Additional parameters passed to
  \code{\link{detrendSTdata}}.}
}
\value{
  A \code{STdata} object, see \code{\link{mesa.data}} for
  an example.
}
\description{
  Creates a \code{STdata} object that can be used as input
  for \code{\link{createSTmodel}}. Names and dates are
  derived from the input data, either using predefined
  fields or \code{rownames} / \code{colnames}; for details
  see the sub-functions linked under the relevant
  Arguments.
}
\examples{
##load the raw data
data(mesa.data.raw)

##extract observations and covariates
obs <- mesa.data.raw$obs
covars <- mesa.data.raw$X

##list with the spatio-temporal covariates
ST.list <- list(lax.conc.1500=mesa.data.raw$lax.conc.1500)

##create STdata object
mesa.data <- createSTdata(obs, covars, SpatioTemporal=ST.list)
print(mesa.data)

##create object with mean 0 spatio temporal covariate
mesa.data.2 <- createSTdata(obs, covars, SpatioTemporal=ST.list,
                            mean.0.ST=TRUE)
print(mesa.data.2)

##create object with mean 0 spatio temporal covariate, and
##trend with two components, and additional dates (every seventh day)
extra.dates <- seq(min(as.Date(rownames(obs))),
                   max(as.Date(rownames(obs))), by=7)
mesa.data.3 <- createSTdata(obs, covars, n.basis=2, extra.dates=extra.dates)
print(mesa.data.3)
}
\author{
  Johan Lindström and Assaf P. Oron
}
\seealso{
  Other STdata functions: \code{\link{c.STmodel}},
  \code{\link{createDataMatrix}},
  \code{\link{createSTmodel}}, \code{\link{detrendSTdata}},
  \code{\link{removeSTcovarMean}},
  \code{\link{updateSTdataTrend}}

  Other STdata methods: \code{\link{plot.STdata}},
  \code{\link{plot.STmodel}}, \code{\link{print.STdata}},
  \code{\link{print.summary.STdata}},
  \code{\link{summary.STdata}}
}

