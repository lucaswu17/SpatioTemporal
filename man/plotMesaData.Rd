\name{plotMesaData}
\encoding{latin1}
\Rdversion{1.1}
\alias{plotMesaData}

\title{
  Illustrates the Observations as Timeseries
}
\description{
  Plots either the time-series of observations or residuals of a
  regression of observations on the temporal trends in
  \code{mesa.data$trend}; for the residuals this plots either the raw
  residuals or the (partial) autocorrelation function.
}
\usage{
plotMesaData(mesa.data, ID, add = FALSE, type = "obs",
             lag.max = NULL)
}
\arguments{
  \item{mesa.data}{
    A data structure containing the observations, either
    \code{\link{mesa.data}} or \cr \code{\link{mesa.data.model}}.
  }
  \item{ID}{
    The location for which we want to plot observations.
    Either a string matching the names in \code{mesa.data$location$ID}
    or an integer; if an integer the functions will plot data from 
    \code{ID=mesa.data.model$location$ID[ID]}.
  }
  \item{add}{
    Add to an existing plot, or start a new plot. Ignored if
    \code{type="acf"} or \cr \code{"pacf"}.
  }
  \item{type}{
    A string indicating which type of plot to construct. Valid options
    are \code{"obs"}, \code{"res"}, \code{"acf"}, and \code{"pacf"}.
  }
  \item{lag.max}{
    Maximum lag for which to calculate the (p)acf, see \code{\link{acf}}
    for details.
  }
}
\details{
  The function plots either timeseries of observations at the location
  specified by \code{ID}, or residuals from the regression of
  observations on the smooth temporal trends in \code{mesa.data$trend}.

  When plotting observations (\code{type="obs"}) this plots both
  observations and the best fit of the trend functions to the
  observations. \code{type="res"} plots the  residuals of observations
  fitted to the temporal trends; this also adds the \code{y=0} line for
  reference. \code{type="acf"} or \code{type="pacf"} plots the
  autocorrelation function or partial autocorrelation function of the
  residuals.
}
\value{
  Does not return anything.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
 \code{\link{plotMonitoringLoc}}, \code{\link{printMesaDataNbrObs}} 
 for additional data summaries.
}
\examples{
##load data
data(mesa.data)

##plot the first site,
par(mfrow=c(4,1),mar=c(2.5,2.5,3,1))
plotMesaData(mesa.data, ID=1, type="obs")
##residuals from the temporal trends,
plotMesaData(mesa.data, ID=1, type="res")
##afc 
plotMesaData(mesa.data, ID=1, type="acf")
##... and pafc for the residuals
plotMesaData(mesa.data, ID=1, type="pacf")

##Same as above but calling the 2nd site by name
par(mfrow=c(4,1),mar=c(2.5,2.5,3,1))
plotMesaData(mesa.data, ID="60370016", type="obs")
plotMesaData(mesa.data, ID="60370016", type="res")
plotMesaData(mesa.data, ID="60370016", type="acf")
plotMesaData(mesa.data, ID="60370016", type="pacf")
}
