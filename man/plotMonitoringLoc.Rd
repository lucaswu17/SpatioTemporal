\name{plotMonitoringLoc}
\encoding{latin1}
\Rdversion{1.1}
\alias{plotMonitoringLoc}

\title{
  Plot of Observation Locations and Dates
}
\description{
  Plots observation locations (or values) against the date of each
  observation. 
}
\usage{
plotMonitoringLoc(mesa.data, main = "", col = NULL, legend.loc =
                  "topleft", legend.names = NULL, add = FALSE,
                  obsplot=FALSE, pch=19, cex=.1, ...)
}
\arguments{
  \item{mesa.data}{
    Data structure holding observations, as well as observation
    locations and dates of observations. See \code{\link{mesa.data}} and
    \code{\link{mesa.data.model}}.
  }
  \item{main}{
    A character string to be used as the main title of the plot.
  }
  \item{col}{
    A vector containing the colors to be used for each monitoring type
    in the plot.
  }
  \item{legend.loc}{
    The location of the legend.  See \code{\link{legend}}.
  }
  \item{legend.names}{
    A vector of character strings to be used in the legend.
  }
  \item{add}{
    If \code{TRUE} adds to an existing plot.
  }
  \item{obsplot}{
    If \code{TRUE} the function will plot a time trend of all raw
    observations, color-coded by type. Otherwise (default), location ID
    is plotted on the y-axis.
  }
  \item{cex, pch}{
    Parameters controlling the plotting symbol; passed on to
    \code{\link{points}}. Default is a small filled circle.
  }
  \item{...}{
    Additional parameters passed to either \code{\link{plot}} or
    \code{\link{points}}, depending on the value of \code{add}.
  }
}
\details{
  This function creates a plot where each of the monitoring locations
  (or observations if \cr \code{obsplot=TRUE}) are on the vertical axis,
  and each of the observation dates are on the horizontal axis.
  If \code{obsplot=TRUE}, the y-axis will show the raw observations
  instead of location ID, and the plot will be a composite time-trend
  of all observations (potentially) color coded by type.

  If \code{mesa.data$location$type} contains location types, the monitoring
  locations are colour-coded by type; if \code{mesa.data} does not
  contain monitoring types, a single type for all locations is assumed.
  The colours to use are specified in \code{col}. By default the
  function uses \code{1:n}, where \code{n} is the number of unique
  monitoring types.   and \code{col} must be length 1.

  The vector \code{legend.names} should have length equal to the number
  of unique location types.  The default legend names are the unique
  values found in \code{mesa.data$location$type}, see
  \code{\link{mesa.data}}. If \code{mesa.data$location} does not contain
  location types, \code{legend.names} must be of length 1 and the
  defaults to \code{"Observations"}. 
}
\value{
  Does not return anything.
}
\author{
  \enc{Johan Lindström}{Johan Lindstrom}
}
\seealso{
  \code{\link{plotMesaData}}, \code{\link{printMesaDataNbrObs}} for further 
  data summaries.
}
\examples{
##load data
data(mesa.data)

##standard plot
plotMonitoringLoc(mesa.data)

##different names/colours/etc
plotMonitoringLoc(mesa.data, main="A nice plot", col=c("green","blue"),
    legend.names=c("Sites of one type", "..and of the other"),
    legend.loc="bottomleft")

##composite time-trend
plotMonitoringLoc(mesa.data, obsplot=TRUE)
}
